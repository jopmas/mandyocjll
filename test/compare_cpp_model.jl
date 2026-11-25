# Runs a case of Mandyoc (folder with param.txt, interfaces.txt, inputs and binary outputs of the C++ code) and compares.
# julia --project=. test/mpirun.jl N compare_cpp_model.jl <folder> <steps> [add_remove]
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/stokes.jl"); include("../src/rheology.jl"); include("../src/reader.jl"); include("../src/particle_control.jl")
using Printf
function read_vec(file)
    open(file) do io
        ntoh(read(io, Int32)); n = ntoh(read(io, Int32)); ntoh.(read!(io, zeros(Float64, n)))
    end
end
function main()
    PETSc.initialize(P); comm = MPI.COMM_WORLD; rank = MPI.Comm_rank(comm)
    dir, steps = ARGS[1], parse(Int, ARGS[2]); add_remove = "add_remove" in ARGS
    d = read_param_file(joinpath(dir, "param.txt"))
    m = ModelParameters(d; fssa_legacy = true)
    layers = read_interfaces(joinpath(dir, "interfaces.txt"), m)
    on(k) = lowercase(get(d, k, "false")) == "true"
    T0 = on("temperature_from_ascii") ? read_ascii_field(joinpath(dir, "input_temperature_0.txt")) : zeros(m.Nx * m.Nz)
    V0 = on("velocity_from_ascii") ? read_ascii_field(joinpath(dir, "input_velocity_0.txt")) : zeros(2m.Nx * m.Nz)
    s = Model(m, layers; temperature = (i, k) -> T0[k*m.Nx+i+1], velocity = (i, k) -> (V0[2(k*m.Nx+i)+1], V0[2(k*m.Nx+i)+2]))
    function cmp(vec, da, dof, name, step; interior = false, scale = identity)
        file = joinpath(dir, "$(name)_$step.bin"); isfile(file) || return NaN
        ref = read_vec(file); si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da)); n2 = [0.0, 0.0]; mx = [0.0, 0.0]
        PETSc.with_local_array!(vec; write = false) do a
            A = reshape(a, dof, mi, mk)
            for k in 1:mk, i in 1:mi, c in 1:dof
                interior && (si + i == m.Nx || sk + k == m.Nz) && continue
                r = ref[dof*((sk+k-1)*m.Nx+si+i-1)+c]; e = scale(A[c, i, k]) - scale(r)
                n2 .+= (e^2, scale(r)^2); mx .= max.(mx, (abs(e), abs(scale(r))))
            end
        end
        n2 = MPI.Allreduce(n2, +, comm)
        sqrt(n2[1] / n2[2])
    end
    function report(step)
        ev = cmp(s.velocity_next, s.da_Veloc, 2, "velocity", step)
        et = cmp(s.temperature, s.da_Thermal, 1, "temperature", step)
        eta = cmp(s.fields.viscosity_factor, s.da_Thermal, 1, "viscosity", step; interior = true, scale = x -> log10(max(x, 1.0)))
        er = cmp(s.fields.strain_rate_field, s.da_Thermal, 1, "strain_rate", step)
        es = cmp(s.fields.strain_field, s.da_Thermal, 1, "strain", step)
        ep = cmp(s.pressure_thermal_mesh, s.da_Thermal, 1, "pressure", step)
        eh = cmp(s.fields.density_field, s.da_Thermal, 1, "density", step)
        rank == 0 && @printf("RESULT np=%d step %d: relative error (L2) vel %.2e | T %.2e | log10 visc %.2e | strain rate %.2e | strain %.2e | pressure %.2e | density %.2e | nonlinear it. %d, dt %.4g yr\n",
            MPI.Comm_size(comm), step, ev, et, eta, er, es, ep, eh, s.nonlinear_iterations, s.time_step_seconds / SECONDS_PER_YEAR)
    end
    e0 = cmp(s.pressure_thermal_mesh, s.da_Thermal, 1, "pressure", -1)
    rank == 0 && @printf("initial lithostatic pressure: %.2e\n", e0)
    initial_velocity!(s); report(0)
    rank == 0 && println("Uzawa ", s.stokes.uzawa_iterations, " denok ", s.stokes.uzawa_residual)
    years(k) = parse(Float64, strip(last(split(read(joinpath(dir, "time_$k.txt"), String), ':'))))
    for k in 1:steps                     # use dt of the C++ code (Calc_dt_calor picks dx/dz with swapped indices)
        s.time_step_seconds = (years(k) - (k == 1 ? 0.0 : years(k - 1))) * SECONDS_PER_YEAR
        step!(s)
        add = add_remove && add_remove_particles!(m, s.particle_swarm, s.da_Veloc, k; cpp_order = true)
        rank == 0 && add_remove && println("particles +", add[1], " -", add[2])
        report(k)
    end
    PETSc.finalize(P)
end
main()
