# Repete um passo térmico gravado pelo Mandyoc C++ (dump) e compara.
# julia --project=. test/mpirun.jl N compare_cpp_thermal.jl <pasta> <Nx> <Nz> <Lx> <depth> <periodic>
include("../src/particle_swarm.jl"); include("../src/thermal.jl")
using Printf
function read_vec(file)          # PetscBinary: int32 classe, int32 n, float64, big-endian
    open(file) do io
        ntoh(read(io, Int32)); n = ntoh(read(io, Int32)); ntoh.(read!(io, zeros(Float64, n)))
    end
end
function main()
    PETSc.initialize(P); comm = MPI.COMM_WORLD
    dir = ARGS[1]; Nx, Nz = parse.(Int, ARGS[2:3]); Lx, depth = parse.(Float64, ARGS[4:5])
    m = ModelParameters(; Nx, Nz, Lx, domain_depth = depth, periodic_boundary = ARGS[6] == "true", rtol = 1e-13,
        use_radiogenic_heating = true, use_adiabatic_heating = true, temperature_bc_left = true, temperature_bc_right = true)
    rd(n) = read_vec(joinpath(dir, n)); at(a) = (i, k) -> a[k*Nx+i+1]; at2(a) = (i, k) -> (a[2(k*Nx+i)+1], a[2(k*Nx+i)+2])
    dv, dT = create_velocity_mesh(m), create_thermal_mesh(m)
    f(a) = fill_field!(LibPETSc.DMCreateGlobalVector(P, dT), dT, at(a))
    fv(a) = fill_field!(LibPETSc.DMCreateGlobalVector(P, dv), dv, at2(a), 2)
    T0, T1 = rd("dump_T0_3.bin"), rd("dump_T1_3.bin")
    T = f(T0); s = ThermalSolver(m, dT, dv, T)
    dt = parse(Float64, read(joinpath(dir, "dump_dt.txt"), String))
    thermal_step!(s, T, fv(rd("dump_V_3.bin")), fv(rd("dump_Vf_3.bin")), f(rd("dump_H_3.bin")), f(rd("dump_kappa_3.bin")), dt)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, dT)); e = 0.0
    PETSc.with_local_array!(T; write = false) do a
        A = reshape(a, mi, mk)
        for k in 1:mk, i in 1:mi; e = max(e, abs(A[i, k] - T1[(sk+k-1)*Nx+si+i])); end
    end
    e = MPI.Allreduce(e, max, comm)
    MPI.Comm_rank(comm) == 0 && @printf("RESULT %s np=%d dif. maxima Julia - C++ = %.3g C (o passo mudou T em ate %.4g C; T max %.6g C; dt %.4g s; its %d)\n",
        e < 1e-6 ? "OK" : "FAIL", MPI.Comm_size(comm), e, maximum(abs, T1 .- T0), maximum(T1), dt, s.iterations)
    PETSc.finalize(P)
end
main()
