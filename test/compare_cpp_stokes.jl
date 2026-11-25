# Repeats a Stokes solve dumped by the C++ Mandyoc and compares.
# julia --project=. test/mpirun.jl N compare_cpp_stokes.jl <folder> <Nx> <Nz> <Lx> <depth> <periodic> <fssa> <denok>
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/stokes.jl")
using Printf
function read_vec(file)
    open(file) do io
        ntoh(read(io, Int32)); n = ntoh(read(io, Int32)); ntoh.(read!(io, zeros(Float64, n)))
    end
end
function main()
    PETSc.initialize(P); comm = MPI.COMM_WORLD; rank = MPI.Comm_rank(comm)
    dir = ARGS[1]; Nx, Nz = parse.(Int, ARGS[2:3]); Lx, depth = parse.(Float64, ARGS[4:5])
    m = ModelParameters(; Nx, Nz, Lx, domain_depth = depth, periodic_boundary = ARGS[6] == "true", solver_relative_tolerance = 1e-13,
        uzawa_residual_tolerance = parse(Float64, ARGS[8]), free_surface_stabilization = ARGS[7] == "true", fssa_legacy = get(ARGS, 10, "true") == "true")
    rd(n) = read_vec(joinpath(dir, "sdump_$(n)_0.bin"))
    at(a) = (i, k) -> a[k*Nx+i+1]; at2(a) = (i, k) -> (a[2(k*Nx+i)+1], a[2(k*Nx+i)+2])
    dv, dT = create_velocity_mesh(m), create_thermal_mesh(m)
    f(n) = fill_field!(LibPETSc.DMCreateGlobalVector(P, dT), dT, at(rd(n)))
    fv(n) = fill_field!(LibPETSc.DMCreateGlobalVector(P, dv), dv, at2(rd(n)), 2)
    mask = rd("mask")
    s = StokesSolver(m, dv, dT, fv("V0"); mask = at2(mask))
    Vn, pr = fv("Vfb"), fv("Pb")
    stokes_solve!(s, Vn, pr, fv("V"), f("geoq"), f("drho"), parse(Float64, read(joinpath(dir, "sdump_dt.txt"), String));
        fssa_density_variation = f(get(ARGS, 9, "drhoA")))
    function cmp(x, name)
        ref = rd(name); si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, dv)); e = 0.0; b = 0.0; lo = Inf; hi = -Inf
        PETSc.with_local_array!(x; write = false) do a
            A = reshape(a, 2, mi, mk)
            for k in 1:mk, i in 1:mi
                name == "P" && (i + si == Nx || k + sk == Nz) && continue
                for c in (name == "P" ? (1:1) : (1:2))
                    r = ref[2((sk+k-1)*Nx+si+i-1)+c]; d = A[c, i, k] - r
                    e = max(e, abs(d)); b = max(b, abs(r)); lo = min(lo, d); hi = max(hi, d)
                end
            end
        end
        e, b, lo, hi = MPI.Allreduce(e, max, comm), MPI.Allreduce(b, max, comm), MPI.Allreduce(lo, min, comm), MPI.Allreduce(hi, max, comm)
        return e / b, (hi - lo) / 2b
    end
    ev, _ = cmp(Vn, "Vf"); ep, ep2 = cmp(pr, "P")
    rank == 0 && @printf("RESULT %s np=%d velocity: relative diff. %.3g | pressure: %.3g (without the constant: %.3g) | Uzawa %d, denok %.3g\n",
        ev < 1e-5 ? "OK" : "FAIL", MPI.Comm_size(comm), ev, ep, ep2, s.uzawa_iterations, s.uzawa_residual)
    PETSc.finalize(P)
end
main()
