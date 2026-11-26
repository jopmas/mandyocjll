# mpiexec -n N julia --project=. test/test_stokes3d.jl   (or: julia --project=. test/mpirun.jl N test_stokes3d.jl)
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/thermal3d.jl")
include("../src/stokes.jl"); include("../src/stokes3d.jl")
using Printf
const comm = MPI.COMM_WORLD
report(ok, s) = MPI.Comm_rank(comm) == 0 && println("RESULT ", ok ? "OK   " : "FAIL ", "np=", MPI.Comm_size(comm), " ", s)
function owned(f, v, da)
    e = 0.0
    PETSc.with_local_array!(v; write = false) do a
        A = owned_view3(da, a, 3)
        for k in axes(A, 4), j in axes(A, 3), i in axes(A, 2); e = max(e, f(view(A, :, i, j, k), i, j, k)); end
    end
    MPI.Allreduce(e, max, comm)
end
function setup(m; v0 = (i, j, k) -> (0.0, 0.0, 0.0), eta = (i, j, k) -> 1.0, drho = (i, j, k) -> 0.0, mask = nothing)
    dv, dT = create_velocity_mesh3d(m), create_thermal_mesh3d(m)
    g(f) = fill_field3!(LibPETSc.DMCreateGlobalVector(P, dT), dT, f)
    V = fill_field3!(LibPETSc.DMCreateGlobalVector(P, dv), dv, v0, 3); Vn = LibPETSc.DMCreateGlobalVector(P, dv)
    pr = LibPETSc.DMCreateGlobalVector(P, dv)
    LibPETSc.VecCopy(P, V, Vn); LibPETSc.VecZeroEntries(P, pr)
    (; s = StokesSolver3D(m, dv, dT, V; mask), V, Vn, pr, eta = g(eta), drho = g(drho), dv)
end
solve!(c, dt = 0.0) = stokes_solve!(c.s, c.Vn, c.pr, c.V, c.eta, c.drho, dt)

function main()
    PETSc.initialize(P)
    base = (; dimensions = 3, Lx = 1.0, Ly = 0.5, domain_depth = 1.0, solver_relative_tolerance = 1e-13,
        uzawa_residual_tolerance = 1e-22, free_surface_stabilization = false)
    # 1. free-slip mode, eta = 1, uniform in y: the 2D solution, with v = 0
    errs = map((6, 11, 21)) do n
        m = ModelParameters(; Nx = n, Ny = 4, Nz = n, base...)
        c = setup(m; drho = (i, j, k) -> cospi(i / (n - 1)) * sinpi(k / (n - 1)))
        solve!(c)
        a = pi^2 / (4pi^4)
        owned(c.Vn, c.dv) do v, i, j, k
            X, Z = i / (n - 1), k / (n - 1)
            max(abs(v[1] + a * sinpi(X) * cospi(Z)), abs(v[2]), abs(v[3] - a * cospi(X) * sinpi(Z))) / a
        end
    end
    o = log2(errs[2] / errs[3])
    report(1.8 < o < 2.2, @sprintf("3D free-slip mode: relative errors %.3g %.3g %.3g, order %.2f", errs..., o))
    # 2. two-layer Couette flow in x (viscosities 1 and 100), exact at the nodes; free slip at y = 0 and Ly
    n = 21; H(k) = k <= 10 ? 0.5 * k / 10 : 0.5 + 0.5 * (k - 10) / 10 / 100
    ex(k) = H(k) / H(20)
    m = ModelParameters(; Nx = n, Ny = 5, Nz = n, base...)
    c = setup(m; v0 = (i, j, k) -> (k == 20 ? 1.0 : i in (0, n - 1) ? ex(k) : 0.0, 0.0, 0.0), eta = (i, j, k) -> k < 10 ? 1.0 : 100.0,
        mask = (i, j, k) -> (k == 0 || k == 20 || i in (0, n - 1)) ? (0.0, 0.0, 0.0) : j in (0, 4) ? (1.0, 0.0, 1.0) : (1.0, 1.0, 1.0))
    solve!(c)
    e = owned((v, i, j, k) -> max(abs(v[1] - ex(k)), abs(v[2]), abs(v[3])), c.Vn, c.dv)
    report(e < 1e-7, @sprintf("3D two-layer Couette: error %.3g (Uzawa %d)", e, c.s.uzawa_iterations))
    PETSc.finalize(P)
end
main()
