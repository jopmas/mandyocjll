# van Keken et al. (1997), case 1a: Rayleigh-Taylor instability, uniform viscosity.
# julia --project=. test/mpirun.jl N van_keken.jl [nx] [ppe in each direction] [nondimensional final time] [max. dt, s]
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/stokes.jl")
using Printf
function main()
    PETSc.initialize(P); comm = MPI.COMM_WORLD; rank = MPI.Comm_rank(comm)
    n = parse(Int, get(ARGS, 1, "81")); pd = parse(Int, get(ARGS, 2, "6")); tend = parse(Float64, get(ARGS, 3, "300")); dtmax = parse(Float64, get(ARGS, 4, "1.6e12"))
    m = ModelParameters(Nx = n, Nz = n, Lx = 91420.0, domain_depth = 100e3, particles_per_element_x = pd, particles_per_element_z = pd,
        solver_relative_tolerance = 1e-7, uzawa_residual_tolerance = 1e-24, free_surface_stabilization = false, use_runge_kutta_4 = true,
        velocity_bc_top_tangential = true, velocity_bc_bottom_tangential = true)
    dv, dT = create_velocity_mesh(m), create_thermal_mesh(m); da = create_advection_mesh(m, dv)
    interface(x) = -80e3 + 2e3 * cospi(x / m.Lx)
    swarm, cell = create_particle_swarm(m, dv; location = :native, init = (x, z, id) -> (z < interface(x) ? 0 : 1, 1e21))
    gt() = LibPETSc.DMCreateGlobalVector(P, dT); gv() = LibPETSc.DMCreateGlobalVector(P, dv)
    rho, eta, drho, T = gt(), gt(), gt(), gt(); V, Vn, pr = gv(), gv(), gv()
    foreach(x -> LibPETSc.VecZeroEntries(P, x), (rho, eta, T, V, Vn, pr))
    s = StokesSolver(m, dv, dT, V); lv = LibPETSc.DMCreateLocalVector(P, da)
    function vrms()
        si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, dv)); e = 0.0; vmax = 0.0
        PETSc.with_local_array!(Vn; write = false) do a
            A = reshape(a, 2, mi, mk)
            for k in 1:mk, i in 1:mi
                w = (si + i - 1 in (0, n - 1) ? 0.5 : 1.0) * (sk + k - 1 in (0, n - 1) ? 0.5 : 1.0)
                e += w * (A[1, i, k]^2 + A[2, i, k]^2); vmax = max(vmax, abs(A[1, i, k]), abs(A[2, i, k]))
            end
        end
        sqrt(MPI.Allreduce(e, +, comm) / (n - 1)^2), MPI.Allreduce(vmax, max, comm)
    end
    U, t0 = 1e-7, 1e12; t = 0.0; best = (0.0, 0.0); step = 0; its = 0; N0 = 0
    out = rank == 0 ? open(joinpath(@__DIR__, "..", "results", "van_keken_$(n)_$(get(ARGS, 4, "1.6e12")).txt"), "w") : nothing
    tic = MPI.Wtime()
    while t <= tend
        swarm_to_mesh!(m, swarm, dT, rho, eta, [-1000.0, 0.0], [1e21, 1e21])
        calc_density_variation!(m, drho, T, rho)
        LibPETSc.VecCopy(P, Vn, V)
        stokes_solve!(s, Vn, pr, V, eta, drho, 0.0); its += s.uzawa_iterations
        vr, vmax = vrms(); vr / U > best[1] && (best = (vr / U, t))
        rank == 0 && println(out, t, " ", vr / U)
        dt = min(0.5element_size_x(m) / vmax, dtmax)
        global_to_local(da, Vn, lv)
        advect_particles!(m, swarm, da, lv, dt)
        t += dt / t0; step += 1
    end
    Np = MPI.Allreduce(Int(LibPETSc.DMSwarmGetLocalSize(P, swarm)), +, comm)
    rank == 0 && (close(out); @printf("RESULT np=%d mesh %d^2 ppe=%d: max. vrms %.6f at t = %.2f | %d steps, %.1f Uzawa/step, particles %d of %d, %.1f s\n",
        MPI.Comm_size(comm), n - 1, pd^2, best..., step, its / step, Np, (n - 1)^2 * pd^2, MPI.Wtime() - tic))
    PETSc.finalize(P)
end
main()
