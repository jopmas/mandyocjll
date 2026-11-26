# mpiexec -n N julia --project=. test/test_thermal3d.jl   (or: julia --project=. test/mpirun.jl N test_thermal3d.jl)
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/thermal3d.jl")
using Printf
const comm = MPI.COMM_WORLD
report(ok, s) = MPI.Comm_rank(comm) == 0 && println("RESULT ", ok ? "OK   " : "FAIL ", "np=", MPI.Comm_size(comm), " ", s)

function run_case(m, T0, vel, steps, dt; H = (i, j, k) -> 0.0, kappa = (i, j, k) -> m.thermal_diffusivity)
    dv, dT = create_velocity_mesh3d(m), create_thermal_mesh3d(m)
    g() = LibPETSc.DMCreateGlobalVector(P, dT)
    T = fill_field3!(g(), dT, T0)
    v = fill_field3!(LibPETSc.DMCreateGlobalVector(P, dv), dv, vel, 3)
    h, kf = fill_field3!(g(), dT, H), fill_field3!(g(), dT, kappa)
    s = ThermalSolver3D(m, dT, dv, T)
    for _ in 1:steps; thermal_step!(s, T, v, v, h, kf, dt); end
    return T, dT, s
end
function max_error(T, da, exact)
    e = 0.0
    PETSc.with_local_array!(T; write = false) do a
        A = owned_view3(da, a)
        for k in axes(A, 3), j in axes(A, 2), i in axes(A, 1); e = max(e, abs(A[i, j, k] - exact(i, j, k))); end
    end
    MPI.Allreduce(e, max, comm)
end

function main()
    PETSc.initialize(P)
    L = 100e3; kap = 1e-6; still = (i, j, k) -> (0.0, 0.0, 0.0)
    base(n, nz; kw...) = ModelParameters(; dimensions = 3, Nx = n, Ny = n, Nz = nz, Lx = L, Ly = L, domain_depth = L,
        solver_relative_tolerance = 1e-13, temperature_bc_left = false, temperature_bc_right = false,
        temperature_bc_front = false, temperature_bc_back = false, kw...)
    # 1. transient diffusion of sin(pi z): exp(-kappa pi^2 t/D^2); order 2
    errs = map((6, 11, 21)) do n
        m = base(5, n); t = 0.1L^2 / kap; steps = 2(n - 1)
        T, d, _ = run_case(m, (i, j, k) -> 1000sinpi(k / (n - 1)), still, steps, t / steps)
        max_error(T, d, (i, j, k) -> 1000sinpi(k / (n - 1)) * exp(-pi^2 * 0.1))
    end
    ord = log2(errs[2] / errs[3])
    report(1.8 < ord < 2.2, @sprintf("3D transient diffusion: errors %.3g %.3g %.3g, order %.2f", errs..., ord))
    # 2. two layers (1e-6 and 0.25e-6): steady state exact at the nodes
    m = base(5, 21; thermal_time_theta = 1.0)
    T, d, _ = run_case(m, (i, j, k) -> 1000 * (1 - k / 20), still, 40, 1e18; kappa = (i, j, k) -> k < 10 ? 1e-6 : 0.25e-6)
    e = max_error(T, d, (i, j, k) -> k <= 10 ? 1000 - 200k / 10 : 800 - 800(k - 10) / 10)
    report(e < 1e-6, @sprintf("3D two layers, steady state: error %.3g", e))
    # 3. radiogenic heating: parabolic profile H D^2/(2 kappa c) z(1-z)
    m = base(5, 21; thermal_time_theta = 1.0, use_radiogenic_heating = true)
    amp = 1e-10 * L^2 / (2kap * m.heat_capacity)
    T, d, _ = run_case(m, (i, j, k) -> 0.0, still, 40, 1e18; H = (i, j, k) -> 1e-10)
    e = max_error(T, d, (i, j, k) -> amp * (k / 20) * (1 - k / 20))
    report(e < 1e-6 * amp, @sprintf("3D radiogenic heating, steady state: error %.3g of %.4g", e, amp / 4))
    # 4. steady advection-diffusion along y, mesh Peclet number 5 (fixed temperature at y = 0 and Ly)
    m = ModelParameters(; dimensions = 3, Nx = 4, Ny = 41, Nz = 4, Lx = L, Ly = L, domain_depth = L, solver_relative_tolerance = 1e-13,
        thermal_time_theta = 1.0, temperature_bc_top = false, temperature_bc_bottom = false, temperature_bc_left = false, temperature_bc_right = false)
    Pe = 400.0; v = Pe * kap / L
    T, d, _ = run_case(m, (i, j, k) -> 1000 * (1 - j / 40), (i, j, k) -> (0.0, v, 0.0), 60, 1e17)
    e = max_error(T, d, (i, j, k) -> 1000 * (1 - exp(Pe * (j / 40 - 1))) / (1 - exp(-Pe)))
    report(e < 1, @sprintf("3D steady advection-diffusion in y, mesh Pe %.0f: max. error %.3g of 1000", Pe / 40 / 2, e))
    PETSc.finalize(P)
end
main()
