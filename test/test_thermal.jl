# mpiexec -n N julia --project=. test/test_thermal.jl
include("../src/particle_swarm.jl"); include("../src/thermal.jl")
using Printf
const comm = MPI.COMM_WORLD
report(ok, s) = MPI.Comm_rank(comm) == 0 && println("RESULT ", ok ? "OK   " : "FAIL ", "np=", MPI.Comm_size(comm), " ", s)

function run_case(m, T0, vel, steps, dt; H = (i, k) -> 0.0, kappa = (i, k) -> m.thermal_diffusivity)
    dv, dT = create_velocity_mesh(m), create_thermal_mesh(m)
    g() = LibPETSc.DMCreateGlobalVector(P, dT)
    T = fill_field!(g(), dT, T0)
    v = fill_field!(LibPETSc.DMCreateGlobalVector(P, dv), dv, vel, 2)
    h, kf = fill_field!(g(), dT, H), fill_field!(g(), dT, kappa)
    s = ThermalSolver(m, dT, dv, T)
    for _ in 1:steps; thermal_step!(s, T, v, v, h, kf, dt); end
    return T, dT, s
end
function max_error(T, da, exact)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da)); e = 0.0
    PETSc.with_local_array!(T; write = false) do a
        A = reshape(a, mi, mk)
        for k in 1:mk, i in 1:mi; e = max(e, abs(A[i, k] - exact(si + i - 1, sk + k - 1))); end
    end
    MPI.Allreduce(e, max, comm)
end

function main()
    PETSc.initialize(P)
    L = 100e3; kap = 1e-6; still = (i, k) -> (0.0, 0.0)
    # 1. transient diffusion: mode sin(pi z) decays as exp(-kappa pi^2 t/D^2); order 2
    errs = map((11, 21, 41)) do n
        m = ModelParameters(Nx = n, Nz = n, Lx = L, domain_depth = L, solver_relative_tolerance = 1e-13)
        t = 0.1L^2 / kap; steps = 2(n - 1)
        T, d, _ = run_case(m, (i, k) -> 1000sinpi(k / (n - 1)), still, steps, t / steps)
        max_error(T, d, (i, k) -> 1000sinpi(k / (n - 1)) * exp(-pi^2 * 0.1))
    end
    ord = log2(errs[2] / errs[3])
    report(1.8 < ord < 2.2, @sprintf("transient diffusion: errors %.3g %.3g %.3g, order %.2f", errs..., ord))
    # 2. two layers (1e-6 and 0.25e-6), steady state exact at the nodes
    m = ModelParameters(Nx = 11, Nz = 21, Lx = L, domain_depth = L, solver_relative_tolerance = 1e-13, thermal_time_theta = 1.0)
    kf(i, k) = k < 10 ? 1e-6 : 0.25e-6
    T, d, _ = run_case(m, (i, k) -> 1000 * (1 - k / 20), still, 40, 1e18; kappa = kf)
    e = max_error(T, d, (i, k) -> k <= 10 ? 1000 - 200k / 10 : 800 - 800(k - 10) / 10)
    report(e < 1e-6, @sprintf("two layers, steady state: error %.3g", e))
    # 3. radiogenic heating: parabolic profile H D^2/(2 kappa c) z(1-z)
    m = ModelParameters(Nx = 11, Nz = 21, Lx = L, domain_depth = L, solver_relative_tolerance = 1e-13, thermal_time_theta = 1.0, use_radiogenic_heating = true)
    amp = 1e-10 * L^2 / (2kap * m.heat_capacity)
    T, d, _ = run_case(m, (i, k) -> 0.0, still, 40, 1e18; H = (i, k) -> 1e-10)
    e = max_error(T, d, (i, k) -> amp * (k / 20) * (1 - k / 20))
    report(e < 1e-6 * amp, @sprintf("radiogenic heating, steady state: error %.3g of %.4g", e, amp / 4))
    # 4. steady advection-diffusion, mesh Peclet number 5: bounded and without oscillations
    m = ModelParameters(Nx = 11, Nz = 41, Lx = L, domain_depth = L, solver_relative_tolerance = 1e-13, thermal_time_theta = 1.0)
    Pe = 400.0; w = Pe * kap / L
    T, d, _ = run_case(m, (i, k) -> 1000 * (1 - k / 40), (i, k) -> (0.0, w), 60, 1e17)
    ex(i, k) = 1000 * (1 - exp(Pe * (k / 40 - 1))) / (1 - exp(-Pe))
    e = max_error(T, d, ex)
    report(e < 1, @sprintf("steady advection-diffusion, mesh Pe %.0f: max. error %.3g of 1000", Pe / 40 / 2, e))
    # 5. periodic: a wave advected around the domain returns to the initial position, with decay
    errs = map((41, 81, 161)) do n
        m = ModelParameters(Nx = n, Nz = 5, Lx = L, domain_depth = L, solver_relative_tolerance = 1e-13, periodic_boundary = true,
            thermal_diffusivity = 1e-6, temperature_bc_top = false, temperature_bc_bottom = false)
        U = 1e-9; steps = 4(n - 1)
        # kappa large enough for mesh Pe < 1: no artificial diffusion
        kp = U * L / 20; t = L / U
        T, d, _ = run_case(m, (i, k) -> 1000sinpi(2i / (n - 1)), (i, k) -> (U, 0.0), steps, t / steps; kappa = (i, k) -> kp)
        max_error(T, d, (i, k) -> 1000sinpi(2i / (n - 1)) * exp(-4pi^2 * kp * t / L^2))
    end
    ord = log2(errs[2] / errs[3])
    report(1.8 < ord < 2.2, @sprintf("periodic, advection+diffusion: errors %.3g %.3g %.3g, order %.2f", errs..., ord))
    PETSc.finalize(P)
end
main()
