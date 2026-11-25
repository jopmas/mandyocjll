# mpiexec -n N julia --project=. test/test_stokes.jl
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/stokes.jl")
using Printf
const comm = MPI.COMM_WORLD
report(ok, s) = MPI.Comm_rank(comm) == 0 && println("RESULT ", ok ? "OK   " : "FAIL ", "np=", MPI.Comm_size(comm), " ", s)
function owned(f, v, da, dof)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da)); e = 0.0
    PETSc.with_local_array!(v; write = false) do a
        A = reshape(a, dof, mi, mk)
        for k in 1:mk, i in 1:mi; e = max(e, f(view(A, :, i, k), si + i - 1, sk + k - 1)); end
    end
    MPI.Allreduce(e, max, comm)
end
function setup(m; v0 = (i, k) -> (0.0, 0.0), eta = (i, k) -> 1.0, drho = (i, k) -> 0.0, mask = nothing)
    dv, dT = create_velocity_mesh(m), create_thermal_mesh(m)
    g(f) = fill_field!(LibPETSc.DMCreateGlobalVector(P, dT), dT, f)
    V = fill_field!(LibPETSc.DMCreateGlobalVector(P, dv), dv, v0, 2); Vn = LibPETSc.DMCreateGlobalVector(P, dv)
    pr = LibPETSc.DMCreateGlobalVector(P, dv)
    LibPETSc.VecCopy(P, V, Vn); LibPETSc.VecZeroEntries(P, pr)
    (; s = StokesSolver(m, dv, dT, V; mask), V, Vn, pr, eta = g(eta), drho = g(drho), dv)
end
solve!(c, dt = 0.0) = stokes_solve!(c.s, c.Vn, c.pr, c.V, c.eta, c.drho, dt)

function main()
    PETSc.initialize(P)
    base = (; Lx = 1.0, domain_depth = 1.0, rtol = 1e-13, uzawa_residual_tolerance = 1e-22, free_surface_stabilization = false)
    # 1. força psi = sin^2(pi x) sin^2(pi z) (velocidade nula nas paredes), eta = 1: f = -eta lap(v)
    #    com força só em z não há solução fechada geral; usa-se o modo free-slip: rho = cos(pi x) sin(pi z)
    #    -> w = cos(pi x) sin(pi z) * k^2/(eta (2k^2)^2) * amplitude, k = pi
    errs = map((11, 21, 41)) do n
        m = ModelParameters(; Nx = n, Nz = n, base...)
        c = setup(m; drho = (i, k) -> cospi(i / (n - 1)) * sinpi(k / (n - 1)))
        solve!(c)
        a = pi^2 / (4pi^4)
        owned(c.Vn, c.dv, 2) do v, i, k
            X, Z = i / (n - 1), k / (n - 1)
            max(abs(v[1] + a * sinpi(X) * cospi(Z)), abs(v[2] - a * cospi(X) * sinpi(Z))) / a
        end
    end
    o = log2(errs[2] / errs[3])
    report(1.8 < o < 2.2, @sprintf("modo free-slip: erros relativos %.3g %.3g %.3g, ordem %.2f", errs..., o))
    # 2. Couette de duas camadas (viscosidade 1 e 100), exato nos nós; contorno periódico
    for periodic in (true, false)
        n = 21; H(k) = k <= 10 ? 0.5 * k / 10 * 1 : 0.5 + 0.5 * (k - 10) / 10 / 100
        ex(k) = H(k) / H(20)
        m = ModelParameters(; Nx = n, Nz = n, base..., periodic_boundary = periodic)
        c = setup(m; v0 = (i, k) -> (k == 20 ? 1.0 : (!periodic && i in (0, n - 1)) ? ex(k) : 0.0, 0.0), eta = (i, k) -> k < 10 ? 1.0 : 100.0,
            mask = (i, k) -> (k == 0 || k == 20 || (!periodic && i in (0, n - 1))) ? (0.0, 0.0) : (1.0, 1.0))
        solve!(c)
        e = owned((v, i, k) -> max(abs(v[1] - ex(k)), abs(v[2])), c.Vn, c.dv, 2)
        report(e < 1e-7, @sprintf("Couette duas camadas, periodico=%s: erro %.3g (Uzawa %d)", periodic, e, c.s.uzawa_iterations))
    end
    # 3. FSSA: perturbação de superfície livre sobre "ar". Sem FSSA a velocidade é w0. Deslocar a
    #    interface de theta*dt*w0 muda a força; a resposta define theta efetivo.
    for legacy in (false, true)
        n = 41; D = 100e3; eta(i, k) = k < 32 ? 1e21 : 1e19
        mk(stab) = ModelParameters(; Nx = n, Nz = n, Lx = D, domain_depth = D, rtol = 1e-13, uzawa_residual_tolerance = 1e-30,
            free_surface_stabilization = stab, fssa_theta = 0.5, fssa_legacy = legacy)
        z0(i, amp) = 0.8D + amp * cospi(i / (n - 1))
        # densidade nodal = fração do volume de controle do nó abaixo da interface
        rho(amp) = (i, k) -> (dz = D / (n - 1); -10 * 3300 * clamp((z0(i, amp) - (k * dz - dz / 2)) / dz, 0, 1))
        w(stab, amp, dt) = (c = setup(mk(stab); eta, drho = rho(amp)); solve!(c, dt); owned((v, i, k) -> (i == 0 && k == 32) ? abs(v[2]) : 0.0, c.Vn, c.dv, 2))
        amp = 200.0; w0 = w(false, amp, 0.0)
        tau = amp / w0; dt = 0.2tau
        theta_eff = (w0 / w(true, amp, dt) - 1) / (dt / tau)
        report(legacy ? true : 0.3 < theta_eff < 0.7, @sprintf("FSSA %s: theta nominal 0.5, theta efetivo %.3f (tau = %.4g anos)", legacy ? "como no C++" : "consistente", theta_eff, tau / 3.156e7))
    end
    PETSc.finalize(P)
end
main()
