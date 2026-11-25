# Tempo por etapa: julia --project=. test/mpirun.jl N bench_stokes.jl [nx] [opções PETSc, ex. -veloc_mat_mumps_icntl_13 1]
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/stokes.jl")
using Printf
import LinearAlgebra
function main()
    PETSc.initialize(P; options = ARGS[2:end]); comm = MPI.COMM_WORLD
    n = parse(Int, get(ARGS, 1, "81"))
    m = ModelParameters(Nx = n, Nz = n, Lx = 91420.0, domain_depth = 100e3, rtol = 1e-7, uzawa_residual_tolerance = 1e-24,
        free_surface_stabilization = false, velocity_bc_top_tangential = true, velocity_bc_bottom_tangential = true)
    dv, dT = create_velocity_mesh(m), create_thermal_mesh(m)
    gt() = LibPETSc.DMCreateGlobalVector(P, dT); gv() = LibPETSc.DMCreateGlobalVector(P, dv)
    V, Vn, pr, b, x = gv(), gv(), gv(), gv(), gv()
    eta = fill_field!(gt(), dT, (i, k) -> 1e21)
    drho = fill_field!(gt(), dT, (i, k) -> k * element_size_z(m) < 20e3 + 2e3cospi(i / (n - 1)) ? 1e4 : 0.0)
    s = StokesSolver(m, dv, dT, V)
    clock(f) = (MPI.Barrier(comm); t = MPI.Wtime(); f(); MPI.Barrier(comm); MPI.Wtime() - t)
    L = LibPETSc
    assemble() = assemble_velocity!(s, V, pr, eta, drho, 0.0)
    assemble(); K = s.velocity_stiffness_matrix; ksp = s.velocity_ksp
    L.KSPSetOperators(P, ksp, K, K); L.KSPSetFromOptions(P, ksp)
    t_first = clock(() -> L.KSPSetUp(P, ksp))
    solve_velocity!(s, s.velocity_rhs_pressure, x)
    t_asm = minimum(clock(assemble) for _ in 1:3)
    t_fac = minimum(begin assemble(); L.KSPSetOperators(P, ksp, K, K); clock(() -> L.KSPSetUp(P, ksp)) end for _ in 1:3)
    t_sol = minimum(clock(() -> solve_velocity!(s, s.velocity_rhs_pressure, x)) for _ in 1:20)
    t_mul = minimum(clock(() -> L.MatMult(P, K, x, b)) for _ in 1:20)
    t_dot = minimum(clock(() -> L.VecDot(P, x, x)) for _ in 1:20)
    MPI.Comm_rank(comm) == 0 && @printf("BENCH np=%d n=%d dof=%d blas=%d | montagem %.1f ms | 1a fatoracao %.1f | fatoracao %.1f | solve %.2f | MatMult %.3f | VecDot %.3f ms\n",
        MPI.Comm_size(comm), n, 2n^2, LinearAlgebra.BLAS.get_num_threads(), 1e3 .* (t_asm, t_first, t_fac, t_sol, t_mul, t_dot)...)
    PETSc.finalize(P)
end
main()
