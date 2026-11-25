# mpiexec -n N julia --project=. test/test_open_boundary.jl <shim|native>
# Uniform flow leaving through x = Lx (case of the `break` in RK4): the particles that leave are
# removed and all the others move exactly U*t.
include("../src/particle_swarm.jl")
function main()
    PETSc.initialize(petsclib); comm = MPI.COMM_WORLD
    location = Symbol(get(ARGS, 1, "shim"))
    m = ModelParameters(Nx = 41, Nz = 13)
    da = create_velocity_mesh(m); da_advection = create_advection_mesh(m, da)
    swarm, cell_dm = create_particle_swarm(m, da; location)
    U = 1e-9; dt = 0.3element_size_x(m) / U; steps = 30
    lv = LibPETSc.DMCreateLocalVector(petsclib, da_advection)
    PETSc.with_local_array!(a -> (a[1:2:end] .= U; a[2:2:end] .= 0.0), lv)
    read_all() = (with_field(copy, swarm, "itag", PInt), with_field(copy, swarm, PIC_COORDINATES, Float64)[1:2:end])
    gather(v) = (c = MPI.Allgather(Cint(length(v)), comm); MPI.Allgatherv!(v, VBuffer(zeros(eltype(v), sum(c)), c), comm))
    I0, X0 = gather.(read_all())
    for _ in 1:steps; advect_particles!(m, swarm, da_advection, lv, dt); end
    I, X = gather.(read_all())
    x0 = Dict(zip(I0, X0)); shift = U * dt * steps
    stay = Set(i for (i, x) in x0 if x + shift <= m.Lx)
    ok = length(I) == length(stay) && Set(I) == stay
    err = maximum(abs(X[p] - x0[I[p]] - shift) for p in eachindex(I))
    MPI.Comm_rank(comm) == 0 && println("RESULT ", ok && err < 1e-3 ? "OK   " : "FAIL ", location, " open boundary np=", MPI.Comm_size(comm),
        " stayed=", length(I), " expected=", length(stay), " left=", length(I0) - length(I), " error_x=", err, " m")
    PETSc.finalize(petsclib)
end
main()
