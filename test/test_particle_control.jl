# mpiexec -n N julia --project=. test/test_particle_control.jl <shim|native>
# Uniform flow entering through x = 0 (the elements become empty): every element keeps a number of particles close
# to the nominal one (at most one is added per step and element), and the result does not depend on the number of ranks.
for f in ("particle_swarm", "thermal", "stokes", "rheology", "particle_control"); include("../src/$f.jl"); end
function main()
    PETSc.initialize(P); comm = MPI.COMM_WORLD; location = Symbol(get(ARGS, 1, "native"))
    m = ModelParameters(Nx = 31, Nz = 13, particles_per_element_x = 6, particles_per_element_z = 6)
    da = create_velocity_mesh(m); adv = create_advection_mesh(m, da)
    swarm, cell = create_particle_swarm(m, da; location, init = (x, z, id) -> (z < -150e3 ? 0 : 1, 1.0))
    U = 1e-9; dt = 0.02element_size_x(m) / U
    lv = LibPETSc.DMCreateLocalVector(P, adv); PETSc.with_local_array!(a -> (a[1:2:end] .= U; a[2:2:end] .= 0), lv)
    added = removed = 0
    for step in 1:150
        advect_particles!(m, swarm, adv, lv, dt)
        a, r = add_remove_particles!(m, swarm, da, step; location); added += a; removed += r
    end
    c = with_field(copy, swarm, PIC_COORDINATES, Float64); l = with_field(copy, swarm, "layer", PInt)
    counts = zeros(Int, m.Nx - 1, m.Nz - 1); wrong = 0; h = UInt64(0)
    for p in eachindex(l)
        counts[trunc(Int, c[2p-1] / element_size_x(m))+1, trunc(Int, (c[2p] + m.domain_depth) / element_size_z(m))+1] += 1
        wrong += l[p] != (c[2p] < -150e3 ? 0 : 1); h += hash((c[2p-1], c[2p], Int(l[p])))   # sum: independent of the order
    end
    counts = MPI.Allreduce(counts, +, comm); wrong = MPI.Allreduce(wrong, +, comm); h = MPI.Allreduce(h, +, comm)
    lo, hi = extrema(counts); ok = lo >= 25 && hi <= 42
    MPI.Comm_rank(comm) == 0 && println("RESULT ", ok ? "OK   " : "FAIL ", location, " np=", MPI.Comm_size(comm), " particles per element ", lo, " to ", hi,
        " (nominal 36), +", added, " -", removed, ", wrong layer ", wrong, ", hash ", string(h; base = 16))
    PETSc.finalize(P)
end
main()
