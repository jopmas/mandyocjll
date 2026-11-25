# mpiexec -n N julia --project=. test/test_particle_swarm.jl <shim|native> <closed|periodic> [passos] [stencil da malha de advecção]
include("../src/particle_swarm.jl")
using Printf

function velocity_at(m, x, z)
    U = 1.0e-9                                   # m/s (~3 cm/ano)
    X = x / m.Lx; Z = (z + m.domain_depth) / m.domain_depth
    m.periodic_boundary || return (U * sinpi(X) * cospi(Z), -U * m.domain_depth / m.Lx * cospi(X) * sinpi(Z))
    return (U * (1 + 0.5Z), 0.3U * sinpi(2X) * sinpi(Z))
end

function main()
    location = Symbol(get(ARGS, 1, "shim")); boundary = get(ARGS, 2, "closed")
    nsteps = parse(Int, get(ARGS, 3, "200")); sw = parse(Int, get(ARGS, 4, "2"))
    PETSc.initialize(petsclib)
    comm = MPI.COMM_WORLD; rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    m = ModelParameters(periodic_boundary = boundary == "periodic", advection_stencil_width = sw)
    da_Veloc = create_velocity_mesh(m); da_advection = create_advection_mesh(m, da_Veloc)
    particle_swarm, cell_dm = create_particle_swarm(m, da_Veloc; location)   # dms, dmcell
    velocity = LibPETSc.DMCreateGlobalVector(petsclib, da_Veloc)             # Veloc
    local_velocity = LibPETSc.DMCreateLocalVector(petsclib, da_advection)        # local_V
    si, sk, _, mi, mk, _ = LibPETSc.DMDAGetCorners(petsclib, da_Veloc)
    PETSc.with_local_array!(velocity) do a
        n = 0
        for k in sk:sk+mk-1, i in si:si+mi-1
            u, w = velocity_at(m, i * element_size_x(m), k * element_size_z(m) - m.domain_depth)
            a[n+=1] = u; a[n+=1] = w
        end
    end
    LibPETSc.DMGlobalToLocal(petsclib, da_advection, velocity, LibPETSc.INSERT_VALUES, local_velocity)
    dt = 0.5 * element_size_x(m) / 1.5e-9
    expected = (m.Nx - 1) * (m.Nz - 1) * particles_per_element(m)
    total() = MPI.Allreduce(Int(LibPETSc.DMSwarmGetLocalSize(petsclib, particle_swarm)), +, comm)
    initial = total(); moved = 0; firstbad = 0
    advect_particles!(m, particle_swarm, da_advection, local_velocity, 0.0)      # compila antes de medir
    MPI.Barrier(comm); t0 = MPI.Wtime()
    for s in 1:nsteps
        moved += advect_particles!(m, particle_swarm, da_advection, local_velocity, dt)
        firstbad == 0 && total() != expected && (firstbad = s)
    end
    MPI.Barrier(comm); elapsed = MPI.Wtime() - t0

    f(name, T) = with_field(copy, particle_swarm, name, T)
    coor, itag, layer, vis = f(PIC_COORDINATES, Float64), f("itag", PInt), f("layer", PInt), f("geoq_fac", Float64)
    badowner = count(p -> owned_cell(CELL_CONTEXT[], coor[2p-1], coor[2p]) < 0, eachindex(itag))
    badowner = MPI.Allreduce(badowner, +, comm)
    gather(v) = (c = MPI.Gather(Cint(length(v)), comm; root = 0);
        MPI.Gatherv!(v, rank == 0 ? VBuffer(zeros(eltype(v), sum(c)), c) : nothing, comm; root = 0))
    C, I, L, G = gather(coor), gather(itag), gather(layer), gather(vis)
    if rank == 0
        o = sortperm(I); ids = I[o] == 1:expected
        payload = all(L .== I .% 7) && all(G .== 1.5 .* I)
        pos = reduce(hcat, [C[[2p-1, 2p]] for p in o])
        file = joinpath(@__DIR__, "..", "results", "positions_$(boundary).bin")
        diff = NaN
        if np == 1 && location == :shim
            write(file, pos)
        elseif isfile(file) && ids
            ref = reshape(reinterpret(Float64, read(file)), 2, :)
            size(ref) == size(pos) && (diff = maximum(abs, ref .- pos))
        end
        ok = ids && payload && (np == 1 || diff < 1e-3) && (location == :native || badowner == 0)
        @printf("RESULT %-4s %-6s %-8s np=%d sw=%d N0=%d N=%d esperado=%d 1o_passo_ruim=%d migracoes=%d ids=%s payload=%s fora_do_dono=%d dpos_max=%.3g m t=%.2fs\n",
            ok ? "OK" : "FAIL", location, boundary, np, sw, initial, length(I), expected, firstbad, moved, ids, payload, badowner, diff, elapsed)
    end
    LibPETSc.VecDestroy(petsclib, velocity); LibPETSc.VecDestroy(petsclib, local_velocity)
    LibPETSc.DMDestroy(petsclib, particle_swarm)
    PETSc.finalize(petsclib)
end
main()
