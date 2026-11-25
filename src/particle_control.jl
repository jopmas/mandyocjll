# Addition and removal of particles (Swarm_add_remove_2d). Include after rheology.jl.
# As in the C++ code, each call adds or removes at most one particle per element.

"Elements whose particles belong to this rank, according to the location rule."
function particle_elements(m, da_Veloc, location)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da_Veloc))
    location == :shim && return si:min(si + mi, m.Nx - 1)-1, sk:min(sk + mk, m.Nz - 1)-1
    # native: the element belongs to the owner of its upper right node (DMLocatePoints_DA_Regular)
    return max(si - 1, 0):si+mi-2, max(sk - 1, 0):sk+mk-2
end

# same sequence as the C++ code with 1 rank; here it does not depend on the rank
function candidate_points(m, step, i, k)
    dx, dz = element_size_x(m), element_size_z(m)
    seed = Ref{Cuint}((UInt64(step) * 1234577 + UInt64(k) * UInt64(m.Nx) + UInt64(i)) % Cuint)   # cell_seed
    uniform() = 2.0 * Float32(ccall(:rand_r, Cint, (Ptr{Cuint},), seed)) / Float64(Float32(typemax(Int32))) - 1.0
    return [(rx = uniform(); rz = uniform();
        (i * dx + (0.5rx + 0.5) * dx, k * dz - m.domain_depth + (0.5rz + 0.5) * dz)) for _ in 1:10]
end

"Indices of the particles of the element or, if it is empty, of the closest neighboring elements."
function nearby_particles(buckets, i, k)
    haskey(buckets, (i, k)) && return buckets[(i, k)]
    for r in 1:5
        near = Int[]
        for kk in k-r:k+r, ii in i-r:i+r
            max(abs(ii - i), abs(kk - k)) == r && append!(near, get(buckets, (ii, kk), Int[]))
        end
        isempty(near) || return near
    end
    return Int[]
end

"""
One pass: returns (added, removed), summed over all ranks. The two particles closest to each other always
tie in the removal: the C++ code breaks the tie by the order in the rank (`cpp_order = true`), and here it is
broken by the position, so the result does not depend on the number of ranks.
"""
function add_remove_particles!(m, particle_swarm, da_Veloc, step; location = :native, cpp_order = false)
    dx, dz, e = element_size_x(m), element_size_z(m), m.boundary_offset
    coor = with_field(copy, particle_swarm, PIC_COORDINATES, Float64)
    n = length(coor) ÷ 2
    buckets = Dict{Tuple{Int, Int}, Vector{Int}}()
    for p in 1:n
        x = clamp(coor[2p-1], e, m.Lx - e); z = clamp(coor[2p], e - m.domain_depth, -e)
        c = (clamp(trunc(Int, x / dx), 0, m.Nx - 2), clamp(trunc(Int, (z + m.domain_depth) / dz), 0, m.Nz - 2))
        push!(get!(() -> Int[], buckets, c), p)
    end
    ppe = particles_per_element(m)                                  # particles_per_ele
    max_particles, min_particles = ppe + ppe ÷ 10 + 2, ppe - ppe ÷ 10 - 2
    xs, zs = particle_elements(m, da_Veloc, location)
    cells = m.periodic_boundary ? collect(keys(buckets)) : vec([(i, k) for i in xs, k in zs])
    union!(cells, keys(buckets))
    sort!(cells; by = reverse)                                      # k outer, i inner
    cpp_order || foreach(b -> sort!(b; by = p -> (coor[2p-1], coor[2p])), values(buckets))
    particles_to_remove = Int[]                                     # p_remove
    new_coords, sources = NTuple{2, Float64}[], Int[]               # p_add_coor and source particle
    for (i, k) in cells
        members = get(buckets, (i, k), Int[]); count = length(members)
        if count > max_particles                                    # remove the one with the closest neighbor
            best, chosen = 1.0e40, 0
            for p in members
                d = 1.0e40
                for q in members
                    q == p && continue
                    d = min(d, (coor[2p-1] - coor[2q-1])^2 + (coor[2p] - coor[2q])^2)
                end
                d < best && ((best, chosen) = (d, p))
            end
            push!(particles_to_remove, chosen); count -= 1
        end
        if count < min_particles                                    # candidate farthest from the existing ones
            near = nearby_particles(buckets, i, k)
            isempty(near) && error("element ($i, $k) and its neighbors have no particles")
            best, chosen, source = 0.0, (0.0, 0.0), 0
            for (c, (x, z)) in enumerate(candidate_points(m, step, i, k))
                d, closest = 1.0e30, 0
                for (o, p) in enumerate(near)
                    d2 = (x - coor[2p-1])^2 + (z - coor[2p])^2 + (cpp_order ? 1e-12 * ((c - 1) + (o - 1) * 10.0) : 0.0)
                    d2 < d && ((d, closest) = (d2, p))
                end
                d > best && ((best, chosen, source) = (d, (x, z), closest))
            end
            push!(new_coords, chosen); push!(sources, source)
        end
    end
    # properties of the source particles, read before the removal
    properties = [(name, T, with_field(a -> a[sources], particle_swarm, name, T)) for (name, T) in
        (("layer", PInt), ("geoq_fac", Float64), ("strain_fac", Float64), ("strain_rate_fac", Float64))]
    for p in sort(particles_to_remove; rev = true)
        LibPETSc.DMSwarmRemovePointAtIndex(petsclib, particle_swarm, PInt(p - 1))
    end
    if !isempty(new_coords)
        old = Int(LibPETSc.DMSwarmGetLocalSize(petsclib, particle_swarm))
        LibPETSc.DMSwarmAddNPoints(petsclib, particle_swarm, PInt(length(new_coords)))
        with_field(particle_swarm, PIC_COORDINATES, Float64) do a
            for (q, (x, z)) in enumerate(new_coords); a[2(old+q)-1] = x; a[2(old+q)] = z; end
        end
        for (name, T, values) in properties
            with_field(a -> (a[old+1:old+length(values)] = values; nothing), particle_swarm, name, T)
        end
        with_field(a -> (a[old+1:end] .= 0; nothing), particle_swarm, "itag", PInt)   # 0: not tracked
    end
    return Tuple(MPI.Allreduce([length(new_coords), length(particles_to_remove)], +, MPI.COMM_WORLD))
end
