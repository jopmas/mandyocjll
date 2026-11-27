# 3D particles and coupled model (DMSwarm_3d.cpp, moveSwarm, Swarm2Mesh_3d, Swarm_add_remove_3d, calc_pressure_3d,
# shift_pressure_3d). Include after rheology.jl. As in 2D, the elements also have the diffusivity of their layer
# (the C++ code does not compute geoq_kappa in 3D) and the shear heating.

trilinear_weights(rx, ry, rz) = ((1 - rx) * (1 - ry) * (1 - rz), rx * (1 - ry) * (1 - rz), (1 - rx) * ry * (1 - rz), rx * ry * (1 - rz),
    (1 - rx) * (1 - ry) * rz, rx * (1 - ry) * rz, (1 - rx) * ry * rz, rx * ry * rz)

"Element (i, j, k) and position inside it, after keeping the point inside the domain (moveSwarm in 3D)."
function locate3(m, x, y, z)
    e = m.boundary_offset
    x = clamp(x, e, m.Lx - e); y = clamp(y, e, m.Ly - e); z = clamp(z, e - m.domain_depth, -e)
    dx, dy, dz = element_size_x(m), element_size_y(m), element_size_z(m)
    i = clamp(trunc(Int, x / dx), 0, m.Nx - 2); j = clamp(trunc(Int, y / dy), 0, m.Ny - 2)
    k = clamp(trunc(Int, (z + m.domain_depth) / dz), 0, m.Nz - 2)
    return i, j, k, x / dx - i, y / dy - j, (z + m.domain_depth) / dz - k, z
end

# ---- swarm (createSwarm_3d) -------------------------------------------------------------
function create_particle_swarm3d(m::ModelParameters, da_Veloc; init)
    si, sj, sk, milocal, mjlocal, mklocal = Int.(LibPETSc.DMDAGetCorners(P, da_Veloc))
    CELL_CONTEXT[] = CellContext(da_Veloc.ptr, si, sk, milocal, mklocal, m, sj, mjlocal)
    cell_dm = LibPETSc.DMShellCreate(P, MPI.COMM_WORLD)
    LibPETSc.DMSetDimension(P, cell_dm, PInt(3))
    shim = Libdl.dlsym(Libdl.dlopen(SHIM_LIBRARY), :mandyoc_set_cell_dm_ops)
    loc = @cfunction(locate_points, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cvoid}))
    nei = @cfunction(get_neighbors, Cint, (Ptr{Cvoid}, Ptr{PInt}, Ptr{Ptr{Cint}}))
    @ccall $shim(cell_dm.ptr::Ptr{Cvoid}, loc::Ptr{Cvoid}, nei::Ptr{Cvoid})::Cint
    swarm = LibPETSc.DMCreate(P, MPI.COMM_WORLD)
    LibPETSc.DMSetType(P, swarm, "swarm")
    LibPETSc.DMSetDimension(P, swarm, PInt(3))
    LibPETSc.DMSwarmSetType(P, swarm, LibPETSc.DMSWARM_PIC)
    LibPETSc.DMSwarmSetCellDM(P, swarm, cell_dm)
    for (name, type) in SWARM_FIELDS
        LibPETSc.DMSwarmRegisterPetscDatatypeField(P, swarm, name, PInt(1), type)
    end
    LibPETSc.DMSwarmFinalizeFieldRegister(P, swarm)
    dx, dy, dz = element_size_x(m), element_size_y(m), element_size_z(m)
    nx, ny, nz = m.particles_per_element_x, m.particles_per_element_y, m.particles_per_element_z
    ppe = particles_per_element(m)
    xs, ys, zs = si:min(si + milocal, m.Nx - 1)-1, sj:min(sj + mjlocal, m.Ny - 1)-1, sk:min(sk + mklocal, m.Nz - 1)-1
    cells = [(i, j, k) for k in zs for j in ys for i in xs]
    LibPETSc.DMSwarmSetLocalSizes(P, swarm, PInt(length(cells) * ppe), PInt(4))
    with_field(swarm, PIC_COORDINATES, Float64) do coor
    with_field(swarm, "itag", PInt) do itag
    with_field(swarm, "layer", PInt) do layer
    with_field(swarm, "geoq_fac", Float64) do viscosity_factor
    with_field(swarm, "strain_fac", Float64) do strain
        p = 0
        for (i, j, k) in cells, cx in 0:nx-1, cy in 0:ny-1, cz in 0:nz-1
            p += 1
            id = ((k * (m.Ny - 1) + j) * (m.Nx - 1) + i) * ppe + (cx * ny + cy) * nz + cz + 1
            r(c) = unit_random(m.random_seed, id, c) * m.particles_perturb_factor
            x = i * dx + (0.5 + cx + r(1)) * dx / nx; y = j * dy + (0.5 + cy + r(2)) * dy / ny
            z = -m.domain_depth + k * dz + (0.5 + cz + r(3)) * dz / nz
            coor[3p-2] = x; coor[3p-1] = y; coor[3p] = z
            itag[p] = id; strain[p] = 0.0
            layer[p], viscosity_factor[p] = init(x, y, z, id)
        end
    end end end end end
    LibPETSc.DMSwarmMigrate(P, swarm, LibPETSc.PETSC_TRUE)
    return swarm, cell_dm
end

# ---- RK4 in 3D -----------------------------------------------------------------------------
"Velocity at (x, y, z), kept inside the domain as the particles (locate3); V from local_view3 of da_advection."
function interpolate_velocity3(m, V, x, y, z)
    i, j, k, rx, ry, rz, _ = locate3(m, x, y, z)
    w8 = trilinear_weights(rx, ry, rz); nodes = element_nodes3(i, j, k)
    return ntuple(c -> sum(w8[q] * V[c, nodes[q]...] for q in 1:8), 3)
end

"Velocity at the time t_n + f*dt of the step: linear between V0 (t_n, Veloc) and V1 (t_n+1, Veloc_fut)."
function velocity_in_time3(m, V0, V1, x, y, z, f)
    a = interpolate_velocity3(m, V0, x, y, z)
    V1 === V0 && return a
    b = interpolate_velocity3(m, V1, x, y, z)
    return (1 - f) .* a .+ f .* b
end

"""
Mean velocity of the step by the fourth-order Runge-Kutta in 3D (x_n+1 = x_n + dt*v), with the velocity varying
linearly in time from V0 (t_n) to V1 (t_n+1). A stage outside the domain is stopped at the boundary (locate3); V0 and
V1 need the ghosts of da_advection (stencil 2) for the stages in the neighbor element.
"""
function runge_kutta_4_velocity3(m, V0, V1, x, y, z, dt)
    A = velocity_in_time3(m, V0, V1, x, y, z, 0.0)
    B = velocity_in_time3(m, V0, V1, x + A[1] * dt / 2, y + A[2] * dt / 2, z + A[3] * dt / 2, 0.5)
    C = velocity_in_time3(m, V0, V1, x + B[1] * dt / 2, y + B[2] * dt / 2, z + B[3] * dt / 2, 0.5)
    D = velocity_in_time3(m, V0, V1, x + C[1] * dt, y + C[2] * dt, z + C[3] * dt, 1.0)
    return (A .+ 2 .* B .+ 2 .* C .+ D) ./ 6
end

# ---- moveSwarm in 3D: T, P, strain, viscosity and (Euler or RK4) position -------------
"""
As update_particles! in 2D: Euler with `local_velocity`, or RK4 with `runge_kutta_4 = (V0, V1)`, the local velocities
at the start and at the end of the step; all from `da_advection`.
"""
function update_particles3!(m, range, layers, swarm, da_advection, local_velocity, da_Thermal, local_temperature,
    local_pressure_thermal_mesh, dt; move = false, runge_kutta_4 = nothing)
    dx, dy, dz = element_size_x(m), element_size_y(m), element_size_z(m)
    e_b = m.boundary_offset
    extra = runge_kutta_4 === nothing ? () : runge_kutta_4
    PETSc.with_local_array!(local_velocity, local_temperature, local_pressure_thermal_mesh, extra...; write = false) do v, t, pr, rk...
        V, T, Pn = local_view3(da_advection, v, 3), local_view3(da_Thermal, t), local_view3(da_Thermal, pr)
        stages = isempty(rk) ? nothing : map(a -> local_view3(da_advection, a, 3), rk)
        with_field(swarm, PIC_COORDINATES, Float64) do coor
        with_field(swarm, "layer", PInt) do layer
        with_field(swarm, "geoq_fac", Float64) do viscosity
        with_field(swarm, "strain_fac", Float64) do strain
        with_field(swarm, "strain_rate_fac", Float64) do strain_rate
            for p in eachindex(layer)
                i, j, k, rx, ry, rz, z = locate3(m, coor[3p-2], coor[3p-1], coor[3p])
                w8 = trilinear_weights(rx, ry, rz); nodes = element_nodes3(i, j, k)
                Tp = sum(w8[n] * T[nodes[n]...] for n in 1:8); Pp = sum(w8[n] * Pn[nodes[n]...] for n in 1:8)
                kx, ky, kz = 2rx - 1, 2ry - 1, 2rz - 1
                e = zeros(6)                                        # strain[0..5]
                n = 0
                for ez in (-1, 1), ey in (-1, 1), ex in (-1, 1)
                    n += 1
                    Nx = ex * (1 + ey * ky) * (1 + ez * kz) / 4 / dx
                    Ny = (1 + ex * kx) * ey * (1 + ez * kz) / 4 / dy
                    Nz = (1 + ex * kx) * (1 + ey * ky) * ez / 4 / dz
                    u, vv, w = V[1, nodes[n]...], V[2, nodes[n]...], V[3, nodes[n]...]
                    e[1] += Nx * u; e[2] += Ny * vv; e[3] += Nz * w
                    e[4] += Ny * u + Nx * vv; e[5] += Nz * vv + Ny * w; e[6] += Nx * w + Nz * u
                end
                mean = (e[1] + e[2] + e[3]) / 3
                e2 = sqrt(((e[1] - mean)^2 + (e[2] - mean)^2 + (e[3] - mean)^2 + (e[4]^2 + e[5]^2 + e[6]^2) / 2) / 2)
                if move && stages === nothing
                    coor[3p-2] += dt * sum(w8[q] * V[1, nodes[q]...] for q in 1:8)
                    coor[3p-1] += dt * sum(w8[q] * V[2, nodes[q]...] for q in 1:8)
                    coor[3p] += dt * sum(w8[q] * V[3, nodes[q]...] for q in 1:8)
                elseif move        # RK4 from the position at t_n (inside the domain); T, P and strain rate stay there
                    x0 = clamp(coor[3p-2], e_b, m.Lx - e_b); y0 = clamp(coor[3p-1], e_b, m.Ly - e_b)
                    u = runge_kutta_4_velocity3(m, stages..., x0, y0, z, dt)
                    coor[3p-2] += dt * u[1]; coor[3p-1] += dt * u[2]; coor[3p] += dt * u[3]
                end
                strain[p] += dt * e2
                strain_rate[p] = e2
                if m.use_strain_healing
                    s = strain[p] - dt * m.strain_healing_rate * exp(-m.strain_healing_activation_factor *
                        (m.strain_healing_reference_temperature + 273) / (Tp + 273))
                    strain[p] = max(s, 0.0)
                end
                viscosity[p] = particle_viscosity(m, range, layers, layer[p] + 1, Tp, Pp, z, e2, strain[p])
            end
        end end end end end
    end
    move && LibPETSc.DMSwarmMigrate(P, swarm, LibPETSc.PETSC_TRUE)
    return (0.0, -m.domain_depth)
end

# ---- Swarm2Mesh_3d ----------------------------------------------------------------------
function particles_to_mesh3!(m, layers, swarm, da_Thermal, f::MeshFields)
    lv() = (v = LibPETSc.DMCreateLocalVector(P, da_Thermal); LibPETSc.VecZeroEntries(P, v); v)
    L = ntuple(_ -> lv(), 8)             # rho, H, strain rate, weight | viscosity, kappa, count | strain
    harmonic = m.use_viscosity_harmonic_mean
    PETSc.with_local_array!(L...) do arrays...
        R, H, S, W, Vi, Ka, N, E = map(a -> local_view3(da_Thermal, a), arrays)
        with_field(swarm, PIC_COORDINATES, Float64) do coor
        with_field(swarm, "layer", PInt) do layer
        with_field(swarm, "geoq_fac", Float64) do visc
        with_field(swarm, "strain_rate_fac", Float64) do rate
        with_field(swarm, "strain_fac", Float64) do strain
            for p in eachindex(layer)
                i, j, k, rx, ry, rz, _ = locate3(m, coor[3p-2], coor[3p-1], coor[3p])
                n = layer[p] + 1
                h = layers.layer_radiogenic_heat[n]
                m.use_shear_heating && (h += m.shear_heating_scale_factor * 4visc[p] * rate[p]^2 / layers.layer_density[n])
                for (w, node) in zip(trilinear_weights(rx, ry, rz), element_nodes3(i, j, k))
                    R[node...] += w * layers.layer_density[n]; H[node...] += w * h; S[node...] += w * rate[p]
                    E[node...] += w * strain[p]; W[node...] += w
                end
                kap = layers.layer_thermal_conductivity[n]
                Vi[i, j, k] += harmonic ? 1 / visc[p] : visc[p]; Ka[i, j, k] += harmonic ? 1 / kap : kap; N[i, j, k] += 1
            end
        end end end end end
    end
    g = LibPETSc.DMCreateGlobalVector(P, da_Thermal); d = LibPETSc.DMCreateGlobalVector(P, da_Thermal)
    for (num, den, out, inv) in ((L[1], L[4], f.density_field, false), (L[2], L[4], f.radiogenic_heat_field, false),
        (L[3], L[4], f.strain_rate_field, false), (L[8], L[4], f.strain_field, false), (L[5], L[7], f.viscosity_factor, harmonic),
        (L[6], L[7], f.thermal_diffusivity_field, harmonic))
        for (a, b) in ((num, g), (den, d))
            LibPETSc.VecZeroEntries(P, b); LibPETSc.DMLocalToGlobal(P, da_Thermal, a, LibPETSc.ADD_VALUES, b)
        end
        PETSc.with_local_array!(out, g, d) do o, a, b
            for n in eachindex(o)        # no particles: keep the previous value
                b[n] > 0 && (o[n] = inv ? b[n] / a[n] : a[n] / b[n])
            end
        end
    end
    foreach(v -> LibPETSc.VecDestroy(P, v), (L..., g, d))
end

"""
Air in 3D: zero temperature where the density is low, and the viscous layer (sticky_blanket_air) in the elements
whose mean density is between 5 and 2160 kg/m3, as in 2D (the 3D C++ code uses the sums of the particles before
they are divided by the weights, bug 30).
"""
function air_corrections3!(m, da_Thermal, f::MeshFields, temperature)
    PETSc.with_local_array!(temperature, f.density_field; write = (true, false)) do T, rho
        for n in eachindex(T); rho[n] < m.air_threshold_density && (T[n] = 0.0); end
    end
    m.sticky_blanket_air || return
    l = LibPETSc.DMCreateLocalVector(P, da_Thermal)
    global_to_local(da_Thermal, f.density_field, l)
    xs, ys, zs = element_ranges3(m, da_Thermal)
    PETSc.with_local_array!(f.viscosity_factor, l; write = (true, false)) do va, ra
        V = owned_view3(da_Thermal, va); R = local_view3(da_Thermal, ra)
        for k in zs, j in ys, i in xs
            mean = sum(R[n...] for n in element_nodes3(i, j, k)) / 8
            5.0 < mean < 2700.0 * 0.8 && (V[i, j, k] = max(V[i, j, k], 50m.min_viscosity))
        end
    end
    LibPETSc.VecDestroy(P, l)
end

# ---- pressure -----------------------------------------------------------------------------
"Element pressure -> nodes (mean of the elements around), with zero mean pressure at the top (shift_pressure_3d)."
function pressure_to_nodes3!(m, da_Veloc, da_Thermal, pressure, local_pressure, pressure_thermal_mesh)
    global_to_local(da_Veloc, pressure, local_pressure)
    top = 0.0; count = 0
    PETSc.with_local_array!(pressure_thermal_mesh, local_pressure; write = (true, false)) do a, pl
        N = owned_view3(da_Thermal, a); E = local_view3(da_Veloc, pl, 3)
        for k in axes(N, 3), j in axes(N, 2), i in axes(N, 1)
            s = 0.0; c = 0
            for kk in (k, k - 1), jj in (j, j - 1), ii in (i, i - 1)
                0 <= ii < m.Nx - 1 && 0 <= jj < m.Ny - 1 && 0 <= kk < m.Nz - 1 && (s += E[1, ii, jj, kk]; c += 1)
            end
            N[i, j, k] = s / c
            k == m.Nz - 1 && (top += N[i, j, k]; count += 1)
        end
    end
    mean = MPI.Allreduce(top, +, MPI.COMM_WORLD) / MPI.Allreduce(count, +, MPI.COMM_WORLD)
    PETSc.with_local_array!(a -> a .-= mean, pressure_thermal_mesh)
end

"Initial lithostatic pressure of each element, from the interfaces at its 4 corners in x and y (calc_pressure_3d)."
function lithostatic_pressure3!(m, layers, da_Veloc, pressure)
    d = layers.interface_depths; n = size(d, 2); rho, g = layers.layer_density, m.gravity_acceleration
    node(i, j) = i + j * m.Nx + 1
    fill_field3!(pressure, da_Veloc, function (i, j, k)
        (i >= m.Nx - 1 || j >= m.Ny - 1 || k >= m.Nz - 1) && return (0.0, 0.0, 0.0)
        z = -(m.Nz - 1 - k) * m.domain_depth / (m.Nz - 1)
        I(q) = (d[node(i, j), q] + d[node(i + 1, j), q] + d[node(i, j + 1), q] + d[node(i + 1, j + 1), q]) / 4
        z > I(n) && return (rho[n+1] * g * -z, 0.0, 0.0)
        p = rho[n+1] * g * -I(n)
        for q in n:-1:2
            z > I(q - 1) && return (p + rho[q] * g * (I(q) - z), 0.0, 0.0)
            p += rho[q] * g * (I(q) - I(q - 1))
        end
        (p + rho[1] * g * (I(1) - z), 0.0, 0.0)
    end, 3)
end

# ---- Swarm_add_remove_3d -----------------------------------------------------------------
function candidate_points3(m, step, i, j, k)
    dx, dy, dz = element_size_x(m), element_size_y(m), element_size_z(m)
    id = UInt64((k * (m.Ny - 1) + j) * (m.Nx - 1) + i)
    u(c) = 2 * (unit_random(UInt64(step), id, c) + 0.5) - 1
    return [(i * dx + (0.5u(3q) + 0.5) * dx, j * dy + (0.5u(3q + 1) + 0.5) * dy, k * dz - m.domain_depth + (0.5u(3q + 2) + 0.5) * dz)
        for q in 1:10]
end

"""
One pass: at most one particle added or removed per element, as in 2D. The choice of the removed particle and of
the candidates does not depend on the number of ranks (the C++ code uses one random sequence per rank).
"""
function add_remove_particles3!(m, swarm, step)
    coor = with_field(copy, swarm, PIC_COORDINATES, Float64)
    n = length(coor) ÷ 3
    buckets = Dict{NTuple{3, Int}, Vector{Int}}()
    for p in 1:n
        i, j, k = locate3(m, coor[3p-2], coor[3p-1], coor[3p])
        push!(get!(() -> Int[], buckets, (i, j, k)), p)
    end
    ppe = particles_per_element(m)
    max_particles, min_particles = ppe + ppe ÷ 10 + 2, ppe - ppe ÷ 10 - 2
    foreach(b -> sort!(b; by = p -> (coor[3p-2], coor[3p-1], coor[3p])), values(buckets))
    d2(p, x, y, z) = (coor[3p-2] - x)^2 + (coor[3p-1] - y)^2 + (coor[3p] - z)^2
    to_remove, new_coords, sources = Int[], NTuple{3, Float64}[], Int[]
    for cell in sort!(collect(keys(buckets)); by = reverse)
        members = buckets[cell]; count = length(members)
        if count > max_particles
            best, chosen = 1.0e40, 0
            for p in members
                d = minimum(d2(q, coor[3p-2], coor[3p-1], coor[3p]) for q in members if q != p)
                d < best && ((best, chosen) = (d, p))
            end
            push!(to_remove, chosen); count -= 1
        end
        if count < min_particles
            best, chosen, source = 0.0, (0.0, 0.0, 0.0), 0
            for (x, y, z) in candidate_points3(m, step, cell...)
                d, closest = findmin(p -> d2(p, x, y, z), members)
                d > best && ((best, chosen, source) = (d, (x, y, z), members[closest]))
            end
            push!(new_coords, chosen); push!(sources, source)
        end
    end
    properties = [(name, T, with_field(a -> a[sources], swarm, name, T)) for (name, T) in
        (("layer", PInt), ("geoq_fac", Float64), ("strain_fac", Float64), ("strain_rate_fac", Float64))]
    for p in sort(to_remove; rev = true)
        LibPETSc.DMSwarmRemovePointAtIndex(P, swarm, PInt(p - 1))
    end
    if !isempty(new_coords)
        old = Int(LibPETSc.DMSwarmGetLocalSize(P, swarm))
        LibPETSc.DMSwarmAddNPoints(P, swarm, PInt(length(new_coords)))
        with_field(swarm, PIC_COORDINATES, Float64) do a
            for (q, (x, y, z)) in enumerate(new_coords); a[3(old+q)-2] = x; a[3(old+q)-1] = y; a[3(old+q)] = z; end
        end
        for (name, T, values) in properties
            with_field(a -> (a[old+1:old+length(values)] = values; nothing), swarm, name, T)
        end
        with_field(a -> (a[old+1:end] .= 0; nothing), swarm, "itag", PInt)
    end
    return Tuple(MPI.Allreduce([length(new_coords), length(to_remove)], +, MPI.COMM_WORLD))
end

# ---- model ------------------------------------------------------------------------------------
function Model3D(m, layers; temperature, velocity, velocity_mask = nothing, initial_strain = (x, y, z, layer) -> 0.0)
    dv, dT = create_velocity_mesh3d(m), create_thermal_mesh3d(m); da = create_advection_mesh3d(m, dv)
    n = size(layers.interface_depths, 2)
    function layer_at(x, y, z, id)
        i = clamp(trunc(Int, x / element_size_x(m)), 0, m.Nx - 2); rx = x / element_size_x(m) - i
        j = clamp(trunc(Int, y / element_size_y(m)), 0, m.Ny - 2); ry = y / element_size_y(m) - j
        d = layers.interface_depths; node(a, b) = a + b * m.Nx + 1
        for q in 1:n
            h = d[node(i, j), q] * (1 - rx) * (1 - ry) + d[node(i + 1, j), q] * rx * (1 - ry) +
                d[node(i, j + 1), q] * (1 - rx) * ry + d[node(i + 1, j + 1), q] * rx * ry
            z < h && return (q - 1, layers.layer_viscosity_factor[q])
        end
        return (n, layers.layer_viscosity_factor[n+1])
    end
    swarm, cell = create_particle_swarm3d(m, dv; init = layer_at)
    with_field(swarm, "strain_fac", Float64) do s
        with_field(swarm, "layer", PInt) do l
        with_field(swarm, PIC_COORDINATES, Float64) do c
            for p in eachindex(s)
                seed = layers.layer_weakening_seed[l[p]+1]
                s[p] = seed >= 0 ? seed : initial_strain(c[3p-2], c[3p-1], c[3p], Int(l[p]))
            end
        end end
    end
    zero_vec(d) = (v = LibPETSc.DMCreateGlobalVector(P, d); LibPETSc.VecZeroEntries(P, v); v)
    T = fill_field3!(zero_vec(dT), dT, temperature); V = fill_field3!(zero_vec(dv), dv, velocity, 3); Vn = zero_vec(dv)
    LibPETSc.VecCopy(P, V, Vn)
    pr = zero_vec(dv); pn = zero_vec(dT)
    stokes = StokesSolver3D(m, dv, dT, V; mask = velocity_mask)
    if n > 0
        lithostatic_pressure3!(m, layers, dv, pr)
        pressure_to_nodes3!(m, dv, dT, pr, stokes.local_pressure, pn)
    end
    # advection mesh with stencil 2 (RK4 stages), as in 2D; the local velocities of the advection and of RK4 use it
    return Model(m, layers, ViscosityRange(m), dv, dT, da, swarm, cell, MeshFields(dT), T, V, Vn, zero_vec(dv), pr, pn,
        LibPETSc.DMCreateLocalVector(P, da), LibPETSc.DMCreateLocalVector(P, da), LibPETSc.DMCreateLocalVector(P, da),
        LibPETSc.DMCreateLocalVector(P, dT), LibPETSc.DMCreateLocalVector(P, dT),
        stokes, ThermalSolver3D(m, dT, dv, T), false, 0.0, 0.0, 0, (0.0, -m.domain_depth))
end

"Time step in 3D: 10% of the time to cross the element along the component of the largest velocity (Calc_dt_calor)."
function time_step3(s::Model)
    c = zeros(3)
    PETSc.with_local_array!(s.velocity_next; write = false) do a
        for q in 1:3; c[q] = maximum(abs, @view a[q:3:end]; init = 0.0); end
    end
    c = MPI.Allreduce(c, max, MPI.COMM_WORLD)
    q = argmax(c); h = (element_size_x(s.m), element_size_y(s.m), element_size_z(s.m))[q]
    return min(0.1h / (c[q] * s.m.time_step_subdivision), s.m.max_time_step * SECONDS_PER_YEAR)
end
