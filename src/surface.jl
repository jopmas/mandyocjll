# Surface tracking and surface processes (sp.cpp and the steps of main.cpp that call it). Include after rheology3d.jl
# and fastscape.jl.
#
# 2D: the C++ code keeps the surface in a DMSwarm (dms_s) split among the ranks and gathers it on rank 0 for each
# process. Here the surface is a plain array replicated on every rank: the advection is done by the rank that owns
# each point and summed, and the rest is computed identically on all ranks. With one rank the steps are those of the
# C++ code; with more ranks this avoids the problems of the split surface listed in docs/bugs_mandyoc_cpp.md.
#
# 3D (not in the C++ code): the surface is the height of a regular grid of points (surface_particles_per_element per
# element in x and y), replicated on every rank and moved by the velocity with a semi-Lagrangian step: the height
# of a point is the height at the position that the velocity brings to the point, plus the vertical displacement.

const SP_MODES = ("none", "diffusion", "sedimentation_only", "sedimentation_rate_limited", "theunissen_sedimentation", "fastscape")
const SP_MODES_3D = ("none", "diffusion", "fastscape")

mutable struct Surface{A <: AbstractArray{Float64}}
    x::Vector{Float64}                   # dms_s coordinates (2D); regular points in x (3D)
    y::Vector{Float64}                   # regular points in y (3D; empty in 2D)
    z::A                                 # height: a vector (2D) or a matrix x by y (3D)
    points_per_element::Int              # dms_s_ppe (surface_particles_per_element)
    processes::Bool                      # sp_surface_processes
    mode::String                         # sp_mode
    diffusivity::Float64                 # sp_d_c (m2/s)
    sea_level::Float64                   # sea_level (m); changed by base_level.txt
    air_to_land::Bool                    # a2l
    continental_slope::Float64           # tangent of continental_slope (degrees in param.txt)
    sediment_strain::Float64             # strain_sed
    sedimentation_rate::Float64          # sedimentation_rate (m2/yr)
    active_sediment_layer::Int           # active_sediment_layer
    sediment_layer_changes::Vector{Tuple{Float64, Int}}     # sediment_layers.txt: (Myr, layer)
    sedimentation_rate_changes::Vector{Tuple{Float64, Float64}}  # sedimentation_rate.txt: (Myr, m2/yr)
    base_level_changes::Vector{Tuple{Float64, Float64}}    # base_level.txt: (Myr, m)
    fastscape::Union{Nothing, FastScapeCoupling}            # sp_mode = fastscape
end
const Surface2D = Surface{Vector{Float64}}
const Surface3D = Surface{Matrix{Float64}}

"Pairs of numbers of a file whose first number is the count (sediment_layers.txt, sedimentation_rate.txt, base_level.txt)."
function read_changes(file, T; header = 1)
    isfile(file) || error("$(basename(file)) not found")
    v = split(read(file, String))
    n = parse(Int, v[1])
    length(v) >= header + 2n || error("$(basename(file)): $n changes announced, $((length(v) - header) ÷ 2) found")
    first_values = v[2:header]
    return [(parse(Float64, v[header+2i-1]), parse(T, v[header+2i])) for i in 1:n], first_values
end

"""
Surface of the scenario, or nothing without sp_surface_tracking (sp_create_surface_swarm_2d and the readers of the
files of the surface processes). The initial surface is the top interface, sampled at points_per_element per element.
"""
function Surface(d::Dict{String, String}, m::ModelParameters, layers::LayerProperties, dir)
    on(k) = lowercase(get(d, k, "false")) == "true"
    num(k, x) = parse(Float64, get(d, k, string(x)))
    tracking, processes, mode = on("sp_surface_tracking"), on("sp_surface_processes"), get(d, "sp_mode", "none")
    mode in SP_MODES || error("unknown sp_mode = $mode; use one of " * join(SP_MODES, ", "))
    processes && mode == "none" && error("sp_surface_processes = True needs an sp_mode other than none")
    !processes && mode != "none" && error("sp_mode = $mode needs sp_surface_processes = True")
    # as reader.cpp: the surface processes switch on the surface tracking
    if processes && !tracking
        tracking = true
        MPI.Comm_rank(MPI.COMM_WORLD) == 0 && println("Activating sp_surface_tracking since sp_surface_processes is active")
    end
    tracking || return nothing
    three = m.dimensions == 3
    three && !(mode in SP_MODES_3D) && error("sp_mode = $mode is not available in 3D; use one of " * join(SP_MODES_3D, ", "))
    ppe = Int(num("surface_particles_per_element", 2))
    dx, e, top = element_size_x(m), m.boundary_offset, layers.interface_depths[:, end]
    x, y, z = Float64[], Float64[], Float64[]
    if three
        # regular points; the height is interpolated (bilinear) from the top interface
        x = collect(range(0, m.Lx, length = (m.Nx - 1) * ppe + 1)); y = collect(range(0, m.Ly, length = (m.Ny - 1) * ppe + 1))
        dy = element_size_y(m); node(a, b) = a + b * m.Nx + 1
        z = [begin
            i = min(trunc(Int, px / dx), m.Nx - 2); j = min(trunc(Int, py / dy), m.Ny - 2)
            rx = px / dx - i; ry = py / dy - j
            (1 - rx) * (1 - ry) * top[node(i, j)] + rx * (1 - ry) * top[node(i + 1, j)] +
                (1 - rx) * ry * top[node(i, j + 1)] + rx * ry * top[node(i + 1, j + 1)]
        end for px in x, py in y]
    else
        for p in 0:m.Nx-1, i in 0:ppe-1
            rx = (i * dx) / ppe; px = p * dx + rx
            px <= m.Lx || continue
            pz = p < m.Nx - 1 ? (1.0 - rx / dx) * top[p+1] + (rx / dx) * top[p+2] : top[p+1]
            push!(x, px == 0 ? px + e : px == m.Lx ? px - e : px); push!(z, pz)
        end
    end
    layer_changes, rate_changes, base_changes, active = Tuple{Float64, Int}[], Tuple{Float64, Float64}[], Tuple{Float64, Float64}[], -1
    if processes        # read whatever the mode, as in the C++ code
        layer_changes, first = read_changes(joinpath(dir, "sediment_layers.txt"), Int; header = 2)
        active = parse(Int, first[1])
    end
    if mode in ("sedimentation_rate_limited", "theunissen_sedimentation")
        rate_changes, _ = read_changes(joinpath(dir, "sedimentation_rate.txt"), Float64)
    end
    if on("variable_base_level")
        base_changes, _ = read_changes(joinpath(dir, "base_level.txt"), Float64)
    end
    # the FastScape grid is the grid of the regular points (in 2D those of interpolate_surface!, from 0 to Lx)
    fastscape = mode == "fastscape" ? FastScapeCoupling(FastScapeParameters(d), z, m.Lx, m.Ly) : nothing
    return Surface(x, y, z, ppe, processes, mode, num("sp_d_c", 0.0), num("sea_level", 0.0), lowercase(get(d, "a2l", "true")) == "true",
        tan(num("continental_slope", 2.0) * pi / 180), num("strain_sed", 0.6), 0.0, active, layer_changes, rate_changes, base_changes,
        fastscape)
end

"""
Changes of the files at the beginning of a step (sp_update_active_sediment_layer, sp_update_base_level,
sp_update_sedimentation_rate): one change of each per step, when the time at the end of the step passes it.
"""
function update_surface_conditions!(sf::Surface, time_years, say)
    passed(c) = !isempty(c) && time_years > 1.0e6 * c[1][1]
    sf.processes && passed(sf.sediment_layer_changes) && (sf.active_sediment_layer = popfirst!(sf.sediment_layer_changes)[2])
    if passed(sf.base_level_changes)
        sf.sea_level = popfirst!(sf.base_level_changes)[2]
        say("base level = $(sf.sea_level) m")
    end
    if sf.mode in ("sedimentation_rate_limited", "theunissen_sedimentation") && passed(sf.sedimentation_rate_changes)
        sf.sedimentation_rate = popfirst!(sf.sedimentation_rate_changes)[2]
        say("sedimentation rate = $(sf.sedimentation_rate) m2/yr, active sediment layer = $(sf.active_sediment_layer)")
    end
end

"Moves the surface with the velocity of the substep (sp_move_surface_swarm_advection_2d)."
# runge_kutta_4 = (V0, V1): local velocities at t_n and t_n+1 on da_advection (RK4 instead of Euler)
function advect_surface!(sf::Surface2D, m, da_advection, local_velocity, dt; runge_kutta_4 = nothing)
    xs, zs, _, xm, zm, _ = LibPETSc.DMDAGetCorners(petsclib, da_advection)
    dx, dz, e = element_size_x(m), element_size_z(m), m.boundary_offset
    shift = zeros(2length(sf.x))
    extra = runge_kutta_4 === nothing ? () : runge_kutta_4
    PETSc.with_local_array!(local_velocity, extra...; write = false) do v, rk...
        V = local_view(da_advection, v, 2)
        stages = map(a -> local_view(da_advection, a, 2), rk)
        for p in eachindex(sf.x)
            cx = sf.x[p] >= m.Lx ? m.Lx - e : sf.x[p] <= 0 ? e : sf.x[p]
            cz = sf.z[p] >= 0 ? -e : sf.z[p] <= -m.domain_depth ? -m.domain_depth + e : sf.z[p]
            i = min(trunc(Int, cx / dx), m.Nx - 2); k = min(trunc(Int, (cz + m.domain_depth) / dz), m.Nz - 2)
            xs <= i < xs + xm && zs <= k < zs + zm || continue          # the rank that owns the element
            if !isempty(stages)
                u, w = runge_kutta_4_velocity(m, stages..., cx, cz, dt)
                shift[2p-1] = dt * u; shift[2p] = dt * w
                continue
            end
            rx = (cx - i * dx) / dx; rz = (cz - (-m.domain_depth + k * dz)) / dz
            shift[2p-1] = dt * linear_interpolation(rx, rz, V[1, i, k], V[1, i+1, k], V[1, i, k+1], V[1, i+1, k+1])
            shift[2p] = dt * linear_interpolation(rx, rz, V[2, i, k], V[2, i+1, k], V[2, i, k+1], V[2, i+1, k+1])
        end
    end
    shift = MPI.Allreduce(shift, +, MPI.COMM_WORLD)
    sf.x .+= @view shift[1:2:end]; sf.z .+= @view shift[2:2:end]
end

"""
Samples the moved surface back on the regular points (sp_surface_swarm_interpolation): linear interpolation between
the closest points on each side. The C++ code skips a point exactly at the sampled position (strict inequalities),
which averages its neighbors; here that point is used.
"""
function interpolate_surface!(sf::Surface2D, m)
    dx, ppe = element_size_x(m), sf.points_per_element
    order = sortperm(sf.x)               # stable: among equal positions, the first one, as in the C++ code
    sx, sz = sf.x[order], sf.z[order]
    x, z = Float64[], Float64[]
    for i in 0:m.Nx-1, j in 0:ppe-1
        rx = (j * dx) / ppe; px = i * dx + rx
        px > m.Lx && break
        r = searchsortedfirst(sx, px; lt = <=)          # first point at the right (cx > px)
        l = searchsortedlast(sx, px)                    # last point at the left or at px (cx <= px)
        l >= 1 && (l = searchsortedfirst(sx, sx[l]))    # first of equal positions
        if l >= 1 && sx[l] == px
            pz = sz[l]
        elseif l >= 1 && r <= length(sx)
            pz = (1.0 - (px - sx[l]) / (sx[r] - sx[l])) * sz[l] + (px - sx[l]) / (sx[r] - sx[l]) * sz[r]
        else
            pz = l >= 1 ? sz[l] : sz[r]                 # beyond the last point: the height of the closest one
        end
        push!(x, px); push!(z, pz)
    end
    sf.x, sf.z = x, z
end

"Mean height of the first quarter of the surface plus the sea level (sp_evaluate_adjusted_mean_elevation_with_sea_level)."
adjusted_mean_elevation(sf::Surface2D) = foldl(+, @view(sf.z[1:length(sf.z)÷4]); init = 0.0) / (length(sf.z) ÷ 4) + sf.sea_level  # in order, as in C++

"Surface processes of the step (sp_evaluate_surface_processes_2d); dt in seconds."
function surface_processes!(sf::Surface, m, dt)
    z, n = sf.z, length(sf.z)
    if sf.mode == "fastscape"
        fastscape_step!(sf.fastscape, sf.z, dt, adjusted_mean_elevation(sf))
    elseif sf isa Surface3D
        sf.mode == "diffusion" && diffuse_surface3!(sf, dt)
    elseif sf.mode == "diffusion"
        dxs = m.Lx / (n - 1)
        sp_dt = dxs * dxs / sf.diffusivity / 4.0
        steps = max(trunc(Int, floor(dt / sp_dt)) + 1, 10)
        r = sf.diffusivity * (dt / steps) / (dxs * dxs)
        z[1] = z[2]; z[n] = z[n-1]
        for _ in 1:steps
            old = copy(z)
            for j in 2:n-1
                z[j] = old[j] + r * (old[j+1] - 2.0 * old[j] + old[j-1])
            end
        end
    elseif sf.mode == "sedimentation_only"
        hsl = adjusted_mean_elevation(sf)
        z .= max.(z, hsl)
    elseif sf.mode == "sedimentation_rate_limited"
        rate_limited_sedimentation!(sf, m, dt)
    elseif sf.mode == "theunissen_sedimentation"
        theunissen_sedimentation!(sf, dt)
    end
end

"Fixed volume per step deposited from each margin, with a continental slope (sp_evaluate_surface_processes_2d_sedimentation_rate_limited)."
function rate_limited_sedimentation!(sf::Surface2D, m, dt)
    x, z, slope = sf.x, sf.z, sf.continental_slope
    hsl = adjusted_mean_elevation(sf)
    volume = sf.sedimentation_rate * dt / SECONDS_PER_YEAR        # m2
    dxs = x[2] - x[1]
    level(xj, brk, left) = left ? (xj < brk ? hsl : hsl - slope * (xj - brk)) : (xj > brk ? hsl : hsl + slope * (xj - brk))
    function fill_volume(brk, left)
        v = 0.0
        for j in eachindex(z)
            h = level(x[j], brk, left)
            z[j] < h && (v += (h - z[j]) * dxs)
        end
        return v
    end
    for left in (true, false)
        brk = left ? 0.0 : m.Lx
        while left ? brk < m.Lx : brk > 0
            fill_volume(brk, left) >= volume && break
            brk += left ? dxs : -dxs
        end
        for j in eachindex(z)
            h = level(x[j], brk, left)
            z[j] < h && (z[j] = (left ? x[j] < brk : x[j] > brk) ? hsl : z[j] + (h - z[j]))
        end
    end
end

"Progradation from both margins (Theunissen et al., 2019; sp_evaluate_surface_processes_2d_theunissen)."
function theunissen_sedimentation!(sf::Surface2D, dt)
    x, z = sf.x, sf.z
    Ld = 35000.0                         # (m) fixed in the C++ code
    hsl = adjusted_mean_elevation(sf)
    dxs = x[2] - x[1]; factor = dxs / Ld
    total = sf.sedimentation_rate * (dt / SECONDS_PER_YEAR)        # m2
    for nodes in (eachindex(z), reverse(eachindex(z)))
        flux = total
        for j in nodes
            h = (flux * factor) / dxs
            room = max(hsl - z[j], 0.0)
            deposit = h <= room ? h : room
            z[j] += deposit
            flux -= deposit * dxs
        end
    end
end

"""
Particles crossed by the surface (sp_update_surface_swarm_particles_properties): air below it becomes the active
sediment layer (a2l; `new_layer`, the magmatic layer after an extraction) and rock above it becomes air. As in the C++ code, the height is the mean of the two surface
points around the particle.
"""
function update_particles_with_surface!(sf::Surface2D, m, layers, swarm; new_layer = sf.active_sediment_layer)
    air = length(layers.layer_density) - 1                       # n_interfaces
    dxs = m.Lx / (length(sf.z) - 1)
    changed = 0
    with_field(swarm, PIC_COORDINATES, Float64) do coor
    with_field(swarm, "layer", PInt) do layer
    with_field(swarm, "geoq_fac", Float64) do factor
    with_field(swarm, "strain_fac", Float64) do strain
        for p in eachindex(layer)
            px, pz = coor[2p-1], coor[2p]
            i = clamp(trunc(Int, px / dxs), 0, length(sf.z) - 2)
            surface = sf.z[i+1] * 0.5 + sf.z[i+2] * 0.5
            if sf.air_to_land && layer[p] == air && pz < surface
                layer[p] = new_layer; factor[p] = layers.layer_viscosity_factor[new_layer+1]
                strain[p] = sf.sediment_strain; changed += 1
            end
            if layer[p] != air && pz > surface
                layer[p] = air; factor[p] = layers.layer_viscosity_factor[air+1]; strain[p] = 0.0; changed += 1
            end
        end
    end end end end
    return MPI.Allreduce(changed, +, MPI.COMM_WORLD)
end

"The steps after the particle control (main.cpp): resampling (2D), then the surface processes and the particles."
function surface_step!(sf::Surface, m, layers, swarm, dt)
    sf isa Surface2D && interpolate_surface!(sf, m)
    sf.processes || return 0
    surface_processes!(sf, m, dt)
    return update_particles_with_surface!(sf, m, layers, swarm)
end

# ---- 3D (not in the C++ code) ---------------------------------------------------------------------------------

"Mean height of the first quarter of the surface in x plus the sea level, as adjusted_mean_elevation in 2D."
adjusted_mean_elevation(sf::Surface3D) = sum(@view sf.z[1:length(sf.x)÷4, :]) / ((length(sf.x) ÷ 4) * length(sf.y)) + sf.sea_level

"Height of the surface at (x, y): bilinear interpolation of the regular points."
function surface_height(sf::Surface3D, x, y)
    nx, ny = size(sf.z); dxs, dys = sf.x[2] - sf.x[1], sf.y[2] - sf.y[1]
    i = clamp(trunc(Int, x / dxs), 0, nx - 2); j = clamp(trunc(Int, y / dys), 0, ny - 2)
    rx = clamp(x / dxs - i, 0.0, 1.0); ry = clamp(y / dys - j, 0.0, 1.0); z = sf.z
    return (1 - rx) * (1 - ry) * z[i+1, j+1] + rx * (1 - ry) * z[i+2, j+1] + (1 - rx) * ry * z[i+1, j+2] + rx * ry * z[i+2, j+2]
end

"""
Velocity (3 x n) at the points (px, py, pz), kept inside the domain as the particles (locate3): each point is
evaluated by the rank that owns its element, and the values are summed over the ranks.
"""
function velocity_at_points3(m, da, local_velocity, px, py, pz)
    xs, ys, zs, xm, ym, zm = Int.(LibPETSc.DMDAGetCorners(petsclib, da))
    out = zeros(3, length(px))
    PETSc.with_local_array!(local_velocity; write = false) do v
        V = local_view3(da, v, 3)
        for p in eachindex(px)
            i, j, k, rx, ry, rz, _ = locate3(m, px[p], py[p], pz[p])
            xs <= i < xs + xm && ys <= j < ys + ym && zs <= k < zs + zm || continue
            w8 = trilinear_weights(rx, ry, rz); nodes = element_nodes3(i, j, k)
            for c in 1:3; out[c, p] = sum(w8[q] * V[c, nodes[q]...] for q in 1:8); end
        end
    end
    return MPI.Allreduce(out, +, MPI.COMM_WORLD)
end

"""
Moves the surface with the velocity of the substep (semi-Lagrangian): the point that the velocity brings to the
regular point (x, y) starts at (x, y) - dt v (two fixed point iterations), and the new height is the height there
plus dt times the vertical velocity there. To first order in dt it is the Lagrangian step of the 2D surface
followed by its resampling.
"""
function advect_surface!(sf::Surface3D, m, da_advection, local_velocity, dt)
    X = vec([x for x in sf.x, _ in sf.y]); Y = vec([y for _ in sf.x, y in sf.y])
    keep(x, y) = (clamp(x, 0.0, m.Lx), clamp(y, 0.0, m.Ly))
    x0, y0 = X, Y; z0 = vec(copy(sf.z))
    for _ in 1:2
        v = velocity_at_points3(m, da_advection, local_velocity, x0, y0, z0)
        d = keep.(X .- dt .* @view(v[1, :]), Y .- dt .* @view(v[2, :]))
        x0, y0 = first.(d), last.(d); z0 = surface_height.(Ref(sf), x0, y0)
    end
    v = velocity_at_points3(m, da_advection, local_velocity, x0, y0, z0)
    vec(sf.z) .= z0 .+ dt .* @view(v[3, :])
end

# Surface in the RK4 step: in 2D the points follow the same RK4 as the particles; in 3D the semi-Lagrangian scheme
# above is kept, over the whole step, with the velocity at the middle of the step (s.velocity_advection).
advect_surface_runge_kutta_4!(sf::Surface2D, s, dt) = advect_surface!(sf, s.m, s.da_advection, s.local_velocity_advection, dt;
    runge_kutta_4 = runge_kutta_4_velocities!(s))
advect_surface_runge_kutta_4!(sf::Surface3D, s, dt) = (global_to_local(s.da_advection, s.velocity_advection,
    s.local_velocity_advection); advect_surface!(sf, s.m, s.da_advection, s.local_velocity_advection, dt))

"Diffusion of the surface (sp_d_c, m2/s): explicit, with the edges kept at the height of their neighbors, as in 2D."
function diffuse_surface3!(sf::Surface3D, dt)
    z = sf.z; nx, ny = size(z); dxs, dys = sf.x[2] - sf.x[1], sf.y[2] - sf.y[1]
    sp_dt = 1 / (4 * sf.diffusivity * (1 / dxs^2 + 1 / dys^2))
    steps = max(trunc(Int, floor(dt / sp_dt)) + 1, 10)
    rx, ry = sf.diffusivity * (dt / steps) / dxs^2, sf.diffusivity * (dt / steps) / dys^2
    z[1, :] .= z[2, :]; z[nx, :] .= z[nx-1, :]; z[:, 1] .= z[:, 2]; z[:, ny] .= z[:, ny-1]
    for _ in 1:steps
        old = copy(z)
        for j in 2:ny-1, i in 2:nx-1
            z[i, j] = old[i, j] + rx * (old[i+1, j] - 2old[i, j] + old[i-1, j]) + ry * (old[i, j+1] - 2old[i, j] + old[i, j-1])
        end
    end
end

"Particles crossed by the surface, as in 2D; the height at the particle is interpolated (bilinear)."
function update_particles_with_surface!(sf::Surface3D, m, layers, swarm; new_layer = sf.active_sediment_layer)
    air = length(layers.layer_density) - 1
    changed = 0
    with_field(swarm, PIC_COORDINATES, Float64) do coor
    with_field(swarm, "layer", PInt) do layer
    with_field(swarm, "geoq_fac", Float64) do factor
    with_field(swarm, "strain_fac", Float64) do strain
        for p in eachindex(layer)
            surface = surface_height(sf, coor[3p-2], coor[3p-1]); pz = coor[3p]
            if sf.air_to_land && layer[p] == air && pz < surface
                layer[p] = new_layer; factor[p] = layers.layer_viscosity_factor[new_layer+1]
                strain[p] = sf.sediment_strain; changed += 1
            end
            if layer[p] != air && pz > surface
                layer[p] = air; factor[p] = layers.layer_viscosity_factor[air+1]; strain[p] = 0.0; changed += 1
            end
        end
    end end end end
    return MPI.Allreduce(changed, +, MPI.COMM_WORLD)
end
