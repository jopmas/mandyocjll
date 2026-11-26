# Rheology (calc_visc.cpp, moveSwarm, Swarm2Mesh_2d, veloc_total, calc_pressure_2d, shift_pressure_2d).
# Include after stokes.jl. Dimensional only; no magmatism.

"Properties of each layer (interfaces.txt), from the deepest to the shallowest."
Base.@kwdef struct LayerProperties
    layer_viscosity_factor::Vector{Float64}        # inter_geoq (C)
    layer_density::Vector{Float64}                 # inter_rho
    layer_radiogenic_heat::Vector{Float64}         # inter_H
    layer_pre_exponential_factor::Vector{Float64}  # inter_A
    layer_power_law_exponent::Vector{Float64}      # inter_n
    layer_activation_energy::Vector{Float64}       # inter_Q
    layer_activation_volume::Vector{Float64}       # inter_V
    layer_cohesion_min::Vector{Float64} = fill(4.0e6, length(layer_density))           # cohesion_min
    layer_cohesion_max::Vector{Float64} = fill(20.0e6, length(layer_density))          # cohesion_max
    layer_friction_angle_min::Vector{Float64} = fill(0.034906, length(layer_density))  # friction_angle_min (rad)
    layer_friction_angle_max::Vector{Float64} = fill(0.261799, length(layer_density))  # friction_angle_max (rad)
    layer_weakening_seed::Vector{Float64} = fill(-1.0, length(layer_density))          # weakening_seed
    layer_thermal_conductivity::Vector{Float64} = fill(1.0e-6, length(layer_density))  # conductivity (it is a diffusivity, m2/s)
    interface_depths::Matrix{Float64} = zeros(0, 0)                                    # interfaces [node in x, interface]
end

"Viscosity range in use (it changes during the initial ramp of Gerya)."
mutable struct ViscosityRange
    computed_min_viscosity::Float64                # visc_MIN_comp
    computed_max_viscosity::Float64                # visc_MAX_comp
end
ViscosityRange(m::ModelParameters) = ViscosityRange(m.min_viscosity, m.max_viscosity)

strain_softening(m, strain, f1, f2) = strain < m.weakening_min ? f1 : strain > m.weakening_max ? f2 :
    f1 - (f1 - f2) * (strain - m.weakening_min) / (m.weakening_max - m.weakening_min)

"Viscosity of a particle (calc_visco_ponto). T in oC, P in Pa, z in m, strain_rate in 1/s."
function particle_viscosity(m::ModelParameters, range::ViscosityRange, l::LayerProperties, layer, T, P, z, strain_rate, strain)
    vf, A, n, Q, V = l.layer_viscosity_factor[layer], l.layer_pre_exponential_factor[layer],
        l.layer_power_law_exponent[layer], l.layer_activation_energy[layer], l.layer_activation_volume[layer]
    ref, rate = m.reference_viscosity, max(strain_rate, 1.0e-36)
    depth = m.use_pressure_in_rheology ? 0.0 : max(-(z + m.air_layer_thickness), 0.0)
    P = m.constant_rheology_pressure >= 0 ? m.constant_rheology_pressure : max(P, 0.0)
    Pe = m.use_pressure_in_rheology ? P : 10.0 * 3300.0 * depth
    TK, Tb = T + 273.0, m.temperature_difference + 273.0
    R = 8.3144
    model = m.rheology_model
    creep() = A^(-1 / n) * rate^((1 - n) / n) * exp((Q + V * Pe) / (n * R * TK))
    eta = if model == 0
        ref
    elseif model == 1
        r = 20.0; q = 225.0 / log(r) - 0.25log(r); g = 15.0 / log(r) - 0.5
        return vf * ref * exp(q / (T / m.temperature_difference + g) - q / (0.5 + g))  # as in the C++ code: no limits
    elseif model in (2, 4)
        ref * exp(120000.0 * (1 / TK - 1 / Tb) / 8.31)
    elseif model == 3
        ref * exp(120000.0 * (1 / TK - 1 / 1573.0) / 8.31)
    elseif model == 5
        ref * 1.0e7 * exp(-TK * 240000.0 / (8.31 * Tb * Tb))
    elseif model == 6
        ref * 1.0e7 * exp(-TK * 240000.0 / (R * Tb * Tb))
    elseif model == 7
        ref * A * exp(-(Q + V * Pe) / (R * TK))
    elseif model == 8
        ref * exp(-Q * TK + V * (-z))
    elseif model == 9 || model == 19
        if !m.use_nonlinear_method                  # the C++ code uses the minimum viscosity; harmless if it is the only one
            m.max_viscosity <= 1.01m.min_viscosity || error("rheology_model $model requires non_linear_method = on (the C++ code would use the minimum viscosity)")
            return m.min_viscosity
        end
        x = model == 9 ? creep() : A^(-1 / n) * rate^((1 - n) / n) * exp((Q + V * P) / (n * R * TK))
        # 19: diffusion creep (dry olivine) limits the dislocation creep when n = 3
        model == 19 && 2.99 < n < 3.01 ? min(x, 5656854249492380.0 * exp((240000.0 + 5.0e-6 * P) / (R * TK)) / 66250.0) : x
    elseif model == 10
        ref * exp(-(6.907755279 * T / 1000.0))
    elseif model == 70
        m.min_viscosity                             # replaced in the plastic block, as in the C++ code
    else
        error("rheology_model $model is unknown")
    end
    m.use_viscosity_factor && (eta *= vf)
    if m.use_nonlinear_method && m.use_plasticity
        c0 = strain_softening(m, strain, l.layer_cohesion_max[layer], l.layer_cohesion_min[layer])
        mu = strain_softening(m, strain, l.layer_friction_angle_max[layer], l.layer_friction_angle_min[layer])
        eta = min(eta, (c0 * cos(mu) + sin(mu) * Pe) / (2 * rate))       # Drucker-Prager
        model == 70 && (eta = 2ref * rate > 1 ? 1 / (2 * rate) : ref)
    end
    eta = clamp(eta, m.min_viscosity, m.max_viscosity)
    # Gerya (p. 215): maps the physical range onto the range in use
    f1 = log(range.computed_max_viscosity / range.computed_min_viscosity)
    return range.computed_min_viscosity * exp(f1 * log(eta / m.min_viscosity) / log(m.max_viscosity / m.min_viscosity))
end

# ---- moveSwarm: T, P, strain, viscosity and (Euler) position ----------------
"""
Updates the strain rate, the accumulated strain and the viscosity of the particles, and moves them
(Euler) if `move`. `local_velocity` is from `da_advection`. Without `move`, dt only accumulates strain.
"""
function update_particles!(m, range, layers, swarm, da_advection, local_velocity, da_Thermal, local_temperature,
    local_pressure_thermal_mesh, dt; move = false)
    dx, dz = element_size_x(m), element_size_z(m); e = m.boundary_offset
    PETSc.with_local_array!(local_velocity, local_temperature, local_pressure_thermal_mesh; write = false) do v, t, pr
        V, T, Pn = local_view(da_advection, v, 2), local_view(da_Thermal, t), local_view(da_Thermal, pr)
        with_field(swarm, PIC_COORDINATES, Float64) do coor
        with_field(swarm, "layer", PInt) do layer
        with_field(swarm, "geoq_fac", Float64) do viscosity
        with_field(swarm, "strain_fac", Float64) do strain
        with_field(swarm, "strain_rate_fac", Float64) do strain_rate
            update_loop!(m, range, layers, V, T, Pn, coor, layer, viscosity, strain, strain_rate, dt, move)
        end end end end end
    end
    move && LibPETSc.DMSwarmMigrate(petsclib, swarm, LibPETSc.PETSC_TRUE)
end

function update_loop!(m, range, layers, V, T, Pn, coor, layer, viscosity, strain, strain_rate, dt, move)
    dx, dz = element_size_x(m), element_size_z(m); e = m.boundary_offset
    image(i, A, d) = m.periodic_boundary ? (i < firstindex(A, d) ? i + m.Nx : i + 1 > lastindex(A, d) ? i - m.Nx : i) : i
    for p in eachindex(layer)
        x = clamp(coor[2p-1], e, m.Lx - e); z = clamp(coor[2p], e - m.domain_depth, -e)
        i = clamp(trunc(Int, x / dx), 0, m.Nx - 2); k = clamp(trunc(Int, (z + m.domain_depth) / dz), 0, m.Nz - 2)
        rx = x / dx - i; rz = (z + m.domain_depth) / dz - k
        j = image(i, V, 2); q = image(i, T, 1)
        Tp = linear_interpolation(rx, rz, T[q, k], T[q+1, k], T[q, k+1], T[q+1, k+1])
        Pp = linear_interpolation(rx, rz, Pn[q, k], Pn[q+1, k], Pn[q, k+1], Pn[q+1, k+1])
        # strain rate at the point
        exx = ezz = g = 0.0
        for (di, dk) in ((0, 0), (1, 0), (0, 1), (1, 1))
            sx, sz = 2di - 1, 2dk - 1
            Nx = sx * (1 + sz * (2rz - 1)) / 2 / dx; Nz = (1 + sx * (2rx - 1)) * sz / 2 / dz
            u, w = V[1, j+di, k+dk], V[2, j+di, k+dk]
            exx += Nx * u; ezz += Nz * w; g += Nx * w + Nz * u
        end
        mean = (exx + ezz) / 3
        e2 = sqrt(((exx - mean)^2 + mean^2 + (ezz - mean)^2 + g^2 / 2) / 2)
        if move
            u = linear_interpolation(rx, rz, V[1, j, k], V[1, j+1, k], V[1, j, k+1], V[1, j+1, k+1])
            w = linear_interpolation(rx, rz, V[2, j, k], V[2, j+1, k], V[2, j, k+1], V[2, j+1, k+1])
            xn = coor[2p-1] + dt * u
            coor[2p-1] = m.periodic_boundary ? mod(xn, m.Lx) : xn
            coor[2p] += dt * w
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
end

# ---- Swarm2Mesh_2d -------------------------------------------------------------
"Mesh fields: nodal density and heat; viscosity and diffusivity of each element."
struct MeshFields{V}
    density_field::V                     # geoq_rho
    radiogenic_heat_field::V             # geoq_H
    strain_field::V                      # geoq_strain
    strain_rate_field::V                 # geoq_strain_rate
    viscosity_factor::V                  # geoq
    thermal_diffusivity_field::V         # geoq_kappa
    density_variation::V                 # dRho
end
MeshFields(da) = MeshFields((v = LibPETSc.DMCreateGlobalVector(P, da); LibPETSc.VecZeroEntries(P, v); v) for _ in 1:7)
MeshFields(g::Base.Generator) = MeshFields(collect(g)...)

function particles_to_mesh!(m, layers, swarm, da_Thermal, f::MeshFields)
    lv() = (v = LibPETSc.DMCreateLocalVector(P, da_Thermal); LibPETSc.VecZeroEntries(P, v); v)
    L = ntuple(_ -> lv(), 8)             # rho, H, strain rate, weight | viscosity, kappa, count | strain
    dx, dz = element_size_x(m), element_size_z(m)
    harmonic = m.use_viscosity_harmonic_mean
    PETSc.with_local_array!(L...) do arrays...
        R, H, S, W, Vi, Ka, N, E = map(a -> local_view(da_Thermal, a), arrays)
        with_field(swarm, PIC_COORDINATES, Float64) do coor
        with_field(swarm, "layer", PInt) do layer
        with_field(swarm, "geoq_fac", Float64) do visc
        with_field(swarm, "strain_rate_fac", Float64) do rate
        with_field(swarm, "strain_fac", Float64) do strain
            for p in eachindex(layer)
                x, z = coor[2p-1], coor[2p] + m.domain_depth
                i = clamp(trunc(Int, x / dx), 0, m.Nx - 2); k = clamp(trunc(Int, z / dz), 0, m.Nz - 2)
                rx = clamp(x / dx - i, 0, 1); rz = clamp(z / dz - k, 0, 1); n = layer[p] + 1
                h = layers.layer_radiogenic_heat[n]
                m.use_shear_heating && (h += m.shear_heating_scale_factor * 4visc[p] * rate[p]^2 / layers.layer_density[n])
                for (ii, kk, w) in ((i, k, (1 - rx) * (1 - rz)), (i + 1, k, rx * (1 - rz)), (i, k + 1, (1 - rx) * rz), (i + 1, k + 1, rx * rz))
                    R[ii, kk] += w * layers.layer_density[n]; H[ii, kk] += w * h; S[ii, kk] += w * rate[p]; E[ii, kk] += w * strain[p]; W[ii, kk] += w
                end
                kap = layers.layer_thermal_conductivity[n]
                Vi[i, k] += harmonic ? 1 / visc[p] : visc[p]; Ka[i, k] += harmonic ? 1 / kap : kap; N[i, k] += 1
            end
        end end end end end
    end
    g = LibPETSc.DMCreateGlobalVector(P, da_Thermal); d = LibPETSc.DMCreateGlobalVector(P, da_Thermal)
    for (num, den, out, inv) in ((L[1], L[4], f.density_field, false), (L[2], L[4], f.radiogenic_heat_field, false),
        (L[3], L[4], f.strain_rate_field, false), (L[8], L[4], f.strain_field, false), (L[5], L[7], f.viscosity_factor, harmonic), (L[6], L[7], f.thermal_diffusivity_field, harmonic))
        for (a, b) in ((num, g), (den, d))
            LibPETSc.VecZeroEntries(P, b); LibPETSc.DMLocalToGlobal(P, da_Thermal, a, LibPETSc.ADD_VALUES, b)
        end
        m.periodic_boundary && (sum_periodic!(m, da_Thermal, g, out === f.viscosity_factor || out === f.thermal_diffusivity_field);
            sum_periodic!(m, da_Thermal, d, out === f.viscosity_factor || out === f.thermal_diffusivity_field))
        PETSc.with_local_array!(out, g, d) do o, a, b
            for n in eachindex(o)        # no particles: keep the previous value
                b[n] > 0 && (o[n] = inv ? b[n] / a[n] : a[n] / b[n])
            end
        end
    end
    foreach(v -> LibPETSc.VecDestroy(P, v), (L..., g, d))
end

"Air: zero temperature where the density is low; viscous layer (sticky_blanket_air) at the interface with the air."
function air_corrections!(m, da_Thermal, f::MeshFields, temperature)
    PETSc.with_local_array!(temperature, f.density_field; write = (true, false)) do T, rho
        for n in eachindex(T); rho[n] < m.air_threshold_density && (T[n] = 0.0); end
    end
    m.sticky_blanket_air || return
    l = LibPETSc.DMCreateLocalVector(P, da_Thermal)
    global_to_local(da_Thermal, f.density_field, l)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da_Thermal))
    xs, zs = element_ranges(m, da_Thermal)
    PETSc.with_local_array!(f.viscosity_factor, l; write = (true, false)) do va, ra
        V = OffsetArray(reshape(va, mi, mk), si:si+mi-1, sk:sk+mk-1); R = local_view(da_Thermal, ra)
        blanket(i, k) = i >= 0 && 5.0 < (R[i, k] + R[i+1, k] + R[i, k+1] + R[i+1, k+1]) / 4 < 2700.0 * 0.8
        for k in zs, i in xs             # this element or the one on the left is at the interface
            (blanket(i, k) || blanket(i - 1, k)) && (V[i, k] = max(V[i, k], 50m.min_viscosity))
        end
    end
    LibPETSc.VecDestroy(P, l)
end

# mean_value_periodic_boundary_2d with esc = 1: nodes 0 and Nx-1 are the same point. Element
# fields are stored at the lower left node and are not added.
function sum_periodic!(m, da, g, element_field)
    element_field && return
    l = LibPETSc.DMCreateLocalVector(P, da)
    global_to_local(da, g, l)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da))
    PETSc.with_local_array!(g, l; write = (true, false)) do ga, la
        G = OffsetArray(reshape(ga, mi, mk), si:si+mi-1, sk:sk+mk-1); A = local_view(da, la)
        for k in sk:sk+mk-1
            si == 0 && (G[0, k] = A[0, k] + A[-1, k])
            si + mi == m.Nx && (G[m.Nx-1, k] = A[m.Nx-1, k] + A[m.Nx, k])
        end
    end
    LibPETSc.VecDestroy(P, l)
end

# ---- pressure -------------------------------------------------------------------
"Element pressure -> nodes, with zero mean pressure at the top (shift_pressure_2d)."
function pressure_to_nodes!(m, da_Veloc, da_Thermal, pressure, local_pressure, pressure_thermal_mesh)
    global_to_local(da_Veloc, pressure, local_pressure)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da_Thermal))
    top = 0.0; count = 0
    PETSc.with_local_array!(pressure_thermal_mesh, local_pressure; write = (true, false)) do a, pl
        N = OffsetArray(reshape(a, mi, mk), si:si+mi-1, sk:sk+mk-1); E = local_view(da_Veloc, pl, 2)
        for k in sk:sk+mk-1, i in si:si+mi-1
            s = 0.0; c = 0
            for (ii, kk) in ((i, k), (i - 1, k), (i, k - 1), (i - 1, k - 1))
                0 <= ii < m.Nx - 1 && 0 <= kk < m.Nz - 1 && (s += E[1, ii, kk]; c += 1)
            end
            N[i, k] = s / c
            k == m.Nz - 1 && (top += N[i, k]; count += 1)
        end
    end
    mean = MPI.Allreduce(top, +, MPI.COMM_WORLD) / MPI.Allreduce(count, +, MPI.COMM_WORLD)
    PETSc.with_local_array!(a -> a .-= mean, pressure_thermal_mesh)
end

"Initial lithostatic pressure of each element, from the interfaces (calc_pressure_2d)."
function lithostatic_pressure!(m, layers, da_Veloc, pressure)
    d = layers.interface_depths; n = size(d, 2); rho, g = layers.layer_density, m.gravity_acceleration
    fill_field!(pressure, da_Veloc, function (i, k)
        (i >= m.Nx - 1 || k >= m.Nz - 1) && return (0.0, 0.0)
        z = -(m.Nz - 1 - k) * m.domain_depth / (m.Nz - 1)
        I(q) = (d[i+1, q] + d[i+2, q]) / 2
        z > I(n) && return (rho[n+1] * g * -z, 0.0)
        p = rho[n+1] * g * -I(n)
        for q in n:-1:2
            z > I(q - 1) && return (p + rho[q] * g * (I(q) - z), 0.0)
            p += rho[q] * g * (I(q) - I(q - 1))
        end
        (p + rho[1] * g * (I(1) - z), 0.0)
    end, 2)
end

# ---- coupled model -------------------------------------------------------------
mutable struct Model{D, V, S, F, TS, SS}
    m::ModelParameters
    layers::LayerProperties
    range::ViscosityRange
    da_Veloc::D                          # da_Veloc
    da_Thermal::D                        # da_Thermal
    da_advection::D
    particle_swarm::S                    # dms
    cell_dm::S                           # dmcell
    fields::F
    temperature::V                       # Temper
    velocity::V                          # Veloc
    velocity_next::V                     # Veloc_fut
    velocity_advection::V                # Veloc_weight
    pressure::V                          # Pressure
    pressure_thermal_mesh::V             # Pressure_aux
    local_velocity_advection::V          # local_V (da_advection)
    local_temperature::V                 # local_Temper
    local_pressure_thermal_mesh::V       # local_P_aux
    stokes::SS
    thermal::TS
    first_velocity_computed::Bool        # Verif_first_veloc
    time_step_seconds::Float64           # dt_calor_sec
    time::Float64                        # tempo (s)
    nonlinear_iterations::Int
end

function Model(m, layers; temperature, velocity, velocity_mask = nothing, location = :native, initial_strain = (x, z, layer) -> 0.0)
    dv, dT = create_velocity_mesh(m), create_thermal_mesh(m); da = create_advection_mesh(m, dv)
    n = size(layers.interface_depths, 2)
    function layer_at(x, z, id)
        i = clamp(trunc(Int, x / element_size_x(m)), 0, m.Nx - 2); rx = x / element_size_x(m) - i
        for q in 1:n
            z < layers.interface_depths[i+1, q] * (1 - rx) + layers.interface_depths[i+2, q] * rx && return (q - 1, layers.layer_viscosity_factor[q])
        end
        return (n, layers.layer_viscosity_factor[n+1])
    end
    swarm, cell = create_particle_swarm(m, dv; location, init = layer_at)
    with_field(swarm, "strain_fac", Float64) do s
        with_field(swarm, "layer", PInt) do l
        with_field(swarm, PIC_COORDINATES, Float64) do c
            for p in eachindex(s)
                seed = layers.layer_weakening_seed[l[p]+1]
                s[p] = seed >= 0 ? seed : initial_strain(c[2p-1], c[2p], Int(l[p]))
            end
        end end
    end
    zero_vec(d) = (v = LibPETSc.DMCreateGlobalVector(P, d); LibPETSc.VecZeroEntries(P, v); v)
    T = fill_field!(zero_vec(dT), dT, temperature); V = fill_field!(zero_vec(dv), dv, velocity, 2); Vn = zero_vec(dv)
    LibPETSc.VecCopy(P, V, Vn)
    pr = zero_vec(dv); pn = zero_vec(dT)
    stokes = StokesSolver(m, dv, dT, V; mask = velocity_mask)
    if n > 0
        lithostatic_pressure!(m, layers, dv, pr)
        pressure_to_nodes!(m, dv, dT, pr, stokes.local_pressure, pn)
    end
    return Model(m, layers, ViscosityRange(m), dv, dT, da, swarm, cell, MeshFields(dT), T, V, Vn, zero_vec(dv), pr, pn,
        LibPETSc.DMCreateLocalVector(P, da), LibPETSc.DMCreateLocalVector(P, dT), LibPETSc.DMCreateLocalVector(P, dT),
        stokes, ThermalSolver(m, dT, dv, T), false, 0.0, 0.0, 0)
end

function particles!(s::Model, dt; move = false)
    global_to_local(s.da_advection, s.velocity_advection, s.local_velocity_advection)
    global_to_local(s.da_Thermal, s.temperature, s.local_temperature)
    global_to_local(s.da_Thermal, s.pressure_thermal_mesh, s.local_pressure_thermal_mesh)
    update_particles!(s.m, s.range, s.layers, s.particle_swarm, s.da_advection, s.local_velocity_advection, s.da_Thermal,
        s.local_temperature, s.local_pressure_thermal_mesh, dt; move)
end

"One assembly and solve (build_veloc + solve_veloc)."
function linear_solve!(s::Model)
    s.first_velocity_computed && LibPETSc.VecCopy(P, s.velocity_next, s.velocity_advection)
    s.first_velocity_computed = true
    particles!(s, 0.0)
    particles_to_mesh!(s.m, s.layers, s.particle_swarm, s.da_Thermal, s.fields)
    air_corrections!(s.m, s.da_Thermal, s.fields, s.temperature)
    calc_density_variation!(s.m, s.fields.density_variation, s.temperature, s.fields.density_field)
    stokes_solve!(s.stokes, s.velocity_next, s.pressure, s.velocity, s.fields.viscosity_factor, s.fields.density_variation, s.time_step_seconds)
    pressure_to_nodes!(s.m, s.da_Veloc, s.da_Thermal, s.pressure, s.stokes.local_pressure, s.pressure_thermal_mesh)
end

# 1 - correlation between two successive velocities (Xi)
function velocity_change(a, b)
    ma, mb, n = LibPETSc.VecSum(P, a), LibPETSc.VecSum(P, b), LibPETSc.VecGetSize(P, a)
    ma /= n; mb /= n
    s = [0.0, 0.0, 0.0]
    PETSc.with_local_array!(a, b; write = false) do x, y
        for i in eachindex(x)
            dx, dy = x[i] - ma, y[i] - mb
            s[1] += dx^2; s[2] += dy^2; s[3] += dx * dy
        end
    end
    s = MPI.Allreduce(s, +, MPI.COMM_WORLD)
    return s[1] == 0 || s[2] == 0 ? 1.0 : 1 - s[3] / sqrt(s[1] * s[2])  # in the C++ code a zero velocity gives NaN and ends the loop
end

"Velocity and pressure, with nonlinear iterations (veloc_total)."
function solve_velocity_total!(s::Model; max_iterations = 700)
    LibPETSc.VecCopy(P, s.velocity_next, s.velocity)
    s.nonlinear_iterations = 0
    if !s.m.use_nonlinear_method
        linear_solve!(s); s.nonlinear_iterations = 1
        return s
    end
    previous = LibPETSc.DMCreateGlobalVector(P, s.da_Veloc)         # Veloc_step1
    xi = Inf
    while s.nonlinear_iterations < max_iterations && xi > s.m.nonlinear_residual_tolerance
        LibPETSc.VecCopy(P, s.velocity_next, previous)
        linear_solve!(s)
        xi = velocity_change(previous, s.velocity_next)
        s.nonlinear_iterations += 1
    end
    LibPETSc.VecDestroy(P, previous)
    xi > s.m.nonlinear_residual_tolerance && MPI.Comm_rank(MPI.COMM_WORLD) == 0 && @warn "the nonlinear iterations did not converge" xi
    return s
end

"Velocity at the beginning or after a change of the boundary velocity, with the viscosity ramp of Gerya (p. 215)."
function ramp_velocity!(s::Model)
    m, r = s.m, s.range
    if m.initial_dynamic_range && m.max_viscosity > m.min_viscosity
        mean = 10.0^(log10(m.min_viscosity) + log10(m.max_viscosity / m.min_viscosity) / 2)
        r.computed_min_viscosity = r.computed_max_viscosity = mean
        solve_velocity_total!(s)
        n = 0
        while r.computed_min_viscosity != m.min_viscosity && r.computed_max_viscosity != m.max_viscosity
            r.computed_min_viscosity = max(mean * 10.0^-n, m.min_viscosity)
            r.computed_max_viscosity = min(mean * 10.0^n, m.max_viscosity)
            solve_velocity_total!(s)
            n += 1
        end
    else
        solve_velocity_total!(s)
    end
end

"Initial solution."
function initial_velocity!(s::Model)
    ramp_velocity!(s)
    LibPETSc.VecCopy(P, s.velocity_next, s.velocity)
    s.time_step_seconds = time_step(s)
    return s
end

const SECONDS_PER_YEAR = 365.0 * 24 * 3600   # seg_per_ano
"Time step (Calc_dt_calor): 10% of the time to cross an element, limited by max_time_step."
function time_step(s::Model)
    u = w = 0.0
    PETSc.with_local_array!(s.velocity_next; write = false) do a
        u = maximum(abs, @view a[1:2:end]; init = 0.0); w = maximum(abs, @view a[2:2:end]; init = 0.0)
    end
    u, w = MPI.Allreduce([u, w], max, MPI.COMM_WORLD)
    h = u >= w ? element_size_x(s.m) : element_size_z(s.m)
    return min(0.1h / (max(u, w) * s.m.time_step_subdivision), s.m.max_time_step * SECONDS_PER_YEAR)
end

"Multiplies the boundary velocity (and the current one) by `scale`, as scale_bcv.txt does (rescaleVeloc)."
function rescale_velocity!(s::Model, scale)
    LibPETSc.VecScale(P, s.velocity_next, scale)
    LibPETSc.VecScale(P, s.stokes.velocity_bc_values, scale)
    s.m.initial_dynamic_range && s.m.max_viscosity > s.m.min_viscosity && ramp_velocity!(s)
end

"One time step: temperature, velocity and advection (4 Euler substeps, as in the C++ code)."
function step!(s::Model; substeps = 4, surface = nothing)
    dt = s.time_step_seconds
    thermal_step!(s.thermal, s.temperature, s.velocity, s.velocity_next, s.fields.radiogenic_heat_field, s.fields.thermal_diffusivity_field, dt)
    solve_velocity_total!(s)
    for c in 0:substeps-1
        f = (0.5 + c) / substeps
        PETSc.with_local_array!(s.velocity_advection, s.velocity, s.velocity_next; write = (true, false, false)) do a, v0, v1
            a .= f .* v1 .+ (1 - f) .* v0
        end
        particles!(s, dt / substeps; move = true)
        surface === nothing || advect_surface!(surface, s.m, s.da_advection, s.local_velocity_advection, dt / substeps)
    end
    s.time += dt
    s.time_step_seconds = time_step(s)
    return s
end
