# Equação térmica 2D (DMT.cpp, DM1_2d.cpp e thermal_Ke.cpp): elementos bilineares,
# esquema theta e difusão artificial ao longo das linhas de corrente.
# Inclua depois de particle_swarm.jl. Fora do escopo: magmatismo, basal_heat, aquecimento
# por cisalhamento e as condições iniciais T_initial.

struct ElementMatrices
    NT::Matrix{Float64}; NT_x::Matrix{Float64}; NT_z::Matrix{Float64}  # NT, NT_x, NT_z [nó, ponto]
    gauss_weight::Vector{Float64}    # prodH*volume; aqui com o jacobiano 1/4, que o Mandyoc omite
    TKe::Matrix{Float64}; TMe::Matrix{Float64}; TFe::Vector{Float64}   # TKe, TMe, TFe
end

function ElementMatrices(m::ModelParameters)                        # montaKeThermal_general_2d
    dx, dz = element_size_x(m), element_size_z(m)                   # dx_const, dz_const
    g = (-sqrt(0.6), 0.0, sqrt(0.6)); w = (5 / 9, 8 / 9, 5 / 9)      # r06, r5p9, r8p9
    NT, NT_x, NT_z = zeros(4, 9), zeros(4, 9), zeros(4, 9); gw = zeros(9)
    q = 0
    for c in 1:3, a in 1:3
        q += 1; kx, kz = g[a], g[c]; gw[q] = w[a] * w[c] * dx * dz / 4
        n = 0
        for ez in (-1, 1), ex in (-1, 1)
            n += 1
            NT[n, q] = (1 + ex * kx) * (1 + ez * kz) / 4
            NT_x[n, q] = ex * (1 + ez * kz) / 2 / dx
            NT_z[n, q] = (1 + ex * kx) * ez / 2 / dz
        end
    end
    TKe = [sum(gw[q] * (NT_x[i, q] * NT_x[j, q] + NT_z[i, q] * NT_z[j, q]) for q in 1:9) for i in 1:4, j in 1:4]
    TMe = [sum(gw[q] * NT[i, q] * NT[j, q] for q in 1:9) for i in 1:4, j in 1:4]
    TFe = [sum(gw[q] * NT[i, q] for q in 1:9) / m.heat_capacity for i in 1:4]
    return ElementMatrices(NT, NT_x, NT_z, gw, TKe, TMe, TFe)
end

# montaKeThermal_simplif_2d: advecção + difusão artificial + difusão
function advection_diffusion!(TCe, e::ElementMatrices, m, u, w, kappa)
    dx, dz = element_size_x(m), element_size_z(m)
    eps_upwind(a) = a < -1 ? -1 - 1 / a : a <= 1 ? 0.0 : 1 - 1 / a   # eps
    kappa_arti = (eps_upwind(u * dx / 2kappa) * u * dx + eps_upwind(w * dz / 2kappa) * w * dz) / 2
    speed2 = u^2 + w^2                                              # mod_vv^2
    fill!(TCe, 0.0)
    @inbounds for q in 1:9, j in 1:4
        advection = e.gauss_weight[q] * (e.NT_x[j, q] * u + e.NT_z[j, q] * w)   # NT_V
        for i in 1:4
            upwind = speed2 > 0 ? kappa_arti * (e.NT_x[i, q] * u + e.NT_z[i, q] * w) / speed2 : 0.0  # NT_u
            TCe[i, j] += (e.NT[i, q] + upwind) * advection
        end
    end
    TCe .+= kappa .* e.TKe
    return TCe
end

mutable struct ThermalSolver{D, V, M, K}
    m::ModelParameters
    element::ElementMatrices
    da_Thermal::D                        # da_Thermal
    da_Veloc::D                          # da_Veloc
    thermal_matrix_lhs::M                # TA
    thermal_rhs::V                       # Tf
    thermal_ksp::K                       # T_ksp
    temperature_bc_values::V             # Temper_0
    temperature_bc_mask::V               # Temper_Cond (1 = livre, 0 = fixo)
    local_temperature::V                 # local_Temper
    local_temperature_bc_mask::V         # local_TC
    local_thermal_rhs::V                 # local_FT
    local_radiogenic_heat_field::V       # local_geoq_H
    local_thermal_diffusivity_field::V   # local_geoq_kappa
    local_velocity::V                    # local_V
    iterations::Int
end

const P = petsclib   # atalho
function local_view(da, a, dof = 1)      # array com índices globais, ghosts incluídos
    gx, gz, _, gm, gn, _ = LibPETSc.DMDAGetGhostCorners(P, da)
    gx, gz, gm, gn = Int(gx), Int(gz), Int(gm), Int(gn)
    return dof == 1 ? OffsetArray(reshape(a, gm, gn), gx:gx+gm-1, gz:gz+gn-1) :
        OffsetArray(reshape(a, dof, gm, gn), 1:dof, gx:gx+gm-1, gz:gz+gn-1)
end
"Preenche um Vec global da malha com f(i, k)."
function fill_field!(vec, da, f, dof = 1)
    si, sk, _, mi, mk, _ = LibPETSc.DMDAGetCorners(P, da)
    PETSc.with_local_array!(vec) do a
        n = 0
        for k in sk:sk+mk-1, i in si:si+mi-1
            v = f(Int(i), Int(k))
            for c in 1:dof
                a[n+=1] = dof == 1 ? v : v[c]
            end
        end
    end
    return vec
end
global_to_local(da, g, l) = LibPETSc.DMGlobalToLocal(P, da, g, LibPETSc.INSERT_VALUES, l)

function ThermalSolver(m::ModelParameters, da_Thermal, da_Veloc, temperature)   # create_thermal
    global_vec() = LibPETSc.DMCreateGlobalVector(P, da_Thermal)
    local_vec() = LibPETSc.DMCreateLocalVector(P, da_Thermal)
    LibPETSc.DMSetMatType(P, da_Thermal, "aij")
    mask = fill_field!(global_vec(), da_Thermal, function (i, k)
        fixed = (k == 0 && m.temperature_bc_bottom) || (k == m.Nz - 1 && m.temperature_bc_top)
        m.periodic_boundary || (fixed |= (i == 0 && m.temperature_bc_left) || (i == m.Nx - 1 && m.temperature_bc_right))
        fixed ? 0.0 : 1.0
    end)
    bc_values = global_vec()
    PETSc.with_local_array!(bc_values, temperature, mask; write = (true, false, false)) do b, T, c
        b .= T .* (1 .- c)
    end
    ksp = LibPETSc.KSPCreate(P, MPI.COMM_WORLD)
    LibPETSc.KSPSetOptionsPrefix(P, ksp, "thermal_")
    s = ThermalSolver(m, ElementMatrices(m), da_Thermal, da_Veloc, LibPETSc.DMCreateMatrix(P, da_Thermal),
        global_vec(), ksp, bc_values, mask, local_vec(), local_vec(), local_vec(), local_vec(), local_vec(),
        LibPETSc.DMCreateLocalVector(P, da_Veloc), 0)
    global_to_local(da_Thermal, mask, s.local_temperature_bc_mask)
    return s
end

# O elemento é do processo dono do nó inferior esquerdo (DMDAGetElementCorners_2d).
function element_ranges(m, da)
    si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, da))
    return si:min(si + mi, m.Nx - 1)-1, sk:min(sk + mk, m.Nz - 1)-1
end

function effective_diffusivity(m, kappa, T, i, k)                   # kappa_eff
    m.high_diffusivity_in_asthenosphere || return kappa
    mean = (T[i, k] + T[i+1, k] + T[i, k+1] + T[i+1, k+1]) / 4      # tt_mean
    mean > 1330.0 && (kappa = m.thermal_diffusivity * (1 + (mean - 1330.0) * 5))
    return min(kappa, 50m.thermal_diffusivity)
end

function mean_velocity(V, i, k)                                     # vv
    u = (V[1, i, k] + V[1, i+1, k] + V[1, i, k+1] + V[1, i+1, k+1]) / 4
    w = (V[2, i, k] + V[2, i+1, k] + V[2, i, k+1] + V[2, i+1, k+1]) / 4
    return u, w
end

function assemble_lhs!(s::ThermalSolver, velocity_next, dt)         # AssembleA_Thermal_2d
    m, A = s.m, s.thermal_matrix_lhs
    LibPETSc.MatZeroEntries(P, A)
    global_to_local(s.da_Veloc, velocity_next, s.local_velocity)
    gx, gz, _, gm, _, _ = Int.(LibPETSc.DMDAGetGhostCorners(P, s.da_Thermal))
    row(i, k) = (i - gx) + (k - gz) * gm                            # índice local com ghosts, base 0
    PETSc.with_local_array!(s.local_velocity, s.local_temperature, s.local_temperature_bc_mask,
        s.local_thermal_diffusivity_field; write = false) do v, t, c, kf
        V = local_view(s.da_Veloc, v, 2)
        T, C, kappa = local_view(s.da_Thermal, t), local_view(s.da_Thermal, c), local_view(s.da_Thermal, kf)
        si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, s.da_Thermal))
        r1, r2 = PInt[0], PInt[0]; one_value = [0.0]
        put(r, c, x) = (r1[1] = r; r2[1] = c; one_value[1] = x;
            LibPETSc.MatSetValuesLocal(P, A, PInt(1), r1, PInt(1), r2, one_value, LibPETSc.ADD_VALUES))
        for k in sk:sk+mk-1, i in si:si+mi-1
            if m.periodic_boundary && i == m.Nx - 1                 # T[Nx-1] = T[0]
                put(row(i, k), row(i, k), 1.0); put(row(i, k), row(m.Nx, k), -1.0)
            else
                put(row(i, k), row(i, k), 1.0 - C[i, k])
            end
        end
        TCe = zeros(4, 4); values = zeros(16); index = zeros(PInt, 4)
        for ek in element_ranges(m, s.da_Thermal)[2], ei in element_ranges(m, s.da_Thermal)[1]
            u, w = mean_velocity(V, ei, ek)
            advection_diffusion!(TCe, s.element, m, u, w, effective_diffusivity(m, kappa[ei, ek], T, ei, ek))
            right = m.periodic_boundary && ei == m.Nx - 2 ? m.Nx : ei + 1
            nodes = ((ei, ek), (right, ek), (ei, ek + 1), (right, ek + 1))
            for a in 1:4
                index[a] = row(nodes[a]...)
                free = C[nodes[a]...] != 0
                for b in 1:4                                        # PETSc lê por linha
                    values[4(a-1)+b] = free ? s.element.TMe[a, b] + m.thermal_time_theta * dt * TCe[a, b] : 0.0
                end
            end
            LibPETSc.MatSetValuesLocal(P, A, PInt(4), index, PInt(4), index, values, LibPETSc.ADD_VALUES)
        end
    end
    LibPETSc.MatAssemblyBegin(P, A, LibPETSc.MAT_FINAL_ASSEMBLY)
    LibPETSc.MatAssemblyEnd(P, A, LibPETSc.MAT_FINAL_ASSEMBLY)
end

function assemble_rhs!(s::ThermalSolver, velocity, dt)              # AssembleF_Thermal_2d
    m, e = s.m, s.element
    one_minus_theta = 1 - m.thermal_time_theta                      # comp_alpha_thermal
    global_to_local(s.da_Veloc, velocity, s.local_velocity)
    LibPETSc.VecZeroEntries(P, s.local_thermal_rhs)
    PETSc.with_local_array!(s.local_thermal_rhs, s.local_velocity, s.local_temperature, s.local_temperature_bc_mask,
        s.local_thermal_diffusivity_field, s.local_radiogenic_heat_field;
        write = (true, false, false, false, false, false)) do f, v, t, c, kf, h
        V = local_view(s.da_Veloc, v, 2); lv(a) = local_view(s.da_Thermal, a)
        F, T, C, kappa, H = lv(f), lv(t), lv(c), lv(kf), lv(h)
        TCe = zeros(4, 4)
        for ek in element_ranges(m, s.da_Thermal)[2], ei in element_ranges(m, s.da_Thermal)[1]
            u, w = mean_velocity(V, ei, ek)
            advection_diffusion!(TCe, e, m, u, w, effective_diffusivity(m, kappa[ei, ek], T, ei, ek))
            old = ((ei, ek), (ei + 1, ek), (ei, ek + 1), (ei + 1, ek + 1))
            right = m.periodic_boundary && ei == m.Nx - 2 ? m.Nx : ei + 1
            new = ((ei, ek), (right, ek), (ei, ek + 1), (right, ek + 1))
            heat = 0.0                                              # H_efetivo
            m.use_radiogenic_heating && (heat += sum(H[n...] for n in old) / 4)
            if m.use_adiabatic_heating
                kelvin = sum(T[n...] for n in old) / 4 + 273.0      # T_mean
                heat -= kelvin * m.thermal_expansion_coefficient * m.gravity_acceleration * w
            end
            for a in 1:4
                C[new[a]...] == 0 && continue
                x = dt * e.TFe[a] * heat
                for b in 1:4
                    x += T[old[b]...] * (e.TMe[a, b] - one_minus_theta * dt * TCe[a, b])
                end
                F[new[a]...] += x
            end
        end
        si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, s.da_Thermal))
        for k in sk:sk+mk-1, i in si:si+mi-1
            if m.periodic_boundary && i == m.Nx - 1
                F[i, k] = 0.0
            elseif C[i, k] == 0
                F[i, k] = T[i, k]
            end
        end
    end
    LibPETSc.VecZeroEntries(P, s.thermal_rhs)
    LibPETSc.DMLocalToGlobal(P, s.da_Thermal, s.local_thermal_rhs, LibPETSc.ADD_VALUES, s.thermal_rhs)
end

"Um passo: temperature (Temper) é atualizada no lugar. Campos por nó/elemento como no Mandyoc."
function thermal_step!(s::ThermalSolver, temperature, velocity, velocity_next,
    radiogenic_heat_field, thermal_diffusivity_field, time_step_seconds)  # build_thermal + solve_thermal
    global_to_local(s.da_Thermal, temperature, s.local_temperature)
    global_to_local(s.da_Thermal, radiogenic_heat_field, s.local_radiogenic_heat_field)
    global_to_local(s.da_Thermal, thermal_diffusivity_field, s.local_thermal_diffusivity_field)
    assemble_lhs!(s, velocity_next, time_step_seconds)
    assemble_rhs!(s, velocity, time_step_seconds)
    A, ksp = s.thermal_matrix_lhs, s.thermal_ksp
    LibPETSc.KSPSetOperators(P, ksp, A, A)
    LibPETSc.KSPSetFromOptions(P, ksp)
    LibPETSc.KSPSetTolerances(P, ksp, s.m.rtol, 1e-50, 1e5, PInt(10000))
    LibPETSc.KSPSolve(P, ksp, s.thermal_rhs, temperature)
    s.iterations = LibPETSc.KSPGetIterationNumber(P, ksp)
    PETSc.with_local_array!(temperature, s.temperature_bc_mask, s.temperature_bc_values; write = (true, false, false)) do T, c, b
        T .= T .* c .+ b
    end
    Int(LibPETSc.KSPGetConvergedReason(P, ksp)) > 0 || error("KSP térmico não convergiu")
    return temperature
end
