# Stokes 2D (DMv.cpp, DM1_v.cpp, calc_drho.cpp): elementos Q1-P0 e gradiente conjugado de
# Uzawa na pressão. Inclua depois de thermal.jl. A pressão do elemento fica no componente 1
# (u) do nó inferior esquerdo da da_Veloc, como no Mandyoc.
# As matrizes têm o jacobiano 1/4 que o Mandyoc omite; a solução é a mesma.

struct StokesElement
    Ke::Matrix{Float64}      # Ke_veloc para viscosidade 1
    VCe::Vector{Float64}     # VCe: gradiente da pressão
    VfMe::Matrix{Float64}    # VfMe: força de volume
end
function StokesElement(m, e::ElementMatrices)  # montaKeVeloc_general_2d, montaCeVeloc2d, montafeVeloc2d
    Ke, VCe, VfMe = zeros(8, 8), zeros(8), zeros(8, 4)
    for q in 1:9
        SN = zeros(3, 8)
        for n in 1:4
            x, z = e.NT_x[n, q], e.NT_z[n, q]
            SN[1, 2n-1] = x; SN[2, 2n] = z; SN[3, 2n-1] = z; SN[3, 2n] = x
        end
        for i in 1:8, j in 1:8
            Ke[i, j] += e.gauss_weight[q] * (2SN[1, i] * SN[1, j] + 2SN[2, i] * SN[2, j] + SN[3, i] * SN[3, j])
        end
        for n in 1:4
            VCe[2n-1] -= e.gauss_weight[q] * e.NT_x[n, q]; VCe[2n] -= e.gauss_weight[q] * e.NT_z[n, q]
            for j in 1:4; VfMe[2n, j] += e.gauss_weight[q] * e.NT[n, q] * e.NT[j, q]; end
        end
    end
    return StokesElement(Ke, VCe, VfMe)
end

mutable struct StokesSolver{D, V, M, K}
    m::ModelParameters
    element::StokesElement
    da_Veloc::D                          # da_Veloc
    da_Thermal::D                        # da_Thermal
    velocity_stiffness_matrix::M         # VA
    pressure_gradient_matrix::M          # VG
    velocity_ksp::K                      # V_ksp
    velocity_rhs_pressure::V             # Vf_P
    velocity_bc_values::V                # Veloc_0
    velocity_bc_mask::V                  # Veloc_Cond (1 = livre, 0 = fixo)
    uzawa_preconditioner::V              # Precon
    rk_vec::V; rk_vec2::V; sk_vec::V; gs_vec::V; uk_vec::V; zk_vec::V; zk_vec2::V
    local_velocity::V                    # local_V
    local_velocity_bc_mask::V            # local_VC
    local_pressure_rhs::V                # local_FP
    local_pressure::V                    # local_P
    local_uzawa_preconditioner::V        # local_Precon
    local_viscosity_factor::V            # local_geoq
    local_density_variation::V           # local_dRho
    local_fssa_density_variation::V      # novo: o C++ usa o dRho da montagem anterior no FSSA
    uzawa_iterations::Int
    uzawa_residual::Float64              # denok
end

"mask = (i, k) -> (u, w) troca as condições de contorno dos parâmetros."
function StokesSolver(m, da_Veloc, da_Thermal, velocity; mask = nothing)   # create_veloc
    gv() = LibPETSc.DMCreateGlobalVector(P, da_Veloc); lv() = LibPETSc.DMCreateLocalVector(P, da_Veloc)
    lt() = LibPETSc.DMCreateLocalVector(P, da_Thermal)
    LibPETSc.DMSetMatType(P, da_Veloc, "aij")
    default = function (i, k)
        u = w = false
        if !m.periodic_boundary
            i == 0 && (u |= m.velocity_bc_left_normal; w |= m.velocity_bc_left_tangential)
            i == m.Nx - 1 && (u |= m.velocity_bc_right_normal; w |= m.velocity_bc_right_tangential)
        end
        k == 0 && (w |= m.velocity_bc_bottom_normal; u |= m.velocity_bc_bottom_tangential)
        k == m.Nz - 1 && (w |= m.velocity_bc_top_normal; u |= m.velocity_bc_top_tangential)
        (Float64(!u), Float64(!w))
    end
    velocity_bc_mask = fill_field!(gv(), da_Veloc, mask === nothing ? default : mask, 2)
    velocity_bc_values = gv()
    LibPETSc.VecCopy(P, velocity, velocity_bc_values)
    ksp = LibPETSc.KSPCreate(P, MPI.COMM_WORLD)
    LibPETSc.KSPSetOptionsPrefix(P, ksp, "veloc_")
    if m.use_direct_solver
        pc = LibPETSc.KSPGetPC(P, ksp)
        LibPETSc.PCSetType(P, pc, "lu"); LibPETSc.PCFactorSetMatSolverType(P, pc, "mumps")
    end
    s = StokesSolver(m, StokesElement(m, ElementMatrices(m)), da_Veloc, da_Thermal,
        LibPETSc.DMCreateMatrix(P, da_Veloc), LibPETSc.DMCreateMatrix(P, da_Veloc), ksp,
        gv(), velocity_bc_values, velocity_bc_mask, gv(), gv(), gv(), gv(), gv(), gv(), gv(), gv(),
        lv(), lv(), lv(), lv(), lv(), lt(), lt(), lt(), 0, 0.0)
    global_to_local(da_Veloc, velocity_bc_mask, s.local_velocity_bc_mask)
    assemble_gradient!(s)
    return s
end

local_row(da) = ((gx, gz, _, gm, _, _) = Int.(LibPETSc.DMDAGetGhostCorners(P, da)); (i, k, c) -> 2((i - gx) + (k - gz) * gm) + c)
finish_matrix(A) = (LibPETSc.MatAssemblyBegin(P, A, LibPETSc.MAT_FINAL_ASSEMBLY); LibPETSc.MatAssemblyEnd(P, A, LibPETSc.MAT_FINAL_ASSEMBLY))
function element_nodes(m, ei, ek)        # nós das equações: na borda periódica, Nx no lugar de Nx-1
    right = m.periodic_boundary && ei == m.Nx - 2 ? m.Nx : ei + 1
    return ((ei, ek), (right, ek), (ei, ek + 1), (right, ek + 1))
end

function assemble_gradient!(s)
    G, row = s.pressure_gradient_matrix, local_row(s.da_Veloc)
    rows, col = zeros(PInt, 8), zeros(PInt, 1)
    xs, zs = element_ranges(s.m, s.da_Veloc)
    for ek in zs, ei in xs
        for (n, (i, k)) in enumerate(element_nodes(s.m, ei, ek)); rows[2n-1] = row(i, k, 0); rows[2n] = row(i, k, 1); end
        col[1] = row(ei, ek, 0)
        LibPETSc.MatSetValuesLocal(P, G, PInt(8), rows, PInt(1), col, s.element.VCe, LibPETSc.ADD_VALUES)
    end
    finish_matrix(G)
end

# Termo do FSSA (Kaus et al., 2010) somado à diagonal dos w do elemento. `density_variation`
# já é a força de volume (-g*rho), então a força muda de theta*dt*w*d(dRho)/dz.
function fssa_term(m, d, dt)
    jump = (-d[1] - d[2] + d[3] + d[4]) / 2                         # rho_mean_bottom
    consistent = m.fssa_theta * dt * jump * element_size_x(m) / 4   # integral nodal: 4 nós, 1/4 cada
    # C++: theta*dt*jump*g*dx/2, num sistema 4x maior
    return m.fssa_legacy ? m.fssa_theta * dt * jump * m.gravity_acceleration * element_size_x(m) / 2 / 4 : consistent
end

function assemble_velocity!(s, velocity, pressure, viscosity_factor, density_variation, dt; fssa_density_variation = density_variation)  # build_veloc
    m, A, el = s.m, s.velocity_stiffness_matrix, s.element
    LibPETSc.MatZeroEntries(P, A)
    global_to_local(s.da_Veloc, velocity, s.local_velocity)
    global_to_local(s.da_Veloc, pressure, s.local_pressure)
    global_to_local(s.da_Thermal, viscosity_factor, s.local_viscosity_factor)
    global_to_local(s.da_Thermal, density_variation, s.local_density_variation)
    global_to_local(s.da_Thermal, fssa_density_variation, s.local_fssa_density_variation)
    LibPETSc.VecZeroEntries(P, s.local_pressure_rhs); LibPETSc.VecZeroEntries(P, s.local_uzawa_preconditioner)
    row = local_row(s.da_Veloc)
    PETSc.with_local_array!(s.local_pressure_rhs, s.local_uzawa_preconditioner, s.local_velocity, s.local_pressure,
        s.local_velocity_bc_mask, s.local_viscosity_factor, s.local_density_variation, s.local_fssa_density_variation;
        write = (true, true, false, false, false, false, false, false)) do f, pc, v, p, c, q, d, d0
        l2(a) = local_view(s.da_Veloc, a, 2); l1(a) = local_view(s.da_Thermal, a)
        F, PC, V, Pr, C, visc, D, D0 = l2(f), l2(pc), l2(v), l2(p), l2(c), l1(q), l1(d), l1(d0)
        si, sk, _, mi, mk, _ = Int.(LibPETSc.DMDAGetCorners(P, s.da_Veloc))
        r, r2, one_value = PInt[0], PInt[0], [0.0]
        put(a, b, x) = (r[1] = a; r2[1] = b; one_value[1] = x;
            LibPETSc.MatSetValuesLocal(P, A, PInt(1), r, PInt(1), r2, one_value, LibPETSc.ADD_VALUES))
        for k in sk:sk+mk-1, i in si:si+mi-1, n in 0:1
            if m.periodic_boundary && i == m.Nx - 1
                put(row(i, k, n), row(i, k, n), 1.0); put(row(i, k, n), row(m.Nx, k, n), -1.0)
            else
                put(row(i, k, n), row(i, k, n), 1.0 - C[n+1, i, k])
            end
        end
        Ke, values, rows, free, fe, de = zeros(8, 8), zeros(64), zeros(PInt, 8), zeros(Bool, 8), zeros(8), zeros(4)
        xs, zs = element_ranges(m, s.da_Veloc)
        for ek in zs, ei in xs
            Ke .= visc[ei, ek] .* el.Ke
            de[1] = D[ei, ek]; de[2] = D[ei+1, ek]; de[3] = D[ei, ek+1]; de[4] = D[ei+1, ek+1]
            if m.free_surface_stabilization
                x = fssa_term(m, (D0[ei, ek], D0[ei+1, ek], D0[ei, ek+1], D0[ei+1, ek+1]), dt)
                for n in (2, 4, 6, 8); Ke[n, n] += x; end
            end
            nodes = element_nodes(m, ei, ek)
            for n in 1:4
                i, k = nodes[n]
                rows[2n-1] = row(i, k, 0); rows[2n] = row(i, k, 1)
                free[2n-1] = C[1, i, k] != 0; free[2n] = C[2, i, k] != 0
            end
            for a in 1:8, b in 1:8; values[8(a-1)+b] = free[a] ? Ke[a, b] : 0.0; end
            LibPETSc.MatSetValuesLocal(P, A, PInt(8), rows, PInt(8), rows, values, LibPETSc.ADD_VALUES)
            for a in 1:8
                PC[1, ei, ek] += el.VCe[a]^2 / Ke[a, a]
                fe[a] = -el.VCe[a] * Pr[1, ei, ek]
                for j in 1:4; fe[a] += el.VfMe[a, j] * de[j]; end
            end
            for n in 1:4, c in 1:2
                free[2n-2+c] && (F[c, nodes[n]...] += fe[2n-2+c])
            end
        end
        for k in sk:sk+mk-1, i in si:si+mi-1, n in 1:2
            if m.periodic_boundary && i == m.Nx - 1
                F[n, i, k] = 0.0
            elseif C[n, i, k] == 0
                F[n, i, k] = V[n, i, k]
            end
        end
    end
    finish_matrix(A)
    for (l, g) in ((s.local_pressure_rhs, s.velocity_rhs_pressure), (s.local_uzawa_preconditioner, s.uzawa_preconditioner))
        LibPETSc.VecZeroEntries(P, g)
        LibPETSc.DMLocalToGlobal(P, s.da_Veloc, l, LibPETSc.ADD_VALUES, g)
    end
    LibPETSc.VecReciprocal(P, s.uzawa_preconditioner)
end

function solve_velocity!(s, b, x)
    LibPETSc.KSPSolve(P, s.velocity_ksp, b, x)
    Int(LibPETSc.KSPGetConvergedReason(P, s.velocity_ksp)) > 0 || error("KSP da velocidade não convergiu")
end

"Resolve o Stokes: velocity_next e pressure são atualizados no lugar (veloc_total, caso linear)."
function stokes_solve!(s::StokesSolver, velocity_next, pressure, velocity, viscosity_factor, density_variation, dt;
    max_iterations = 400, fssa_density_variation = density_variation)
    assemble_velocity!(s, velocity, pressure, viscosity_factor, density_variation, dt; fssa_density_variation)
    L, K, G = LibPETSc, s.velocity_stiffness_matrix, s.pressure_gradient_matrix
    M, mask = s.uzawa_preconditioner, s.velocity_bc_mask
    L.KSPSetOperators(P, s.velocity_ksp, K, K)
    L.KSPSetFromOptions(P, s.velocity_ksp)
    L.KSPSetInitialGuessNonzero(P, s.velocity_ksp, L.PETSC_TRUE)
    L.KSPSetTolerances(P, s.velocity_ksp, s.m.rtol, 1e-50, 1e5, PInt(10000))
    solve_velocity!(s, s.velocity_rhs_pressure, velocity_next)             # K V0 = F
    L.VecPointwiseMult(P, velocity_next, velocity_next, mask)
    L.VecAXPY(P, velocity_next, 1.0, s.velocity_bc_values)
    L.MatMultTranspose(P, G, velocity_next, s.rk_vec2)                     # r0 = G^T V0
    L.VecPointwiseMult(P, s.zk_vec2, M, s.rk_vec2)
    # 16 = 4^2: o G do Mandyoc é 4x maior, e assim a tolerância tem o mesmo sentido
    residual() = 16 * Float64(L.VecDot(P, s.rk_vec2, s.rk_vec2))           # denok
    s.uzawa_residual = residual(); s.uzawa_iterations = 0
    for k in 1:max_iterations-1
        s.uzawa_residual > s.m.uzawa_residual_tolerance || break
        if k == 1
            L.VecCopy(P, s.zk_vec2, s.sk_vec)
        else
            L.VecCopy(P, s.zk_vec2, s.zk_vec)
            L.VecPointwiseMult(P, s.zk_vec2, M, s.rk_vec2)
            beta = L.VecDot(P, s.zk_vec2, s.rk_vec2) / L.VecDot(P, s.zk_vec, s.rk_vec)   # betak
            L.VecAYPX(P, s.sk_vec, beta, s.zk_vec2)
        end
        L.MatMult(P, G, s.sk_vec, s.gs_vec)
        L.VecPointwiseMult(P, s.gs_vec, s.gs_vec, mask)
        solve_velocity!(s, s.gs_vec, s.uk_vec)                             # K uk = G sk
        L.VecPointwiseMult(P, s.uk_vec, s.uk_vec, mask)
        alpha = L.VecDot(P, s.zk_vec2, s.rk_vec2) / L.VecDot(P, s.gs_vec, s.uk_vec)      # alphak
        L.VecAXPY(P, pressure, alpha, s.sk_vec)
        L.VecAXPY(P, velocity_next, -alpha, s.uk_vec)
        L.VecCopy(P, s.rk_vec2, s.rk_vec)
        L.MatMultTranspose(P, G, s.uk_vec, s.rk_vec2)
        L.VecAYPX(P, s.rk_vec2, -alpha, s.rk_vec)
        s.uzawa_iterations = k; s.uzawa_residual = residual()
    end
    s.uzawa_residual > s.m.uzawa_residual_tolerance && MPI.Comm_rank(MPI.COMM_WORLD) == 0 &&
        @warn "Uzawa parou sem atingir a tolerância" s.uzawa_residual
    return velocity_next
end

"density_variation (dRho) = a*T - g*rho, com a = alfa*g*rho_manto (calc_drho)."
function calc_density_variation!(m, density_variation, temperature, density_field)
    a = m.thermal_expansion_coefficient * m.gravity_acceleration * m.mantle_density
    PETSc.with_local_array!(density_variation, temperature, density_field; write = (true, false, false)) do d, T, rho
        d .= a .* T .- m.gravity_acceleration .* rho
    end
end

"Partículas -> malha (Swarm2Mesh_2d): densidade nos nós, viscosidade média por elemento."
function swarm_to_mesh!(m, particle_swarm, da_Thermal, density_field, viscosity_factor, layer_density, layer_viscosity)
    lv() = LibPETSc.DMCreateLocalVector(P, da_Thermal); gv() = LibPETSc.DMCreateGlobalVector(P, da_Thermal)
    l = (lv(), lv(), lv(), lv()); g = (gv(), gv())                  # rho, peso, visc, contagem
    dx, dz = element_size_x(m), element_size_z(m)
    PETSc.with_local_array!(l...) do ra, ca, va, na
        R, W, V, N = local_view(da_Thermal, ra), local_view(da_Thermal, ca), local_view(da_Thermal, va), local_view(da_Thermal, na)
        fill!(R, 0); fill!(W, 0); fill!(V, 0); fill!(N, 0)
        with_field(particle_swarm, PIC_COORDINATES, Float64) do coor
        with_field(particle_swarm, "layer", PInt) do layer
            for p in eachindex(layer)
                x, z = coor[2p-1], coor[2p] + m.domain_depth
                i = clamp(trunc(Int, x / dx), 0, m.Nx - 2); k = clamp(trunc(Int, z / dz), 0, m.Nz - 2)
                rx = x / dx - i; rz = z / dz - k
                rho = layer_density[layer[p]+1]; eta = layer_viscosity[layer[p]+1]
                for (ii, kk, w) in ((i, k, (1 - rx) * (1 - rz)), (i + 1, k, rx * (1 - rz)), (i, k + 1, (1 - rx) * rz), (i + 1, k + 1, rx * rz))
                    R[ii, kk] += w * rho; W[ii, kk] += w
                end
                V[i, k] += m.harmonic_viscosity_mean ? 1 / eta : eta; N[i, k] += 1
            end
        end end
    end
    for (num, den, out) in ((l[1], l[2], density_field), (l[3], l[4], viscosity_factor))
        for (a, b) in ((num, g[1]), (den, g[2]))
            LibPETSc.VecZeroEntries(P, b); LibPETSc.DMLocalToGlobal(P, da_Thermal, a, LibPETSc.ADD_VALUES, b)
        end
        harmonic = out === viscosity_factor && m.harmonic_viscosity_mean
        PETSc.with_local_array!(out, g[1], g[2]) do o, a, b
            for n in eachindex(o)                                   # sem partículas: mantém o valor anterior
                b[n] > 0 && (o[n] = harmonic ? b[n] / a[n] : a[n] / b[n])
            end
        end
    end
    foreach(v -> LibPETSc.VecDestroy(P, v), (l..., g...))
end
