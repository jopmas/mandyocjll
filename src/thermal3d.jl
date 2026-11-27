# 3D mesh and thermal equation (DMT.cpp, DM1_3d.cpp and thermal_Ke.cpp in 3D): trilinear hexahedra, theta scheme
# and the streamline artificial diffusion of the 2D code. Include after thermal.jl.
# The C++ code uses the diffusivity of the mantle in 3D; here, as in 2D, that of each element (the k row of
# interfaces.txt and high_kappa_in_asthenosphere).

# ---- meshes (create_veloc and create_thermal in 3D) ---------------------------------
function create_mesh3d(m::ModelParameters, dof, width)
    m.periodic_boundary && error("periodic_boundary is not available in 3D (as in the C++ code)")
    none = LibPETSc.DM_BOUNDARY_NONE
    da = LibPETSc.DMDACreate3d(P, MPI.COMM_WORLD, none, none, none, LibPETSc.DMDA_STENCIL_BOX,
        PInt.((m.Nx, m.Ny, m.Nz, m.Px, m.Py, m.Pz, dof, width))..., C_NULL, C_NULL, C_NULL)
    LibPETSc.DMSetUp(P, da)
    LibPETSc.DMDASetUniformCoordinates(P, da, 0.0, m.Lx, 0.0, m.Ly, -m.domain_depth, 0.0)
    return da
end
create_velocity_mesh3d(m) = create_mesh3d(m, 3, 1)                  # da_Veloc
create_thermal_mesh3d(m) = create_mesh3d(m, 1, 1)                   # da_Thermal
# Same partition as da_Veloc, with ghosts for the RK4 stages in the neighbor element (as create_advection_mesh).
function create_advection_mesh3d(m, da_Veloc)
    da = create_mesh3d(m, 3, m.advection_stencil_width)
    @assert LibPETSc.DMDAGetCorners(P, da) == LibPETSc.DMDAGetCorners(P, da_Veloc)
    return da
end

"Array with global indices [dof,] i, j, k, ghosts included."
function local_view3(da, a, dof = 1)
    gx, gy, gz, gm, gn, gp = Int.(LibPETSc.DMDAGetGhostCorners(P, da))
    return dof == 1 ? OffsetArray(reshape(a, gm, gn, gp), gx:gx+gm-1, gy:gy+gn-1, gz:gz+gp-1) :
        OffsetArray(reshape(a, dof, gm, gn, gp), 1:dof, gx:gx+gm-1, gy:gy+gn-1, gz:gz+gp-1)
end

"Owned part of a global Vec with global indices [dof,] i, j, k."
function owned_view3(da, a, dof = 1)
    sx, sy, sz, mx, my, mz = Int.(LibPETSc.DMDAGetCorners(P, da))
    return dof == 1 ? OffsetArray(reshape(a, mx, my, mz), sx:sx+mx-1, sy:sy+my-1, sz:sz+mz-1) :
        OffsetArray(reshape(a, dof, mx, my, mz), 1:dof, sx:sx+mx-1, sy:sy+my-1, sz:sz+mz-1)
end

"Fills a global Vec of the 3D mesh with f(i, j, k) (a tuple if dof > 1)."
function fill_field3!(vec, da, f, dof = 1)
    PETSc.with_local_array!(vec) do a
        A = owned_view3(da, a, dof)
        for k in axes(A, dof == 1 ? 3 : 4), j in axes(A, dof == 1 ? 2 : 3), i in axes(A, dof == 1 ? 1 : 2)
            v = f(i, j, k)
            dof == 1 ? (A[i, j, k] = v) : (for c in 1:dof; A[c, i, j, k] = v[c]; end)
        end
    end
    return vec
end

"Elements of the rank: those of its lower left front node (DMDAGetElementCorners_3d)."
function element_ranges3(m, da)
    sx, sy, sz, mx, my, mz = Int.(LibPETSc.DMDAGetCorners(P, da))
    return sx:min(sx + mx, m.Nx - 1)-1, sy:min(sy + my, m.Ny - 1)-1, sz:min(sz + mz, m.Nz - 1)-1
end

# the 8 nodes of an element, x first, then y, then z (ind[0..7] of the C++ code)
element_nodes3(i, j, k) = ((i, j, k), (i + 1, j, k), (i, j + 1, k), (i + 1, j + 1, k),
    (i, j, k + 1), (i + 1, j, k + 1), (i, j + 1, k + 1), (i + 1, j + 1, k + 1))

# ---- element matrices (montaKeThermal_general_3d) ------------------------------------
struct ElementMatrices3D
    NT::Matrix{Float64}; NT_x::Matrix{Float64}; NT_y::Matrix{Float64}; NT_z::Matrix{Float64}   # [node, point]
    gauss_weight::Vector{Float64}        # prodH*volume, with the 1/8 Jacobian that Mandyoc leaves out
    TKe::Matrix{Float64}; TMe::Matrix{Float64}; TFe::Vector{Float64}
end

function ElementMatrices3D(m::ModelParameters)
    dx, dy, dz = element_size_x(m), element_size_y(m), element_size_z(m)
    g = (-sqrt(0.6), 0.0, sqrt(0.6)); w = (5 / 9, 8 / 9, 5 / 9)
    NT, NT_x, NT_y, NT_z = zeros(8, 27), zeros(8, 27), zeros(8, 27), zeros(8, 27); gw = zeros(27)
    q = 0
    for c in 1:3, b in 1:3, a in 1:3
        q += 1; kx, ky, kz = g[a], g[b], g[c]; gw[q] = w[a] * w[b] * w[c] * dx * dy * dz / 8
        n = 0
        for ez in (-1, 1), ey in (-1, 1), ex in (-1, 1)
            n += 1
            NT[n, q] = (1 + ex * kx) * (1 + ey * ky) * (1 + ez * kz) / 8
            NT_x[n, q] = ex * (1 + ey * ky) * (1 + ez * kz) / 4 / dx
            NT_y[n, q] = (1 + ex * kx) * ey * (1 + ez * kz) / 4 / dy
            NT_z[n, q] = (1 + ex * kx) * (1 + ey * ky) * ez / 4 / dz
        end
    end
    TKe = [sum(gw[q] * (NT_x[i, q] * NT_x[j, q] + NT_y[i, q] * NT_y[j, q] + NT_z[i, q] * NT_z[j, q]) for q in 1:27) for i in 1:8, j in 1:8]
    TMe = [sum(gw[q] * NT[i, q] * NT[j, q] for q in 1:27) for i in 1:8, j in 1:8]
    TFe = [sum(gw[q] * NT[i, q] for q in 1:27) / m.heat_capacity for i in 1:8]
    return ElementMatrices3D(NT, NT_x, NT_y, NT_z, gw, TKe, TMe, TFe)
end

# montaKeThermal_simplif_3d: advection + artificial diffusion + diffusion
function advection_diffusion3!(TCe, e::ElementMatrices3D, m, u, v, w, kappa)
    dx, dy, dz = element_size_x(m), element_size_y(m), element_size_z(m)
    eps_upwind(a) = a < -1 ? -1 - 1 / a : a <= 1 ? 0.0 : 1 - 1 / a
    kappa_arti = (eps_upwind(u * dx / 2kappa) * u * dx + eps_upwind(v * dy / 2kappa) * v * dy + eps_upwind(w * dz / 2kappa) * w * dz) / 2
    speed2 = u^2 + v^2 + w^2
    fill!(TCe, 0.0)
    @inbounds for q in 1:27, j in 1:8
        advection = e.gauss_weight[q] * (e.NT_x[j, q] * u + e.NT_y[j, q] * v + e.NT_z[j, q] * w)
        for i in 1:8
            upwind = speed2 > 0 ? kappa_arti * (e.NT_x[i, q] * u + e.NT_y[i, q] * v + e.NT_z[i, q] * w) / speed2 : 0.0
            TCe[i, j] += (e.NT[i, q] + upwind) * advection
        end
    end
    TCe .+= kappa .* e.TKe
    return TCe
end

function effective_diffusivity3(m, kappa, T, i, j, k)
    m.high_diffusivity_in_asthenosphere || return kappa
    mean = sum(T[n...] for n in element_nodes3(i, j, k)) / 8
    mean > 1330.0 && (kappa = m.thermal_diffusivity * (1 + (mean - 1330.0) * 5))
    return min(kappa, 50m.thermal_diffusivity)
end

mean_velocity3(V, i, j, k) = ntuple(c -> sum(V[c, n...] for n in element_nodes3(i, j, k)) / 8, 3)

# ---- solver ---------------------------------------------------------------------------
mutable struct ThermalSolver3D{D, V, M, K}
    m::ModelParameters
    element::ElementMatrices3D
    da_Thermal::D
    da_Veloc::D
    thermal_matrix_lhs::M                # TA
    thermal_rhs::V                       # Tf
    thermal_ksp::K                       # T_ksp
    temperature_bc_values::V             # Temper_0
    temperature_bc_mask::V               # Temper_Cond (1 = free, 0 = fixed)
    local_temperature::V
    local_temperature_bc_mask::V
    local_thermal_rhs::V
    local_radiogenic_heat_field::V
    local_thermal_diffusivity_field::V
    local_velocity::V
    iterations::Int
end

function ThermalSolver3D(m::ModelParameters, da_Thermal, da_Veloc, temperature)
    global_vec() = LibPETSc.DMCreateGlobalVector(P, da_Thermal)
    local_vec() = LibPETSc.DMCreateLocalVector(P, da_Thermal)
    LibPETSc.DMSetMatType(P, da_Thermal, "aij")
    mask = fill_field3!(global_vec(), da_Thermal, function (i, j, k)
        fixed = (k == 0 && m.temperature_bc_bottom) || (k == m.Nz - 1 && m.temperature_bc_top) ||
            (i == 0 && m.temperature_bc_left) || (i == m.Nx - 1 && m.temperature_bc_right) ||
            (j == 0 && m.temperature_bc_front) || (j == m.Ny - 1 && m.temperature_bc_back)
        fixed ? 0.0 : 1.0
    end)
    bc_values = global_vec()
    PETSc.with_local_array!(bc_values, temperature, mask; write = (true, false, false)) do b, T, c
        b .= T .* (1 .- c)
    end
    ksp = LibPETSc.KSPCreate(P, MPI.COMM_WORLD)
    LibPETSc.KSPSetOptionsPrefix(P, ksp, "thermal_")
    s = ThermalSolver3D(m, ElementMatrices3D(m), da_Thermal, da_Veloc, LibPETSc.DMCreateMatrix(P, da_Thermal),
        global_vec(), ksp, bc_values, mask, local_vec(), local_vec(), local_vec(), local_vec(), local_vec(),
        LibPETSc.DMCreateLocalVector(P, da_Veloc), 0)
    global_to_local(da_Thermal, mask, s.local_temperature_bc_mask)
    return s
end

function local_row3(da)
    gx, gy, gz, gm, gn, _ = Int.(LibPETSc.DMDAGetGhostCorners(P, da))
    return (i, j, k) -> (i - gx) + (j - gy) * gm + (k - gz) * gm * gn
end

function assemble_lhs3!(s::ThermalSolver3D, velocity_next, dt)      # AssembleA_Thermal_3d
    m, A = s.m, s.thermal_matrix_lhs
    LibPETSc.MatZeroEntries(P, A)
    global_to_local(s.da_Veloc, velocity_next, s.local_velocity)
    row = local_row3(s.da_Thermal)
    PETSc.with_local_array!(s.local_velocity, s.local_temperature, s.local_temperature_bc_mask,
        s.local_thermal_diffusivity_field; write = false) do v, t, c, kf
        V = local_view3(s.da_Veloc, v, 3)
        T, C, kappa = local_view3(s.da_Thermal, t), local_view3(s.da_Thermal, c), local_view3(s.da_Thermal, kf)
        sx, sy, sz, mx, my, mz = Int.(LibPETSc.DMDAGetCorners(P, s.da_Thermal))
        r1 = PInt[0]; one_value = [0.0]
        for k in sz:sz+mz-1, j in sy:sy+my-1, i in sx:sx+mx-1
            r1[1] = row(i, j, k); one_value[1] = 1.0 - C[i, j, k]
            LibPETSc.MatSetValuesLocal(P, A, PInt(1), r1, PInt(1), r1, one_value, LibPETSc.ADD_VALUES)
        end
        TCe = zeros(8, 8); values = zeros(64); index = zeros(PInt, 8)
        xs, ys, zs = element_ranges3(m, s.da_Thermal)
        for ek in zs, ej in ys, ei in xs
            u, vy, w = mean_velocity3(V, ei, ej, ek)
            advection_diffusion3!(TCe, s.element, m, u, vy, w, effective_diffusivity3(m, kappa[ei, ej, ek], T, ei, ej, ek))
            nodes = element_nodes3(ei, ej, ek)
            for a in 1:8
                index[a] = row(nodes[a]...)
                free = C[nodes[a]...] != 0
                for b in 1:8
                    values[8(a-1)+b] = free ? s.element.TMe[a, b] + m.thermal_time_theta * dt * TCe[a, b] : 0.0
                end
            end
            LibPETSc.MatSetValuesLocal(P, A, PInt(8), index, PInt(8), index, values, LibPETSc.ADD_VALUES)
        end
    end
    LibPETSc.MatAssemblyBegin(P, A, LibPETSc.MAT_FINAL_ASSEMBLY)
    LibPETSc.MatAssemblyEnd(P, A, LibPETSc.MAT_FINAL_ASSEMBLY)
end

function assemble_rhs3!(s::ThermalSolver3D, velocity, dt)           # AssembleF_Thermal_3d
    m, e = s.m, s.element
    one_minus_theta = 1 - m.thermal_time_theta
    global_to_local(s.da_Veloc, velocity, s.local_velocity)
    LibPETSc.VecZeroEntries(P, s.local_thermal_rhs)
    PETSc.with_local_array!(s.local_thermal_rhs, s.local_velocity, s.local_temperature, s.local_temperature_bc_mask,
        s.local_thermal_diffusivity_field, s.local_radiogenic_heat_field;
        write = (true, false, false, false, false, false)) do f, v, t, c, kf, h
        V = local_view3(s.da_Veloc, v, 3); lv(a) = local_view3(s.da_Thermal, a)
        F, T, C, kappa, H = lv(f), lv(t), lv(c), lv(kf), lv(h)
        TCe = zeros(8, 8)
        xs, ys, zs = element_ranges3(m, s.da_Thermal)
        for ek in zs, ej in ys, ei in xs
            u, vy, w = mean_velocity3(V, ei, ej, ek)
            advection_diffusion3!(TCe, e, m, u, vy, w, effective_diffusivity3(m, kappa[ei, ej, ek], T, ei, ej, ek))
            nodes = element_nodes3(ei, ej, ek)
            heat = 0.0                                              # H_efetivo
            m.use_radiogenic_heating && (heat += sum(H[n...] for n in nodes) / 8)
            if m.use_adiabatic_heating
                kelvin = sum(T[n...] for n in nodes) / 8 + 273.0
                heat -= kelvin * m.thermal_expansion_coefficient * m.gravity_acceleration * w
            end
            for a in 1:8
                C[nodes[a]...] == 0 && continue
                x = dt * e.TFe[a] * heat
                for b in 1:8
                    x += T[nodes[b]...] * (e.TMe[a, b] - one_minus_theta * dt * TCe[a, b])
                end
                F[nodes[a]...] += x
            end
        end
        sx, sy, sz, mx, my, mz = Int.(LibPETSc.DMDAGetCorners(P, s.da_Thermal))
        for k in sz:sz+mz-1, j in sy:sy+my-1, i in sx:sx+mx-1
            C[i, j, k] == 0 && (F[i, j, k] = T[i, j, k])
        end
    end
    LibPETSc.VecZeroEntries(P, s.thermal_rhs)
    LibPETSc.DMLocalToGlobal(P, s.da_Thermal, s.local_thermal_rhs, LibPETSc.ADD_VALUES, s.thermal_rhs)
end

"Heat flow at the base in 3D (only in 2D in the C++ code): the bottom plane is the plane above plus the step."
function basal_heat_flow3!(s::ThermalSolver3D, temperature)
    m = s.m
    step = m.basal_heat_flow * element_size_z(m) / (m.thermal_diffusivity * m.mantle_density * m.heat_capacity)
    sz, mz = Int.(LibPETSc.DMDAGetCorners(P, s.da_Thermal))[[3, 6]]
    sz == 0 || return temperature
    mz >= 2 || error("basal_heat needs at least two planes of nodes per rank in z")
    PETSc.with_local_array!(temperature) do t
        T = owned_view3(s.da_Thermal, t)
        T[:, :, 0] .= T[:, :, 1] .+ step
    end
    return temperature
end

function thermal_step!(s::ThermalSolver3D, temperature, velocity, velocity_next,
    radiogenic_heat_field, thermal_diffusivity_field, time_step_seconds; melt_increment_field = nothing)
    melt_increment_field === nothing || error("magmatism is not available in 3D")
    global_to_local(s.da_Thermal, temperature, s.local_temperature)
    global_to_local(s.da_Thermal, radiogenic_heat_field, s.local_radiogenic_heat_field)
    global_to_local(s.da_Thermal, thermal_diffusivity_field, s.local_thermal_diffusivity_field)
    assemble_lhs3!(s, velocity_next, time_step_seconds)
    assemble_rhs3!(s, velocity, time_step_seconds)
    A, ksp = s.thermal_matrix_lhs, s.thermal_ksp
    LibPETSc.KSPSetOperators(P, ksp, A, A)
    LibPETSc.KSPSetFromOptions(P, ksp)
    LibPETSc.KSPSetTolerances(P, ksp, s.m.solver_relative_tolerance, 1e-50, 1e5, PInt(10000))
    LibPETSc.KSPSolve(P, ksp, s.thermal_rhs, temperature)
    s.iterations = LibPETSc.KSPGetIterationNumber(P, ksp)
    PETSc.with_local_array!(temperature, s.temperature_bc_mask, s.temperature_bc_values; write = (true, false, false)) do T, c, b
        T .= T .* c .+ b
    end
    s.m.basal_heat_flow > 0 && basal_heat_flow3!(s, temperature)
    Int(LibPETSc.KSPGetConvergedReason(P, ksp)) > 0 || error("the thermal KSP did not converge")
    return temperature
end
