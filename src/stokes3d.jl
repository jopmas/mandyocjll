# 3D Stokes (DMv.cpp, DM1_v.cpp in 3D): Q1-P0 hexahedra and the Uzawa conjugate gradients of the 2D code.
# Include after stokes.jl and thermal3d.jl. The pressure of an element is in component 1 (u) of its lower left
# front node of da_Veloc, as in Mandyoc. The matrices have the 1/8 Jacobian that Mandyoc leaves out.

struct StokesElement3D
    Ke_veloc::Matrix{Float64}            # Ke_veloc for unit viscosity (24 x 24)
    VCe::Vector{Float64}                 # VCe: pressure gradient
    VfMe::Matrix{Float64}                # VfMe: body force (24 x 8)
end
function StokesElement3D(e::ElementMatrices3D)    # montaKeVeloc_general_3d, montaCeVeloc3d, montafeVeloc3d
    Ke, VCe, VfMe = zeros(24, 24), zeros(24), zeros(24, 8)
    for q in 1:27
        SN = zeros(6, 24)
        for n in 1:8
            x, y, z = e.NT_x[n, q], e.NT_y[n, q], e.NT_z[n, q]; a = 3(n - 1)
            SN[1, a+1] = x; SN[2, a+2] = y; SN[3, a+3] = z
            SN[4, a+1] = y; SN[4, a+2] = x
            SN[5, a+2] = z; SN[5, a+3] = y
            SN[6, a+1] = z; SN[6, a+3] = x
        end
        for i in 1:24, j in 1:24
            Ke[i, j] += e.gauss_weight[q] * (2 * (SN[1, i] * SN[1, j] + SN[2, i] * SN[2, j] + SN[3, i] * SN[3, j]) +
                SN[4, i] * SN[4, j] + SN[5, i] * SN[5, j] + SN[6, i] * SN[6, j])
        end
        for n in 1:8
            a = 3(n - 1)
            VCe[a+1] -= e.gauss_weight[q] * e.NT_x[n, q]; VCe[a+2] -= e.gauss_weight[q] * e.NT_y[n, q]
            VCe[a+3] -= e.gauss_weight[q] * e.NT_z[n, q]
            for j in 1:8; VfMe[a+3, j] += e.gauss_weight[q] * e.NT[n, q] * e.NT[j, q]; end
        end
    end
    return StokesElement3D(Ke, VCe, VfMe)
end

mutable struct StokesSolver3D{D, V, M, K}
    m::ModelParameters
    element::StokesElement3D
    da_Veloc::D
    da_Thermal::D
    velocity_stiffness_matrix::M         # VA
    pressure_gradient_matrix::M          # VG
    velocity_ksp::K                      # V_ksp
    velocity_rhs_pressure::V             # Vf_P
    velocity_bc_values::V                # Veloc_0
    velocity_bc_mask::V                  # Veloc_Cond (1 = free, 0 = fixed)
    uzawa_preconditioner::V              # Precon
    rk_vec::V; rk_vec2::V; sk_vec::V; gs_vec::V; uk_vec::V; zk_vec::V; zk_vec2::V
    local_velocity::V
    local_velocity_bc_mask::V
    local_pressure_rhs::V
    local_pressure::V
    local_uzawa_preconditioner::V
    local_viscosity_factor::V
    local_density_variation::V
    local_fssa_density_variation::V
    uzawa_iterations::Int
    uzawa_residual::Float64              # denok
end

"""
Boundary mask of the parameters (create_veloc in 3D): normal and tangential (both components) per face. As in the
C++ code, `back` is y = 0 and `front` is y = Ly for the velocity (the opposite of the temperature).
"""
function default_velocity_mask3(m)
    return function (i, j, k)
        u = v = w = false
        i == 0 && (u |= m.velocity_bc_left_normal; m.velocity_bc_left_tangential && (v = w = true))
        i == m.Nx - 1 && (u |= m.velocity_bc_right_normal; m.velocity_bc_right_tangential && (v = w = true))
        j == 0 && (v |= m.velocity_bc_back_normal; m.velocity_bc_back_tangential && (u = w = true))
        j == m.Ny - 1 && (v |= m.velocity_bc_front_normal; m.velocity_bc_front_tangential && (u = w = true))
        k == 0 && (w |= m.velocity_bc_bottom_normal; m.velocity_bc_bottom_tangential && (u = v = true))
        k == m.Nz - 1 && (w |= m.velocity_bc_top_normal; m.velocity_bc_top_tangential && (u = v = true))
        (Float64(!u), Float64(!v), Float64(!w))
    end
end

"mask = (i, j, k) -> (u, v, w) replaces the boundary conditions of the parameters."
function StokesSolver3D(m, da_Veloc, da_Thermal, velocity; mask = nothing)
    gv() = LibPETSc.DMCreateGlobalVector(P, da_Veloc); lv() = LibPETSc.DMCreateLocalVector(P, da_Veloc)
    lt() = LibPETSc.DMCreateLocalVector(P, da_Thermal)
    LibPETSc.DMSetMatType(P, da_Veloc, "aij")
    velocity_bc_mask = fill_field3!(gv(), da_Veloc, mask === nothing ? default_velocity_mask3(m) : mask, 3)
    velocity_bc_values = gv()
    LibPETSc.VecCopy(P, velocity, velocity_bc_values)
    ksp = LibPETSc.KSPCreate(P, MPI.COMM_WORLD)
    LibPETSc.KSPSetOptionsPrefix(P, ksp, "veloc_")
    if m.use_direct_solver
        pc = LibPETSc.KSPGetPC(P, ksp)
        LibPETSc.PCSetType(P, pc, "lu"); LibPETSc.PCFactorSetMatSolverType(P, pc, "mumps")
    end
    s = StokesSolver3D(m, StokesElement3D(ElementMatrices3D(m)), da_Veloc, da_Thermal,
        LibPETSc.DMCreateMatrix(P, da_Veloc), LibPETSc.DMCreateMatrix(P, da_Veloc), ksp,
        gv(), velocity_bc_values, velocity_bc_mask, gv(), gv(), gv(), gv(), gv(), gv(), gv(), gv(),
        lv(), lv(), lv(), lv(), lv(), lt(), lt(), lt(), 0, 0.0)
    global_to_local(da_Veloc, velocity_bc_mask, s.local_velocity_bc_mask)
    assemble_gradient3!(s)
    return s
end

function local_row3v(da)
    gx, gy, gz, gm, gn, _ = Int.(LibPETSc.DMDAGetGhostCorners(P, da))
    return (i, j, k, c) -> 3((i - gx) + (j - gy) * gm + (k - gz) * gm * gn) + c
end

function assemble_gradient3!(s)
    G, row = s.pressure_gradient_matrix, local_row3v(s.da_Veloc)
    rows, col = zeros(PInt, 24), zeros(PInt, 1)
    xs, ys, zs = element_ranges3(s.m, s.da_Veloc)
    for ek in zs, ej in ys, ei in xs
        for (n, (i, j, k)) in enumerate(element_nodes3(ei, ej, ek)), c in 0:2; rows[3(n-1)+c+1] = row(i, j, k, c); end
        col[1] = row(ei, ej, ek, 0)
        LibPETSc.MatSetValuesLocal(P, G, PInt(24), rows, PInt(1), col, s.element.VCe, LibPETSc.ADD_VALUES)
    end
    finish_matrix(G)
end

"""
FSSA term in 3D: theta*dt*d(dRho)/dz*w with nodal integration (dx*dy*dz/8 per node). In the C++ code the term is
overwritten before it is used, so the 3D FSSA has no effect (bug 29): with fssa_legacy = True there is no term.
"""
function fssa_term3(m, d, dt)
    m.fssa_legacy && return 0.0
    jump = (d[5] + d[6] + d[7] + d[8] - d[1] - d[2] - d[3] - d[4]) / 4
    return m.fssa_theta * dt * jump * element_size_x(m) * element_size_y(m) / 8
end

function assemble_velocity3!(s, velocity, pressure, viscosity_factor, density_variation, dt; fssa_density_variation = density_variation)
    m, A, el = s.m, s.velocity_stiffness_matrix, s.element
    LibPETSc.MatZeroEntries(P, A)
    global_to_local(s.da_Veloc, velocity, s.local_velocity)
    global_to_local(s.da_Veloc, pressure, s.local_pressure)
    global_to_local(s.da_Thermal, viscosity_factor, s.local_viscosity_factor)
    global_to_local(s.da_Thermal, density_variation, s.local_density_variation)
    global_to_local(s.da_Thermal, fssa_density_variation, s.local_fssa_density_variation)
    LibPETSc.VecZeroEntries(P, s.local_pressure_rhs); LibPETSc.VecZeroEntries(P, s.local_uzawa_preconditioner)
    row = local_row3v(s.da_Veloc)
    PETSc.with_local_array!(s.local_pressure_rhs, s.local_uzawa_preconditioner, s.local_velocity, s.local_pressure,
        s.local_velocity_bc_mask, s.local_viscosity_factor, s.local_density_variation, s.local_fssa_density_variation;
        write = (true, true, false, false, false, false, false, false)) do f, pc, v, p, c, q, d, d0
        l3(a) = local_view3(s.da_Veloc, a, 3); l1(a) = local_view3(s.da_Thermal, a)
        F, PC, V, Pr, C, visc, D, D0 = l3(f), l3(pc), l3(v), l3(p), l3(c), l1(q), l1(d), l1(d0)
        sx, sy, sz, mx, my, mz = Int.(LibPETSc.DMDAGetCorners(P, s.da_Veloc))
        r, one_value = PInt[0], [0.0]
        for k in sz:sz+mz-1, j in sy:sy+my-1, i in sx:sx+mx-1, n in 0:2
            r[1] = row(i, j, k, n); one_value[1] = 1.0 - C[n+1, i, j, k]
            LibPETSc.MatSetValuesLocal(P, A, PInt(1), r, PInt(1), r, one_value, LibPETSc.ADD_VALUES)
        end
        Ke, values, rows, free, fe, de = zeros(24, 24), zeros(576), zeros(PInt, 24), zeros(Bool, 24), zeros(24), zeros(8)
        xs, ys, zs = element_ranges3(m, s.da_Veloc)
        for ek in zs, ej in ys, ei in xs
            nodes = element_nodes3(ei, ej, ek)
            Ke .= visc[ei, ej, ek] .* el.Ke_veloc
            for n in 1:8; de[n] = D[nodes[n]...]; end
            if m.free_surface_stabilization
                x = fssa_term3(m, ntuple(n -> D0[nodes[n]...], 8), dt)
                for n in 1:8; Ke[3n, 3n] += x; end
            end
            for n in 1:8, c in 1:3
                i, j, k = nodes[n]
                rows[3(n-1)+c] = row(i, j, k, c - 1); free[3(n-1)+c] = C[c, i, j, k] != 0
            end
            for a in 1:24, b in 1:24; values[24(a-1)+b] = free[a] ? Ke[a, b] : 0.0; end
            LibPETSc.MatSetValuesLocal(P, A, PInt(24), rows, PInt(24), rows, values, LibPETSc.ADD_VALUES)
            for a in 1:24
                PC[1, ei, ej, ek] += el.VCe[a]^2 / Ke[a, a]
                fe[a] = -el.VCe[a] * Pr[1, ei, ej, ek]
                for b in 1:8; fe[a] += el.VfMe[a, b] * de[b]; end
            end
            for n in 1:8, c in 1:3
                free[3(n-1)+c] && (F[c, nodes[n]...] += fe[3(n-1)+c])
            end
        end
        for k in sz:sz+mz-1, j in sy:sy+my-1, i in sx:sx+mx-1, n in 1:3
            C[n, i, j, k] == 0 && (F[n, i, j, k] = V[n, i, j, k])
        end
    end
    finish_matrix(A)
    for (l, g) in ((s.local_pressure_rhs, s.velocity_rhs_pressure), (s.local_uzawa_preconditioner, s.uzawa_preconditioner))
        LibPETSc.VecZeroEntries(P, g)
        LibPETSc.DMLocalToGlobal(P, s.da_Veloc, l, LibPETSc.ADD_VALUES, g)
    end
    LibPETSc.VecReciprocal(P, s.uzawa_preconditioner)
end

"Solves Stokes in 3D: velocity_next and pressure are updated in place (veloc_total, linear case)."
function stokes_solve!(s::StokesSolver3D, velocity_next, pressure, velocity, viscosity_factor, density_variation, dt;
    max_iterations = 400, fssa_density_variation = density_variation)
    assemble_velocity3!(s, velocity, pressure, viscosity_factor, density_variation, dt; fssa_density_variation)
    # 64 = 8^2: G is 8 times larger in Mandyoc (no 1/8 Jacobian), so the tolerance has the same meaning
    uzawa!(s, velocity_next, pressure, 64.0, max_iterations)
end
