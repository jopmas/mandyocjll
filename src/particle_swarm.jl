using PETSc, MPI, OffsetArrays
using PETSc: LibPETSc
include("parameters.jl")

const petsclib = first(l for l in PETSc.petsclibs if PETSc.scalartype(l) === Float64)
const PInt = PETSc.inttype(petsclib)          # PetscInt: 32 or 64 bits, depending on PETSc
const PIC_COORDINATES = "DMSwarmPIC_coor"     # DMSwarmPICField_coor
const SHIM_LIBRARY = joinpath(@__DIR__, "..", "csrc",
    "libmandyoc_shim.$(Sys.isapple() ? "dylib" : "so")")

# ---- velocity mesh (create_veloc) ------------------------------------
function create_mesh(m::ModelParameters, dof, width)
    bx = m.periodic_boundary ? LibPETSc.DM_BOUNDARY_PERIODIC : LibPETSc.DM_BOUNDARY_NONE
    da = LibPETSc.DMDACreate2d(petsclib, MPI.COMM_WORLD, bx, LibPETSc.DM_BOUNDARY_NONE,
        LibPETSc.DMDA_STENCIL_BOX, PInt.((m.Nx, m.Nz, m.Px, m.Pz, dof, width))..., C_NULL, C_NULL)
    LibPETSc.DMSetUp(petsclib, da)
    LibPETSc.DMDASetUniformCoordinates(petsclib, da, 0.0, m.Lx, -m.domain_depth, 0.0, 0.0, 0.0)
    return da
end
create_velocity_mesh(m) = create_mesh(m, 2, stencil_width(m))       # da_Veloc
create_thermal_mesh(m) = create_mesh(m, 1, stencil_width(m))        # da_Thermal
# Same partition as da_Veloc, with ghosts for the RK4 stages in the neighboring element.
function create_advection_mesh(m, da_Veloc)
    da = create_mesh(m, 2, m.advection_stencil_width)
    @assert LibPETSc.DMDAGetCorners(petsclib, da) == LibPETSc.DMDAGetCorners(petsclib, da_Veloc)
    return da
end

# ---- callbacks of cell_dm (DMLocatePoints_DMDARegular_2d etc.) --------------
struct SFNode
    rank::PInt
    index::PInt
end
mutable struct CellContext
    da_ptr::Ptr{Cvoid}
    si::Int; sk::Int; milocal::Int; mklocal::Int
    m::ModelParameters
    sj::Int; mjlocal::Int                # 3D (DMLocatePoints_DMDARegular_3d)
end
CellContext(da_ptr, si, sk, milocal, mklocal, m) = CellContext(da_ptr, si, sk, milocal, mklocal, m, 0, 1)
const CELL_CONTEXT = Ref{CellContext}()

const POINT_NOT_FOUND = -367   # DMLOCATEPOINT_POINT_NOT_FOUND
# rule of Mandyoc: the rank that owns the lower left node owns the element
function owned_cell(c::CellContext, x, z)
    m = c.m
    (x < 0 || x > m.Lx || z < -m.domain_depth || z > 0) && return POINT_NOT_FOUND
    mi = trunc(Int, x / element_size_x(m)) - c.si
    mk = trunc(Int, (z + m.domain_depth) / element_size_z(m)) - c.sk
    (0 <= mi < c.milocal && 0 <= mk < c.mklocal) || return POINT_NOT_FOUND
    return mi + mk * c.milocal
end
function owned_cell(c::CellContext, x, y, z)
    m = c.m
    (x < 0 || x > m.Lx || y < 0 || y > m.Ly || z < -m.domain_depth || z > 0) && return POINT_NOT_FOUND
    mi = trunc(Int, x / element_size_x(m)) - c.si
    mj = trunc(Int, y / element_size_y(m)) - c.sj
    mk = trunc(Int, (z + m.domain_depth) / element_size_z(m)) - c.sk
    (0 <= mi < c.milocal && 0 <= mj < c.mjlocal && 0 <= mk < c.mklocal) || return POINT_NOT_FOUND
    return mi + (mj + mk * c.mjlocal) * c.milocal
end

const petsc_path = petsclib.petsc_library
function locate_points(dm::Ptr{Cvoid}, pos::Ptr{Cvoid}, ltype::Cint, cell_sf::Ptr{Cvoid})::Cint
    n = Ref{PInt}(0); a = Ref{Ptr{Float64}}(C_NULL)
    e = @ccall $(dlsym_petsc(:VecGetLocalSize))(pos::Ptr{Cvoid}, n::Ptr{PInt})::Cint
    e == 0 || return e
    e = @ccall $(dlsym_petsc(:VecGetArrayRead))(pos::Ptr{Cvoid}, a::Ptr{Ptr{Float64}})::Cint
    e == 0 || return e
    c = CELL_CONTEXT[]; dim = c.m.dimensions
    npoints = n[] ÷ dim
    coor = unsafe_wrap(Array, a[], n[])
    cells = dim == 3 ? [SFNode(0, owned_cell(c, coor[3p-2], coor[3p-1], coor[3p])) for p in 1:npoints] :
        [SFNode(0, owned_cell(c, coor[2p-1], coor[2p])) for p in 1:npoints]
    e = @ccall $(dlsym_petsc(:VecRestoreArrayRead))(pos::Ptr{Cvoid}, a::Ptr{Ptr{Float64}})::Cint
    e == 0 || return e
    # 0 = PETSC_COPY_VALUES: the array belongs to Julia, and PETSc copies it
    GC.@preserve cells @ccall $(dlsym_petsc(:PetscSFSetGraph))(cell_sf::Ptr{Cvoid}, npoints::PInt,
        npoints::PInt, C_NULL::Ptr{PInt}, 0::Cint, cells::Ptr{SFNode}, 0::Cint)::Cint
end
function get_neighbors(dm::Ptr{Cvoid}, n::Ptr{PInt}, ranks::Ptr{Ptr{Cint}})::Cint
    @ccall $(dlsym_petsc(:DMGetNeighbors))(CELL_CONTEXT[].da_ptr::Ptr{Cvoid}, n::Ptr{PInt}, ranks::Ptr{Ptr{Cint}})::Cint
end
using Libdl
const PETSC_HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
function dlsym_petsc(s)
    PETSC_HANDLE[] == C_NULL && (PETSC_HANDLE[] = Libdl.dlopen(petsc_path))
    Libdl.dlsym(PETSC_HANDLE[], s)
end

# ---- swarm (createSwarm_2d) --------------------------------------------------
const SWARM_FIELDS = (("itag", LibPETSc.PETSC_INT), ("layer", LibPETSc.PETSC_INT),
    ("geoq_fac", LibPETSc.PETSC_DOUBLE), ("strain_fac", LibPETSc.PETSC_DOUBLE),
    ("strain_rate_fac", LibPETSc.PETSC_DOUBLE), ("cont", LibPETSc.PETSC_INT))
# registered only with the magmatism, as in the C++ code: depletion (X), melt fraction (Phi) and its last increment (dPhi)
const MELT_FIELDS = ("X", "Phi", "dPhi")
melt_field_names(m) = m.use_magmatism ? MELT_FIELDS : ()

"Calls `f` with the field as a Julia Array and returns it to PETSc, even if `f` throws."
function with_field(f, swarm, name, ::Type{T}) where {T}
    n = LibPETSc.DMSwarmGetLocalSize(petsclib, swarm)
    bs, type, ptr = LibPETSc.DMSwarmGetField(petsclib, swarm, name)
    try
        return f(unsafe_wrap(Array, Ptr{T}(ptr), bs * n))
    finally
        LibPETSc.DMSwarmRestoreField(petsclib, swarm, name, bs, type)
    end
end

# pseudorandom number that depends only on the particle: the initial position does not
# change with the number of ranks (rand_r in Mandyoc depends on the rank)
unit_random(seed, id, k) = (hash((seed, id, k)) >> 11) * 2.0^-53 - 0.5

function create_particle_swarm(m::ModelParameters, da_Veloc; location = :shim, init = (x, z, id) -> (id % 7, 1.5id))
    si, sk, _, milocal, mklocal, _ = LibPETSc.DMDAGetCorners(petsclib, da_Veloc)
    CELL_CONTEXT[] = CellContext(da_Veloc.ptr, si, sk, milocal, mklocal, m)
    cell_dm = da_Veloc                                              # dmcell
    if location == :shim
        cell_dm = LibPETSc.DMShellCreate(petsclib, MPI.COMM_WORLD)  # dmcell
        LibPETSc.DMSetDimension(petsclib, cell_dm, PInt(2))
        shim = Libdl.dlsym(Libdl.dlopen(SHIM_LIBRARY), :mandyoc_set_cell_dm_ops)
        loc = @cfunction(locate_points, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cvoid}))
        nei = @cfunction(get_neighbors, Cint, (Ptr{Cvoid}, Ptr{PInt}, Ptr{Ptr{Cint}}))
        @ccall $shim(cell_dm.ptr::Ptr{Cvoid}, loc::Ptr{Cvoid}, nei::Ptr{Cvoid})::Cint
    end
    particle_swarm = LibPETSc.DMCreate(petsclib, MPI.COMM_WORLD)    # dms
    LibPETSc.DMSetType(petsclib, particle_swarm, "swarm")
    LibPETSc.DMSetDimension(petsclib, particle_swarm, PInt(2))
    LibPETSc.DMSwarmSetType(petsclib, particle_swarm, LibPETSc.DMSWARM_PIC)
    LibPETSc.DMSwarmSetCellDM(petsclib, particle_swarm, cell_dm)
    for (name, type) in SWARM_FIELDS
        LibPETSc.DMSwarmRegisterPetscDatatypeField(petsclib, particle_swarm, name, PInt(1), type)
    end
    for name in melt_field_names(m)
        LibPETSc.DMSwarmRegisterPetscDatatypeField(petsclib, particle_swarm, name, PInt(1), LibPETSc.PETSC_DOUBLE)
    end
    LibPETSc.DMSwarmFinalizeFieldRegister(petsclib, particle_swarm)

    dx, dz = element_size_x(m), element_size_z(m)                   # dx_const, dz_const
    ppe = particles_per_element(m)                                  # particles_per_ele
    cells = [(i, k) for k in sk:min(sk + mklocal, m.Nz - 1)-1 for i in si:min(si + milocal, m.Nx - 1)-1]
    LibPETSc.DMSwarmSetLocalSizes(petsclib, particle_swarm, PInt(length(cells) * ppe), PInt(4))
    with_field(particle_swarm, PIC_COORDINATES, Float64) do coor    # array
    with_field(particle_swarm, "itag", PInt) do itag               # iarray
    with_field(particle_swarm, "layer", PInt) do layer             # layer_array
    with_field(particle_swarm, "geoq_fac", Float64) do viscosity_factor  # rarray (geoq)
    with_field(particle_swarm, "strain_fac", Float64) do strain
        p = 0
        for (i, k) in cells, cx in 0:m.particles_per_element_x-1, cz in 0:m.particles_per_element_z-1
            p += 1
            id = (k * (m.Nx - 1) + i) * ppe + cx * m.particles_per_element_z + cz + 1
            rx = unit_random(m.random_seed, id, 1); rz = unit_random(m.random_seed, id, 2)
            coor[2p-1] = i * dx + (0.5 + cx + rx * m.particles_perturb_factor) * dx / m.particles_per_element_x
            coor[2p] = -m.domain_depth + k * dz + (0.5 + cz + rz * m.particles_perturb_factor) * dz / m.particles_per_element_z
            itag[p] = id; strain[p] = 0.0
            layer[p], viscosity_factor[p] = init(coor[2p-1], coor[2p], id)
        end
    end end end end end
    for (name, value) in zip(melt_field_names(m), (1.0, 0.0, 0.0))
        with_field(a -> (a .= value; nothing), particle_swarm, name, Float64)
    end
    LibPETSc.DMSwarmMigrate(petsclib, particle_swarm, LibPETSc.PETSC_TRUE)
    return particle_swarm, cell_dm
end

# ---- advection (moveSwarm, kinematic part) ----------------------------------
linear_interpolation(rx, rz, v0, v1, v2, v3) =
    (1 - rx) * (1 - rz) * v0 + rx * (1 - rz) * v1 + (1 - rx) * rz * v2 + rx * rz * v3

function interpolate_velocity(m, V, x, z)                           # V: VV
    dx, dz = element_size_x(m), element_size_z(m)
    i = clamp(trunc(Int, x / dx), 0, m.Nx - 2)
    k = clamp(trunc(Int, (z + m.domain_depth) / dz), 0, m.Nz - 2)
    rx = x / dx - i; rz = (z + m.domain_depth) / dz - k
    if m.periodic_boundary                  # use the periodic image available on this rank
        i < firstindex(V, 2) && (i += m.Nx)
        i + 1 > lastindex(V, 2) && (i -= m.Nx)
    end
    u = linear_interpolation(rx, rz, V[1, i, k], V[1, i+1, k], V[1, i, k+1], V[1, i+1, k+1])
    w = linear_interpolation(rx, rz, V[2, i, k], V[2, i+1, k], V[2, i, k+1], V[2, i+1, k+1])
    return u, w
end

# Mandyoc does a `break` here in RK4, which skips the remaining particles of the rank; here only
# the position of the stage is wrapped (periodic) or stopped at the boundary (closed)
function keep_inside(m, x, z)
    e = m.boundary_offset                                           # epsilon_x
    x = m.periodic_boundary ? mod(x, m.Lx) : clamp(x, e, m.Lx - e)
    return x, clamp(z, e - m.domain_depth, -e)
end

"Velocity at the time t_n + f*dt of the step: linear between V0 (t_n, Veloc) and V1 (t_n+1, Veloc_fut)."
function velocity_in_time(m, V0, V1, x, z, f)
    u0, w0 = interpolate_velocity(m, V0, x, z)
    V1 === V0 && return u0, w0
    u1, w1 = interpolate_velocity(m, V1, x, z)
    return (1 - f) * u0 + f * u1, (1 - f) * w0 + f * w1
end

"""
Mean velocity of the step by the fourth-order Runge-Kutta (x_n+1 = x_n + dt*v), with the velocity varying linearly
in time from V0 (t_n) to V1 (t_n+1); V1 = V0 for a velocity constant in the step. V0 and V1 need ghosts for the
stages in the neighbor element (da_advection, stencil 2).
"""
function runge_kutta_4_velocity(m, V0, V1, x, z, dt)
    uA, wA = velocity_in_time(m, V0, V1, x, z, 0.0)
    uB, wB = velocity_in_time(m, V0, V1, keep_inside(m, x + uA * dt / 2, z + wA * dt / 2)..., 0.5)
    uC, wC = velocity_in_time(m, V0, V1, keep_inside(m, x + uB * dt / 2, z + wB * dt / 2)..., 0.5)
    uD, wD = velocity_in_time(m, V0, V1, keep_inside(m, x + uC * dt, z + wC * dt)..., 1.0)
    return (uA + 2uB + 2uC + uD) / 6, (wA + 2wB + 2wC + wD) / 6
end

function advect_coordinates!(m, coor, V, dt)
    for p in 1:length(coor)÷2
        x, z = coor[2p-1], coor[2p]
        uA, wA = m.use_runge_kutta_4 ? runge_kutta_4_velocity(m, V, V, x, z, dt) : interpolate_velocity(m, V, x, z)
        # final position: particles crossing an open boundary leave, and DMSwarmMigrate removes them
        x += dt * uA
        coor[2p-1] = m.periodic_boundary ? mod(x, m.Lx) : x
        coor[2p] = z + dt * wA
    end
end

function advect_particles!(m, particle_swarm, da_advection, local_velocity, dt)  # dms, local_V
    gx, gz, _, gm, gn, _ = LibPETSc.DMDAGetGhostCorners(petsclib, da_advection)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    with_field(particle_swarm, "cont", PInt) do cont
        fill!(cont, rank)                                           # rank before the migration
    end
    PETSc.with_local_array!(local_velocity; write = false) do a
        V = OffsetArray(reshape(a, 2, Int(gm), Int(gn)), 1:2, gx:gx+gm-1, gz:gz+gn-1)
        with_field(coor -> advect_coordinates!(m, coor, V, dt), particle_swarm, PIC_COORDINATES, Float64)
    end
    LibPETSc.DMSwarmMigrate(petsclib, particle_swarm, LibPETSc.PETSC_TRUE)
    moved = with_field(c -> count(!=(rank), c), particle_swarm, "cont", PInt)
    return MPI.Allreduce(moved, +, MPI.COMM_WORLD)
end
