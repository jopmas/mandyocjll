# Output in zarr (format v2, conventions of xarray), with the same names, dimensions, units and types as
# the .nc files of tapIOca (meshes2netcdf.jl and particles2netcdf.jl): <variable>.zarr instead of <variable>.nc.
# Each step is a new chunk, and the variables are split among the ranks, which write in parallel.
# Include after rheology.jl.
using Zlib_jll

# ---- minimal zarr v2 -------------------------------------------------------------
zarr_type(::Type{Float32}) = "<f4"
zarr_type(::Type{Float64}) = "<f8"
zarr_type(::Type{Int8}) = "|i1"
zarr_type(::Type{Int32}) = "<i4"
json(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"") * "\""
json(x::Real) = string(x)
json(x::Nothing) = "null"
json(x::Union{Tuple, AbstractVector}) = "[" * join(json.(x), ", ") * "]"
json(x::AbstractDict) = "{" * join((json(string(k)) * ": " * json(x[k]) for k in sort(collect(keys(x)))), ", ") * "}"

function write_atomic(path, data)        # a reader never sees a partially written file during the simulation
    mkpath(dirname(path)); tmp = path * ".tmp"
    write(tmp, data); mv(tmp, path; force = true)
end

function zlib_compress(data::Vector{UInt8}, level)
    out = Vector{UInt8}(undef, ccall((:compressBound, libz), Culong, (Culong,), length(data)))
    n = Ref{Culong}(length(out))
    r = ccall((:compress2, libz), Cint, (Ptr{UInt8}, Ptr{Culong}, Ptr{UInt8}, Culong, Cint), out, n, data, length(data), level)
    r == 0 || error("compress2: $r")
    return resize!(out, n[])
end

struct ZarrArray
    name::String
    type::DataType
    dims::Vector{String}                 # order of xarray (C), time first
    chunks::Vector{Int}
    attributes::Dict{String, Any}
end
zarray_json(a::ZarrArray, shape, level) = json(Dict("zarr_format" => 2, "shape" => shape, "chunks" => a.chunks, "dtype" => zarr_type(a.type),
    "compressor" => Dict("id" => "zlib", "level" => level), "fill_value" => nothing,     # xarray masks a fill_value of zarr v2
    "order" => "C", "filters" => nothing))
zattrs_json(a::ZarrArray) = json(merge(a.attributes, Dict{String, Any}("_ARRAY_DIMENSIONS" => a.dims)))

"Writes a chunk; in C order the last index varies fastest, so `data` has x first (order of Julia)."
write_chunk(store, a::ZarrArray, index, data, level) = write_atomic(joinpath(store, a.name, join(index, ".")),
    zlib_compress(collect(reinterpret(UInt8, vec(convert(Array{a.type}, data)))), level))

"Metadata of the store (.zgroup, .zarray, .zattrs and the consolidated .zmetadata)."
function write_metadata(store, arrays, shapes, level, attributes = Dict{String, Any}())
    meta = Dict(".zgroup" => "{\"zarr_format\": 2}", ".zattrs" => json(attributes))
    for a in arrays
        meta["$(a.name)/.zarray"] = zarray_json(a, shapes[a.name], level)
        meta["$(a.name)/.zattrs"] = zattrs_json(a)
    end
    for (k, v) in meta; write_atomic(joinpath(store, k), v); end
    write_atomic(joinpath(store, ".zmetadata"),
        "{\"zarr_consolidated_format\": 1, \"metadata\": {" * join((json(k) * ": " * meta[k] for k in sort(collect(keys(meta)))), ", ") * "}}")
end

# ---- variables of Mandyoc --------------------------------------------------------
# name => (type, unit), as in meshes2netcdf.jl
const OUTPUT_VARIABLES = ["density" => (Float32, "kg/m3"), "viscosity" => (Float64, "Pa.s"), "pressure" => (Float32, "Pa"),
    "strain" => (Float64, "dimensionless"), "strain_rate" => (Float64, "1/s"), "temperature" => (Float32, "deg C"),
    "velocity" => (Float64, "m/s"), "heat" => (Float64, "W/kg"), "thermal_diffusivity" => (Float64, "m2/s"),
    "Phi" => (Float64, "dimensionless"), "dPhi" => (Float64, "dimensionless"), "X_depletion" => (Float64, "dimensionless")]
const MELT_OUTPUTS = ("Phi", "dPhi", "X_depletion")    # written with the magmatism
const TIME_CHUNK = 1024                  # time and step are in a single chunk, rewritten at each output

mutable struct OutputWriter
    path::String
    level::Int                           # zlib level (0 to 9)
    variables::Vector{String}
    export_lithology::Bool               # export_lithology
    write_particles::Bool                # print_step_files
    times::Vector{Float64}               # Myr
    steps::Vector{Int}
    seconds::Float64                     # time spent writing
end
function OutputWriter(path; level = 1, export_thermal_diffusivity = false, export_lithology = false, write_particles = false,
    magmatism = false)
    names = [k for (k, _) in OUTPUT_VARIABLES if (k != "thermal_diffusivity" || export_thermal_diffusivity) && (magmatism || !(k in MELT_OUTPUTS))]
    return OutputWriter(path, level, names, export_lithology, write_particles, Float64[], Int[], 0.0)
end

"Stops if the folder already has the output of a previous run."
function check_output_folder(path)
    stale = [joinpath(path, "$v.zarr") for v in vcat(first.(OUTPUT_VARIABLES), ["lithology", "particles_trajectories", "surface"])]
    any(isdir, stale) && error("there are already zarr outputs in $path; move or delete them before running again")
end

"Gathers a field of the DMDA on rank `root` (natural order, x first)."
function gather_field(vec, da, dof, Nx, Nz, root)
    comm = MPI.COMM_WORLD
    corners = MPI.Gather(Int32[Int32(c) for c in LibPETSc.DMDAGetCorners(P, da)[[1, 2, 4, 5]]], comm; root)
    values = PETSc.with_local_array!(copy, vec; write = false)
    MPI.Comm_rank(comm) == root || (MPI.Gatherv!(values, nothing, comm; root); return nothing)
    c = reshape(corners, 4, :); counts = [Int32(dof * c[3, r] * c[4, r]) for r in axes(c, 2)]
    buffer = zeros(sum(counts)); MPI.Gatherv!(values, VBuffer(buffer, counts), comm; root)
    out = zeros(dof, Nx, Nz); o = 0
    for r in axes(c, 2)
        n = counts[r]; block = reshape(view(buffer, o+1:o+n), dof, Int(c[3, r]), Int(c[4, r])); o += n
        out[:, c[1, r]+1:c[1, r]+c[3, r], c[2, r]+1:c[2, r]+c[4, r]] .= block
    end
    return out
end

time_arrays() = [ZarrArray("time", Float32, ["time"], [TIME_CHUNK], Dict{String, Any}("units" => "Myr", "long_name" => "time", "axis" => "T")),
    ZarrArray("step", Int32, ["time"], [TIME_CHUNK], Dict{String, Any}("units" => "", "long_name" => "step"))]
axis_array(name, n) = ZarrArray(name, Float32, [name], [n], Dict{String, Any}("units" => "m", "long_name" => name, "axis" => uppercase(name)))

function write_times(w, store)
    t, s = time_arrays()
    for c in 0:(length(w.times)-1)÷TIME_CHUNK
        r = c*TIME_CHUNK+1:min((c + 1) * TIME_CHUNK, length(w.times))
        pad(v) = vcat(v[r], zeros(eltype(v), TIME_CHUNK - length(r)))
        write_chunk(store, t, (c,), pad(w.times), w.level); write_chunk(store, s, (c,), pad(w.steps), w.level)
    end
end

# regular mesh: one chunk per step; `fields` = [(name, type, unit, matrix [x, z])]
function write_mesh_store(w, variable, Nx, Nz, Lx, Lz, fields)
    store = joinpath(w.path, "$variable.zarr"); nt = length(w.times)
    x, z = axis_array("x", Nx), axis_array("z", Nz)
    arrays = vcat(time_arrays(), [x, z], [ZarrArray(n, T, ["time", "z", "x"], [1, Nz, Nx],
        Dict{String, Any}("long_name" => n, "units" => u)) for (n, T, u, _) in fields])
    if nt == 1
        write_chunk(store, x, (0,), range(0, Lx, length = Nx), w.level); write_chunk(store, z, (0,), range(-Lz, 0, length = Nz), w.level)
    end
    for (a, (_, _, _, data)) in zip(arrays[5:end], fields); write_chunk(store, a, (nt - 1, 0, 0), data, w.level); end
    write_times(w, store)
    shapes = Dict{String, Vector{Int}}("time" => [nt], "step" => [nt], "x" => [Nx], "z" => [Nz])
    for (n, _, _, _) in fields; shapes[n] = [nt, Nz, Nx]; end
    write_metadata(store, arrays, shapes, w.level)
end

"Gathers a field of the 3D DMDA on rank `root` (natural order, x first): [dof, x, y, z]."
function gather_field3(vec, da, dof, Nx, Ny, Nz, root)
    comm = MPI.COMM_WORLD
    corners = MPI.Gather(Int32[Int32(c) for c in LibPETSc.DMDAGetCorners(P, da)], comm; root)
    values = PETSc.with_local_array!(copy, vec; write = false)
    MPI.Comm_rank(comm) == root || (MPI.Gatherv!(values, nothing, comm; root); return nothing)
    c = reshape(corners, 6, :); counts = [Int32(dof * c[4, r] * c[5, r] * c[6, r]) for r in axes(c, 2)]
    buffer = zeros(sum(counts)); MPI.Gatherv!(values, VBuffer(buffer, counts), comm; root)
    out = zeros(dof, Nx, Ny, Nz); o = 0
    for r in axes(c, 2)
        n = counts[r]; block = reshape(view(buffer, o+1:o+n), dof, Int(c[4, r]), Int(c[5, r]), Int(c[6, r])); o += n
        out[:, c[1, r]+1:c[1, r]+c[4, r], c[2, r]+1:c[2, r]+c[5, r], c[3, r]+1:c[3, r]+c[6, r]] .= block
    end
    return out
end

# regular 3D mesh: one chunk per step, dims (time, z, y, x); `fields` = [(name, type, unit, array [x, y, z])]
function write_mesh_store3(w, variable, m, fields)
    store = joinpath(w.path, "$variable.zarr"); nt = length(w.times)
    x, y, z = axis_array("x", m.Nx), axis_array("y", m.Ny), axis_array("z", m.Nz)
    arrays = vcat(time_arrays(), [x, y, z], [ZarrArray(n, T, ["time", "z", "y", "x"], [1, m.Nz, m.Ny, m.Nx],
        Dict{String, Any}("long_name" => n, "units" => u)) for (n, T, u, _) in fields])
    if nt == 1
        write_chunk(store, x, (0,), range(0, m.Lx, length = m.Nx), w.level); write_chunk(store, y, (0,), range(0, m.Ly, length = m.Ny), w.level)
        write_chunk(store, z, (0,), range(-m.domain_depth, 0, length = m.Nz), w.level)
    end
    for (a, (_, _, _, data)) in zip(arrays[6:end], fields); write_chunk(store, a, (nt - 1, 0, 0, 0), data, w.level); end
    write_times(w, store)
    shapes = Dict{String, Vector{Int}}("time" => [nt], "step" => [nt], "x" => [m.Nx], "y" => [m.Ny], "z" => [m.Nz])
    for (n, _, _, _) in fields; shapes[n] = [nt, m.Nz, m.Ny, m.Nx]; end
    write_metadata(store, arrays, shapes, w.level)
end

"Lithology on a 5 times finer mesh, as tapIOca does with the litho_* files."
function lithology_mesh(m, particle_swarm, root)
    n, k = 5(m.Nx - 1) + 1, 5(m.Nz - 1) + 1
    local_mesh = fill(Int8(-1), n, k)
    with_field(particle_swarm, PIC_COORDINATES, Float64) do c
        with_field(particle_swarm, "layer", PInt) do layer
            for p in eachindex(layer)
                i = trunc(Int, c[2p-1] * 5 / element_size_x(m)); j = trunc(Int, -c[2p] * 5 / element_size_z(m))
                0 <= i < n && 0 <= j < k && (local_mesh[i+1, k-j] = layer[p])
            end
        end
    end
    mesh = MPI.Reduce(local_mesh, max, MPI.COMM_WORLD; root)
    mesh === nothing && return nothing
    while any(<(0), mesh)                # cells without particles: copy from a neighbor
        old = copy(mesh); changed = false
        for j in 1:k, i in 1:n
            old[i, j] < 0 || continue
            for (a, b) in ((i, j + 1), (i, j - 1), (i - 1, j), (i + 1, j))
                if 1 <= a <= n && 1 <= b <= k && old[a, b] >= 0
                    mesh[i, j] = old[a, b]; changed = true
                    break
                end
            end
        end
        changed || break
    end
    return mesh
end

"Trajectories of the first particle of each element (those the C++ code writes to step_*)."
function write_trajectories(w, m, particle_swarm, root)
    comm = MPI.COMM_WORLD; ppe = particles_per_element(m)
    tag = with_field(copy, particle_swarm, "itag", PInt); coor = with_field(copy, particle_swarm, PIC_COORDINATES, Float64)
    layer = with_field(copy, particle_swarm, "layer", PInt); strain = with_field(copy, particle_swarm, "strain_fac", Float64)
    keep = findall(t -> t > 0 && (t - 1) % ppe == 0, tag)
    local_values = zeros(5, length(keep))
    for (n, p) in enumerate(keep)
        local_values[:, n] .= ((tag[p] - 1) ÷ ppe, coor[2p-1], coor[2p], layer[p], strain[p])
    end
    counts = MPI.Gather(Int32(length(local_values)), comm; root)
    MPI.Comm_rank(comm) == root || (MPI.Gatherv!(local_values, nothing, comm; root); return)
    all_values = zeros(5, sum(counts) ÷ 5); MPI.Gatherv!(local_values, VBuffer(all_values, counts), comm; root)
    ne = (m.Nx - 1) * (m.Nz - 1); nt = length(w.times); store = joinpath(w.path, "particles_trajectories.zarr")
    x, z, l, s = fill(-999.0, ne), fill(-999.0, ne), fill(-1, ne), fill(-999.0, ne)   # removed particles, as in particles2netcdf.jl
    for c in eachcol(all_values); i = Int(c[1]) + 1; x[i] = c[2]; z[i] = c[3]; l[i] = Int(c[4]); s[i] = c[5]; end
    variable(n, T, a) = ZarrArray(n, T, ["time", "id"], [1, ne], Dict{String, Any}("long_name" => n, a...))
    id = ZarrArray("id", Int32, ["id"], [ne], Dict{String, Any}("long_name" => "particle_id", "cf_role" => "trajectory_id"))
    arrays = vcat(time_arrays(), [id, variable("x", Float32, ("units" => "m", "axis" => "X")), variable("z", Float32, ("units" => "m", "axis" => "Z")),
        variable("lithology", Int8, ("coordinates" => "time z x",)), variable("strain", Float32, ("units" => "dimensionless",))])
    nt == 1 && write_chunk(store, id, (0,), 10000 .+ (0:ne-1), w.level)
    for (a, data) in zip(arrays[4:7], (x, z, l, s)); write_chunk(store, a, (nt - 1, 0), data, w.level); end
    write_times(w, store)
    shapes = Dict{String, Vector{Int}}("time" => [nt], "step" => [nt], "id" => [ne])
    for a in arrays[4:7]; shapes[a.name] = [nt, ne]; end
    write_metadata(store, arrays, shapes, w.level, Dict{String, Any}("Conventions" => "CF-1.8", "featureType" => "trajectory",
        "description" => "particle trajectories", "reference_timestep" => "0"))
end

"Surface (surface.zarr): the height of the regular points, as surface.nc of tapIOca (x from 0 to Lx)."
function write_surface_store(w, surface, Lx)
    store = joinpath(w.path, "surface.zarr"); nt = length(w.times); n = length(surface.z)
    x = axis_array("x", n)
    height = ZarrArray("surface", Float32, ["time", "x"], [1, n], Dict{String, Any}("long_name" => "surface", "units" => "m"))
    arrays = vcat(time_arrays(), [x, height])
    nt == 1 && write_chunk(store, x, (0,), range(0, Lx, length = n), w.level)
    write_chunk(store, height, (nt - 1, 0), surface.z, w.level)
    write_times(w, store)
    write_metadata(store, arrays, Dict{String, Vector{Int}}("time" => [nt], "step" => [nt], "x" => [n], "surface" => [nt, n]), w.level)
end

"Appends the current state of the model to the zarr stores."
function write_output!(w::OutputWriter, s::Model, step; surface = nothing)
    start = MPI.Wtime(); comm = MPI.COMM_WORLD; m = s.m; rank, np = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    push!(w.times, s.time / SECONDS_PER_YEAR / 1e6); push!(w.steps, step)
    MPI.Barrier(comm)
    f = s.fields
    sources = Dict("density" => (f.density_field, s.da_Thermal, 1), "viscosity" => (f.viscosity_factor, s.da_Thermal, 1),
        "pressure" => (s.pressure_thermal_mesh, s.da_Thermal, 1), "strain" => (f.strain_field, s.da_Thermal, 1),
        "strain_rate" => (f.strain_rate_field, s.da_Thermal, 1), "temperature" => (s.temperature, s.da_Thermal, 1),
        "velocity" => (s.velocity_next, s.da_Veloc, 2), "heat" => (f.radiogenic_heat_field, s.da_Thermal, 1),
        "thermal_diffusivity" => (f.thermal_diffusivity_field, s.da_Thermal, 1), "Phi" => (f.melt_fraction_field, s.da_Thermal, 1),
        "dPhi" => (f.melt_increment_field, s.da_Thermal, 1), "X_depletion" => (f.depletion_field, s.da_Thermal, 1))
    info = Dict(OUTPUT_VARIABLES)
    if m.dimensions == 3
        for (n, name) in enumerate(w.variables)
            vec, da, dof = sources[name]; root = (n - 1) % np; T, u = info[name]
            dof = name == "velocity" ? 3 : dof
            a = gather_field3(vec, da, dof, m.Nx, m.Ny, m.Nz, root)
            a === nothing && continue
            a[abs.(a) .< 1e-50] .= 0
            fields = dof == 3 ? [("vx", T, u, a[1, :, :, :]), ("vy", T, u, a[2, :, :, :]), ("vz", T, u, a[3, :, :, :])] :
                [(name, T, u, a[1, :, :, :])]
            write_mesh_store3(w, name, m, fields)
        end
        MPI.Barrier(comm); w.seconds += MPI.Wtime() - start
        return
    end
    for (n, name) in enumerate(w.variables)                          # each variable is written by one rank
        vec, da, dof = sources[name]; root = (n - 1) % np; T, u = info[name]
        a = gather_field(vec, da, dof, m.Nx, m.Nz, root)
        a === nothing && continue
        a[abs.(a) .< 1e-50] .= 0
        fields = dof == 2 ? [("vx", T, u, a[1, :, :]), ("vz", T, u, a[2, :, :])] : [(name, T, u, a[1, :, :])]
        write_mesh_store(w, name, m.Nx, m.Nz, m.Lx, m.domain_depth, fields)
    end
    if w.export_lithology
        root = length(w.variables) % np; mesh = lithology_mesh(m, s.particle_swarm, root)
        mesh === nothing || write_mesh_store(w, "lithology", size(mesh)..., m.Lx, m.domain_depth, [("lithology", Int8, "ID", mesh)])
    end
    w.write_particles && write_trajectories(w, m, s.particle_swarm, (length(w.variables) + 1) % np)
    surface !== nothing && rank == (length(w.variables) + 2) % np && write_surface_store(w, surface, m.Lx)
    MPI.Barrier(comm)
    w.seconds += MPI.Wtime() - start
end
