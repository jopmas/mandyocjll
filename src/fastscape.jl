# Coupling with FastScape (FastScapeLib, Braun & Willett, 2013; Yuan et al., 2019), sp_mode = fastscape.
# Include before surface.jl. The library (Fastscapelib_jll) keeps a single grid in global (Fortran module) variables,
# so FastScape runs only on rank 0 and the new surface is broadcast to the other ranks.
#
# The FastScape grid is the surface grid: the regular points of the surface (surface_particles_per_element per
# element in x, and in y in 3D). In 2D the profile is extended in y to a strip of 4 identical rows, spaced by
# fastscape_strip_width, with periodic boundaries in y: each row is a river of that width, so the drainage area of
# a point is the length of the profile upstream times fastscape_strip_width (with 3 rows the periodic boundary of
# FastScape mixes the rows; test/test_surface.jl checks the steady profile).
#
# At each step the tectonic change of the surface (its advection by the velocity of mandyocjll since the end of the
# last FastScape call) is given to FastScape as an uplift rate, and FastScape erodes and deposits during dt, in
# substeps of at most fastscape_dt_max years. The uplift is applied at each substep, interleaved with the erosion.
using Fastscapelib_jll: libfastscapelib_fortran

# ---- FastScape API (FastScape_api.f90: arguments by reference, arrays with x varying fastest) -------------------
fs_init() = ccall((:fastscape_init_, libfastscapelib_fortran), Cvoid, ())
fs_setup() = ccall((:fastscape_setup_, libfastscapelib_fortran), Cvoid, ())
fs_execute_step() = ccall((:fastscape_execute_step_, libfastscapelib_fortran), Cvoid, ())
fs_set_nx_ny(nx, ny) = ccall((:fastscape_set_nx_ny_, libfastscapelib_fortran), Cvoid, (Ref{Cint}, Ref{Cint}), nx, ny)
fs_set_xl_yl(xl, yl) = ccall((:fastscape_set_xl_yl_, libfastscapelib_fortran), Cvoid, (Ref{Float64}, Ref{Float64}), xl, yl)
fs_set_dt(dt) = ccall((:fastscape_set_dt_, libfastscapelib_fortran), Cvoid, (Ref{Float64},), dt)
fs_set_bc(code) = ccall((:fastscape_set_bc_, libfastscapelib_fortran), Cvoid, (Ref{Cint},), code)
fs_init_h(h) = ccall((:fastscape_init_h_, libfastscapelib_fortran), Cvoid, (Ptr{Float64},), h)
fs_set_u(u) = ccall((:fastscape_set_u_, libfastscapelib_fortran), Cvoid, (Ptr{Float64},), u)
fs_copy_h!(h) = ccall((:fastscape_copy_h_, libfastscapelib_fortran), Cvoid, (Ptr{Float64},), h)
fs_set_erosional_parameters(kf, kfsed, m, n, kd, kdsed, g1, g2, p) = ccall((:fastscape_set_erosional_parameters_, libfastscapelib_fortran),
    Cvoid, (Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
    kf, kfsed, m, n, kd, kdsed, g1, g2, p)
fs_set_marine_parameters(sl, p1, p2, z1, z2, r, l, kds1, kds2) = ccall((:fastscape_set_marine_parameters_, libfastscapelib_fortran),
    Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
    sl, p1, p2, z1, z2, r, l, kds1, kds2)

"Parameters of FastScape read from param.txt (keys fastscape_*); units of FastScape: m and yr."
Base.@kwdef struct FastScapeParameters
    kf::Float64 = 1.0e-5                 # fastscape_kf: stream power coefficient of the bedrock (m^(1-2m)/yr)
    kfsed::Float64 = kf                  # fastscape_kfsed: the same for the sediments
    m::Float64 = 0.4                     # fastscape_m: exponent of the drainage area
    n::Float64 = 1.0                     # fastscape_n: exponent of the slope
    kd::Float64 = 1.0e-2                 # fastscape_kd: hillslope diffusivity of the bedrock (m2/yr)
    kdsed::Float64 = kd                  # fastscape_kdsed: the same for the sediments
    g::Float64 = 1.0                     # fastscape_g: deposition coefficient of the bedrock (dimensionless)
    gsed::Float64 = g                    # fastscape_gsed: the same for the sediments
    p::Float64 = -2.0                    # fastscape_p: exponent of the multiple flow direction; < -1.5: single direction
    boundary::String = "0101"            # fastscape_boundary: bottom (y = 0), right, top and left; 1 = fixed base level
    dt_max::Float64 = 5000.0             # fastscape_dt_max: longest substep (yr)
    strip_width::Float64 = 10000.0       # fastscape_strip_width (2D): width of the river of each row (m)
    marine::Bool = false                 # fastscape_marine: marine transport and deposition below the sea level
    porosity_silt::Float64 = 0.63        # fastscape_porosity_silt
    porosity_sand::Float64 = 0.49        # fastscape_porosity_sand
    porosity_depth_silt::Float64 = 1960.0 # fastscape_porosity_depth_silt: e-folding depth (m)
    porosity_depth_sand::Float64 = 3700.0 # fastscape_porosity_depth_sand (m)
    sand_ratio::Float64 = 0.5            # fastscape_sand_ratio: fraction of sand in the flux from the continent
    mixed_layer::Float64 = 100.0         # fastscape_mixed_layer: thickness of the mixed layer at the sea floor (m)
    kd_silt::Float64 = 300.0             # fastscape_kd_silt: marine diffusivity of the silt (m2/yr)
    kd_sand::Float64 = 500.0             # fastscape_kd_sand (m2/yr)
end

function FastScapeParameters(d::Dict{String, String})
    num(k, x) = parse(Float64, get(d, "fastscape_$k", string(x)))
    kf, kd, g = num("kf", 1.0e-5), num("kd", 1.0e-2), num("g", 1.0)
    boundary = strip(get(d, "fastscape_boundary", "0101"))
    occursin(r"^[01]{4}$", boundary) || error("fastscape_boundary = $boundary: use 4 digits 0 or 1 (bottom, right, top, left)")
    return FastScapeParameters(; kf, kfsed = num("kfsed", kf), m = num("m", 0.4), n = num("n", 1.0), kd, kdsed = num("kdsed", kd),
        g, gsed = num("gsed", g), p = num("p", -2.0), boundary, dt_max = num("dt_max", 5000.0), strip_width = num("strip_width", 10000.0),
        marine = lowercase(get(d, "fastscape_marine", "false")) == "true",
        porosity_silt = num("porosity_silt", 0.63), porosity_sand = num("porosity_sand", 0.49),
        porosity_depth_silt = num("porosity_depth_silt", 1960.0), porosity_depth_sand = num("porosity_depth_sand", 3700.0),
        sand_ratio = num("sand_ratio", 0.5), mixed_layer = num("mixed_layer", 100.0), kd_silt = num("kd_silt", 300.0),
        kd_sand = num("kd_sand", 500.0))
end

"State of the coupling: the surface at the end of the last FastScape call (on all ranks) and the grid of FastScape."
mutable struct FastScapeCoupling
    p::FastScapeParameters
    nx::Int; ny::Int                     # grid of FastScape (ny = 4 in 2D)
    rows::Int                            # rows of the surface (1 in 2D)
    previous::Vector{Float64}            # surface after the last call, x varying fastest
end

"""
Creates the grid of FastScape on rank 0 with the initial surface `z` (a vector in 2D, nx x ny in 3D), spanning
`Lx` (and `Ly` in 3D).
"""
function FastScapeCoupling(p::FastScapeParameters, z, Lx, Ly)
    three = z isa AbstractMatrix
    nx = size(z, 1); rows = three ? size(z, 2) : 1; ny = three ? rows : 4
    three && ny < 4 && error("sp_mode = fastscape needs at least 4 surface points in y (ny and surface_particles_per_element)")
    yl = three ? Ly : (ny - 1) * p.strip_width
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        fs_init(); fs_set_nx_ny(nx, ny); fs_setup(); fs_set_xl_yl(Lx, yl)
        fs_init_h(three ? vec(copy(z)) : repeat(vec(z), ny))
        fs_set_erosional_parameters(fill(p.kf, nx * ny), p.kfsed, p.m, p.n, fill(p.kd, nx * ny), p.kdsed, p.g, p.gsed, p.p)
        # in 2D the strip is periodic in y, whatever the digits of y in fastscape_boundary
        fs_set_bc(parse(Int, three ? p.boundary : "0" * p.boundary[2] * "0" * p.boundary[4]))
    end
    return FastScapeCoupling(p, nx, ny, rows, vec(copy(z)))
end

"""
Erosion and deposition during dt (s) on the surface `z` (moved by the tectonics since the last call), with the
sea level `sea_level` (m, for fastscape_marine). Changes z in place, on all ranks.
"""
function fastscape_step!(fs::FastScapeCoupling, z, dt, sea_level)
    years = dt / SECONDS_PER_YEAR
    n = max(1, ceil(Int, years / fs.p.dt_max - 1e-9))
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        uplift = (vec(z) .- fs.previous) ./ years                   # tectonic change, m/yr
        fs_set_u(fs.rows == fs.ny ? uplift : repeat(uplift, fs.ny))
        p = fs.p
        p.marine && fs_set_marine_parameters(sea_level, p.porosity_silt, p.porosity_sand, p.porosity_depth_silt,
            p.porosity_depth_sand, p.sand_ratio, p.mixed_layer, p.kd_silt, p.kd_sand)
        fs_set_dt(years / n)
        for _ in 1:n; fs_execute_step(); end
        h = zeros(fs.nx, fs.ny); fs_copy_h!(h)
        new = fs.rows == fs.ny ? vec(h) : vec(sum(h; dims = 2) ./ fs.ny)   # 2D: mean of the rows of the strip
    else
        new = zeros(length(z))
    end
    MPI.Bcast!(new, 0, MPI.COMM_WORLD)
    vec(z) .= new                        # vec of a Matrix shares its memory
    fs.previous .= new
    return z
end
