# mpiexec -n N julia --project=. test/test_surface.jl [FOLDER]   (or: julia --project=. test/mpirun.jl N test_surface.jl)
# Surface tracking in 2D and 3D with a prescribed velocity, and the coupling with FastScape (steady river profile).
# With FOLDER, the river profiles are written to FOLDER/river_profile_{2d,3d}.csv (figure of docs/apostila).
using Mandyocjll, MPI, Printf
const M = Mandyocjll
const comm = MPI.COMM_WORLD
report(ok, s) = MPI.Comm_rank(comm) == 0 && println("RESULT ", ok ? "OK   " : "FAIL ", "np=", MPI.Comm_size(comm), " ", s)

# velocity (U, V, W(x, y)): independent of z, so the height is h0 at the start of the path plus the integral of W
const L, D, U, V, W0, H0 = 100e3, 50e3, 2e-10, 1e-10, 5e-11, -10e3
W(x, y) = W0 * sinpi(2x / L) * (1 + 0.5cospi(y / L))
h0(x, y) = H0 + 300sinpi(x / L) * cospi(y / L)
function exact(x, y, t; n = 200)                 # midpoint rule along the characteristic
    s = (0.5:n) .* (t / n)
    return h0(x - U * t, y - V * t) + sum(W(x - U * (t - si), y - V * (t - si)) for si in s) * t / n
end

surface(x, y, z, mode = "none", fs = nothing) = M.Surface(x, y, z, 2, mode != "none", mode, 0.0, 0.0, true, 0.0, 0.6, 0.0, 0,
    Tuple{Float64, Int}[], Tuple{Float64, Float64}[], Tuple{Float64, Float64}[], fs)

function velocity_vector(da, f, dof)
    v = M.LibPETSc.DMCreateGlobalVector(M.P, da)
    loc = M.LibPETSc.DMCreateLocalVector(M.P, da)
    M.global_to_local(da, dof == 3 ? M.fill_field3!(v, da, f, 3) : M.fill_field!(v, da, f, 2), loc)
    return loc
end

# returns the largest error in the interior (away from the inflow edges, where the C++ rule keeps the edge points)
function advection_error(three, steps)
    t = 2.0e13; dt = t / steps                   # 0.63 Myr: displacements of 4 km (U) and 2 km (V), uplift of 3 km
    if three
        m = M.ModelParameters(; dimensions = 3, Nx = 21, Ny = 11, Nz = 6, Lx = L, Ly = L, domain_depth = D)
        da = M.create_velocity_mesh3d(m)
        dx, dy, dz = M.element_size_x(m), M.element_size_y(m), M.element_size_z(m)
        lv = velocity_vector(da, (i, j, k) -> (U, V, W(i * dx, j * dy)), 3)
        x = collect(range(0, L, length = 41)); y = collect(range(0, L, length = 21))
        sf = surface(x, y, [h0(a, b) for a in x, b in y])
        for _ in 1:steps; M.advect_surface!(sf, m, da, lv, dt); end
        return maximum(abs(sf.z[i, j] - exact(x[i], y[j], t)) for i in eachindex(x), j in eachindex(y) if x[i] > 10e3 && y[j] > 5e3)
    else
        m = M.ModelParameters(; Nx = 21, Nz = 6, Lx = L, domain_depth = D)
        da = M.create_advection_mesh(m, M.create_velocity_mesh(m))
        dx = M.element_size_x(m)
        lv = velocity_vector(da, (i, k) -> (U, W(i * dx, 0.0)), 2)
        x = collect(range(0, L, length = 41))
        sf = surface(x, Float64[], [h0(a, 0.0) for a in x])
        for _ in 1:steps
            M.advect_surface!(sf, m, da, lv, dt); M.interpolate_surface!(sf, m)
        end
        return maximum(abs(sf.z[i] - exact(x[i], 0.0, t; n = 200)) for i in eachindex(x) if x[i] > 10e3)
    end
end

# steady state of the stream power law with n = 1 and a fixed base level only at x = 0: the drainage area of the
# node i (from 0) is (nx - i) dx w, so h_i - h_0 = sum_{j <= i} (U/K) dx ((nx - j) dx w)^-m
function river_profile(three)
    nx, w, Lr, K, m_, uplift = 51, 20e3, 50e3, 2e-5, 0.4, 1e-3        # uplift in m/yr
    p = M.FastScapeParameters(; kf = K, m = m_, n = 1.0, kd = 0.0, g = 0.0, boundary = "0001", dt_max = 1e4, strip_width = w)
    x = collect(range(0, Lr, length = nx)); dx = x[2]
    z = three ? zeros(nx, 5) : zeros(nx)
    fs = M.FastScapeCoupling(p, z, Lr, three ? 4w : 0.0)
    dt = 2e4 * M.SECONDS_PER_YEAR
    for _ in 1:1000                              # 20 Myr; the response time is about 1 Myr
        z[2:end, :] .+= uplift * 2e4             # the tectonic change of the step; the base level (x = 0) stays
        M.fastscape_step!(fs, z, dt, -1e9)
    end
    prof = three ? z[:, 3] : z
    expected = [sum((uplift / K) * dx * ((nx - j) * dx * w)^(-m_) for j in 1:i; init = 0.0) for i in 0:nx-1]
    spread = three ? maximum(abs.(z .- z[:, 1:1])) : 0.0
    if !isempty(ARGS) && MPI.Comm_rank(comm) == 0
        mkpath(ARGS[1])
        open(joinpath(ARGS[1], "river_profile_$(three ? "3d" : "2d").csv"), "w") do f
            println(f, "x,fastscape,exact")
            for i in eachindex(x); println(f, x[i], ",", prof[i] - prof[1], ",", expected[i]); end
        end
    end
    return maximum(abs.((prof .- prof[1]) .- expected)) / maximum(expected), spread, maximum(expected)
end

function main()
    M.PETSc.initialize(M.P)
    e2 = advection_error(false, 40); e2b = advection_error(false, 80)
    # the error is dominated by the interpolation in space (W on the mesh, the height between the points)
    report(e2 < 45 && e2b < 45, @sprintf("2D surface, prescribed velocity: largest error %.2f m (40 steps), %.2f m (80 steps); uplift 3 km", e2, e2b))
    e3 = advection_error(true, 40); e3b = advection_error(true, 80)
    report(e3 < 45 && e3b < 45, @sprintf("3D surface, prescribed velocity: largest error %.2f m (40 steps), %.2f m (80 steps); uplift 3 km", e3, e3b))
    for three in (false, true)                   # FastScape runs on rank 0, the other ranks receive the surface
        r, spread, relief = river_profile(three)
        report(r < 1e-3 && spread < 1e-6, @sprintf("FastScape %s, steady river profile: relative error %.2e, relief %.0f m, spread in y %.1e m",
            three ? "3D grid" : "2D strip", r, relief, spread))
    end
    M.PETSc.finalize(M.P)
end
main()
