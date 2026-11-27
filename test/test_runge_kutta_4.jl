# julia --project=. test/test_runge_kutta_4.jl
# Order of the RK4 of the particles (runge_kutta_4_velocity in 2D, runge_kutta_4_velocity3 in 3D), without PETSc:
# a rigid rotation whose angular velocity grows linearly in time, ω(t) = ω0 (1 + β t) (V0 at t = 0, V1 at t = dt). The
# field is linear, so the interpolation is exact and the error of one step is that of the time integration:
# O(dt^5) for RK4 (ratio ~32 when dt is halved) and O(dt^2) for one Euler step.
using Mandyocjll, OffsetArrays, Printf
const M = Mandyocjll

function rotation_error(m, dims, dt; rk4 = true, β = 1.0e-12, ω = 1.0e-12)
    xc, yc = m.Lx / 2, (dims == 3 ? m.Ly / 2 : 0.0); zc = -m.domain_depth / 2
    dx, dz = M.element_size_x(m), M.element_size_z(m)
    if dims == 2      # rotation in the x-z plane
        V0 = OffsetArray(zeros(2, m.Nx, m.Nz), 1:2, 0:m.Nx-1, 0:m.Nz-1)
        for k in 0:m.Nz-1, i in 0:m.Nx-1
            x, z = i * dx, k * dz - m.domain_depth
            V0[1, i, k] = -ω * (z - zc); V0[2, i, k] = ω * (x - xc)
        end
    else              # rotation in the x-y plane, plus a uniform vertical velocity
        dy = M.element_size_y(m)
        V0 = OffsetArray(zeros(3, m.Nx, m.Ny, m.Nz), 1:3, 0:m.Nx-1, 0:m.Ny-1, 0:m.Nz-1)
        for k in 0:m.Nz-1, j in 0:m.Ny-1, i in 0:m.Nx-1
            x, y = i * dx, j * dy
            V0[1, i, j, k] = -ω * (y - yc); V0[2, i, j, k] = ω * (x - xc); V0[3, i, j, k] = ω * 1.0e3
        end
    end
    V1 = (1 + β * dt) .* V0
    r = 0.2 * min(m.Lx, m.domain_depth); θ = ω * (dt + β * dt^2 / 2)  # exact angle of the step
    if dims == 2
        x, z = xc + r, zc
        u, w = rk4 ? M.runge_kutta_4_velocity(m, V0, V1, x, z, dt) : M.interpolate_velocity(m, V0, x, z)
        return hypot(x + dt * u - (xc + r * cos(θ)), z + dt * w - (zc + r * sin(θ)))
    else
        x, y, z = xc + r, yc, zc
        v = rk4 ? M.runge_kutta_4_velocity3(m, V0, V1, x, y, z, dt) : M.interpolate_velocity3(m, V0, x, y, z)
        exact = (xc + r * cos(θ), yc + r * sin(θ), zc + ω * 1.0e3 * (dt + β * dt^2 / 2))
        return sqrt(sum(((x, y, z) .+ dt .* v .- exact) .^ 2))
    end
end

function main()
    ok = true
    for (dims, m) in ((2, M.ModelParameters(Nx = 41, Nz = 41, Lx = 1.0e6, domain_depth = 1.0e6)),
                      (3, M.ModelParameters(dimensions = 3, Nx = 21, Ny = 21, Nz = 21, Lx = 1.0e6, Ly = 1.0e6, domain_depth = 1.0e6)))
        dts = [0.08, 0.04, 0.02, 0.01] ./ 1.0e-12              # rotation angle of 0.08 to 0.01 rad
        e4 = [rotation_error(m, dims, dt) for dt in dts]
        e1 = [rotation_error(m, dims, dt; rk4 = false) for dt in dts]
        o4 = log2.(e4[1:end-1] ./ e4[2:end]); o1 = log2.(e1[1:end-1] ./ e1[2:end])
        good = all(4.7 .< o4 .< 5.3) && all(1.7 .< o1 .< 2.3)
        ok &= good
        @printf("RESULT %-4s %dD rotation, velocity varying in time: one-step error order RK4 %s (expected 5), Euler %s (expected 2); RK4 error %.2e m at dt = %.3g yr\n",
            good ? "OK" : "FAIL", dims, join((@sprintf("%.2f", o) for o in o4), ", "), join((@sprintf("%.2f", o) for o in o1), ", "),
            e4[end], dts[end] / M.SECONDS_PER_YEAR)
    end
    ok || exit(1)
end
main()
