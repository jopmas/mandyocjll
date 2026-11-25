# julia --project=. test/test_rheology.jl <src of the C++ Mandyoc> <include of PETSc> [include of MPI]
# Compares the viscosity of each model with the function calc_visco_ponto of the C++ code.
include("../src/particle_swarm.jl"); include("../src/thermal.jl"); include("../src/stokes.jl"); include("../src/rheology.jl")
using Printf
src, inc = abspath(ARGS[1]), ARGS[2:end]
lib = joinpath(mktempdir(), "libvisc.so")
undefined_flags = Sys.isapple() ? ["-undefined", "dynamic_lookup"] : String[]
define = "-DCALC_VISC=\"$src/calc_visc.cpp\""
run(`c++ -O1 -w -shared -fPIC $undefined_flags $define $("-I" .* inc) -o $lib $(joinpath(@__DIR__, "viscosity_harness.cpp"))`)
h = Libdl.dlopen(lib); set_p = Libdl.dlsym(h, :set_params); visc_p = Libdl.dlsym(h, :visc)
seed = UInt64(1); rnd(a, b) = (global seed = hash(seed); a + (b - a) * (seed >> 11) * 2.0^-53)
worst = Dict{Int, Float64}(); total = 0
for model in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 70), trial in 1:4000
    nl = model in (9, 19) || rnd(0, 1) > 0.5; pl = rnd(0, 1) > 0.3; pr = rnd(0, 1) > 0.3; vf = rnd(0, 1) > 0.2
    lo, hi = 10.0^rnd(17, 19), 10.0^rnd(22, 26); t = rnd(0, 1)
    m = ModelParameters(rheology_model = model, use_nonlinear_method = nl, use_plasticity = pl, use_pressure_in_rheology = pr,
        use_viscosity_factor = vf, min_viscosity = lo, max_viscosity = hi, reference_viscosity = 10.0^rnd(19, 23),
        temperature_difference = rnd(1000, 1600), air_layer_thickness = pr ? 0.0 : rnd(0, 4e4), constant_rheology_pressure = rnd(0, 1) > 0.9 ? 1e8 : -1.0)
    range = ViscosityRange(lo * 10^rnd(0, 1), hi / 10^rnd(0, 1))
    n = rnd(0, 1) > 0.5 ? 3.0 : rnd(1, 4.5)
    l = LayerProperties(layer_viscosity_factor = [10.0^rnd(-1, 1.5)], layer_density = [3300.0], layer_radiogenic_heat = [0.0],
        layer_pre_exponential_factor = [model == 7 ? 10.0^rnd(-2, 14) : 10.0^rnd(-28, -13)], layer_power_law_exponent = [n], layer_activation_energy = [model == 8 ? rnd(0, 0.01) : rnd(0, 5.5e5)],
        layer_activation_volume = [model == 8 ? rnd(0, 1e-5) : rnd(0, 25e-6)], layer_cohesion_min = [rnd(1e6, 8e6)], layer_cohesion_max = [rnd(1e7, 3e7)],
        layer_friction_angle_min = [rnd(0.01, 0.1)], layer_friction_angle_max = [rnd(0.15, 0.6)])
    T, Pr, z, e2, s = rnd(0, 1500), rnd(-1e8, 1e10), rnd(-3e5, 0), rnd(0, 1) > 0.95 ? 0.0 : 10.0^rnd(-20, -12), rnd(0, 2)
    d = [hi, lo, range.computed_max_viscosity, range.computed_min_viscosity, m.reference_viscosity, m.temperature_difference, m.air_layer_thickness,
        m.constant_rheology_pressure, m.weakening_min, m.weakening_max, l.layer_cohesion_min[1], l.layer_cohesion_max[1], l.layer_friction_angle_min[1], l.layer_friction_angle_max[1]]
    ccall(set_p, Cvoid, (Cint, Cint, Cint, Cint, Cint, Ptr{Float64}), model, vf, nl, pl, pr, d)
    c = ccall(visc_p, Float64, (Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64), T, Pr, z, l.layer_viscosity_factor[1], e2, s, l.layer_pre_exponential_factor[1], n,
        l.layer_activation_energy[1], l.layer_activation_volume[1])
    j = particle_viscosity(m, range, l, 1, T, Pr, z, e2, s)
    worst[model] = max(get(worst, model, 0.0), abs(j - c) / abs(c)); global total += 1
end
for k in sort(collect(keys(worst))); @printf("RESULT %s rheology %2d: max. relative difference Julia - C++ %.2g\n", worst[k] < 1e-11 ? "OK  " : "FAIL", k, worst[k]); end
println(total, " points")
