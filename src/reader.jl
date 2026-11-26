# Reader of param.txt and interfaces.txt of Mandyoc (same keys). Only what has been ported.
function read_param_file(file)
    d = Dict{String, String}()
    for line in eachline(file)
        l = strip(first(split(line, '#')))
        occursin('=', l) || continue
        k, v = strip.(split(l, '='; limit = 2)); d[k] = first(split(v))
    end
    return d
end
function ModelParameters(d::Dict{String, String}; overrides...)
    f(k, default) = haskey(d, k) ? parse(Float64, d[k]) : default
    b(k, default) = haskey(d, k) ? lowercase(d[k]) in ("true", "on", "fixed") : default
    Nx, Nz = parse(Int, d["nx"]), parse(Int, d["nz"]); Lx, depth = f("lx", 0.0), f("lz", 0.0)
    ppe = Int(f("particles_per_element", 81.0))
    nz = trunc(Int, sqrt(ppe * (depth / (Nz - 1)) / (Lx / (Nx - 1)))); nx = trunc(Int, ppe / nz)
    g(k, n) = (v = Int(f(k, 0.0)); v > 0 ? v : n)
    def = ModelParameters()
    min_viscosity = f("viscosity_min", 1e17); max_viscosity = f("viscosity_max", 1e25)
    max_viscosity == min_viscosity && (max_viscosity = 1.0001min_viscosity)      # as in the C++ code
    return ModelParameters(; Nx, Nz, Lx, domain_depth = depth, particles_per_element_x = g("particles_per_element_x", nx),
        particles_per_element_z = g("particles_per_element_z", nz), particles_perturb_factor = f("particles_perturb_factor", 0.5),
        use_runge_kutta_4 = get(d, "RK4", "Euler") == "Runge-Kutta", periodic_boundary = b("periodic_boundary", false),
        thermal_diffusivity = f("thermal_diffusivity_coefficient", def.thermal_diffusivity),
        thermal_expansion_coefficient = f("thermal_expansion_coefficient", def.thermal_expansion_coefficient),
        gravity_acceleration = f("gravity_acceleration", 10.0), heat_capacity = f("heat_capacity", 1250.0),
        mantle_density = f("density_mantle", 3300.0), temperature_difference = f("temperature_difference", 1500.0),
        use_radiogenic_heating = b("radiogenic_component", false), use_adiabatic_heating = b("adiabatic_component", false),
        use_shear_heating = b("shear_heating", false), shear_heating_scale_factor = f("Xi_shear", 1.0), high_diffusivity_in_asthenosphere = b("high_kappa_in_asthenosphere", false),
        temperature_bc_top = b("top_temperature", true), temperature_bc_bottom = b("bot_temperature", true),
        temperature_bc_left = b("left_temperature", true), temperature_bc_right = b("right_temperature", true),
        solver_relative_tolerance = f("rtol", 1e-5), uzawa_residual_tolerance = f("denok", 1e-4),
        nonlinear_residual_tolerance = f("Xi_min", 1e-14), use_nonlinear_method = b("non_linear_method", false),
        use_plasticity = b("plasticity", true), rheology_model = Int(f("rheology_model", 0.0)),
        reference_viscosity = f("viscosity_reference", 1e21), max_viscosity, min_viscosity,
        use_viscosity_factor = b("geoq", true), viscosity_factor_scale = f("geoq_fac", 1.0),
        use_pressure_in_rheology = get(d, "viscosity_dependence", "depth") == "pressure",
        constant_rheology_pressure = f("pressure_const", -1.0), weakening_min = f("weakening_min", 0.05), weakening_max = f("weakening_max", 1.05),
        use_viscosity_harmonic_mean = get(d, "viscosity_mean_method", "harmonic") == "harmonic",
        sticky_blanket_air = b("sticky_blanket_air", false), air_threshold_density = f("air_threshold_density", 100.0),
        use_direct_solver = get(d, "solver", "direct") == "direct", free_surface_stabilization = b("free_surface_stab", true),
        fssa_theta = f("theta_FSSA", 0.5), initial_dynamic_range = b("initial_dynamic_range", false),
        use_strain_healing = b("strain_healing", false), strain_healing_rate = f("A_healing", def.strain_healing_rate),
        strain_healing_activation_factor = f("k_healing", def.strain_healing_activation_factor),
        strain_healing_reference_temperature = f("T0_healing", def.strain_healing_reference_temperature),
        max_time_step = f("dt_max", 1e6), time_step_subdivision = f("sub_division_time_step", 1.0),
        velocity_bc_top_normal = b("top_normal_velocity", true), velocity_bc_top_tangential = b("top_tangential_velocity", false),
        velocity_bc_bottom_normal = b("bot_normal_velocity", true), velocity_bc_bottom_tangential = b("bot_tangential_velocity", false),
        velocity_bc_left_normal = b("left_normal_velocity", true), velocity_bc_left_tangential = b("left_tangential_velocity", false),
        velocity_bc_right_normal = b("right_normal_velocity", true), velocity_bc_right_tangential = b("right_tangential_velocity", false), overrides...)
end

function read_interfaces(file, m::ModelParameters)
    rows = Dict{String, Vector{Float64}}(); data = Vector{Float64}[]
    for line in eachline(file)
        t = split(line); isempty(t) && continue
        isletter(first(t[1])) ? (rows[t[1]] = parse.(Float64, t[2:end])) : push!(data, parse.(Float64, t))
    end
    depths = permutedims(reduce(hcat, data))
    size(depths, 1) == m.Nx || error("interfaces.txt has $(size(depths, 1)) rows, but nx = $(m.Nx)")
    n = length(rows["rho"]); get_row(k, default) = get(rows, k, fill(default, n))
    return LayerProperties(; layer_viscosity_factor = rows["C"], layer_density = rows["rho"], layer_radiogenic_heat = rows["H"],
        layer_pre_exponential_factor = rows["A"], layer_power_law_exponent = rows["n"], layer_activation_energy = rows["Q"],
        layer_activation_volume = rows["V"], layer_cohesion_min = get_row("cohesion_min", 4.0e6), layer_cohesion_max = get_row("cohesion_max", 20.0e6),
        layer_friction_angle_min = haskey(rows, "friction_angle_min") ? deg2rad.(rows["friction_angle_min"]) : fill(0.034906, n),
        layer_friction_angle_max = haskey(rows, "friction_angle_max") ? deg2rad.(rows["friction_angle_max"]) : fill(0.261799, n),
        layer_weakening_seed = get_row("weakening_seed", -1.0), layer_thermal_conductivity = get_row("k", m.thermal_diffusivity),
        interface_depths = depths)
end

"ASCII field of Mandyoc (4 header lines; lines starting with P are skipped)."
read_ascii_field(file) = [parse(Float64, l) for l in Iterators.drop(eachline(file), 4) if !isempty(strip(l)) && !startswith(l, "P")]

"Initial temperature of T_initial 0 and 1 (Thermal_init_2d)."
function initial_temperature(m::ModelParameters, T_initial, i, k)
    T_initial in (0, 1) || error("T_initial = $T_initial has not been ported: use temperature_from_ascii = True")
    d = m.temperature_difference
    k == m.Nz - 1 && return 0.0
    k == 0 && return d
    return clamp(d * (m.Nz - 1 - k) / (m.Nz - 1) + 100cos(i * 3.14159 / (m.Nx - 1)), 0.0, d)
end

"scale_bcv.txt: number of changes, then time (Myr) and factor of each one."
function read_scale_bcv(file)
    x = parse.(Float64, split(read(file, String)))
    n = Int(x[1])
    return [(x[2q], x[2q+1]) for q in 1:n]
end
