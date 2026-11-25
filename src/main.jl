# Runs a scenario of Mandyoc; called by bin/mandyocjll on each rank (included by Mandyocjll.jl).
using Printf

# options that change the result and have not been ported yet: stop instead of ignoring them
function check_unsupported(d)
    on(k) = lowercase(get(d, k, "false")) in ("true", "on")
    bad = [k for k in ("sp_surface_processes", "magmatism", "nondimensionalization", "multi_velocity", "plot_sediment") if on(k)]
    get(d, "sp_mode", "none") == "none" || push!(bad, "sp_mode")
    get(d, "RK4", "Euler") == "Euler" || push!(bad, "RK4")
    parse(Float64, get(d, "basal_heat", "-1")) > 0 && push!(bad, "basal_heat")
    parse(Int, get(d, "dimensions", "2")) == 2 || push!(bad, "dimensions")
    isempty(bad) || error("options of param.txt not ported yet: " * join(bad, ", "))
    on("interfaces_from_ascii") || error("use interfaces_from_ascii = True")
end

"-seed and -strain_seed: initial strain of the given layers (parse_options and reader in the C++ code)."
function apply_seeds!(layers, seeds, strains, say)
    isempty(seeds) && (isempty(strains) || error("specify the seed layers with -seed (required by -strain_seed)"); return)
    if isempty(strains)
        say("Using the default value 0.5 for -strain_seed (all seed layers)")
        strains = fill(0.5, length(seeds))
    end
    length(seeds) == length(strains) || error("-seed and -strain_seed must have the same number of values")
    n = length(layers.layer_density)
    for (l, s) in zip(seeds, strains)
        0 <= l < n || error("-seed: layer $l does not exist (interfaces.txt has layers 0 to $(n - 1))")
        if layers.layer_weakening_seed[l+1] >= 0
            say("seed layer: $l - using the weakening_seed of interfaces.txt ($(layers.layer_weakening_seed[l+1])) instead of $s")
        else
            layers.layer_weakening_seed[l+1] = s
            say("seed layer: $l - strain: $s")
        end
    end
end

function run_scenario(dir, out; seeds = Int[], strain_seeds = Float64[], h_air = nothing, Px = -1, Pz = -1, petsc_options = String[])
    PETSc.initialize(P; options = petsc_options); comm = MPI.COMM_WORLD; rank, np = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    say(args...) = rank == 0 && (println(args...); flush(stdout))
    rank == 0 && check_output_folder(out)
    MPI.Barrier(comm)
    d = read_param_file(joinpath(dir, "param.txt")); check_unsupported(d)
    overrides = h_air === nothing ? (;) : (; air_layer_thickness = h_air)
    m = ModelParameters(d; fssa_legacy = lowercase(get(d, "fssa_legacy", "false")) == "true", Px, Pz, overrides...)
    m.use_pressure_in_rheology || h_air !== nothing || error("specify the thickness of the air layer with -h_air (depth dependent rheology)")
    layers = read_interfaces(joinpath(dir, "interfaces.txt"), m)
    say("mandyocjll: $np MPI rank(s), mesh $(m.Nx) x $(m.Nz), $(particles_per_element(m)) particles per element")
    apply_seeds!(layers, seeds, strain_seeds, say)
    on(k) = lowercase(get(d, k, "false")) == "true"; num(k, x) = parse(Float64, get(d, k, string(x)))
    on("sp_surface_tracking") && say("warning: sp_surface_tracking has not been ported yet; the surface will not be written")
    T_initial = Int(num("T_initial", 0))
    T0 = on("temperature_from_ascii") ? read_ascii_field(joinpath(dir, "input_temperature_0.txt")) : [initial_temperature(m, T_initial, i, k) for k in 0:m.Nz-1 for i in 0:m.Nx-1]
    mask = on("bcv_extern") ? read_ascii_field(joinpath(dir, "bcv_0.txt")) : nothing       # 1 = free, 0 = fixed
    velocity_mask = mask === nothing ? nothing : (i, k) -> (mask[2(k*m.Nx+i)+1], mask[2(k*m.Nx+i)+2])
    velocity_changes = on("variable_bcv") ? read_scale_bcv(joinpath(dir, "scale_bcv.txt")) : Tuple{Float64, Float64}[]
    V0 = on("velocity_from_ascii") ? read_ascii_field(joinpath(dir, "input_velocity_0.txt")) : zeros(2m.Nx * m.Nz)
    random_strain = num("random_initial_strain", 0.0)
    initial_strain(x, z, l) = random_strain * (unit_random(m.random_seed, hash((x, z)), 3) + 0.5)
    location = m.periodic_boundary ? :shim : :native                 # the periodic DMDA has a different node spacing
    s = Model(m, layers; temperature = (i, k) -> T0[k*m.Nx+i+1], velocity = (i, k) -> (V0[2(k*m.Nx+i)+1], V0[2(k*m.Nx+i)+2]), location, initial_strain, velocity_mask)
    if rank == 0 && abspath(out) != abspath(dir)                     # tapIOca reads param.txt from the output folder
        mkpath(out); cp(joinpath(dir, "param.txt"), joinpath(out, "param.txt"); force = true)
    end
    MPI.Barrier(comm)
    w = OutputWriter(out; export_thermal_diffusivity = on("export_thermal_diffusivity"), export_lithology = on("export_lithology"),
        write_particles = on("print_step_files"))
    step_max, time_max = Int(num("step_max", 1e9)), num("time_max", 1e30) * SECONDS_PER_YEAR
    step_print = Int(num("step_print", 10)); initial_print_step = Int(num("initial_print_step", 0))
    initial_print_max_time = num("initial_print_max_time", 1e6) * SECONDS_PER_YEAR
    start = MPI.Wtime()
    initial_velocity!(s)
    say(@sprintf("step 0: %d nonlinear iterations, %.1f s", s.nonlinear_iterations, MPI.Wtime() - start))
    write_output!(w, s, 0)
    for step in 1:step_max
        s.time + s.time_step_seconds <= time_max || break
        dt = s.time_step_seconds
        # as in the C++ code, the time compared is that of the end of the step, and one change is applied per step
        if !isempty(velocity_changes) && s.time + dt > velocity_changes[1][1] * 1e6 * SECONDS_PER_YEAR
            (t, scale) = popfirst!(velocity_changes)
            rescale_velocity!(s, scale)
            say("velocity multiplied by $scale ($t Myr)")
        end
        step!(s)
        added, removed = add_remove_particles!(m, s.particle_swarm, s.da_Veloc, step; location)
        every = initial_print_step > 0 && s.time < initial_print_max_time ? initial_print_step : step_print
        step % every == 0 && write_output!(w, s, step)
        n = MPI.Allreduce(Int(LibPETSc.DMSwarmGetLocalSize(P, s.particle_swarm)), +, comm)
        say(@sprintf("step %d: t = %.4g yr, dt = %.4g yr, %d nonlinear iterations, %d particles (+%d -%d)", step, s.time / SECONDS_PER_YEAR,
            dt / SECONDS_PER_YEAR, s.nonlinear_iterations, n, added, removed))
    end
    say(@sprintf("done: %.1f s, %.2f s of them writing the outputs (%d times)", MPI.Wtime() - start, w.seconds, length(w.times)))
    PETSc.finalize(P)
end
