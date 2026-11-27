# Runs a scenario of Mandyoc; called by bin/mandyocjll on each rank (included by Mandyocjll.jl).
using Printf

# options that change the result and have not been ported yet: stop instead of ignoring them
function check_unsupported(d)
    on(k) = lowercase(get(d, k, "false")) in ("true", "on")
    bad = String[]
    if parse(Int, get(d, "dimensions", "2")) == 3     # options the 3D C++ code does not have either
        get(d, "RK4", "Euler") == "Euler" || push!(bad, "RK4 in 3D")
        on("magmatism") && push!(bad, "magmatism in 3D")
        on("periodic_boundary") && push!(bad, "periodic_boundary in 3D")
    end
    isempty(bad) || error("options of param.txt not ported yet: " * join(bad, ", "))
    # as in the C++ code, where "variable" (its default) is no longer available and stops the run
    get(d, "viscosity_per_element", "variable") == "constant" ||
        error("use viscosity_per_element = constant (variable, the default of Mandyoc, is no longer available)")
    on("interfaces_from_ascii") || error("use interfaces_from_ascii = True")
end

"""
nondimensionalization = True: the C++ code solves the equations in units of the depth, the maximum viscosity,
the gravity, the mantle density and the diffusivity, and mixes these units with years in the time step, time_max
and the times of the input files (bug 27). The change of units does not change the solution, so here the equations
are solved in SI units, as without the option, and only the scales are written (dimensional_convertion.txt).
The tolerance denok, which depends on the units, is converted (run_scenario), and so is the extra g of the FSSA
term of fssa_legacy (bug 7).
"""
function nondimensional_scales(d, m::ModelParameters)
    h0, visc0 = m.domain_depth, parse(Float64, get(d, "viscosity_max", "1e25"))   # visc_MAX before the 1.0001 factor
    g0, rho0, kappa0 = m.gravity_acceleration, m.mantle_density, m.thermal_diffusivity
    return (; h0, visc0, g0, rho0, kappa0, veloc0 = rho0 * g0 * h0^2 / visc0)
end

function write_dimensional_scales(sc, folders, say)
    text = @sprintf("depth %g\nviscosity %g\ngravity %g\ndensity %g\ndiffusivity %g\n", sc.h0, sc.visc0, sc.g0, sc.rho0, sc.kappa0)
    foreach(f -> write(joinpath(f, "dimensional_convertion.txt"), text), unique(abspath.(folders)))
    say("nondimensionalization = True: solved in SI units (same solution, outputs in SI units); scales in dimensional_convertion.txt")
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

function run_scenario(dir, out; seeds = Int[], strain_seeds = Float64[], h_air = nothing, Px = -1, Py = -1, Pz = -1, petsc_options = String[])
    PETSc.initialize(P; options = petsc_options); comm = MPI.COMM_WORLD; rank, np = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    say(args...) = rank == 0 && (println(args...); flush(stdout))
    rank == 0 && check_output_folder(out)
    MPI.Barrier(comm)
    d = read_param_file(joinpath(dir, "param.txt")); check_unsupported(d)
    overrides = h_air === nothing ? (;) : (; air_layer_thickness = h_air)
    m = ModelParameters(d; fssa_legacy = lowercase(get(d, "fssa_legacy", "false")) == "true", Px, Py, Pz, overrides...)
    if lowercase(get(d, "nondimensionalization", "false")) == "true"
        # denok bounds the squared divergence of the velocity (G^T v, in m2/s): in the C++ code it is in units of
        # veloc0*h0, so it is converted to SI units
        sc = nondimensional_scales(d, m)
        m = ModelParameters(d; fssa_legacy = m.fssa_legacy, Px, Py, Pz, overrides...,
            uzawa_residual_tolerance = m.uzawa_residual_tolerance * (sc.veloc0 * sc.h0)^2, fssa_legacy_gravity_unit = 1 / sc.g0)
    end
    # denok bounds the squared divergence (G^T v) in SI units: a large value stops the Uzawa iterations at once (bug 28)
    m.uzawa_residual_tolerance > 1e-10 && say(@sprintf("warning: denok = %g (SI units) may leave the velocity compressible; tapIOca uses 1e-14", m.uzawa_residual_tolerance))
    m.use_pressure_in_rheology || h_air !== nothing || error("specify the thickness of the air layer with -h_air (depth dependent rheology)")
    layers = read_interfaces(joinpath(dir, "interfaces.txt"), m)
    three = m.dimensions == 3
    mesh = three ? "$(m.Nx) x $(m.Ny) x $(m.Nz)" : "$(m.Nx) x $(m.Nz)"
    say("mandyocjll: $np MPI rank(s), mesh $mesh, $(particles_per_element(m)) particles per element")
    apply_seeds!(layers, seeds, strain_seeds, say)
    on(k) = lowercase(get(d, k, "false")) == "true"; num(k, x) = parse(Float64, get(d, k, string(x)))
    T_initial = Int(num("T_initial", 0))
    three && !on("temperature_from_ascii") && error("in 3D, use temperature_from_ascii = True")
    T0 = on("temperature_from_ascii") ? read_ascii_field(joinpath(dir, "input_temperature_0.txt")) : [initial_temperature(m, T_initial, i, k) for k in 0:m.Nz-1 for i in 0:m.Nx-1]
    mask = on("bcv_extern") ? read_ascii_field(joinpath(dir, "bcv_0.txt")) : nothing       # 1 = free, 0 = fixed
    # nodal fields of the ascii files: x first, then y (3D), then z
    node(i, k) = k * m.Nx + i; node(i, j, k) = (k * m.Ny + j) * m.Nx + i
    velocity_mask = mask === nothing ? nothing : three ? (i, j, k) -> Tuple(mask[3node(i, j, k)+c] for c in 1:3) :
        (i, k) -> (mask[2node(i, k)+1], mask[2node(i, k)+2])
    velocity_changes = on("variable_bcv") ? read_scale_bcv(joinpath(dir, "scale_bcv.txt")) : Tuple{Float64, Float64}[]
    V0 = on("velocity_from_ascii") ? read_ascii_field(joinpath(dir, "input_velocity_0.txt")) : zeros((three ? 3m.Ny : 2) * m.Nx * m.Nz)
    velocity_files = on("multi_velocity") ? read_multi_velocity(joinpath(dir, "multi_veloc.txt")) : Float64[]   # mv_time
    velocity_file = 0                                                # i_veloc - 1
    random_strain = num("random_initial_strain", 0.0)
    initial_strain(x, z, l) = random_strain * (unit_random(m.random_seed, hash((x, z)), 3) + 0.5)
    location = m.periodic_boundary ? :shim : :native                 # the periodic DMDA has a different node spacing
    s = three ? Model3D(m, layers; temperature = (i, j, k) -> T0[node(i, j, k)+1], velocity = (i, j, k) -> Tuple(V0[3node(i, j, k)+c] for c in 1:3),
            velocity_mask, initial_strain = (x, y, z, l) -> random_strain * (unit_random(m.random_seed, hash((x, y, z)), 3) + 0.5)) :
        Model(m, layers; temperature = (i, k) -> T0[node(i, k)+1], velocity = (i, k) -> (V0[2node(i, k)+1], V0[2node(i, k)+2]), location, initial_strain, velocity_mask)
    if three && (on("export_lithology") || on("print_step_files"))
        say("warning: export_lithology and print_step_files are not available in 3D (ignored)")
        d["export_lithology"] = d["print_step_files"] = "False"
    end
    surface = Surface(d, m, layers, dir)
    surface === nothing || say("surface: " * join(size(surface.z), " x ") * " points" * (surface.processes ? ", sp_mode = $(surface.mode)" : ""))
    extraction = m.use_magmatism && m.use_melt_extraction && surface !== nothing ? MeltExtraction() : nothing
    m.use_magmatism && m.use_melt_extraction && surface === nothing &&
        say("magmatism_extraction = on needs sp_surface_tracking = True: the melt is not extracted (as in the C++ code)")
    if rank == 0 && abspath(out) != abspath(dir)                     # tapIOca reads param.txt from the output folder
        mkpath(out); cp(joinpath(dir, "param.txt"), joinpath(out, "param.txt"); force = true)
    end
    rank == 0 && on("nondimensionalization") && write_dimensional_scales(nondimensional_scales(d, m), (dir, out), say)
    MPI.Barrier(comm)
    w = OutputWriter(out; export_thermal_diffusivity = on("export_thermal_diffusivity"), export_lithology = on("export_lithology"),
        write_particles = on("print_step_files"), magmatism = m.use_magmatism)
    step_max, time_max = Int(num("step_max", 1e9)), num("time_max", 1e30) * SECONDS_PER_YEAR
    step_print = Int(num("step_print", 10)); initial_print_step = Int(num("initial_print_step", 0))
    initial_print_max_time = num("initial_print_max_time", 1e6) * SECONDS_PER_YEAR
    start = MPI.Wtime()
    warned_above_top = false
    initial_velocity!(s)
    say(@sprintf("step 0: %d nonlinear iterations, %.1f s", s.nonlinear_iterations, MPI.Wtime() - start))
    write_output!(w, s, 0; surface)
    for step in 1:step_max
        s.time + s.time_step_seconds <= time_max || break
        dt = s.time_step_seconds
        # as in the C++ code, the time compared is that of the end of the step, and one change is applied per step
        if !isempty(velocity_changes) && s.time + dt > velocity_changes[1][1] * 1e6 * SECONDS_PER_YEAR
            (t, scale) = popfirst!(velocity_changes)
            rescale_velocity!(s, scale)
            say("velocity multiplied by $scale ($t Myr)")
        end
        # multi_veloc_change, after rescaleVeloc: the next input_velocity_<n>.txt (as in the C++ code, it replaces
        # the boundary values, including a factor of scale_bcv.txt applied before)
        if !isempty(velocity_files) && s.time + dt > velocity_files[1] * 1e6 * SECONDS_PER_YEAR
            t = popfirst!(velocity_files)
            field = nothing
            if on("velocity_from_ascii")
                velocity_file += 1
                say("Reading input_velocity_$velocity_file.txt ($t Myr)")
                Vn = read_ascii_field(joinpath(dir, "input_velocity_$velocity_file.txt"))
                field = three ? (i, j, k) -> Tuple(Vn[3node(i, j, k)+c] for c in 1:3) : (i, k) -> (Vn[2node(i, k)+1], Vn[2node(i, k)+2])
            end
            replace_velocity!(s, field)
        end
        surface === nothing || update_surface_conditions!(surface, (s.time + dt) / SECONDS_PER_YEAR, say)
        step!(s; surface)
        added, removed = three ? add_remove_particles3!(m, s.particle_swarm, step) : add_remove_particles!(m, s.particle_swarm, s.da_Veloc, step; location)
        surface === nothing || surface_step!(surface, m, s.layers, s.particle_swarm, dt)
        every = initial_print_step > 0 && s.time < initial_print_max_time ? initial_print_step : step_print
        if step % every == 0
            # as in the C++ code, the melt is extracted only at the output steps
            extraction === nothing || extract_melt!(extraction, s, surface, say)
        end
        if surface !== nothing && !warned_above_top && maximum(surface.z) > 0
            # the C++ code drops the points of its surface swarm that leave the mesh; here they are kept
            say(@sprintf("warning: the surface is above the top of the model (z = %.4g m)", maximum(surface.z))); warned_above_top = true
        end
        step % every == 0 && write_output!(w, s, step; surface)
        n = MPI.Allreduce(Int(LibPETSc.DMSwarmGetLocalSize(P, s.particle_swarm)), +, comm)
        say(@sprintf("step %d: t = %.4g yr, dt = %.4g yr, %d nonlinear iterations, %d particles (+%d -%d)", step, s.time / SECONDS_PER_YEAR,
            dt / SECONDS_PER_YEAR, s.nonlinear_iterations, n, added, removed))
    end
    say(@sprintf("done: %.1f s, %.2f s of them writing the outputs (%d times)", MPI.Wtime() - start, w.seconds, length(w.times)))
    PETSc.finalize(P)
end
