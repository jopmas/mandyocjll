# Builds bin/mandyocjll.so (or .dylib): Julia with mandyocjll already compiled. Run by `make`.
using PackageCompiler
include(joinpath(@__DIR__, "..", "bin", "image.jl"))

# small scenario, run with 1 and 2 ranks to record what has to be compiled
function write_scenario(folder)
    nx, nz, lx, lz = 21, 11, 400e3, 200e3
    write(joinpath(folder, "param.txt"), """
        nx = $nx
        nz = $nz
        lx = $lx
        lz = $lz
        particles_per_element = 20
        step_max = 3
        step_print = 1
        time_max = 1.0e9
        dt_max = 1.0e5
        viscosity_max = 1.0e24
        viscosity_min = 1.0e19
        viscosity_per_element = constant
        viscosity_dependence = pressure
        viscosity_mean_method = harmonic
        non_linear_method = on
        plasticity = on
        rheology_model = 9
        initial_dynamic_range = True
        sticky_blanket_air = True
        Xi_min = 1.0e-3
        denok = 1.0e-13
        rtol = 1.0e-7
        geoq = on
        interfaces_from_ascii = True
        temperature_from_ascii = True
        velocity_from_ascii = True
        export_lithology = True
        print_step_files = True
        radiogenic_component = on
        adiabatic_component = on
        random_initial_strain = 0.1
        """)
    open(joinpath(folder, "interfaces.txt"), "w") do io
        println(io, "C 1.0 1.0 1.0\nrho 3378.0 2800.0 1.0\nH 0.0 2.0e-10 0.0\nA 1.393e-14 8.574e-28 1.0e-18\nn 3.0 4.0 1.0\nQ 429.0e3 222.0e3 0.0\nV 15.0e-6 0.0 0.0")
        foreach(_ -> println(io, "-80000.0 -40000.0"), 1:nx)
    end
    header = "# 1\n# 2\n# 3\n# 4\n"
    open(joinpath(folder, "input_temperature_0.txt"), "w") do io
        print(io, header)
        for k in 0:nz-1, i in 0:nx-1
            println(io, clamp(1300 * ((nz - 1 - k) * lz / (nz - 1) - 40e3) / 120e3, 0, 1300))
        end
    end
    open(joinpath(folder, "input_velocity_0.txt"), "w") do io
        print(io, header)
        for k in 0:nz-1, i in 0:nx-1
            println(io, 3e-10 * (i / (nx - 1) - 0.5) * (k == 0 || k == nz - 1 ? 1 : 1), "\n", k == 0 ? 3e-10 * lz / lx : 0.0)
        end
    end
end

work = mktempdir()
write_scenario(work)
rm(STAMP; force = true)                  # the launcher must not use an old image here
for ranks in (1, 2)
    @info "Recording the precompilation statements with $ranks rank(s)"
    run(addenv(`$(Base.julia_cmd()) --project=$ROOT --startup-file=no $(joinpath(ROOT, "bin", "mandyocjll")) -n $ranks -i $work -o $(joinpath(work, "out$ranks")) -seed 0 -strain_seed 0.1`,
        "MANDYOCJLL_TRACE" => joinpath(work, "trace")))
end
traces = filter(f -> startswith(basename(f), "trace"), readdir(work; join = true))
temporary = joinpath(ROOT, "bin", "new_image." * (Sys.isapple() ? "dylib" : "so"))
create_sysimage(["Mandyocjll", "MPI"]; project = ROOT, sysimage_path = temporary, precompile_statements_file = traces,
    include_transitive_dependencies = false)   # the JLLs would load their own MPI, METIS... besides the ones in use
mv(temporary, IMAGE; force = true)
write(STAMP, joinpath(Sys.BINDIR, "julia"), "\n", image_checksum())
@info "Built $IMAGE"
