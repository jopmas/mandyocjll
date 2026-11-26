# Writes LocalPreferences.toml so that mandyocjll uses a PETSc and the MPI it was built with:
#
#     julia tools/use_petsc.jl PETSC_DIR [PETSC_ARCH] [OUTPUT]     this PETSc (PETSC_ARCH "" for a prefix install)
#     julia tools/use_petsc.jl --auto [OUTPUT]                      the first compatible one of:
#         1. PETSC_DIR and PETSC_ARCH of the environment;
#         2. Linux (clusters): a PETSc built from source in ~/opt/petsc (each folder with lib/libpetsc, newest
#            first), then the PETSc of the conda environment tapioca (TAPIOCA_ENV overrides the name);
#            macOS: the conda environment tapioca first, since PETSc built from source on the Mac often fails.
#     If none is compatible, nothing is written and PETSc_jll and MPICH_jll are used.
#
# `make` runs it (with --auto when LocalPreferences.toml does not exist; with PETSC_DIR=... PETSC_ARCH=... for a given
# PETSc), then rebuilds the shim and the image, which depend on it. A PETSc is compatible when it has real double
# scalars, MUMPS and OpenMPI or MPICH, libpetsc loads, and libmpi and mpiexec are found.

using Libdl

fail(message) = (println(stderr, "use_petsc: ", message); exit(1))
const EXT = Sys.isapple() ? "dylib" : "so"

"Settings of the PETSc in prefix (PETSC_DIR/PETSC_ARCH), or the reason why it can not be used."
function inspect(prefix)
    lib = joinpath(prefix, "lib", "libpetsc.$EXT")
    isfile(lib) || return "no lib/libpetsc.$EXT"
    conf = joinpath(prefix, "include", "petscconf.h")
    isfile(conf) || return "no include/petscconf.h"
    defines = Set(m[1] for m in eachmatch(r"(?m)^#define (\w+) 1$", read(conf, String)))
    has(name) = name in defines
    has("PETSC_USE_REAL_DOUBLE") || return "not double precision"
    has("PETSC_USE_COMPLEX") && return "complex scalars"
    has("PETSC_HAVE_MUMPS") || return "no MUMPS"
    abi = has("PETSC_HAVE_OPENMPI") ? "OpenMPI" : has("PETSC_HAVE_MPICH") ? "MPICH" : return "MPI neither OpenMPI nor MPICH"
    libmpi = find_libmpi(lib, prefix)
    libmpi === nothing && return "libmpi of libpetsc not found"
    candidates = [joinpath(dirname(dirname(libmpi)), "bin", "mpiexec"), joinpath(prefix, "bin", "mpiexec")]
    i = findfirst(isfile, candidates)
    i === nothing && return "mpiexec not found"
    handle = Libdl.dlopen(realpath(lib); throw_error = false)
    handle === nothing && return "libpetsc does not load"
    Libdl.dlclose(handle)
    return (; prefix, lib = realpath(lib), abi, libmpi, mpiexec = candidates[i],
        int = has("PETSC_USE_64BIT_INDICES") ? "Int64" : "Int32")
end

"libmpi that libpetsc is linked to (otool on macOS, ldd on Linux)."
function find_libmpi(lib, prefix)
    text = try
        Sys.isapple() ? read(`otool -L $lib`, String) : read(`ldd $lib`, String)
    catch
        return nothing
    end
    for line in split(text, '\n')
        m = match(r"(\S*/libmpi\.[0-9.]*(?:so|dylib)[0-9.]*)", line)
        m === nothing && continue
        path = startswith(m[1], "@rpath/") ? joinpath(prefix, "lib", basename(m[1])) : m[1]  # conda: rpath is lib/
        isfile(path) && return realpath(path)
    end
    return nothing
end

"Folders of the conda environment `name` (conda env list, then the usual places)."
function conda_prefixes(name)
    found = String[]
    conda = Sys.which("conda")
    if conda !== nothing
        try
            for line in eachline(`$conda env list`)
                t = split(line)
                !isempty(t) && t[1] == name && push!(found, t[end])
            end
        catch
        end
    end
    home = homedir()
    append!(found, joinpath.([joinpath(home, "miniconda3"), joinpath(home, "miniforge3"), joinpath(home, "anaconda3"),
        "/opt/miniconda3", "/opt/homebrew/Caskroom/miniforge/base"], "envs", name))
    return unique(filter(isdir, found))
end

"Candidates in the order of priority: (description, prefix)."
function candidates()
    list = Tuple{String, String}[]
    if haskey(ENV, "PETSC_DIR")
        push!(list, ("PETSC_DIR/PETSC_ARCH of the environment", joinpath(ENV["PETSC_DIR"], get(ENV, "PETSC_ARCH", ""))))
    end
    built_list = Tuple{String, String}[]
    built = joinpath(homedir(), "opt", "petsc")
    if isdir(built)
        archs = filter(d -> isfile(joinpath(built, d, "lib", "libpetsc.$EXT")), readdir(built))
        sort!(archs; by = d -> mtime(joinpath(built, d, "lib", "libpetsc.$EXT")), rev = true)
        built_list = [("PETSc built in ~/opt/petsc", joinpath(built, d)) for d in archs]
    end
    env = get(ENV, "TAPIOCA_ENV", "tapioca")
    conda_list = [("PETSc of the conda environment $env", p) for p in conda_prefixes(env)]
    return vcat(list, Sys.isapple() ? [conda_list; built_list] : [built_list; conda_list])
end

function write_preferences(output, p)
    write(output, """
    # Written by tools/use_petsc.jl for $(p.prefix)
    # Specific to this computer; not in git. Delete it to go back to the PETSc_jll and MPICH_jll packages.
    [MPIPreferences]
    _format = "1.0"
    abi = "$(p.abi)"
    binary = "system"
    cclibs = []
    libmpi = "$(p.libmpi)"
    mpiexec = "$(p.mpiexec)"
    preloads = []
    preloads_env_switch = ""

    [PETSc]
    library_path = "$(p.lib)"
    PetscScalar = "Float64"
    PetscInt = "$(p.int)"
    """)
    println("$output: PETSc of $(p.prefix) ($(p.int)), $(p.abi) ($(p.libmpi))")
end

const DEFAULT_OUTPUT = joinpath(dirname(@__DIR__), "LocalPreferences.toml")
function main()
if get(ARGS, 1, "") == "--auto"
    length(ARGS) <= 2 || fail("usage: julia tools/use_petsc.jl --auto [OUTPUT]")
    output = length(ARGS) == 2 ? abspath(ARGS[2]) : DEFAULT_OUTPUT
    for (description, prefix) in candidates()
        p = inspect(prefix)
        if p isa String
            println("skipped: $description, $prefix: $p")
        else
            println("using: $description, $prefix")
            write_preferences(output, p)
            exit(0)
        end
    end
    println("no compatible PETSc found: using the PETSc_jll and MPICH_jll packages (slow in parallel on macOS)")
    isfile(output) && rm(output)
else
    length(ARGS) in 1:3 || fail("usage: julia tools/use_petsc.jl PETSC_DIR [PETSC_ARCH] [OUTPUT] | --auto [OUTPUT]")
    prefix = joinpath(abspath(expanduser(ARGS[1])), get(ARGS, 2, ""))
    p = inspect(prefix)
    p isa String && fail("the PETSc of $prefix can not be used: $p")
    write_preferences(length(ARGS) == 3 ? abspath(ARGS[3]) : DEFAULT_OUTPUT, p)
end
end

main()
