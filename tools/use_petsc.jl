# Writes LocalPreferences.toml so that mandyocjll uses a given PETSc and the MPI it was built with:
#
#     julia tools/use_petsc.jl PETSC_DIR [PETSC_ARCH] [OUTPUT]
#     make PETSC_DIR=$HOME/opt/petsc PETSC_ARCH=optimized-v3.24.1-openmpi
#
# PETSC_DIR/PETSC_ARCH is a PETSc built from source (PETSC_ARCH empty for a prefix install, e.g. a conda environment).
# The integer size, the scalar type and the MPI (OpenMPI or MPICH) come from petscconf.h; libmpi is the one that
# libpetsc is linked to, and mpiexec is the one next to it. Run `make` afterwards (the shim and the image depend on it).

const USAGE = "usage: julia tools/use_petsc.jl PETSC_DIR [PETSC_ARCH] [OUTPUT]"
fail(message) = (println(stderr, "use_petsc: ", message); exit(1))
length(ARGS) in 1:3 || fail(USAGE)
const DIR = abspath(expanduser(ARGS[1]))
const ARCH = get(ARGS, 2, "")
const OUTPUT = length(ARGS) == 3 ? abspath(ARGS[3]) : joinpath(dirname(@__DIR__), "LocalPreferences.toml")
const PREFIX = isempty(ARCH) ? DIR : joinpath(DIR, ARCH)

const EXT = Sys.isapple() ? "dylib" : "so"
const LIBPETSC = joinpath(PREFIX, "lib", "libpetsc.$EXT")
isfile(LIBPETSC) || fail("$LIBPETSC not found (check PETSC_DIR and PETSC_ARCH)")
const CONF = joinpath(PREFIX, "include", "petscconf.h")
isfile(CONF) || fail("$CONF not found")
const DEFINES = Set(m[1] for m in eachmatch(r"(?m)^#define (\w+) 1$", read(CONF, String)))
has(name) = name in DEFINES

has("PETSC_USE_REAL_DOUBLE") || fail("mandyocjll needs a PETSc with double precision")
const SCALAR = has("PETSC_USE_COMPLEX") ? "ComplexF64" : "Float64"
SCALAR == "Float64" || fail("mandyocjll needs a PETSc with real scalars (this one is complex)")
const INT = has("PETSC_USE_64BIT_INDICES") ? "Int64" : "Int32"
has("PETSC_HAVE_MUMPS") || println(stderr, "use_petsc: warning: this PETSc has no MUMPS; set solver = iterative in param.txt")
const ABI = has("PETSC_HAVE_OPENMPI") ? "OpenMPI" : has("PETSC_HAVE_MPICH") ? "MPICH" :
    fail("the MPI of this PETSc is neither OpenMPI nor MPICH")

"libmpi that libpetsc is linked to (otool on macOS, ldd on Linux)."
function find_libmpi()
    text = Sys.isapple() ? read(`otool -L $LIBPETSC`, String) : read(`ldd $LIBPETSC`, String)
    for line in split(text, '\n')
        m = match(r"(\S*/libmpi\.[0-9.]*(?:so|dylib)[0-9.]*)", line)
        m === nothing && continue
        path = m[1]
        if startswith(path, "@rpath/")      # macOS: the rpath of a conda environment is its lib folder
            path = joinpath(PREFIX, "lib", basename(path))
        end
        isfile(path) && return realpath(path)
    end
    fail("libpetsc is not linked to libmpi (or it was not found): $LIBPETSC")
end
const LIBMPI = find_libmpi()
const MPIEXEC = let candidates = [joinpath(dirname(dirname(LIBMPI)), "bin", "mpiexec"), joinpath(PREFIX, "bin", "mpiexec")]
    i = findfirst(isfile, candidates)
    i === nothing ? fail("mpiexec not found in $(join(dirname.(candidates), " or "))") : candidates[i]
end

write(OUTPUT, """
# Written by tools/use_petsc.jl for $PREFIX
# Specific to this computer; not in git. Delete it to go back to the PETSc_jll and MPICH_jll packages.
[MPIPreferences]
_format = "1.0"
abi = "$ABI"
binary = "system"
cclibs = []
libmpi = "$LIBMPI"
mpiexec = "$MPIEXEC"
preloads = []
preloads_env_switch = ""

[PETSc]
library_path = "$(realpath(LIBPETSC))"
PetscScalar = "$SCALAR"
PetscInt = "$INT"
""")
println("$OUTPUT: PETSc of $PREFIX ($SCALAR, $INT), $ABI ($LIBMPI)")
println("Now run `make` (or `make shim image`).")
