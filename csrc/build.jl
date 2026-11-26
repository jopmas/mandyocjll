# Builds the shim against the PETSc in use (PETSc_jll or the one in LocalPreferences.toml).
# Usage: julia --project=. csrc/build.jl   (run again after changing PETSc)
using PETSc, MPI
info = PETSc.library_info()
if info.source == :preferences                      # <prefix>/lib/libpetsc*, headers in <prefix>/include
    prefix = dirname(dirname(info.path))
    # a PETSc built from source has the headers of the arch in PETSC_DIR/PETSC_ARCH/include and the others in
    # PETSC_DIR/include; mpi.h is next to the libmpi in use (e.g. a module of the cluster)
    libmpi = MPI.API.libmpi
    mpi_include = isabspath(libmpi) ? joinpath(dirname(dirname(realpath(libmpi))), "include") : ""
    includes = filter(isdir, unique([joinpath(prefix, "include"), joinpath(dirname(prefix), "include"), mpi_include]))
else
    loaded(name) = only(m for (k, m) in Base.loaded_modules if k.name == name)
    includes = [joinpath(loaded("PETSc_jll").artifact_dir, "lib", "petsc", "double_real_Int64", "include"),
        joinpath(loaded("MPICH_jll").artifact_dir, "include")]
end
shim_library = joinpath(@__DIR__, "libmandyoc_shim.$(Sys.isapple() ? "dylib" : "so")")
undefined_flags = Sys.isapple() ? `-undefined dynamic_lookup` : ``
run(`cc -O2 -shared -fPIC $undefined_flags $("-I" .* includes) -o $shim_library $(joinpath(@__DIR__, "cell_dm_shim.c"))`)
println(shim_library)
