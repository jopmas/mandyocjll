# Builds the shim against the PETSc in use (PETSc_jll or the one in LocalPreferences.toml).
# Usage: julia --project=. csrc/build.jl   (run again after changing PETSc)
using PETSc, MPI
info = PETSc.library_info()
if info.source == :preferences                      # <prefix>/lib/libpetsc*, headers in <prefix>/include
    prefix = dirname(dirname(info.path))
    includes = [joinpath(prefix, "include")]
else
    loaded(name) = only(m for (k, m) in Base.loaded_modules if k.name == name)
    includes = [joinpath(loaded("PETSc_jll").artifact_dir, "lib", "petsc", "double_real_Int64", "include"),
        joinpath(loaded("MPICH_jll").artifact_dir, "include")]
end
shim_library = joinpath(@__DIR__, "libmandyoc_shim.$(Sys.isapple() ? "dylib" : "so")")
undefined_flags = Sys.isapple() ? `-undefined dynamic_lookup` : ``
run(`cc -O2 -shared -fPIC $undefined_flags $("-I" .* includes) -o $shim_library $(joinpath(@__DIR__, "cell_dm_shim.c"))`)
println(shim_library)
