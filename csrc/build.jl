# Compila o shim contra os headers do PETSc_jll. Uso: julia --project=. csrc/build.jl
using PETSc, MPI
loaded(name) = only(m for (k, m) in Base.loaded_modules if k.name == name)
petsc_include = joinpath(loaded("PETSc_jll").artifact_dir, "lib", "petsc", "double_real_Int64", "include")
mpi_include = joinpath(loaded("MPICH_jll").artifact_dir, "include")
shim_library = joinpath(@__DIR__, "libmandyoc_shim.$(Sys.isapple() ? "dylib" : "so")")
undefined_flags = Sys.isapple() ? `-undefined dynamic_lookup` : ``
run(`cc -O2 -shared -fPIC $undefined_flags -I$petsc_include -I$mpi_include -o $shim_library $(joinpath(@__DIR__, "cell_dm_shim.c"))`)
println(shim_library)
