# Development

## Layout

| Path | |
|---|---|
| `src/Mandyocjll.jl` | the module |
| `src/parameters.jl` | `ModelParameters` and `LayerProperties` (all the settings of a run) |
| `src/reader.jl` | `param.txt`, `interfaces.txt` and the ASCII fields |
| `src/particle_swarm.jl` | 2D DMDAs, DMSwarm, particle location (shim or native) and advection |
| `src/particle_control.jl` | adding and removing particles (2D) |
| `src/thermal.jl`, `src/thermal3d.jl` | heat equation in 2D and 3D (the 3D file also has the 3D DMDAs) |
| `src/stokes.jl`, `src/stokes3d.jl` | Stokes (Q1-P0, FSSA) in 2D and 3D; the Uzawa iterations (`uzawa!`) are shared |
| `src/rheology.jl` | viscosity, plasticity, strain, particles to mesh, pressure, nonlinear loop, viscosity ramp, the `Model` |
| `src/rheology3d.jl` | the same in 3D: particles, particles to mesh, pressure, addition and removal, `Model3D` |
| `src/surface.jl` | surface tracking and surface processes (2D) |
| `src/magmatism.jl` | melt extraction and extrusion (the melting is in `src/rheology.jl`) |
| `src/output.jl` | zarr writer |
| `src/main.jl` | `run_scenario`: reads a scenario and runs it |
| `bin/mandyocjll` | command line, MPI launch, compiled image |
| `bin/image.jl`, `build/` | compiled image (PackageCompiler) |
| `csrc/` | C shim that fills `cell_dm->ops` for the DMSwarm location (option `shim`) |
| `tools/zarr2nc.py`, `tools/use_petsc.jl` | zarr to NetCDF; choice of PETSc and MPI |
| `docs/apostila/` | handbook (Portuguese): `build_apostila.py` builds the PDF, `scripts/make_figures.py` the figures |
| `test/` | tests and comparisons with the C++ code |
| `examples/` | examples of Mandyoc, also as SVG scenarios |

## Conventions

- Names follow `docs/rename_table.csv`: each variable or function that comes from Mandyoc has a descriptive name and
  the old one in a comment (`particle_swarm # dms`).
- Comments and printed messages are in English.
- An option of `param.txt` that changes the result and has not been ported must stop the run (`check_unsupported`
  in `src/main.jl`), never be ignored.
- A bug of the C++ code is fixed here, documented in `docs/bugs_mandyoc_cpp.md` (with a suggested fix for the C++
  code) and, when useful, reproducible with an option (`fssa_legacy`, `cpp_order`). Rebuild the PDF with
  `python3 docs/build_pdf.py` (needs xelatex).

## After a change

1. Run the tests of [validation.md](validation.md) that cover the change.
2. Run an example (`examples/*/run.sh`) with 1 and 4 ranks; the results must not depend on the number of ranks.
3. `make image`, or `mandyocjll` will run without the image (with a warning).

PETSc is called through the low-level `LibPETSc` interface of PETSc.jl 0.5.6. Two pitfalls found: `MatStencil` is in
the wrong order in that version (use `MatSetValuesLocal`), and `DMLOCATEPOINT_POINT_NOT_FOUND` is -367.

`test/bench_stokes.jl` measures the time of each stage of a Stokes solve.
