<div class="banner">
        <h1>Mandyocjll</h1>
        <p> <strong>MANtle DYnamics simulatOr Code in Julia</strong></p>
    </div>

## About

Mandyocjll is a port of [Mandyoc](https://github.com/ggciag/mandyoc) to Julia: a finite element code, written on top
of the PETSc library (through PETSc.jl), to simulate the thermomechanical evolution of the lithosphere and the
Earth's mantle in 2D and 3D. It reads the same input files as Mandyoc (`param.txt`, `interfaces.txt`,
`input_*_0.txt`), solves the same equations with the same numerical methods, and writes the outputs as zarr stores,
with the names and units of the NetCDF files of [tapIOca](https://github.com/jopmas/tapioca).

- Particles in cell (PETSc DMSwarm) on a regular mesh (DMDA), split among MPI ranks.
- Stokes flow with Q1-P0 elements and Uzawa iterations, free surface stabilization (FSSA), nonlinear (Picard)
  iterations and the viscosity ramp of Gerya.
- 13 rheologies, Drucker-Prager plasticity with strain softening, strain healing.
- Heat equation with radiogenic, adiabatic and shear heating, diffusivity per layer and basal heat flow.
- Surface tracking and surface processes (diffusion and three sedimentation modes), magmatism (melting, latent heat,
  extraction), time-dependent boundary velocities, nondimensional inputs.
- New: erosion and sedimentation with [FastScape](https://github.com/fastscape-lem/fastscapelib-fortran)
  (`sp_mode = fastscape`, in 2D and 3D) and the surface tracking in 3D.
- Results independent of the number of MPI ranks.

Each part was validated against the C++ code with the same inputs (see [docs/validation.md](docs/validation.md)).
The bugs found in the C++ code during the port are listed, with suggested fixes, in
[docs/bugs_mandyoc_cpp.pdf](docs/bugs_mandyoc_cpp.pdf) (in Portuguese).

## Quick start

```bash
git clone git@github.com:jopmas/mandyocjll.git ~/opt/mandyocjll
cd ~/opt/mandyocjll
make                                      # dependencies, PETSc/MPI choice, C shim and compiled image (a few minutes)

cd examples/vanKeken1997_case1a
bash run.sh 4                             # writes the inputs (Python) and runs with 4 MPI ranks
python ../plot_example.py .
```

Run `mandyocjll` in the folder of a scenario, as the C++ `mandyoc`:

```bash
~/opt/mandyocjll/bin/mandyocjll -n 4 -seed 0,5,8 -strain_seed 0.0,1.0,1.0
```

`mandyocjll -h` lists the options. Each variable is written to `<variable>.zarr`; `tools/zarr2nc.py <folder>`
converts them to the `.nc` files read by tapIOca. With tapIOca (branch `feat/svg-scenario-inputs`), `--code
mandyocjll` (or `code = mandyocjll` in the `setup` shape of the SVG) creates run scripts (`run_linux.sh`,
`run_mac.sh`, `run_hypatia.sh`) that call mandyocjll.

## Documentation

| | |
|---|---|
| [docs/](docs/README.md) | installation, running, parameters, outputs, numerical methods, validation and development |
| [docs/apostila/apostila_mandyocjll.pdf](docs/apostila/apostila_mandyocjll.pdf) | handbook (in Portuguese): the model, the code, a tutorial and the validation, with figures |
| [examples/](examples/README.md) | examples of Mandyoc; [examples/examples_svg/](examples/examples_svg/README.md): the same examples as SVG scenarios of tapIOca |
| [docs/rename_table.csv](docs/rename_table.csv) | names of Mandyoc and the descriptive names used here |

## Layout

| Path | |
|---|---|
| `src/parameters.jl`, `src/reader.jl` | settings of a run; reading of `param.txt`, `interfaces.txt` and the ASCII fields |
| `src/particle_swarm.jl`, `src/particle_control.jl` | meshes, particles, location, advection, addition and removal |
| `src/thermal.jl`, `src/thermal3d.jl` | heat equation in 2D and 3D |
| `src/stokes.jl`, `src/stokes3d.jl` | Stokes (Q1-P0, Uzawa, FSSA) in 2D and 3D |
| `src/rheology.jl`, `src/rheology3d.jl` | viscosity, particles to mesh, pressure, nonlinear loop and the coupled model |
| `src/surface.jl`, `src/fastscape.jl`, `src/magmatism.jl` | surface tracking and processes (2D and 3D); coupling with FastScape; melt extraction |
| `src/output.jl` | zarr writer |
| `src/main.jl`, `bin/mandyocjll` | `run_scenario` and the command line |
| `csrc/` | small C shim that connects the particle location in Julia to PETSc |
| `tools/` | `use_petsc.jl` (PETSc/MPI choice) and `zarr2nc.py` |
| `test/` | tests and comparisons with the C++ code |

## Tests

```bash
julia --project=. test/mpirun.jl 4 test_thermal.jl        # also test_stokes, test_thermal3d, test_stokes3d,
julia --project=. test/mpirun.jl 4 test_stokes3d.jl       # test_surface, test_open_boundary (shim), test_particle_control (native)
```

Each test prints `RESULT OK` or `RESULT FAIL`; the full list is in [docs/validation.md](docs/validation.md).

## Not ported yet

Restart, `T_initial` other than 0 and 1, and in 3D (as in the C++ code)
magmatism and the periodic boundary; in 3D the surface processes are only `diffusion` and `fastscape`.

## License

BSD 3-Clause, as Mandyoc (see [LICENSE](LICENSE)). The examples in `examples/` come from Mandyoc (see
`examples/LICENSE_MANDYOC`). `sp_mode = fastscape` calls FastScapeLib (`Fastscapelib_jll`), distributed under the
GNU GPL v3, which is installed as a dependency and is not part of this code.
