# Running a scenario

`bin/mandyocjll` is the executable. Run it in the folder of the scenario, as the C++ `mandyoc`:

```bash
~/opt/mandyocjll/bin/mandyocjll -n 4 -seed 0,5,8 -strain_seed 0.0,1.0,1.0
julia ~/opt/mandyocjll/bin/mandyocjll -n 4 -seed 0,5,8 -strain_seed 0.0,1.0,1.0
mpiexec -n 4 ~/opt/mandyocjll/bin/mandyocjll -seed 0,5,8 -strain_seed 0.0,1.0,1.0
```

The three forms are equivalent. With `-n`, mandyocjll starts the ranks itself with the `mpiexec` of the MPI it uses;
with `mpiexec`, it must be the `mpiexec` of that MPI (`LocalPreferences.toml`), and mandyocjll stops with a message if
it is not. To call it as `mandyocjll`, link it into a folder of the `PATH`:

```bash
ln -s ~/opt/mandyocjll/bin/mandyocjll ~/.local/bin/
```

## Options

| Option | |
|---|---|
| `-n RANKS` | number of MPI ranks (default: 1); not with `mpiexec` |
| `-i FOLDER` | folder of the scenario (default: current folder) |
| `-o FOLDER` | folder of the outputs (default: folder of the scenario); `param.txt` is copied to it |
| `-seed L1,L2,...` | layers (0 = bottom layer) with initial strain |
| `-strain_seed S1,S2,...` | initial strain of these layers (default: 0.5 for each) |
| `-h_air H` | thickness of the air layer [m], needed by the depth dependent rheology (`viscosity_dependence = depth`) |
| `-Px N`, `-Pz N` | number of ranks in x and z (default: chosen by PETSc) |
| `-h`, `--help` | help |

Any other option goes to PETSc (e.g. `-log_view`, `-veloc_ksp_view`). `-strain_seed_constant`, accepted by the C++
code where it has no effect, is ignored. A layer with `weakening_seed` in `interfaces.txt` uses that value instead of
`-strain_seed`. A layer of `-seed` that does not exist stops the run (the C++ code writes out of the array, see the
bugs document).

## What is printed

```
mandyocjll: 4 MPI rank(s), mesh 401 x 76, 36 particles per element
step 0: 8 nonlinear iterations, 7.1 s
step 1: t = 1e+04 yr, dt = 1e+04 yr, 9 nonlinear iterations, 1080000 particles (+0 -0)
velocity multiplied by -1.0 (20.0 Myr)
...
done: 478.3 s, 0.90 s of them writing the outputs (4 times)
```

For each step: the time at its end, the time step, the nonlinear iterations of Stokes, the number of particles and
those added and removed by the particle control. Warnings are printed when Uzawa or the nonlinear iterations stop
before the tolerance. Errors (missing files, options not ported, an output folder that already has zarr outputs) stop
all the ranks with a message and exit code 1.

## Stopping conditions and outputs

The run stops at `step_max` or when the next step would pass `time_max`. Outputs are written at step 0 and every
`step_print` steps (every `initial_print_step` steps while the time is below `initial_print_max_time`). mandyocjll
refuses to write to a folder that already has its zarr outputs: move or delete them first.

## With tapIOca

`tapioca-svg` writes the run scripts of a scenario drawn in an SVG. With `code = mandyocjll` in the `setup` shape (or
`--code mandyocjll`), `run_linux.sh`, `run_mac.sh` and `run_hypatia.sh` call mandyocjll:

```bash
tapioca-svg scenario.svg -o scenario --code mandyocjll --ncores 8 --mandyoc-options="-seed 0,2 -strain_seed 0.0,1.0"
cd scenario
bash run_mac.sh            # or run_linux.sh; on Hypatia: sbatch run_hypatia.sh
```

The script runs `$MANDYOCJLL -n <ncores>` (default `$HOME/opt/mandyocjll/bin/mandyocjll`), stops if the run fails and
then creates the `.nc` files with `tools/zarr2nc.py`, the snapshots and the video (scripts of post_processing) and
the zip of the scenario. tapioca-svg warns when the scenario uses options that mandyocjll does not have yet. This
needs the branch `feat/svg-scenario-inputs` of tapIOca. On Hypatia, mandyocjll must be installed and built (`make`)
on the cluster; that script has not been tested there yet.

## Running the examples

See [`../examples/README.md`](../examples/README.md) and [`../examples/examples_svg/README.md`](../examples/examples_svg/README.md).
