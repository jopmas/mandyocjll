# Tutorial: from a scenario to the figures

The whole path of a scenario: build the inputs, run mandyocjll and process the outputs. The commands below were run
for the handbook (`docs/apostila`), on a Mac mini with the conda environment `tapioca`.

## 1. Build a scenario with tapIOca

The simplest way is to draw the scenario in an SVG and let `tapioca-svg` write the inputs and the run script. The
examples of Mandyoc are already drawn in `examples/examples_svg`. The continental rift, with a 10 km mesh and 20 Myr:

```bash
conda activate tapioca
cd ~/opt/mandyocjll/examples/examples_svg/continental_rift
tapioca-svg continental_rift.svg -o inputs --dx 10 --dz 10 --time-max 20 \
    --code mandyocjll --system mac --ncores 4
```

`--code mandyocjll` (or `code = mandyocjll` in the `setup` shape of the SVG) makes the run script call mandyocjll;
`--system` chooses `run_mac.sh`, `run_linux.sh` or `run_hypatia.sh`. The folder `inputs` gets `param.txt`,
`interfaces.txt`, `input_temperature_0.txt`, `input_velocity_0.txt`, figures of the scenario
(`numerical_setup.png`, `initial_temperature.png`) and the run script. Edit `param.txt` before running if needed
(for example `step_print`).

## 2. Build a scenario without tapIOca

The inputs are plain text files; each example in `examples/` has a `make_inputs.py` that writes them
(`examples/continental_rift/make_inputs.py` is a good template):

- `interfaces.txt`: one line per property of the layers (`C`, `rho`, `H`, `A`, `n`, `Q`, `V` and, optional, `k`,
  `cohesion_min`, `cohesion_max`, `friction_angle_min`, `friction_angle_max`, `weakening_seed`), with one value per
  layer from the deepest to the shallowest; then one line per node of the top with the depth (m, negative) of each
  interface. In 3D there is one line per node `(i, j)`, with `i` (x) varying fastest.
- `input_temperature_0.txt`, `input_velocity_0.txt`: 4 header lines, then one value per node (the velocity: one
  per component), with x varying fastest, then y (3D), then z from the bottom to the top.
- `param.txt`: the keys of Mandyoc; [parameters.md](parameters.md) lists those read by mandyocjll.

## 3. Run

```bash
cd inputs
bash run_mac.sh
```

The script runs mandyocjll (the log goes to `FD.out`) and, if the run ends without error, creates the `.nc` files
(`tools/zarr2nc.py`), the snapshots and the video (scripts of `post_processing`) and zips the scenario. Without the
script, run mandyocjll in the folder of the scenario:

```bash
~/opt/mandyocjll/bin/mandyocjll -n 4                       # add -seed/-strain_seed if the scenario needs them
```

The rift above (161 × 31 nodes, 2,000 steps of 10 kyr) took 211 s with 4 ranks. On Hypatia, submit the script
written with `--system hypatia` (`sbatch run_hypatia.sh`).

## 4. Process the outputs

The zarr stores can be read with xarray, also while the run goes on:

```python
import xarray as xr
import matplotlib.pyplot as plt

T = xr.open_zarr("temperature.zarr").temperature        # dimensions (time, z, x)
v = xr.open_zarr("velocity.zarr")                       # vx, vz
lith = xr.open_zarr("lithology.zarr").lithology         # mesh 5 times finer (export_lithology = True)

print(T.time.values)                                    # times of the outputs (Myr)
T.sel(time=10, method="nearest").plot(cmap="inferno")   # temperature near 10 Myr
plt.savefig("temperature_10Myr.png")

profile = T.isel(time=-1).sel(x=800e3, method="nearest") # vertical profile at the center, last output
vz_top = v.vz.isel(z=-1)                                # vertical velocity at the top, all outputs
lith.isel(time=-1).plot(cmap="tab20")                   # lithology of the last output
```

To use the tools of tapIOca, which read NetCDF:

```bash
python ~/opt/mandyocjll/tools/zarr2nc.py .               # all the variables
python ~/opt/mandyocjll/tools/zarr2nc.py . temperature   # only the temperature
```

[outputs.md](outputs.md) describes every store.

## 5. A 3D scenario

A 3D scenario has `dimensions = 3` and the keys `ny` and `ly` in `param.txt`; the faces in y have the keys `front_*`
and `back_*`. The example `examples/continental_rift_non_dimensional/model2` is a 3D rift:

```bash
cd ~/opt/mandyocjll/examples/continental_rift_non_dimensional/model2
python make_inputs.py
~/opt/mandyocjll/bin/mandyocjll -n 4 -seed 0,2 -strain_seed 0.0,1.0
```

The outputs have the extra dimension `y` (`T.sel(y=30e3, method="nearest")` gives a section). Magmatism, the surface
tracking and the periodic boundary are not available in 3D, as in the C++ code.
