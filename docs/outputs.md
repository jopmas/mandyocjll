# Outputs

Each variable is one zarr store (format v2, zlib) in the output folder, with the names, dimensions, units and types
of the NetCDF files that tapIOca writes from the text outputs of Mandyoc (`meshes2netcdf.jl`, `particles2netcdf.jl`):

| Store | Variables | Units | Written |
|---|---|---|---|
| `density.zarr` | `density` | kg/m3 | always |
| `viscosity.zarr` | `viscosity` | Pa.s | always |
| `pressure.zarr` | `pressure` | Pa | always |
| `strain.zarr` | `strain` | | always |
| `strain_rate.zarr` | `strain_rate` | 1/s | always |
| `temperature.zarr` | `temperature` | deg C | always |
| `velocity.zarr` | `vx`, `vz` | m/s | always |
| `heat.zarr` | `heat` | W/kg | always |
| `thermal_diffusivity.zarr` | `thermal_diffusivity` | m2/s | `export_thermal_diffusivity = True` |
| `lithology.zarr` | `lithology` (layer ID, `int8`) on a mesh 5 times finer | | `export_lithology = True` |
| `surface.zarr` | `surface` (height of the surface at `x`, measured from the top of the model, negative below it) | m | `sp_surface_tracking = True` or `sp_surface_processes = True` |
| `particles_trajectories.zarr` | `x`, `z`, `lithology`, `strain` of the first particle of each element, by `id` | m | `print_step_files = True` |

The mesh stores have the coordinates `time` (Myr), `step`, `x` and `z` (m, `z` from `-lz` to 0) and the dimensions
`(time, z, x)`, one chunk per output. The trajectories have the dimensions `(time, id)`; removed particles have
`-999` (`lithology`: -1), as in `particles2netcdf.jl`. Each output appends to the stores, so they can be read while
the run goes on.

## Reading

```python
import xarray as xr

T = xr.open_zarr("temperature.zarr").temperature
T.sel(time=10, method="nearest").plot()

ds = xr.merge([xr.open_zarr(f"{v}.zarr") for v in ("temperature", "velocity")], compat="override")
```

`examples/plot_example.py FOLDER [STEP]` plots density, temperature, viscosity, strain rate and velocity of one output.

## NetCDF

The post-processing scripts of tapIOca and post_processing read `.nc` files today. `tools/zarr2nc.py` writes them
next to the stores, with the same contents; it needs only numpy, xarray and h5netcdf:

```bash
python ~/opt/mandyocjll/tools/zarr2nc.py OUTPUT_FOLDER            # all the stores
python ~/opt/mandyocjll/tools/zarr2nc.py OUTPUT_FOLDER temperature
```

The run scripts of tapIOca with `code = mandyocjll` call it before the snapshots.

## Other files

`param.txt` is copied to the output folder (tapIOca reads it there). Nothing else is written: no `time_*.txt`,
`bc_velocity_*.txt` nor `step_*` files.
