"""Plots the zarr outputs of a scenario: python plot_example.py FOLDER [STEP]

Temperature, viscosity, strain rate and velocity of the last output (or of STEP), saved in FOLDER/plot_<step>.png.
Needs numpy, matplotlib, xarray and zarr (conda environment of tapIOca).
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import xarray as xr

matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = Path(sys.argv[1])
ds = xr.merge([xr.open_zarr(folder / f"{v}.zarr") for v in ("density", "temperature", "viscosity", "strain_rate", "velocity")], compat="override")
i = -1 if len(sys.argv) < 3 else int(np.flatnonzero(ds.step.values == int(sys.argv[2]))[0])
d = ds.isel(time=i).load()
scale, unit = (1.0e3, "km") if float(d.x.max()) > 1.0e4 else (1.0, "m")
x, z = d.x / scale, d.z / scale
speed = np.hypot(d.vx, d.vz) * 100 * 365.25 * 24 * 3600
with np.errstate(divide="ignore"):
    fields = [("Density [kg/m$^3$]", d.density, "viridis"), ("Temperature [$^o$C]", d.temperature, "inferno"),
              ("log$_{10}$ viscosity [Pa s]", np.log10(d.viscosity.where(d.viscosity > 0)), "magma"),
              ("log$_{10}$ strain rate [1/s]", np.log10(d.strain_rate.where(d.strain_rate > 0)), "cividis"),
              ("Speed [cm/yr]", speed, "plasma")]
fig, axes = plt.subplots(len(fields), 1, figsize=(10, 3 * len(fields)), sharex=True, constrained_layout=True)
for ax, (title, values, cmap) in zip(axes, fields):
    im = ax.pcolormesh(x, z, values, cmap=cmap, shading="nearest")
    fig.colorbar(im, ax=ax, label=title)
    ax.set_aspect("equal")
    ax.set_ylabel(f"z [{unit}]")
k = max(1, d.sizes["x"] // 25)
axes[-1].quiver(x[::k], z[::k], d.vx[::k, ::k], d.vz[::k, ::k], color="w", width=0.002)
axes[-1].set_xlabel(f"x [{unit}]")
fig.suptitle(f"{folder.resolve().name}: {float(d.time):.3g} Myr (step {int(d.step)})")
out = folder / f"plot_{int(d.step)}.png"
fig.savefig(out, dpi=120)
print(out)
