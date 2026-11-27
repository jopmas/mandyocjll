"""Figures of the mandyocjll handbook (docs/apostila), from the test runs described in its chapter on validation.

    python docs/apostila/scripts/make_figures.py RUNS_FOLDER [FIGURE ...]

RUNS_FOLDER has one subfolder per test (the names below); each figure needs only its own subfolders and is skipped,
with a message, if they are missing. The runs were made with the inputs of docs/validation.md (Julia: *.zarr;
C++: the text files of Mandyoc, <variable>_<step>.txt).
"""
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from cmcrameri import cm

R = Path(sys.argv[1])
OUT = Path(__file__).resolve().parents[1] / "figures"
OUT.mkdir(exist_ok=True)
JL, CPP, REF, INK2, GRID = "#2a78d6", "#eb6834", "#52514e", "#52514e", "#e4e3de"   # categorical slots 1 and 2
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "figure.dpi": 100,
                     "savefig.dpi": 200, "savefig.bbox": "tight", "lines.linewidth": 2,
                     "figure.constrained_layout.use": True})


def asc(f):
    return pd.read_csv(f, header=None, skiprows=2, comment="P")[0].to_numpy(float)


def zarr(folder, name, var=None):
    return xr.open_zarr(R / folder / f"{name}.zarr")[var or name].to_numpy().astype(float)


def needs(*folders):
    missing = [f for f in folders if not (R / f).exists()]
    if missing:
        print("skipped (missing:", ", ".join(missing) + ")")
    return not missing


def save(fig, name):
    fig.savefig(OUT / name); plt.close(fig); print(OUT / name)


def field_map(ax, a, extent, cmap, label, vmin=None, vmax=None):
    im = ax.imshow(a, origin="lower", extent=extent, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.grid(False)
    plt.colorbar(im, ax=ax, label=label, fraction=0.046, pad=0.02)


def diff_map(ax, d, extent, label):
    m = np.nanmax(np.abs(d)) or 1.0
    field_map(ax, d, extent, cm.vik, label, -m, m)


# ---- 1. two layers of different conductivity: steady state against the analytical solution ----
def two_layers():
    f = "ex2/two_layers_with_variable_conductivity"
    if not needs(f): return
    T = zarr(f, "temperature"); z = xr.open_zarr(R / f / "temperature.zarr").z.to_numpy() / 1e3
    Tb, Tt = T[0, 0, 5], T[0, -1, 5]; Ti = (1.0e-6 * Tb + 0.25e-6 * Tt) / (1.25e-6)
    exact = np.where(z < -20, Tb + (Ti - Tb) * (z + 40) / 20, Ti + (Tt - Ti) * (z + 20) / 20)
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.4, 2.8), gridspec_kw={"width_ratios": [1.3, 1]})
    a.plot(exact, z, color=CPP, label="analytical"); a.plot(T[-1, :, 5], z, "o", color=JL, ms=4, label="mandyocjll (50 Myr)")
    a.plot(T[0, :, 5], z, color=REF, lw=1, ls=":", label="initial")
    a.set_xlabel("temperature (°C)"); a.set_ylabel("z (km)"); a.legend(loc="lower left")
    b.plot(T[-1, :, 5] - exact, z, color=JL); b.axhline(-20, color=REF, lw=0.8, ls="--")
    b.set_xlabel("mandyocjll − analytical (°C)"); b.text(0.02, 0.52, "interface", transform=b.transAxes, color=INK2)
    save(fig, "fig_two_layers.png")


# ---- 2. punch: velocity of mandyocjll and difference to the C++ code ----
def punch():
    if not needs("punch_full", "cpp_punch"): return
    v = xr.open_zarr(R / "punch_full/velocity.zarr").isel(time=0)
    nz, nx = v.vx.shape
    u, w = v.vx.to_numpy(), v.vz.to_numpy()
    c = asc(R / "cpp_punch/velocity_0.txt").reshape(nz, nx, 2)
    s = np.hypot(u, w); sc = np.hypot(c[..., 0], c[..., 1])
    ext = [0, float(v.x[-1]), float(v.z[0]), float(v.z[-1])]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 2.8))
    field_map(a, s / s.max(), ext, cm.lajolla, "|v| / max |v|"); a.set_title("mandyocjll", loc="left")
    diff_map(b, (s - sc) / s.max(), ext, "(mandyocjll − C++) / max |v|"); b.set_title("difference to the C++ code", loc="left")
    for ax in (a, b): ax.set_xlabel("x"); ax.set_ylabel("z")
    save(fig, "fig_punch.png")


# ---- 3. variable velocity field: horizontal velocity at the surface across a reversal of scale_bcv.txt ----
def variable_velocity():
    if not needs("jl_vvf", "cpp_vvf"): return
    v = xr.open_zarr(R / "jl_vvf/velocity.zarr"); nz, nx = v.vx.shape[1:]; x = v.x.to_numpy() / 1e3
    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    for k in range(v.sizes["time"]):
        c = asc(R / f"cpp_vvf/velocity_{k}.txt").reshape(nz, nx, 2)[-1, :, 0]
        yr = 3.15576e7 * 100
        ax.plot(x, c * yr, color=CPP, lw=3, alpha=0.5, label="C++" if k == 0 else None)
        ax.plot(x, v.vx[k, -1].to_numpy() * yr, color=JL, lw=1.2, label="mandyocjll" if k == 0 else None)
        ax.annotate(f"step {k}", (x[-1], v.vx[k, -1, -1] * yr), xytext=(4, 0), textcoords="offset points", color=INK2, va="center")
    ax.set_xlabel("x (km)"); ax.set_ylabel("vx at the top (cm/yr)"); ax.legend(loc="upper left")
    save(fig, "fig_variable_velocity.png")


# ---- 4. strain healing and shear heating: strain and heat of mandyocjll and differences ----
def healing_shear():
    if not needs("heal_jl", "heal_cpp", "xs2_jl", "xs2_cpp"): return
    s = xr.open_zarr(R / "heal_jl/strain.zarr"); nz, nx = s.strain.shape[1:]; k = s.sizes["time"] - 1
    ext = [0, float(s.x[-1]) / 1e3, float(s.z[0]) / 1e3, float(s.z[-1]) / 1e3]
    a1 = s.strain[k].to_numpy(); c1 = asc(R / f"heal_cpp/strain_{k}.txt").reshape(nz, nx)
    h = zarr("xs2_jl", "heat")[k]; c2 = asc(R / f"xs2_cpp/heat_{k}.txt").reshape(nz, nx)
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 4.0))
    field_map(ax[0, 0], a1, ext, cm.lajolla, "strain"); ax[0, 0].set_title("strain_healing: strain (mandyocjll)", loc="left")
    diff_map(ax[0, 1], (a1 - c1) / np.abs(c1).max(), ext, "relative difference"); ax[0, 1].set_title("mandyocjll − C++", loc="left")
    field_map(ax[1, 0], h, ext, cm.lajolla, "heat (W/kg)"); ax[1, 0].set_title("Xi_shear = 5: heat (mandyocjll)", loc="left")
    diff_map(ax[1, 1], (h - c2) / np.abs(c2).max(), ext, "relative difference"); ax[1, 1].set_title("mandyocjll − C++", loc="left")
    for a in ax.flat: a.set_xlabel("x (km)"); a.set_ylabel("z (km)")
    save(fig, "fig_healing_shear.png")


# ---- 5. surface processes: change of the surface in the four modes ----
def surface_processes():
    modes = ["diffusion", "sedimentation_only", "sedimentation_rate_limited", "theunissen_sedimentation"]
    if not needs(*(f"sp/{m}_{c}" for m in modes for c in ("jl", "cpp"))): return
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 4.0), sharex=True)
    for ax, m in zip(axs.flat, modes):
        s = xr.open_zarr(R / f"sp/{m}_jl/surface.zarr"); x = s.x.to_numpy() / 1e3; z0 = s.surface[0].to_numpy()
        k = s.sizes["time"] - 1
        c = np.loadtxt(R / f"sp/{m}_cpp/surface_{k}.txt")[:, 1]; c0 = np.loadtxt(R / f"sp/{m}_cpp/surface_0.txt")[:, 1]
        ax.plot(x, c - c0, color=CPP, lw=3, alpha=0.5, label="C++")
        ax.plot(x, s.surface[k].to_numpy() - z0, color=JL, lw=1.2, label="mandyocjll")
        ax.set_title(m, loc="left"); ax.set_ylabel("Δ surface (m)"); ax.ticklabel_format(axis="y", useOffset=False)
    axs[0, 0].legend(loc="upper right")
    for ax in axs[1]: ax.set_xlabel("x (km)")
    save(fig, "fig_surface_processes.png")


# ---- 6. magmatism: melt fraction, extrusion and the dependence of dPhi on the nonlinear iterations (bug 22) ----
def magmatism():
    if not needs("mg/e970_jl", "mg/e970_cpp"): return
    p = xr.open_zarr(R / "mg/e970_jl/Phi.zarr"); nz, nx = p.Phi.shape[1:]
    ext = [0, float(p.x[-1]) / 1e3, float(p.z[0]) / 1e3, float(p.z[-1]) / 1e3]
    fig = plt.figure(figsize=(7.4, 4.2)); g = fig.add_gridspec(2, 2)
    a = fig.add_subplot(g[0, :]); field_map(a, p.Phi[1].to_numpy(), ext, cm.lajolla, "Phi")
    a.set_title("melt fraction after step 1 (magmatism_Ts0 = 970 °C)", loc="left"); a.set_xlabel("x (km)"); a.set_ylabel("z (km)")
    a.set_ylim(-250, 0)
    b = fig.add_subplot(g[1, 0]); s = xr.open_zarr(R / "mg/e970_jl/surface.zarr"); x = s.x.to_numpy() / 1e3
    c = np.loadtxt(R / "mg/e970_cpp/surface_1.txt")[:, 1]
    b.plot(x, c / 1e3, color=CPP, lw=3, alpha=0.5, label="C++"); b.plot(x, s.surface[1].to_numpy() / 1e3, color=JL, lw=1.2, label="mandyocjll")
    b.set_xlim(0, 80); b.set_xlabel("x (km)"); b.set_ylabel("surface (km)"); b.legend(loc="upper right")
    b.set_title("extrusion at x = 0 (bug 25)", loc="left")
    c = fig.add_subplot(g[1, 1])
    runs = [("1e-3", 2), ("1e-6", 10), ("1e-9", 34)]
    if all((R / f"mg/it{t}").exists() for t, _ in runs):
        d = [zarr(f"mg/it{t}", "dPhi")[0].max() for t, _ in runs]
        c.bar(range(3), d, color=JL, width=0.6); c.set_yscale("log"); c.set_xticks(range(3), [f"Xi_min {t}\n{n} iterations" for t, n in runs])
        c.set_ylabel("max dPhi (step 0)"); c.set_title("latent heat source vs iterations (bug 22)", loc="left")
    save(fig, "fig_magmatism.png")


# ---- 7. denok and nondimensionalization: divergence of the velocity (bug 28) ----
def denok():
    if not needs("nd/dim2_cpp", "nd/dim3_cpp", "nd/nd3_cpp", "nd/nd4_jl"): return
    h0, visc0, g0, rho0 = 3e5, 1e25, 10, 3300; v0 = rho0 * g0 * h0**2 / visc0
    dx, dz = 1600e3 / 160, 300e3 / 30
    def ratio(v):
        u, w = v[..., 0], v[..., 1]
        dudx = ((u[:-1, 1:] - u[:-1, :-1]) + (u[1:, 1:] - u[1:, :-1])) / 2 / dx
        dwdz = ((w[1:, :-1] - w[:-1, :-1]) + (w[1:, 1:] - w[:-1, 1:])) / 2 / dz
        dudz = ((u[1:, :-1] - u[:-1, :-1]) + (u[1:, 1:] - u[:-1, 1:])) / 2 / dz
        dwdx = ((w[:-1, 1:] - w[:-1, :-1]) + (w[1:, 1:] - w[1:, :-1])) / 2 / dx
        return np.linalg.norm(dudx + dwdz) / np.linalg.norm(np.sqrt(dudx**2 + dwdz**2 + dudz**2 + dwdx**2))
    rows = [("C++, SI, denok = 1e-5", asc(R / "nd/dim2_cpp/velocity_0.txt").reshape(31, 161, 2)),
            ("C++, SI, denok = 7.9e-14", asc(R / "nd/dim3_cpp/velocity_0.txt").reshape(31, 161, 2)),
            ("C++, nondimensional, denok = 1e-5", asc(R / "nd/nd3_cpp/velocity_0.txt").reshape(31, 161, 2) * v0)]
    v = xr.open_zarr(R / "nd/nd4_jl/velocity.zarr").isel(time=0)
    rows.append(("mandyocjll, nondimensional (converted)", np.stack([v.vx.to_numpy(), v.vz.to_numpy()], -1)))
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.0), gridspec_kw={"width_ratios": [1.1, 1]})
    vals = [ratio(r) for _, r in rows]
    a.barh(range(4)[::-1], vals, color=[CPP, CPP, CPP, JL], height=0.6); a.set_xscale("log")
    a.set_yticks(range(4)[::-1], [n for n, _ in rows]); a.set_xlabel("‖div v‖ / ‖grad v‖ (step 0)")
    x = np.linspace(0, 1600, 161)
    for (n, r), col, ls in zip(rows, [CPP, CPP, CPP, JL], [":", "-", "--", "-"]):
        b.plot(x, r[20, :, 1] * 3.15576e9, color=col, ls=ls, lw=1.4 if col == JL else 2.2, alpha=1 if col == JL else 0.6, label=n)
    b.set_xlabel("x (km)"); b.set_ylabel("vz at z = −100 km (cm/yr)"); b.legend(loc="lower center", fontsize=6.5)
    save(fig, "fig_denok_nondimensional.png")


# ---- 8. 3D: temperature and velocity in a section, and difference to the C++ code ----
def three_d():
    if not needs("d3/j1", "d3/c1"): return
    t = xr.open_zarr(R / "d3/j1/temperature.zarr"); v = xr.open_zarr(R / "d3/j1/velocity.zarr")
    nz, ny, nx = t.temperature.shape[1:]; k = t.sizes["time"] - 1; j = ny // 2
    ext = [0, float(t.x[-1]) / 1e3, float(t.z[0]) / 1e3, float(t.z[-1]) / 1e3]
    T = t.temperature[k, :, j].to_numpy(); Tc = asc(R / f"d3/c1/temperature_{k}.txt").reshape(nz, ny, nx)[:, j]
    w = v.vz[k, :, j].to_numpy(); wc = asc(R / f"d3/c1/velocity_{k}.txt").reshape(nz, ny, nx, 3)[:, j, :, 2]
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 3.8))
    field_map(ax[0, 0], T, ext, cm.lajolla, "°C"); ax[0, 0].set_title(f"temperature, section y = {float(t.y[j]) / 1e3:.0f} km", loc="left")
    diff_map(ax[0, 1], T - Tc, ext, "°C"); ax[0, 1].set_title("mandyocjll − C++", loc="left")
    field_map(ax[1, 0], w * 3.15576e9, ext, cm.vik, "cm/yr", -np.abs(w).max() * 3.15576e9, np.abs(w).max() * 3.15576e9)
    ax[1, 0].set_title("vertical velocity", loc="left")
    diff_map(ax[1, 1], (w - wc) / np.abs(wc).max(), ext, "relative"); ax[1, 1].set_title("mandyocjll − C++", loc="left")
    for a in ax.flat: a.set_xlabel("x (km)"); a.set_ylabel("z (km)")
    save(fig, "fig_3d.png")


# ---- 9. tutorial: lithology of the continental rift at three times ----
def tutorial():
    run = "tut20/inputs" if (R / "tut20/inputs/lithology.zarr").exists() else "tut/inputs"
    if not needs(f"{run}/lithology.zarr"): return
    l = xr.open_zarr(R / run / "lithology.zarr")
    ks = [0, l.sizes["time"] // 2, l.sizes["time"] - 1]
    fig, axs = plt.subplots(3, 1, figsize=(7.0, 4.6), sharex=True)
    ext = [0, float(l.x[-1]) / 1e3, float(l.z[0]) / 1e3, float(l.z[-1]) / 1e3]
    for ax, k in zip(axs, ks):
        ax.imshow(l.lithology[k].to_numpy(), origin="lower", extent=ext, aspect="auto", cmap=cm.batlowS, interpolation="nearest")
        ax.grid(False); ax.set_ylim(-200, 0); ax.set_xlim(400, 1200); ax.set_ylabel("z (km)")
        ax.set_title(f"{float(l.time[k]):.2f} Myr", loc="left")
    axs[-1].set_xlabel("x (km)")
    save(fig, "fig_tutorial_lithology.png")


# ---- 10. FastScape: steady river profile (test/test_surface.jl) and the surface of a rift with each sp_mode ----
def fastscape():
    runs = [("none", "no surface processes"), ("diffusion", "sp_mode = diffusion"), ("fastscape", "sp_mode = fastscape (marine)")]
    if not needs("surface_test", *(f"fastscape/{r}" for r, _ in runs)): return
    fig = plt.figure(figsize=(7.4, 5.6))
    upper, lower = fig.subfigures(2, 1, height_ratios=[1, 1.15])
    a, e = upper.subplots(1, 2, gridspec_kw={"width_ratios": [2, 1]})
    p2 = pd.read_csv(R / "surface_test/river_profile_2d.csv"); p3 = pd.read_csv(R / "surface_test/river_profile_3d.csv")
    a.plot(p2.x / 1e3, p2.exact, color=CPP, lw=3, alpha=0.5, label="exact (discrete)")
    a.plot(p2.x / 1e3, p2.fastscape, color=JL, lw=1.2, label="2D strip")
    a.plot(p3.x[::5] / 1e3, p3.fastscape[::5], "o", color=JL, ms=3.5, label="3D grid")
    a.set_xlabel("distance from the base level (km)"); a.set_ylabel("height (m)"); a.legend(loc="lower right")
    a.set_title("steady river profile, uplift 1 mm/yr", loc="left")
    e.semilogy(p2.x / 1e3, np.abs(p2.fastscape - p2.exact) / p2.exact.max() + 1e-17, color=JL, lw=1.2)
    e.set_xlabel("distance (km)"); e.set_ylabel("relative error"); e.set_title("error", loc="left"); e.set_ylim(1e-17, 1e-12)
    times = [0, 5, 10, 15, 20]
    seq = [cm.batlow(i / (len(times) - 1) * 0.85) for i in range(len(times))]
    for n, ((r, label), ax) in enumerate(zip(runs, lower.subplots(1, 3, sharey=True))):
        s = xr.open_zarr(R / f"fastscape/{r}/surface.zarr"); z0 = float(s.surface[0, 0])
        for t, c in zip(times, seq):
            ax.plot(s.x / 1e3, (s.surface.sel(time=t, method="nearest") - z0) / 1e3, color=c, lw=1.1, label=f"{t} Myr")
        ax.set_xlim(500, 1100); ax.set_ylim(-1.2, 2.0); ax.set_title(label, loc="left", fontsize=8)
        ax.set_xlabel("x (km)")
        if n == 0: ax.set_ylabel("surface (km)"); ax.legend(loc="upper left", fontsize=7, ncol=2)
    save(fig, "fig_fastscape.png")


FIGURES = dict(two_layers=two_layers, punch=punch, variable_velocity=variable_velocity, healing_shear=healing_shear,
               surface_processes=surface_processes, magmatism=magmatism, denok=denok, three_d=three_d, tutorial=tutorial,
               fastscape=fastscape)
for name in sys.argv[2:] or FIGURES:
    print(f"-- {name}"); FIGURES[name]()
