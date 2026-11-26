'''Zarr (zlib 1) size of one output of craton_colision_lwr, from its ASCII files (same dtypes/chunks as
src/zarr_output.cpp). Values come with ASCII precision, so the f8 fields are a lower bound.'''
import glob, sys, time, zlib
import numpy as np
import pandas as pd

n = sys.argv[1]
nx, nz = 2501, 351
DT = {"density": "f4", "pressure": "f4", "temperature": "f4", "viscosity": "f8", "strain": "f8",
      "strain_rate": "f8", "heat": "f8", "velocity": "f8"}


def field(name):
    vals = pd.read_csv(f"{name}_{n}.txt", header=None, comment=None, skiprows=0, names=["v"], dtype=str)["v"]
    vals = pd.to_numeric(vals.str.strip(), errors="coerce").dropna().to_numpy()
    return vals


def zsize(a):
    return len(zlib.compress(np.ascontiguousarray(a).tobytes(), 1))


out = {}
for name, dt in DT.items():
    v = field(name)
    raw = v.astype(dt)
    t0 = time.time()
    if name == "velocity":
        s = zsize(raw[0::2]) + zsize(raw[1::2])
    else:
        s = zsize(raw)
    out[name] = (v.size, s, time.time() - t0)

cols = ["x", "z", "id", "layer", "strain"]
parts = pd.concat([pd.read_csv(f, sep=r"\s+", header=None, names=cols) for f in glob.glob(f"step_{n}_*.txt")])
p = parts.sort_values("id")
s = zsize(p.x.to_numpy("f4")) + zsize(p.z.to_numpy("f4")) + zsize(p.layer.to_numpy("i1")) + \
    zsize(p.strain.to_numpy("f4"))
out["particles"] = (len(p), s, 0)

fnx, fnz = (nx - 1) * 5 + 1, (nz - 1) * 5 + 1
lit = np.full((fnz, fnx), -1, dtype=np.int8)
for f in glob.glob(f"litho_{n}_*.txt"):
    d = pd.read_csv(f, sep=r"\s+", header=None).to_numpy()
    i, j, l = d[:, 0], d[:, 1], d[:, 2]
    ok = (i >= 0) & (i < fnx) & (j >= 0) & (j < fnz)
    np.maximum.at(lit, (j[ok], i[ok]), l[ok].astype(np.int8))
out["lithology"] = (lit.size, zsize(lit[::-1]), 0)

sv = np.loadtxt(f"surface_{n}.txt", dtype="f4") if glob.glob(f"surface_{n}.txt") else np.zeros(1)
out["surface"] = (sv.size, zsize(sv.ravel().astype("f4")), 0)

for k, (cnt, s, t) in out.items():
    print(f"{k} {cnt} {s}")
