"""Converts the .zarr outputs of mandyocjll to the .nc files read by tapIOca (same names, dimensions and types).

    python tools/zarr2nc.py <output folder> [variable ...]

It needs only numpy, xarray and h5netcdf: the zarr v2 (zlib) stores are read directly.
"""
import json
import sys
import zlib
from itertools import product
from pathlib import Path

import numpy as np
import xarray as xr


def read_array(path):
    meta = json.loads((path / ".zarray").read_text())
    attrs = json.loads((path / ".zattrs").read_text())
    shape, chunks = meta["shape"], meta["chunks"]
    out = np.zeros(shape, dtype=meta["dtype"])
    for index in product(*(range(-(-n // c)) for n, c in zip(shape, chunks))):
        f = path / ".".join(map(str, index))
        if not f.exists():
            continue
        raw = f.read_bytes()
        if meta["compressor"]:
            raw = zlib.decompress(raw)
        block = np.frombuffer(raw, meta["dtype"]).reshape(chunks)
        where = tuple(slice(i * c, min((i + 1) * c, n)) for i, c, n in zip(index, chunks, shape))
        out[where] = block[tuple(slice(0, w.stop - w.start) for w in where)]
    dims = attrs.pop("_ARRAY_DIMENSIONS")
    return dims, out, attrs


def open_store(store):
    arrays = {p.parent.name: read_array(p.parent) for p in sorted(store.glob("*/.zarray"))}
    coords = {n: a for n, a in arrays.items() if a[0] == [n]}
    data = {n: a for n, a in arrays.items() if n not in coords}
    return xr.Dataset(data, coords, json.loads((store / ".zattrs").read_text()))


if __name__ == "__main__":
    folder = Path(sys.argv[1])
    names = sys.argv[2:] or [p.stem for p in sorted(folder.glob("*.zarr"))]
    for name in names:
        ds = open_store(folder / f"{name}.zarr")
        encoding = {v: {"zlib": True, "complevel": 5, "shuffle": True, "_FillValue": None} for v in ds.variables}
        ds.to_netcdf(folder / f"{name}.nc", engine="h5netcdf", encoding=encoding)
        print(f"{name}.nc", dict(ds.sizes))
