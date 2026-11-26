# Installation

## Requirements

- Julia 1.12 or newer ([juliaup](https://github.com/JuliaLang/juliaup) is the simplest way to install it). Tested
  with 1.12.6 and 1.13.0.
- PETSc with MUMPS and an MPI that matches it. Without configuration, the packages bring their own (PETSc_jll and
  MPICH_jll), which work but are slow in parallel on macOS (see below).
- Python with numpy, xarray and zarr to create inputs and read the outputs: the conda environment `tapioca` of
  tapIOca has all of them, and also PETSc 3.24.5 and OpenMPI.

## Build

```bash
git clone git@github.com:jopmas/mandyocjll.git ~/opt/mandyocjll
cd ~/opt/mandyocjll
make
```

`make` installs the Julia packages (`Pkg.instantiate`), compiles a small C shim (`csrc/`) and builds the compiled image
`bin/mandyocjll.dylib` (`.so` on Linux) with PackageCompiler. It takes a few minutes and about 350 MB. Without the
image, mandyocjll still runs, but each run spends 10 to 15 s compiling before the first step; with it, a run starts
in about 5 s.

| Target | |
|---|---|
| `make` | shim and image |
| `make shim` | only the C shim (`csrc/libmandyoc_shim`) |
| `make image` | only the image |
| `make clean` | removes what was built |

The image is tied to the code, the packages and the PETSc/MPI preferences (a checksum is written to
`bin/mandyocjll.stamp`). After a change of `src/`, of `Project.toml`/`Manifest.toml` or of `LocalPreferences.toml`,
`mandyocjll` warns that the image is out of date and runs without it; run `make image` again.

## Using the PETSc and the MPI of conda

The MPICH that comes with MPI.jl on macOS uses TCP only: a 1 MB message between two processes takes about 1260 µs,
against 24 µs with the OpenMPI of conda, and the MUMPS solve with 4 ranks becomes several times slower. To use the
PETSc and OpenMPI of a conda environment, create `LocalPreferences.toml` in the folder of mandyocjll (`<env>` is the
folder of the environment, e.g. `/opt/miniconda3/envs/tapioca`):

```toml
[MPIPreferences]
_format = "1.0"
abi = "OpenMPI"
binary = "system"
cclibs = []
libmpi = "<env>/lib/libmpi.dylib"
mpiexec = "<env>/bin/mpiexec"
preloads = []
preloads_env_switch = ""

[PETSc]
library_path = "<env>/lib/libpetsc.dylib"
PetscScalar = "Float64"
PetscInt = "Int32"
```

`PetscInt` must be the integer size of that PETSc (32 bits for the conda packages). The file is specific to each
computer and is in `.gitignore`; without it, PETSc_jll and MPICH_jll are used. After creating, changing or deleting
it, run `make` again (the shim and the image depend on it).

`tools/use_petsc.jl` writes this file from a PETSc folder (the integer size, the scalar type and the MPI come from its
`petscconf.h`, and `libmpi` is the one `libpetsc` is linked to). When `LocalPreferences.toml` does not exist, `make`
chooses the first compatible PETSc (real double, MUMPS, OpenMPI or MPICH, loads, `libmpi` and `mpiexec` found):

1. `PETSC_DIR`/`PETSC_ARCH` of the environment;
2. on Linux (clusters), a PETSc built in `~/opt/petsc` (newest first), then the conda environment `tapioca`;
   on macOS, the conda environment `tapioca` first, since PETSc built from source on the Mac often fails.

Without a compatible PETSc, PETSc_jll and MPICH_jll are used. `make petsc` chooses again, and a given PETSc wins:

```bash
make PETSC_DIR=/opt/miniconda3/envs/tapioca                                 # PETSc of a conda environment
make PETSC_DIR=$HOME/opt/petsc PETSC_ARCH=optimized-v3.24.1-openmpi           # PETSc built from source (Hypatia)
```

On a cluster, use the PETSc built with the MPI module there, and run `make` on the cluster.

## Threads

Each MPI rank should use one OpenMP thread. The BLAS of the conda PETSc opens one thread per core in every rank,
and the threads of the ranks compete for the cores: a run with 4 ranks was 9 times slower (Mandyoc C++ with the same
PETSc: 5.6 times). `mandyocjll` sets `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` when it runs with more than one
rank, unless they are already set, and the run scripts of tapIOca export `OMP_NUM_THREADS=1` for both codes.
