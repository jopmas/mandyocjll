# mandyocjll documentation

mandyocjll is a port of [Mandyoc](https://github.com/ggciag/mandyoc) (C++/PETSc) to Julia. It reads the same input
files (`param.txt`, `interfaces.txt`, `input_*_0.txt`), runs the same 2D and 3D thermomechanical model (particles in cell,
Q1-P0 Stokes with Uzawa, nonlinear viscoplastic rheology, heat equation) and writes the outputs in zarr, with the
names and units of the NetCDF files of [tapIOca](https://github.com/jopmas/tapioca).

| Page | Contents |
|---|---|
| [installation.md](installation.md) | Julia, PETSc, MPI, `make` and the compiled image |
| [running.md](running.md) | the `mandyocjll` command, its options, MPI, tapIOca run scripts |
| [parameters.md](parameters.md) | what is read from `param.txt` and `interfaces.txt`, and what is not ported yet |
| [outputs.md](outputs.md) | the zarr stores, reading them in Python, conversion to NetCDF |
| [numerics.md](numerics.md) | the numerical methods and where mandyocjll differs from the C++ code |
| [validation.md](validation.md) | tests and comparisons with the C++ code |
| [development.md](development.md) | layout of the code, names, tests, rebuilding the image |
| [bugs_mandyoc_cpp.pdf](bugs_mandyoc_cpp.pdf) | bugs found in the C++ code during the port (in Portuguese) |
| [rename_table.csv](rename_table.csv) | old Mandyoc names and the names used in mandyocjll |

The examples of Mandyoc are in [`../examples`](../examples), and the same examples drawn as SVG scenarios of tapIOca
are in [`../examples/examples_svg`](../examples/examples_svg).

Quick start:

```bash
cd ~/opt/mandyocjll
make                                          # once: dependencies and compiled image (a few minutes)
cd examples/vanKeken1997_case1a
bash run.sh 4                                 # inputs (Python) and a run with 4 MPI ranks
python ../plot_example.py .
```
