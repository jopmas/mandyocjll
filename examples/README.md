# Examples

The examples of Mandyoc (`mandyoc/examples`) for mandyocjll. `make_inputs.py` of each folder has the code cells of the
notebook of Mandyoc (BSD 3-Clause, see `LICENSE_MANDYOC`) that write the input files; `run.sh` creates the inputs, if
needed, and runs mandyocjll in the folder of the example:

```bash
conda activate tapioca                 # numpy, scipy, matplotlib, xarray, zarr
cd vanKeken1997_case1a
bash run.sh 4                          # number of MPI ranks (default: 1)
python ../plot_example.py .            # plot_<step>.png of the last output (or: python ../plot_example.py . STEP)
```

`PYTHON` and `MANDYOCJLL` override the programs. Run the long examples in a copy of the folder, or in a cluster.

| Example | Runs in mandyocjll | Options |
|---|---|---|
| `vanKeken1997_case1a` | yes | |
| `Crameri2012_case2` | yes | |
| `continental_rift` | yes | `-seed 0,2 -strain_seed 0.0,1.0` |
| `slab_detachment` | yes | `-seed 0 -strain_seed 0.0` |
| `two_layers_with_variable_conductivity` | yes | |
| `punch` | yes (slow: up to 700 nonlinear iterations) | |
| `variable_velocity_field` | yes | `-seed 0,2 -strain_seed 0.0,1.0` |
| `continental_rift_with_surface_processes` | yes | `-seed 0,2 -strain_seed 0.0,1.0` |
| `continental_rift_non_dimensional/model1` | yes (solved in SI units; denok of the example is converted) | `-seed 0,2 -strain_seed 0.0,1.0` |
| `continental_rift_non_dimensional/model2`, `model3` | no: 3D | `-seed 0,2 -strain_seed 0.0,1.0` |

The examples that do not run yet are kept as targets of the port: mandyocjll stops and names the missing options.

`examples_svg` has the same examples drawn as SVG scenarios of tapIOca.
