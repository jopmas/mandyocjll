# Examples from SVGs

The examples of Mandyoc drawn as SVG scenarios of tapIOca, which creates the input files and the run script. They are
the SVGs of `examples/classical_examples` of tapIOca with `code = mandyocjll` in the `setup` shape, so the run script
calls mandyocjll. The differences from the notebooks of Mandyoc are in the comments at the top of each SVG.

```bash
conda activate tapioca
bash run_example.sh                          # lists the examples
bash run_example.sh vanKeken1997_case1a 4    # inputs and outputs in vanKeken1997_case1a/inputs, 4 MPI ranks
TAPIOCA_SVG_OPTIONS="--dx 10 --dz 10" bash run_example.sh continental_rift 4    # coarser
```

Or step by step:

```bash
cd continental_rift
tapioca-svg continental_rift.svg -o inputs --ncores 4     # --system linux, mac or hypatia; --code mandyoc for the C++ code
cd inputs
bash run_mac.sh                                           # or run_linux.sh; sbatch run_hypatia.sh
```

Without `--no-post-processing`, the script also creates the `.nc` files, the snapshots and the video. The weak seeds
are layers of the drawing, so `-seed` is not needed.

Needs the setting `code` of tapIOca (branch `feat/svg-scenario-inputs`). `continental_rift_non_dimensional` does not run in
mandyocjll yet (tapioca-svg warns); `model2` and `model3` are the xz
section of 3D models.
