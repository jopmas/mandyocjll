# Validation

All results below were obtained on macOS (arm64) with PETSc 3.24.5 and OpenMPI of the conda environment `tapioca`.
The C++ code is `ggciag/mandyoc` (commit `956af67`) or `mandyoc_jpm` (commit `64dc7d1`), built in temporary copies.

## Tests of the parts

Run in the folder of mandyocjll (`test/mpirun.jl N` starts N ranks with the MPI of MPI.jl):

```bash
julia --project=. test/mpirun.jl 1 shim closed      # particles: reference, run first
julia --project=. test/mpirun.jl 4 shim closed
julia --project=. test/mpirun.jl 4 native closed
julia --project=. test/mpirun.jl 1 shim periodic    # reference of the periodic case
julia --project=. test/mpirun.jl 4 shim periodic
julia --project=. test/mpirun.jl 4 test_thermal.jl
julia --project=. test/mpirun.jl 4 test_stokes.jl
julia --project=. test/mpirun.jl 4 test_open_boundary.jl shim
julia --project=. test/mpirun.jl 4 test_particle_control.jl native
julia --project=. test/test_rheology.jl <mandyoc>/src <petsc>/include
julia --project=. test/mpirun.jl 1 van_keken.jl 81 6 300
```

Each test prints `RESULT OK` or `RESULT FAIL`:

| Test | What it checks | Result |
|---|---|---|
| `test_particle_swarm.jl` (default of `mpirun.jl`) | migration of 243,000 particles between ranks, closed and periodic, with the shim and the native location | same particles and positions with 1, 2, 4 and 6 ranks |
| `test_thermal.jl` | transient diffusion (order 2), two layers and radiogenic heating (steady state), advection-diffusion, periodic | errors 1e-10 or analytical order |
| `test_stokes.jl` | two-layer Couette flow, closed and periodic; effective weight of the FSSA | 1e-14; FSSA as intended and as in C++ |
| `test_open_boundary.jl` | particles leaving at an open boundary | exact count, positions 1e-10 m |
| `test_particle_control.jl` | particles per element and layers after adding and removing | within the limits, same hash for any number of ranks |
| `test_rheology.jl` | viscosity of the 13 rheologies against the C++ function (compiled from the Mandyoc source) | 1e-14 |
| `van_keken.jl` | van Keken et al. (1997), case 1a | vrms peak 0.003103 (C++: 0.003111) |

## Comparisons with the C++ code

Thermal step of a rift (201 × 51 nodes, with advection, radiogenic and adiabatic heating): maximum difference
1e-9 °C for changes of 7 °C. One Stokes solve (van Keken, 81 × 81, `fssa_legacy = True`): 1e-13 without FSSA, 1e-11
with it. The full model (`test/compare_cpp_model.jl`, relative L2 norm):

| Case | Step | Velocity | log10 viscosity | Temperature |
|---|:-:|:-:|:-:|:-:|
| viscoplastic rift (201 × 51) | 3 | 2.8e-8 | 1.3e-8 | 8.1e-12 |
| slab, power law (101 × 61) | 4 | 4.3e-7 | 6.1e-9 | - |
| slab, with particle control | 6 | 3.1e-7 | strain 6.9e-6 | - |
| rift, with particle control | 30 | 8.0e-7 | strain 1.4e-6 | - |

The rift gives the same result with 1 and 4 ranks.

## Examples of Mandyoc

Short runs of the examples in `examples/`:

| Example | Comparison | Result |
|---|---|---|
| `two_layers_with_variable_conductivity` | analytical steady state after 10,000 steps | max. error 0.87 °C of 700 °C |
| `variable_velocity_field` | C++ with the same inputs, 3 steps across a reversal of `scale_bcv.txt` (4 ranks, `fssa_legacy`, no perturbation of the particles) | velocity 1.4e-6, temperature 2.2e-6 |
| `punch` | C++ and the reference output of the tests of Mandyoc (`tests/data/punch`) | velocity 6.3e-4 (C++), 9.6e-4 (reference); C++ vs reference: 4.1e-4 |
| `continental_rift` from its SVG (161 × 31) with `strain_healing` | C++, 4 steps, 2 ranks | velocity 4e-6, strain 2e-6 |
| `continental_rift` from its SVG (161 × 31) with `shear_heating`, `Xi_shear = 5` | C++, 4 steps, 2 ranks | heat 1.6e-5, temperature 2.1e-6 |
| `vanKeken1997_case1a`, `Crameri2012_case2`, `continental_rift`, `slab_detachment` | run, plots | no comparison yet beyond the tests above |
| `continental_rift_with_surface_processes`, `continental_rift_non_dimensional` | stop with the list of options not ported | |

The punch does not converge in either code: the nonlinear iterations stop at 700 with `Xi` around 1e-8 (C++: 8e-9),
above `Xi_min = 1e-12`, so the differences are at the level of the last iterations.

## Outputs

The zarr stores were read with the official `zarr` 3.4.2 and xarray, and `tools/zarr2nc.py` gives the same values
as the NetCDF files written by tapIOca from the text outputs. The compiled image gives the same results as the
uncompiled code (1e-11).
