# Input files and parameters

mandyocjll reads the input files of Mandyoc. Keys of `param.txt` that it does not use are ignored; keys of features
that change the result and have not been ported stop the run with the list of them. The C++ code stops at keys it
does not know, such as `fssa_legacy`: remove them before running the same folder with Mandyoc.

## Files

| File | When | |
|---|---|---|
| `param.txt` | always | `key = value`; `#` starts a comment |
| `interfaces.txt` | always (`interfaces_from_ascii = True` is required) | layer properties and the depth of each interface at each x node |
| `input_temperature_0.txt` | `temperature_from_ascii = True` | otherwise `T_initial` 0 or 1 |
| `input_velocity_0.txt` | `velocity_from_ascii = True` | initial velocity and the values of the fixed boundaries |
| `bcv_0.txt` | `bcv_extern = True` | mask of the velocity: 1 = free, 0 = fixed (per node and component) |
| `scale_bcv.txt` | `variable_bcv = True` | number of changes, then `time [Myr] factor` per line: the boundary velocity is multiplied by the factor once the end of a step passes that time |

The ASCII fields have 4 header lines and one value per line, x varying fastest (velocity: `vx`, `vz` per node).

## interfaces.txt

Rows `C`, `rho`, `H`, `A`, `n`, `Q`, `V` (required) and `cohesion_min`, `cohesion_max`, `friction_angle_min`,
`friction_angle_max` (degrees), `weakening_seed` and `k` (thermal diffusivity per layer), with one column per layer
from the bottom. Then one line per x node with the depth of each interface. The number of lines must be `nx`.

## param.txt

Read by mandyocjll (C++ names; the Julia names are in `src/parameters.jl` and in `rename_table.csv`):

| Group | Keys |
|---|---|
| Mesh | `nx`, `nz`, `lx`, `lz`, `periodic_boundary` |
| Particles | `particles_per_element`, `particles_perturb_factor`, `random_initial_strain` |
| Time | `step_max`, `time_max`, `dt_max`, `sub_division_time_step`, `step_print`, `initial_print_step`, `initial_print_max_time` |
| Stokes solver | `solver` (direct: MUMPS; iterative), `rtol`, `denok`, `Xi_min`, `non_linear_method`, `initial_dynamic_range`, `free_surface_stab`, `theta_FSSA`, `fssa_legacy` |
| Rheology | `rheology_model` (0-10, 19, 70), `viscosity_reference`, `viscosity_max`, `viscosity_min`, `viscosity_mean_method`, `viscosity_dependence`, `pressure_const`, `plasticity`, `weakening_min`, `weakening_max`, `geoq`, `geoq_fac`, `viscosity_per_element` (must be `constant`), `sticky_blanket_air`, `air_threshold_density`, `strain_healing`, `A_healing`, `k_healing`, `T0_healing` |
| Thermal | `temperature_difference`, `thermal_diffusivity_coefficient`, `thermal_expansion_coefficient`, `gravity_acceleration`, `density_mantle`, `heat_capacity`, `radiogenic_component`, `adiabatic_component`, `shear_heating`, `Xi_shear`, `high_kappa_in_asthenosphere`, `T_initial` (0, 1) |
| Boundaries | `top_/bot_/left_/right_temperature`, `top_/bot_/left_/right_normal_velocity`, `..._tangential_velocity` |
| Inputs | `interfaces_from_ascii`, `temperature_from_ascii`, `velocity_from_ascii`, `bcv_extern`, `variable_bcv` |
| Outputs | `print_step_files` (particle trajectories), `export_lithology`, `export_thermal_diffusivity` |

`fssa_legacy` is new: `True` reproduces the free surface stabilization term of the C++ code (see
[numerics.md](numerics.md)). With `viscosity_max = viscosity_min`, the maximum is raised by 0.01%, as in the C++ code.
`viscosity_per_element` must be `constant`: `variable`, the default of the C++ code, is no longer available there
(it stops the run), so a `param.txt` without this key stops both codes. Rheologies 9 and 19 with `non_linear_method = off` stop the run, unless the viscosity is constant
(`viscosity_max = viscosity_min`): the C++ code silently uses the minimum viscosity (bug 12).

## Not ported yet

These options stop the run (`options of param.txt not ported yet: ...`):

- surface processes: `sp_surface_processes = True`, `sp_mode` other than `none`, `plot_sediment`;
- `magmatism`;
- `nondimensionalization`;
- `multi_velocity`;
- `RK4 = Runge-Kutta` (the RK4 of the C++ code is marked as not working; see the bugs document);
- `basal_heat > 0`;
- `dimensions = 3`.

Ignored on purpose: `binary_output` and `checkered` (outputs are zarr), `external_heat`, `surface_velocity` and
`multigrid` (read but not used by the C++ code). The RK4 will be ported when it works in the C++ code.

`sp_surface_tracking = True` only prints a warning: the surface is not tracked nor written. Also missing:
`T_initial` other than 0 and 1 (use `temperature_from_ascii`), `binary_output`, restart, and the per-step
`step_*` text files (the trajectories are in `particles_trajectories.zarr`).
