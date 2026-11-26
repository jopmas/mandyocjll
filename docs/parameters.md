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
| Thermal | `temperature_difference`, `thermal_diffusivity_coefficient`, `thermal_expansion_coefficient`, `gravity_acceleration`, `density_mantle`, `heat_capacity`, `radiogenic_component`, `adiabatic_component`, `shear_heating`, `Xi_shear`, `basal_heat`, `high_kappa_in_asthenosphere`, `T_initial` (0, 1) |
| Boundaries | `top_/bot_/left_/right_temperature`, `top_/bot_/left_/right_normal_velocity`, `..._tangential_velocity` |
| Inputs | `interfaces_from_ascii`, `temperature_from_ascii`, `velocity_from_ascii`, `bcv_extern`, `variable_bcv`, `multi_velocity` (`multi_veloc.txt` and `input_velocity_<n>.txt`) |
| Surface | `sp_surface_tracking`, `sp_surface_processes`, `sp_mode` (`diffusion`, `sedimentation_only`, `sedimentation_rate_limited`, `theunissen_sedimentation`), `sp_d_c`, `surface_particles_per_element`, `sea_level`, `a2l`, `sp_continental_slope`, `sp_sediment_strain`, `variable_base_level` |
| Magmatism | `magmatism`, `magmatism_extraction` (needs the surface tracking), `magmatic_layer`, `magmatism_Ts0`, `magmatism_dTdP`, `magmatism_dTdX`, `magmatism_dS` |
| Outputs | `print_step_files` (particle trajectories), `export_lithology`, `export_thermal_diffusivity` |

`fssa_legacy` is new: `True` reproduces the free surface stabilization term of the C++ code (see
[numerics.md](numerics.md)). With `viscosity_max = viscosity_min`, the maximum is raised by 0.01%, as in the C++ code.
`viscosity_per_element` must be `constant`: `variable`, the default of the C++ code, is no longer available there
(it stops the run), so a `param.txt` without this key stops both codes. Rheologies 9 and 19 with `non_linear_method = off` stop the run, unless the viscosity is constant
(`viscosity_max = viscosity_min`): the C++ code silently uses the minimum viscosity (bug 12).
`nondimensionalization = True` is accepted, but the equations are solved in SI units (the same solution; the C++
code mixes units in the time step, `time_max` and the input files with times, bug 27): the outputs are in SI units,
the scales are written to `dimensional_convertion.txt`, and `denok`, given in the units of the C++ run, is converted.
`denok` is an absolute bound (bug 28): in SI units use about 1e-14, as tapIOca does (a warning is printed above 1e-10).
`magmatism_dS` is used both in the melting and in the latent heat; the C++ code uses 400 in the latent heat (bug 23).

## Not ported yet

These options stop the run (`options of param.txt not ported yet: ...`):

- `RK4 = Runge-Kutta` (the RK4 of the C++ code is marked as not working; see the bugs document);
- `dimensions = 3`.

Ignored on purpose: `binary_output`, `checkered` and `plot_sediment` (outputs are zarr), `external_heat`, `surface_velocity` and
`multigrid` and `aggradation_rate` (read but not used by the C++ code). The RK4 will be ported when it works in the C++ code.

Not ported:
`T_initial` other than 0 and 1 (use `temperature_from_ascii`), `binary_output`, restart, and the per-step
`step_*` text files (the trajectories are in `particles_trajectories.zarr`).
