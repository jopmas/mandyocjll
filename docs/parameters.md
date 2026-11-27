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
| Mesh | `dimensions` (2, 3), `nx`, `ny` (3D), `nz`, `lx`, `ly` (3D), `lz`, `periodic_boundary` (2D) |
| Particles | `particles_per_element`, `particles_perturb_factor`, `random_initial_strain` |
| Time | `step_max`, `time_max`, `dt_max`, `sub_division_time_step`, `step_print`, `initial_print_step`, `initial_print_max_time` |
| Stokes solver | `solver` (direct: MUMPS; iterative), `rtol`, `denok`, `Xi_min`, `non_linear_method`, `initial_dynamic_range`, `free_surface_stab`, `theta_FSSA`, `fssa_legacy` |
| Rheology | `rheology_model` (0-10, 19, 70), `viscosity_reference`, `viscosity_max`, `viscosity_min`, `viscosity_mean_method`, `viscosity_dependence`, `pressure_const`, `plasticity`, `weakening_min`, `weakening_max`, `geoq`, `geoq_fac`, `viscosity_per_element` (must be `constant`), `sticky_blanket_air`, `air_threshold_density`, `strain_healing`, `A_healing`, `k_healing`, `T0_healing` |
| Thermal | `temperature_difference`, `thermal_diffusivity_coefficient`, `thermal_expansion_coefficient`, `gravity_acceleration`, `density_mantle`, `heat_capacity`, `radiogenic_component`, `adiabatic_component`, `shear_heating`, `Xi_shear`, `basal_heat`, `high_kappa_in_asthenosphere`, `T_initial` (0, 1) |
| Boundaries | `top_/bot_/left_/right_temperature`, `top_/bot_/left_/right_normal_velocity`, `..._tangential_velocity`; in 3D also `front_/back_temperature` and `front_/back_normal_velocity`, `..._tangential_velocity` (as in the C++ code, `front` is y = 0 for the temperature and y = Ly for the velocity, bug 30) |
| Inputs | `interfaces_from_ascii`, `temperature_from_ascii`, `velocity_from_ascii`, `bcv_extern`, `variable_bcv`, `multi_velocity` (`multi_veloc.txt` and `input_velocity_<n>.txt`) |
| Surface | `sp_surface_tracking`, `sp_surface_processes`, `sp_mode` (`diffusion`, `sedimentation_only`, `sedimentation_rate_limited`, `theunissen_sedimentation`, `fastscape`; in 3D only `diffusion` and `fastscape`), `sp_d_c`, `surface_particles_per_element`, `sea_level`, `a2l`, `sp_continental_slope`, `sp_sediment_strain`, `variable_base_level` |
| FastScape (`sp_mode = fastscape`, new) | `fastscape_kf`, `fastscape_kfsed`, `fastscape_m`, `fastscape_n`, `fastscape_kd`, `fastscape_kdsed`, `fastscape_g`, `fastscape_gsed`, `fastscape_p`, `fastscape_boundary`, `fastscape_dt_max`, `fastscape_strip_width` (2D), `fastscape_marine` and its parameters (table below) |
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

## Erosion and sedimentation with FastScape

`sp_mode = fastscape` (not in the C++ code) couples the surface with FastScapeLib (Braun & Willett, 2013; Yuan et
al., 2019): stream power erosion, hillslope diffusion and, with `fastscape_marine = True`, marine transport and
deposition. It needs `sp_surface_processes = True` and `sediment_layers.txt`, as the other modes: the eroded rock
becomes air and the air below the new surface becomes the active sediment layer. The units are those of FastScape
(m and years). See [numerics.md](numerics.md#erosion-and-sedimentation-with-fastscape) for the coupling.

| Key | Default | Meaning |
|---|---|---|
| `fastscape_kf` | 1e-5 | stream power coefficient of the bedrock (m^(1-2m)/yr) |
| `fastscape_kfsed` | `fastscape_kf` | the same for the sediments |
| `fastscape_m`, `fastscape_n` | 0.4, 1 | exponents of the drainage area and of the slope |
| `fastscape_kd` | 1e-2 | hillslope diffusivity of the bedrock (m2/yr) |
| `fastscape_kdsed` | `fastscape_kd` | the same for the sediments |
| `fastscape_g`, `fastscape_gsed` | 1, `fastscape_g` | fluvial deposition coefficients of the bedrock and of the sediments (0: no deposition) |
| `fastscape_p` | -2 | exponent of the multiple flow direction; below -1.5, single flow direction |
| `fastscape_boundary` | `0101` | the sides y = 0, x = Lx, y = Ly and x = 0 (order of FastScape): 1 is a fixed base level, 0 a free side (periodic when both opposite sides are 0). In 2D only the x sides count |
| `fastscape_dt_max` | 5000 | longest substep of FastScape (yr) |
| `fastscape_strip_width` | 10000 | 2D: width of the river of the profile (m); the drainage area of a point is the length of the profile upstream times this width |
| `fastscape_marine` | False | marine transport below the sea level (`sea_level` above the mean height of the first quarter of the surface in x, as in the sedimentation modes) |
| `fastscape_porosity_silt`, `fastscape_porosity_sand` | 0.63, 0.49 | surface porosities |
| `fastscape_porosity_depth_silt`, `fastscape_porosity_depth_sand` | 1960, 3700 | e-folding depths of the porosity (m) |
| `fastscape_sand_ratio` | 0.5 | fraction of sand in the flux from the continent |
| `fastscape_mixed_layer` | 100 | thickness of the mixed layer at the sea floor (m) |
| `fastscape_kd_silt`, `fastscape_kd_sand` | 300, 500 | marine diffusivities (m2/yr) |

## Not ported yet

These options stop the run (`options of param.txt not ported yet: ...`):

- in 3D (`dimensions = 3`), as in the C++ code: `magmatism` and `periodic_boundary`; of the
  surface processes, only `sp_mode = diffusion` and `fastscape` (the surface tracking in 3D is new).
  `temperature_from_ascii = True` is required, and `export_lithology` and `print_step_files` are ignored.

Ignored on purpose: `binary_output`, `checkered` and `plot_sediment` (outputs are zarr), `external_heat`, `surface_velocity` and
`multigrid` and `aggradation_rate` (read but not used by the C++ code).

`RK4 = Runge-Kutta` (2D and 3D): one fourth-order Runge-Kutta step per time step instead of the 4 Euler substeps. The velocity
varies linearly in time from `Veloc` (t_n) to `Veloc_fut` (t_n+1), the stages that leave the domain are wrapped
(periodic) or stopped at the boundary, and they read the velocity of the advection mesh (stencil 2). T, P and the
strain rate are taken at the position at t_n, with the velocity at the middle of the step. In 2D the surface is
advected in the same way; in 3D the surface keeps its semi-Lagrangian scheme, over the whole step, with the velocity
at the middle of the step. Same scheme as `RK4` of mandyoc_jpm (branches `rk4` and `zarr_outputs`).

Not ported:
`T_initial` other than 0 and 1 (use `temperature_from_ascii`), `binary_output`, restart, and the per-step
`step_*` text files (the trajectories are in `particles_trajectories.zarr`).
