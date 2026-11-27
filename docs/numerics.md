# Numerical methods

mandyocjll follows the algorithms of Mandyoc (Sacek et al., 2022, and the C++ code) step by step. This page lists
them briefly and then what is different. The bugs of the C++ code found during the port, and how they are handled
here, are in [bugs_mandyoc_cpp.pdf](bugs_mandyoc_cpp.pdf).

## Model

- **Mesh**: two PETSc DMDAs on the same nodes, `da_Thermal` (1 dof) and `da_Veloc` (2 dofs), with the stencils of
  Mandyoc. A third DMDA, used only to advect the particles, has stencil width 2.
- **Particles**: a DMSwarm of type PIC on the cells of `da_Veloc`, with layer, strain and the fields of the rheology.
  The initial positions are regular in each element with a random perturbation (`particles_perturb_factor`); the
  random numbers depend on the particle ID and not on the rank. Particles move with the bilinear velocity of their
  element (Euler, as in Mandyoc); those that leave the domain at open boundaries are removed, and the particle
  control adds or removes particles so that each element keeps between the limits of Mandyoc.
- **Particles to mesh**: nodal density and element viscosity factor (mean of the particles of the element,
  arithmetic or harmonic), as `Swarm2Mesh_2d`.
- **Heat equation**: bilinear elements, Crank-Nicolson (`alpha_thermal` = 0.5), SUPG-like upwinding of Mandyoc,
  radiogenic, adiabatic and shear heating, diffusivity per layer, fixed or insulating boundaries, periodic boundary.
- **Stokes**: Q1-P0 elements with the penalty-free Uzawa scheme of Mandyoc (conjugate gradients on the pressure), the
  velocity system solved by MUMPS (`solver = direct`) or by a Krylov method, free surface stabilization (FSSA),
  boundary conditions per side and component or from `bcv_0.txt`, periodic boundary.
- **Rheology**: the 13 rheologies of Mandyoc, Drucker-Prager plasticity with strain softening, strain healing,
  viscosity limits, sticky air (`sticky_blanket_air`), the viscosity ramp of Gerya (`initial_dynamic_range`) at the
  start and after each change of `scale_bcv.txt`, and the nonlinear (Picard) iterations up to 700 or `Xi_min`.
- **Surface**: the surface of `sp_surface_tracking`, advected with the particles and resampled every step, and the
  surface processes of `sp_mode` (diffusion and the sedimentation modes, with `sediment_layers.txt`,
  `sedimentation_rate.txt` and `base_level.txt`), which turn air particles below the surface into sediment and rock
  above it into air (`sp.cpp`). New: the coupling with FastScape (`sp_mode = fastscape`) and the surface in 3D
  (below).
- **Magmatism**: melt fraction, depletion and last increment of each particle, from a solidus that rises with the
  pressure and the depletion; the viscosity factor of the melt (Lu & Huismans, 2021); the latent heat in the heat
  equation; and, at the output steps, the extraction of the new melt to the surface as a Gaussian (sigma = 10 km)
  above the shallowest melting particle, as `magmatic_layer` (`magmatic.cpp`). As in the C++ code, the melt is
  computed at each particle update, with the temperature of the mesh (bug 22).
- **Inputs over time**: `scale_bcv.txt` (velocity factors) and `multi_velocity` (a new `input_velocity_<n>.txt`,
  also as the boundary values, at each time of `multi_veloc.txt`), both followed by the viscosity ramp if
  `initial_dynamic_range`; `basal_heat`, the heat flow at the base, imposed after each thermal solve on the bottom
  row with the conductivity of the mantle.
- **3D**: the same model with trilinear hexahedra (27 Gauss points), Q1-P0 Stokes with the same Uzawa iterations,
  particles located with the same rule (the rank of the lower left front node owns the element) and the Euler
  advection. As in 2D, the elements have the diffusivity of their layer and the shear heating, which the 3D C++
  code does not have; the 3D FSSA has the consistent weight (it has no effect in the C++ code, bug 29).
- **Time step**: 10% of the time to cross an element at the largest velocity, limited by `dt_max`.

## Differences from the C++ code

| | C++ | mandyocjll |
|---|---|---|
| FSSA weight | the term is not weighted by `theta_FSSA` as intended (5 times larger with g = 10) and uses the `dRho` of the previous assembly (bugs 7 and 8) | consistent weight and the current `dRho`; `fssa_legacy = True` reproduces the C++ term |
| Initial positions of the particles | depend on the number of ranks (`rand_r` per rank) | independent of the number of ranks |
| Removal of particles | depends on the order of the particles in the rank (bug 15) | sorted by position: independent of the number of ranks |
| RK4 | `break` stops the advection of the other particles; reads outside the ghosts (bugs 1, 2) | not available yet (the run stops) |
| Advection mesh | stencil of `da_Veloc` | separate DMDA with stencil 2 |
| Time step with dx ≠ dz | uses the wrong element size (bug 11) | corrected |
| Nonlinear loop with zero initial velocity | stops at the first iteration (`Xi = nan`, bug 10) | continues |
| Rheologies 9 and 19 without nonlinear iterations | minimum viscosity (bug 12) | error, unless the viscosity is constant |
| `-seed` with a layer that does not exist | writes out of the array (bug 16) | error |
| Uzawa or nonlinear loop without convergence | silent (bug 9) | warning |
| Threads | as many OpenMP threads per rank as the BLAS opens | one per rank when `-n` > 1 |
| Surface | a swarm split among the ranks: steps in the initial surface, lost and duplicated points and particles not updated with more than one rank (bugs 18-21) | an array replicated on every rank: independent of the number of ranks |
| Melt extraction | at x = 0 if no particle melted in the last update (bug 25); a tie goes to the highest rank | as the C++ code, but a tie (depths within 1 cm) goes to the smallest x: independent of the number of ranks |
| Latent heat | `dS = 400` whatever `magmatism_dS` is (bug 23) | `magmatism_dS` |
| `nondimensionalization` | equations in units of the depth, maximum viscosity, gravity, mantle density and diffusivity; the time step, `time_max` and the times of the input files mix these units with years (bug 27) | equations in SI units (same solution), outputs in SI, scales in `dimensional_convertion.txt` |
| 3D | no diffusivity per layer, shear heating or basal_heat; FSSA without effect (bug 29); particle seeds per rank | as in 2D: diffusivity per layer, shear heating, basal_heat, consistent FSSA; particles independent of the number of ranks |
| Outputs | text files per step, organized by tapIOca | zarr stores |

With these differences, the results of the two codes are not bit-identical. See [validation.md](validation.md) for
the size of the differences.

## Surface in 3D

The C++ code has no surface in 3D. Here the surface is the height of a regular grid of points
(`surface_particles_per_element` per element in x and in y), replicated on every rank. At each advection substep
it moves with the velocity of the substep by a semi-Lagrangian step: the point that the velocity brings to the grid
point (x, y) starts at (x, y) - dt v, found with two fixed point iterations, and the new height is the height there
(bilinear interpolation) plus dt times the vertical velocity there. To first order in dt this is the 2D scheme (the
points move with the velocity and the surface is resampled on the regular points). The velocity of each point is
computed by the rank that owns its element and summed over the ranks, so the surface does not depend on the number
of ranks. The particles above the surface become air, and the air below it the active sediment layer, with the
height interpolated at the particle (in 2D the C++ code uses the mean of the two neighbors). `sp_mode = diffusion`
is explicit in x and y, with the edges kept at the height of their neighbors, as in 2D.

## Erosion and sedimentation with FastScape

`sp_mode = fastscape` calls FastScapeLib (`Fastscapelib_jll`, the Fortran library of FastScape) on rank 0; the new
surface is broadcast to the other ranks. The grid of FastScape is the grid of the surface points, so no
interpolation is needed: in 3D the regular points; in 2D the profile of regular points (from 0 to `lx`) is extended
to a strip of 4 identical rows spaced by `fastscape_strip_width`, periodic in y, so that each row is a river of
that width (with 3 rows the periodic boundary of FastScape mixes the rows).

At each step, after the advection of the surface:

1. the change of the surface since the end of the last call of FastScape (the tectonics, and the extrusion of the
   melt) is given to FastScape as an uplift rate, `(z - z_previous) / dt`;
2. FastScape runs during dt in `ceil(dt / fastscape_dt_max)` substeps; each substep applies the uplift and then the
   stream power law, the hillslope diffusion and, with `fastscape_marine`, the marine transport;
3. the new height (in 2D the mean of the rows) replaces the surface, and the particles crossed by it change layer.

FastScape keeps its own basement, so it knows where the sediments are (`fastscape_kfsed`, `fastscape_kdsed`,
`fastscape_gsed`). The sides with a fixed base level (`fastscape_boundary`) are not eroded, but they move with the
tectonics. `test/test_surface.jl` checks the steady river profile with uniform uplift against the exact discrete
solution, in the 2D strip and on a 3D grid.
