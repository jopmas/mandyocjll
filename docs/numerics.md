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
  above it into air (`sp.cpp`).
- **Magmatism**: melt fraction, depletion and last increment of each particle, from a solidus that rises with the
  pressure and the depletion; the viscosity factor of the melt (Lu & Huismans, 2021); the latent heat in the heat
  equation; and, at the output steps, the extraction of the new melt to the surface as a Gaussian (sigma = 10 km)
  above the shallowest melting particle, as `magmatic_layer` (`magmatic.cpp`). As in the C++ code, the melt is
  computed at each particle update, with the temperature of the mesh (bug 22).
- **Inputs over time**: `scale_bcv.txt` (velocity factors) and `multi_velocity` (a new `input_velocity_<n>.txt`,
  also as the boundary values, at each time of `multi_veloc.txt`), both followed by the viscosity ramp if
  `initial_dynamic_range`; `basal_heat`, the heat flow at the base, imposed after each thermal solve on the bottom
  row with the conductivity of the mantle.
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
| Outputs | text files per step, organized by tapIOca | zarr stores |

With these differences, the results of the two codes are not bit-identical. See [validation.md](validation.md) for
the size of the differences.
