"""Writes the input files of the example slab_detachment of Mandyoc.

Code cells of examples/slab_detachment/create_input_files_slab.ipynb of Mandyoc (https://github.com/ggciag/mandyoc),
Copyright (c) 2020, MANDYOC Developers, BSD 3-Clause License (see LICENSE_MANDYOC).
"""
import matplotlib

matplotlib.use("Agg")  # the figures of the notebook are not shown

import numpy as np
import matplotlib.pyplot as plt


# total model horizontal extent (m)
Lx = 1000 * 1.0e3
# total model vertical extent (m)
Lz = 660 * 1.0e3
# number of points in horizontal direction
Nx = 1001
# number of points in vertical direction
Nz = 601
# thickness of lithosphere (m)
H_litho = 80 * 1.0e3
# width of slab (m)
W_slab = 80 * 1.0e3
# thickness of slab (m)
H_slab = 250 * 1.0e3

x = np.linspace(0, Lx, Nx)
z = np.linspace(Lz, 0, Nz)
X, Z = np.meshgrid(x, z)


##############################################################################
# Interfaces (bottom first)
##############################################################################
center_of_slab = int((Nx - 1) / 2)
n_slab = int(W_slab / x[1]) # x[1] = dx

slab_range = range(center_of_slab - n_slab // 2, center_of_slab + n_slab // 2)

litho = np.ones(Nx) * H_litho
litho[slab_range] = litho[slab_range] + H_slab

interfaces = {
    "litho": litho
}

# Create the interfaces file
with open("interfaces.txt", "w") as f:
    layer_properties = f"""
        C   1.0         1.0
        rho 3150.0      3300.0
        H   0.0         0.0
        A   1.0e-21     1.96438e-47
        n   1.0         4.0
        Q   0.0e3       0.0e3
        V   0.0e-6      0.0
    """

    for line in layer_properties.split("\n"):
        line = line.strip()
        if len(line):
            f.write(" ".join(line.split()) + "\n")

    # layer interfaces
    data = -1 * np.array(tuple(interfaces.values())).T
    np.savetxt(f, data, fmt="%.1f")

# Plot interfaces
##############################################################################
fig, ax = plt.subplots(figsize=(16, 8))

for label, layer in interfaces.items():
    ax.plot(x, layer, label=f"{label}")

ax.set_xticks(np.arange(0, Lx + 1, 100 * 1.0e3))
ax.set_yticks(np.arange(-Lz, 0 + 1, 10 * 1.0e3))

ax.set_xlim([0, Lx])
ax.set_ylim([Lz, 0])

plt.legend()
plt.savefig("interfaces.png")
plt.close()

##############################################################################
# Parameters file
##############################################################################
params = f"""
nx = {Nx}
nz = {Nz}
lx = {Lx}
lz = {Lz}

# Simulation options
multigrid                           = 1
solver                              = direct        # default is direct [direct/iterative]
denok                               = 1.0e-15       # default is 1.0E-4
particles_per_element               = 40          # default is 81
particles_perturb_factor            = 0.7           # default is 0.5 [values are between 0 and 1]
rtol                                = 1.0e-7        # the absolute size of the residual norm (relevant only for iterative methods), default is 1.0E-5
RK4                                 = Euler         # default is Euler [Euler/Runge-Kutta]
Xi_min                              = 1.0e-7       # default is 1.0E-14
random_initial_strain               = 0.3           # default is 0.0
pressure_const                      = -1.0          # default is -1.0 (not used) - useful only in horizontal 2D models
initial_dynamic_range               = True         # default is False [True/False]
periodic_boundary                   = False         # default is False [True/False]
high_kappa_in_asthenosphere         = False         # default is False [True/False]
K_fluvial                           = 2.0e-7        # default is 2.0E-7
m_fluvial                           = 1.0           # default is 1.0
sea_level                           = 0.0           # default is 0.0
basal_heat                          = 0.0          # default is -1.0

# Surface processes
sp_surface_tracking                 = False         # default is False [True/False]
sp_surface_processes                = False         # default is False [True/False]
sp_d_c                              = 1.0          # default is 0.0
plot_sediment                       = False         # default is False [True/False]
a2l                                 = True          # default is True [True/False]

free_surface_stab                   = True          # default is True [True/False]
theta_FSSA                          = 0.5           # default is 0.5 (only relevant when free_surface_stab = True)

# Time constrains
step_max                            = 5000          # Maximum time-step of the simulation
time_max                            = 25.0e6        # Maximum time of the simulation [years]
dt_max                              = 1.0e6        # Maximum time between steps of the simulation [years]
step_print                          = 10            # Make file every <step_print>
sub_division_time_step              = 0.5           # default is 1.0
initial_print_step                  = 1             # default is 0
initial_print_max_time              = 1.0e6         # default is 1.0E6 [years]

# Viscosity
viscosity_reference                 = 3.0e22        # Reference viscosity [Pa.s]
viscosity_max                       = 1.0e25        # Maximum viscosity [Pa.s]
viscosity_min                       = 1.0e21        # Minimum viscosity [Pa.s]
viscosity_per_element               = constant      # default is variable [constant/variable]
viscosity_mean_method               = harmonic      # default is harmonic [harmonic/arithmetic]
viscosity_dependence                = pressure      # default is depth [pressure/depth]

# External ASCII inputs/outputs
interfaces_from_ascii               = True          # default is False [True/False]
n_interfaces                        = {len(interfaces.keys())}           # Number of interfaces int the interfaces.txt file
variable_bcv                        = False         # default is False [True/False]
temperature_from_ascii              = True          # default is False [True/False]
velocity_from_ascii                 = False         # default is False [True/False]
binary_output                       = False         # default is False [True/False]
sticky_blanket_air                  = False         # default is False [True/False]


print_step_files                    = True          # default is True [True/False]
checkered                           = False         # Print one element in the print_step_files (default is False [True/False])

sp_mode                             = none

geoq                                = on
geoq_fac                            = 100.0

# Physical parameters
temperature_difference              = 0.0
thermal_expansion_coefficient       = 3.28e-5
thermal_diffusivity_coefficient     = 1.0e-6
gravity_acceleration                = 9.81
density_mantle                      = 3300.0
external_heat                       = 0.0e-12
heat_capacity                       = 1250.0

non_linear_method                   = on
adiabatic_component                 = off
radiogenic_component                = off
plasticity                          = off

# Velocity boundary conditions
top_normal_velocity                 = fixed
top_tangential_velocity             = free
bot_normal_velocity                 = fixed
bot_tangential_velocity             = free
left_normal_velocity                = fixed
left_tangential_velocity            = fixed
right_normal_velocity               = fixed
right_tangential_velocity           = fixed

surface_velocity                    = 0.0e-2
multi_velocity                      = False         # default is False [True/False]

# Temperature boundary conditions
top_temperature                     = fixed
bot_temperature                     = fixed
left_temperature                    = fixed
right_temperature                   = fixed

rheology_model                      = 9
T_initial                           = 3
"""
# Create the parameter file
with open("param.txt", "w") as f:
    for line in params.split("\n"):
        line = line.strip()
        f.write(line + "\n")

##############################################################################
# Initial temperature field
##############################################################################

T = np.zeros_like(z) # Temperature
T = np.zeros_like(X) * T[:, None]

# Save the initial temperature file
np.savetxt("input_temperature_0.txt", np.reshape(T, (Nx * Nz)), header="T1\nT2\nT3\nT4")

# Plot temperature field and thermal profile
##############################################################################

fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
ax0.pcolormesh(X / 1.0e3, Z / 1.0e3, T)
ax0.set_ylim(Lz / 1.0e3, 0)
ax0.set_xlabel("km")
ax0.set_ylabel("km")
ax1.set_xlabel("$^\circ$C")

for label, layer in interfaces.items():
    ax0.plot(x / 1.0e3, layer/ 1.0e3, color="C0", label=f"{label}")
    ax1.axhline(y=interfaces[label][0] / 1.0e3, color=f"C0", label=f"{label}")

ax1.plot(T[:, 0], z / 1.0e3, "-k")

ax1.set_ylim(Lz / 1.0e3, 0)
plt.legend()
plt.savefig("temperature_field.png")
plt.close()
