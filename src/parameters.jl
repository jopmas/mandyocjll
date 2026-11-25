# Nomes novos: rename_table.csv do mandyoc_jpm. O nome antigo fica no comentário
# (nomes sem comentário são iguais aos do Mandyoc ou novos neste projeto).
Base.@kwdef struct ModelParameters
    Nx::Int = 101                              # Nx
    Nz::Int = 31                               # Nz
    Lx::Float64 = 1000.0e3                     # Lx
    domain_depth::Float64 = 300.0e3            # depth
    particles_per_element_x::Int = 9           # nx_ppe
    particles_per_element_z::Int = 9           # nz_ppe
    particles_perturb_factor::Float64 = 0.5    # particles_perturb_factor
    random_seed::UInt64 = 1                    # seed
    use_runge_kutta_4::Bool = true             # RK4
    periodic_boundary::Bool = false            # periodic_boundary
    boundary_offset::Float64 = 1.0e-7          # epsilon_x
    advection_stencil_width::Int = 2           # novo: malha só da advecção (ver bug do RK4)
    thermal_time_theta::Float64 = 0.5          # alpha_thermal
    thermal_diffusivity::Float64 = 1.0e-6      # kappa
    thermal_expansion_coefficient::Float64 = 3.28e-5  # alpha_exp_thermo
    gravity_acceleration::Float64 = 10.0       # gravity
    heat_capacity::Float64 = 1250.0            # c_heat_capacity
    use_radiogenic_heating::Bool = false       # WITH_RADIOGENIC_H
    use_adiabatic_heating::Bool = false        # WITH_ADIABATIC_H
    high_diffusivity_in_asthenosphere::Bool = false   # high_kappa_in_asthenosphere
    temperature_bc_top::Bool = true            # bcT_top   (true = fixed)
    temperature_bc_bottom::Bool = true         # bcT_bot
    temperature_bc_left::Bool = false          # bcT_left
    temperature_bc_right::Bool = false         # bcT_right
    rtol::Float64 = 1.0e-5                     # rtol
    mantle_density::Float64 = 3300.0           # RHOM
    uzawa_residual_tolerance::Float64 = 1.0e-4 # denok_min
    use_direct_solver::Bool = true             # direct_solver
    free_surface_stabilization::Bool = true    # free_surface_stab
    fssa_theta::Float64 = 0.5                  # theta_FSSA
    fssa_legacy::Bool = false                  # novo: true repete o termo do C++ (ver docs)
    harmonic_viscosity_mean::Bool = false      # visc_harmonic_mean
    velocity_bc_top_normal::Bool = true        # bcv_top_normal   (true = fixed)
    velocity_bc_top_tangential::Bool = false   # bcv_top_slip
    velocity_bc_bottom_normal::Bool = true     # bcv_bot_normal
    velocity_bc_bottom_tangential::Bool = false # bcv_bot_slip
    velocity_bc_left_normal::Bool = true       # bcv_left_normal
    velocity_bc_left_tangential::Bool = false  # bcv_left_slip
    velocity_bc_right_normal::Bool = true      # bcv_right_normal
    velocity_bc_right_tangential::Bool = false # bcv_right_slip
end
element_size_x(m::ModelParameters) = m.Lx / (m.Nx - 1)            # dx_const
element_size_z(m::ModelParameters) = m.domain_depth / (m.Nz - 1)  # dz_const
particles_per_element(m::ModelParameters) =                       # particles_per_ele
    m.particles_per_element_x * m.particles_per_element_z
stencil_width(m::ModelParameters) = m.periodic_boundary ? 2 : 1   # stencil_width (create_veloc)
