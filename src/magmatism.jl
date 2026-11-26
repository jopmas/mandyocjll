# Extraction of the melt to the surface (magmatic.cpp). The melting itself is in update_loop! (rheology.jl),
# the latent heat in assemble_rhs! (thermal.jl). Include after surface.jl.
using Printf

"Volume of melt already extracted (previous_magmatic_volume)."
mutable struct MeltExtraction
    previous_melt_volume::Float64        # previous_magmatic_volume (m2)
end
MeltExtraction() = MeltExtraction(0.0)

"""
calc_magmatic_extraction and magmatic_extrusive: the melt created since the last extraction (sum of Phi on the
mesh nodes times the element area) is added to the surface as a Gaussian (sigma = 10 km) centered on the
shallowest particle that melted in the last particle update, and the air below the new surface becomes
magmatic_layer. Returns the number of particles that changed layer.
"""
function extract_melt!(ex::MeltExtraction, s::Model, sf::Surface, say)
    m, comm = s.m, MPI.COMM_WORLD
    total = LibPETSc.VecSum(P, s.fields.melt_fraction_field) * element_size_x(m) * element_size_z(m)
    volume = total - ex.previous_melt_volume
    ex.previous_melt_volume = total
    say(@sprintf("Volume of magmatism created: %g m2", volume))
    volume > 0 || return 0
    x, z = s.shallowest_melt
    z_high = MPI.Allreduce(z, max, comm)
    # a tie (MELT_TIE_DEPTH) is won by the smallest x, also between ranks (the C++ code takes the highest rank)
    x_high = MPI.Allreduce(z >= z_high - MELT_TIE_DEPTH ? x : Inf, min, comm)
    say(@sprintf("x_magma_high: %f\nz_magma_high: %f", x_high, z_high))
    sigma = 10000.0                      # (m) fixed in the C++ code
    amplitude = volume / (sigma * sqrt(2 * 3.14159))
    dxs = sf.x[2] - sf.x[1]              # dx_sed: the first point is shifted by epsilon_x
    for j in eachindex(sf.z)
        d = (j - 1) * dxs - x_high
        sf.z[j] += amplitude * exp(-d * d / (2 * sigma * sigma))
    end
    return update_particles_with_surface!(sf, m, s.layers, s.particle_swarm; new_layer = m.magmatic_layer)
end
