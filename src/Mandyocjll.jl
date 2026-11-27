module Mandyocjll
for f in ("particle_swarm", "thermal", "thermal3d", "stokes", "stokes3d", "rheology", "rheology3d", "fastscape", "surface", "magmatism", "reader", "particle_control", "output", "main")
    include("$f.jl")
end
end
