module Mandyocjll
for f in ("particle_swarm", "thermal", "stokes", "rheology", "surface", "magmatism", "reader", "particle_control", "output", "main")
    include("$f.jl")
end
end
