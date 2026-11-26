module Mandyocjll
for f in ("particle_swarm", "thermal", "stokes", "rheology", "surface", "reader", "particle_control", "output", "main")
    include("$f.jl")
end
end
