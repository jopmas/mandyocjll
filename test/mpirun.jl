# julia --project=. test/mpirun.jl N [argumentos do teste]   (usa o MPI do MPI.jl, não o do conda)
using MPI
script, rest = endswith(get(ARGS, 2, ""), ".jl") ? (ARGS[2], ARGS[3:end]) : ("test_particle_swarm.jl", ARGS[2:end])
run(`$(mpiexec()) -n $(ARGS[1]) $(Base.julia_cmd()) --project=. $(joinpath(@__DIR__, script)) $(rest)`)
