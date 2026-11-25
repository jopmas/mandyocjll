# mandyocjll: teste de viabilidade do Mandyoc em Julia

Partículas (DMDA + DMSwarm), equação térmica e Stokes linear em 2D. Ainda não há reologia.

- `src/parameters.jl`, `src/particle_swarm.jl`: nomes do `rename_table.csv` (mandyoc_jpm),
  com o nome antigo no comentário. A coluna `decisao` da tabela ainda estava vazia.
- `src/thermal.jl`: térmica 2D (sem magmatismo, basal_heat, shear heating nem T_initial).
- `src/stokes.jl`: Stokes 2D Q1-P0 com Uzawa, FSSA e partículas -> malha.
- `docs/bugs_mandyoc_cpp.md` e `.pdf`: bugs achados no C++.
- `csrc/`: shim em C que preenche `cell_dm->ops` (opção `shim`); a opção `native` não usa C.
- `results/summary.txt`: resultados.

## Rodar

    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    julia --project=. csrc/build.jl                     # só para a opção shim
    julia --project=. test/mpirun.jl 1 shim closed      # gera a referência; rodar antes
    julia --project=. test/mpirun.jl 4 native closed
    julia --project=. test/mpirun.jl 4 shim periodic    # referência: 1 shim periodic

    julia --project=. test/mpirun.jl 4 test_thermal.jl
    julia --project=. test/mpirun.jl 4 test_stokes.jl
    julia --project=. test/mpirun.jl 1 van_keken.jl 81 6 300
    julia --project=. test/mpirun.jl 4 test_open_boundary.jl shim

Use `test/mpirun.jl`, que chama o MPI do MPI.jl; o `mpiexec` do conda é de outro MPI.

## Diferenças em relação ao Mandyoc

- RK4: em vez de `break` quando um estágio sai do domínio, a posição do estágio é trazida de volta.
- Malha só da advecção, com stencil 2; `da_Veloc` e `da_Thermal` têm o stencil do Mandyoc.
- Posição inicial independente do número de processos (o `rand_r` do Mandyoc depende do rank).
- FSSA com peso consistente (`fssa_legacy = true` repete o C++) e `dRho` atual.

## Desempenho em paralelo (Mac)

O MPICH que vem com o MPI.jl no macOS é compilado só com TCP (`ch4:ofi:tcp`): 1 MB entre dois
processos leva cerca de 1260 µs, contra 24 µs no OpenMPI do conda. O solve do MUMPS com 4 processos
fica várias vezes mais lento. Em 81 x 81 nós, nem o C++ ganha com mais processos: use 1.

Para usar outro PETSc e outro MPI (os dois têm de combinar), crie `LocalPreferences.toml`:

    [MPIPreferences]
    _format = "1.0"
    abi = "OpenMPI"
    binary = "system"
    cclibs = []
    libmpi = "<env>/lib/libmpi.dylib"
    mpiexec = "<env>/bin/mpiexec"
    preloads = []
    preloads_env_switch = ""

    [PETSc]
    library_path = "<env>/lib/libpetsc.dylib"
    PetscScalar = "Float64"
    PetscInt = "Int32"

e adicione `MPIPreferences` ao `Project.toml`. Testado com PETSc 3.24.5 do conda (PetscInt de 32 bits),
só com `location = :native`: o shim precisa ser compilado contra o PETSc usado. Mesmo assim, o Julia
pode baixar variantes de algumas dependências ao carregar.
`test/bench_stokes.jl` mede o tempo por etapa.
