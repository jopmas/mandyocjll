# mandyocjll: teste de viabilidade do Mandyoc em Julia

Partículas (DMDA + DMSwarm), equação térmica, Stokes não linear e reologia em 2D e 3D.

- `src/parameters.jl`, `src/particle_swarm.jl`: nomes do `rename_table.csv` (mandyoc_jpm),
  com o nome antigo no comentário. A coluna `decisao` da tabela ainda está vazia.
- `src/thermal.jl`, `src/thermal3d.jl`: térmica 2D e 3D.
- `src/stokes.jl`, `src/stokes3d.jl`: Stokes 2D e 3D Q1-P0 com Uzawa e FSSA.
- `src/rheology3d.jl`: partículas, partículas -> malha, pressão e modelo acoplado em 3D.
- `src/rheology.jl`: 13 reologias, plasticidade, deformação, laço não linear, rampa de Gerya, ar e laço de tempo.
- `src/reader.jl`: lê `param.txt` e `interfaces.txt` do Mandyoc.
- `src/particle_control.jl`: adição e remoção de partículas.
- `src/output.jl`: saída em zarr v2, com os nomes, dimensões, unidades e tipos dos `.nc` do tapIOca.
- `bin/mandyocjll`, `src/main.jl`: rodam um cenário; `tools/zarr2nc.py` converte os `.zarr` nos `.nc` do tapIOca.
- `docs/bugs_mandyoc_cpp.md` e `.pdf`: bugs achados no C++ (`python3 docs/build_pdf.py` gera o PDF).
- `docs/rename_table.csv`: cópia da planilha de nomes.
- `csrc/`: shim em C que preenche `cell_dm->ops` (opção `shim`); a opção `native` não usa C.
- `results/summary.txt`: resultados dos primeiros testes.
- `examples/`: exemplos do Mandyoc; `docs/`: documentação.

## Rodar

    julia --project=. -e 'using Pkg; Pkg.instantiate()'
    julia --project=. csrc/build.jl                     # só para a opção shim
    julia --project=. test/mpirun.jl 1 shim closed      # gera a referência; rodar antes
    julia --project=. test/mpirun.jl 4 native closed
    julia --project=. test/mpirun.jl 4 shim periodic    # referência: 1 shim periodic

    julia --project=. test/mpirun.jl 4 test_thermal.jl
    julia --project=. test/mpirun.jl 4 test_stokes.jl
    julia --project=. test/mpirun.jl 1 van_keken.jl 81 6 300
    julia --project=. test/test_rheology.jl <mandyoc>/src <petsc>/include
    julia --project=. test/mpirun.jl 4 compare_cpp_model.jl <pasta do caso, com a saída binária do C++> <passos> add_remove
    julia --project=. test/mpirun.jl 4 test_open_boundary.jl shim
    julia --project=. test/mpirun.jl 4 test_particle_control.jl native

Os testes usam `test/mpirun.jl`, que chama o MPI em uso pelo MPI.jl.

## Diferenças em relação ao Mandyoc

- RK4: em vez de `break` quando um estágio sai do domínio, a posição do estágio é trazida de volta.
- Malha só da advecção, com stencil 2; `da_Veloc` e `da_Thermal` têm o stencil do Mandyoc.
- Posição inicial independente do número de processos (o `rand_r` do Mandyoc depende do rank).
- FSSA com peso consistente (`fssa_legacy = true` repete o C++) e `dRho` atual.

## Documentação, exemplos e compilação

A documentação (em inglês) está em [`docs/`](docs/README.md): instalação, uso, parâmetros, saídas, métodos numéricos,
validação e desenvolvimento. Os exemplos do Mandyoc estão em [`examples/`](examples/README.md), e os mesmos exemplos
desenhados em SVG para o tapIOca estão em [`examples/examples_svg/`](examples/examples_svg/README.md).

`make` instala os pacotes, compila o shim e gera a imagem compilada (`bin/mandyocjll.dylib`, PackageCompiler, alguns
minutos). Com ela, uma rodada começa em ~5 s em vez de 13–15 s. Depois de mudar `src/`, os pacotes ou o
`LocalPreferences.toml`, rode `make image` de novo; até lá, o `mandyocjll` avisa e roda sem a imagem.

## Rodar um cenário

`bin/mandyocjll` é o executável. Rode na pasta do cenário, como o `mandyoc`:

    ~/opt/mandyocjll/bin/mandyocjll -n 4 -seed 0,5,8 -strain_seed 0.0,1.0,1.0
    julia ~/opt/mandyocjll/bin/mandyocjll -n 4 -seed 0,5,8 -strain_seed 0.0,1.0,1.0
    mpiexec -n 4 ~/opt/mandyocjll/bin/mandyocjll -seed 0,5,8 -strain_seed 0.0,1.0,1.0

As três formas são equivalentes; na última, o `mpiexec` tem de ser o do MPI em uso (`LocalPreferences.toml`).
`mandyocjll -h` lista as opções. Cada variável vira `<variável>.zarr`; `tools/zarr2nc.py <pasta>` gera os `.nc` que
o tapIOca lê hoje. Com mais de um processo, o `mandyocjll` usa 1 thread de OpenMP por processo
(`OMP_NUM_THREADS=1`): com as threads padrão do BLAS do conda, 4 processos ficavam 9 vezes mais lentos.
Opções do `param.txt` ainda não portadas interrompem a rodada (ver `docs/parameters.md`).

No tapIOca (branch `feat/svg-scenario-inputs`), `code = mandyocjll` na forma `setup` do SVG, ou `--code mandyocjll`,
gera `run_linux.sh`, `run_mac.sh` ou `run_hypatia.sh` que chamam o mandyocjll.

## Ainda não portado

RK4, restart, `T_initial` diferente de 0 e 1, e em 3D (como no C++) magmatismo, superfície e fronteira periódica.

## Desempenho em paralelo (Mac)

O MPICH que vem com o MPI.jl no macOS é compilado só com TCP (`ch4:ofi:tcp`): 1 MB entre dois
processos leva cerca de 1260 µs, contra 24 µs no OpenMPI do conda. O solve do MUMPS com 4 processos
fica várias vezes mais lento. Em 81 x 81 nós, o solve direto
não acelera com mais processos, nem no C++; aqui, 1 processo é o mais rápido.

Para usar outro PETSc e outro MPI (os dois têm de combinar), crie `LocalPreferences.toml` (exemplo com `<env>` = pasta do ambiente conda):

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

`MPIPreferences` já está nos `[extras]` do `Project.toml`. O arquivo é específico de cada máquina e está no
`.gitignore`; sem ele, valem o MPICH e o PETSc_jll. Depois de criar, mudar ou apagar o arquivo,
recompile o shim: `julia --project=. csrc/build.jl`. Rode também `Pkg.instantiate()` e `Pkg.precompile()` (baixa ~190 MB de
variantes das dependências): sem elas, cada início do Julia refaz a precompilação (~16 s) e mostra um
erro sobre `SCALAPACK32`, inofensivo. Testado com PETSc 3.24.5 do ambiente conda `tapioca` (PetscInt de 32 bits).
`test/bench_stokes.jl` mede o tempo por etapa.
