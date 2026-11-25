---
title: "Bugs encontrados no Mandyoc (C++)"
subtitle: "Levantados ao portar o código para Julia (mandyocjll)"
date: "25 de novembro de 2026"
lang: pt-BR
geometry: margin=2.2cm
fontsize: 10pt
colorlinks: true
---

Código analisado: `mandyoc_jpm`, commit `64dc7d1`. Foi compilado com PETSc 3.24.5 no macOS (arm64), em cópias temporárias. **Nenhum repositório do Mandyoc foi alterado.** As linhas citadas são desse commit.

**Confirmado** significa que o problema foi reproduzido executando o C++. **Por leitura** significa que ele só foi visto no código.

| Nº | Problema | Quando ocorre | Evidência | Situação no Julia |
|:-:|:-----------------|:------------|:---------|:----------------|
| 1 | `break` no RK4 para a advecção das outras partículas | RK4, 2D, fluxo pelas bordas | Confirmado | Resolvido e testado |
| 2 | RK4 lê fora dos ghosts | RK4, 2D, em paralelo | Confirmado (SEGV) | Resolvido e testado |
| 3 | `MPI_Bcast` de 8 bytes em variável de 4 | Em paralelo | Confirmado (ASan) | Não se aplica ainda |
| 4 | Contorno periódico aborta com 1 processo em x | Periódico | Confirmado | Resolvido e testado |
| 5 | RK4 não move a superfície | RK4 com rastreamento da superfície | Por leitura | Ainda não implementado |
| 6 | `local_divV` é perdido e vira `NULL` | Magmatismo | Por leitura | Ainda não implementado |
| 7 | Peso do FSSA não é `theta_FSSA` (5 vezes maior com g = 10) | FSSA ativo (padrão), 2D | Confirmado | Resolvido e testado |
| 8 | FSSA usa o `dRho` da montagem anterior | FSSA ativo (padrão), 2D | Confirmado | Resolvido e testado |
| 9 | Uzawa para no limite de iterações sem aviso | Sem convergência | Confirmado | Emite aviso |

O RK4 já é marcado como `not working yet!` em `header.h:50`, e todos os testes usam Euler. Os itens 1, 2 e 5 só afetam quem ativar `RK4 = Runge-Kutta`.

\newpage

# 1. O `break` no RK4 interrompe o laço das partículas

**Onde:** `src/DMSwarm_move.cpp`, função `moveSwarm`, ramo `RK4==1` (três ocorrências).

Se um estágio intermediário sai do domínio, o `break` encerra o laço sobre as partículas, e não só o tratamento daquela. As partículas seguintes **daquele processo** não são movidas nessa passagem. O resultado depende do número de processos.

**Medição:** `continental_rift` reamostrado para 201 × 51 nós, 25 passos. Uma passagem é uma chamada de `moveSwarm`.

| Processos | Passagens afetadas | Partículas não movidas (pior passagem) | Média |
|:-:|:-:|:-:|:-:|
| 1 | 21 de 201 | 74% | 7,3% |
| 2 | 21 de 195 | 37% | 3,7% |
| 4 | 21 de 213 | 18% | 1,7% |

**Sugestão:** trocar `break` por `continue`, para pular só a partícula que saiu:

```cpp
if (xB>Lx || xB<0 || zB>0 || zB<-depth){
    array[2*p  ] = xB;
    array[2*p+1] = zB;
    continue;   // era break
}
```

Com isso (e a correção 2), o mesmo caso terminou com 399.863 partículas em 1, 2 e 4 processos. Antes, eram 400.151, 400.155 e 400.057. Com contorno periódico o `continue` não basta: é preciso aplicar o ajuste de `Lx` à posição antes de pular. Isso não foi testado no C++.

**Julia:** resolvido. A posição de cada estágio volta ao domínio (periódico) ou para na borda (fechado); quem cruza uma borda aberta sai, sem afetar as demais. Teste com fluxo uniforme saindo do domínio, em 1, 2 e 4 processos: saíram exatamente as 8748 partículas esperadas e as outras 30.132 andaram a distância exata (erro de 9e-11 m). Com parede fechada e periódico, em 1, 2, 4 e 6 processos, nenhuma partícula se perdeu e as posições são idênticas às da rodada serial.

# 2. RK4 lê velocidades fora da região de ghosts

**Onde:** `src/DMv.cpp` (`stencil_width = 1`) e `src/DMSwarm_move.cpp` (estágios B, C e D).

Um estágio pode cair no elemento vizinho, que pertence a outro processo. Esse elemento usa um nó além da única camada de ghosts, e o acesso `VV[k][i]` sai do array.

**Medição:** van Keken com 2 e 4 processos termina em **SEGV**; com Euler, termina normalmente. No caso de 2 processos, o registro logo antes da falha mostra o nó `k=42` exigido (`k=41` e `k+1`), com ghosts até 41. No rift, houve 405 acessos inválidos com 2 processos e 1691 com 4, **sem nenhum aviso**.

**Não basta** usar `stencil_width = 2` em `da_Veloc`: a execução aborta. Provável causa (por leitura): `DMDAGetElementCorners_2d` calcula o primeiro elemento supondo stencil 1 sem contorno periódico (`si+1+periodic_boundary`).

**Sugestão:** uma DMDA só para a advecção, com a mesma partição e stencil 2:

```cpp
DMDAGetInfo(da_Veloc,NULL,&M,&N,NULL,&m,&n,NULL,&dof,NULL,&bx,&bz,NULL,NULL);
DMDAGetOwnershipRanges(da_Veloc,&lx,&lz,NULL);
DMDACreate2d(PETSC_COMM_WORLD,bx,bz,DMDA_STENCIL_BOX,M,N,m,n,dof,2,
             lx,lz,&da_advection);
DMSetUp(da_advection);
DMCreateLocalVector(da_advection,&local_V_advection);
// em moveSwarm (2D), usar da_advection e local_V_advection para VV
DMGlobalToLocal(da_advection,Veloc_weight,INSERT_VALUES,local_V_advection);
```

Testado numa cópia do C++: os acessos inválidos foram a zero em 1, 2 e 4 processos. Vale para passos em que a partícula anda menos de um elemento.

**Julia:** resolvido com a mesma solução (`create_advection_mesh`), testado com 1, 2, 4 e 6 processos.

# 3. `MPI_Bcast` com tipo errado

**Onde:** `src/reader.cpp:495`.

`particles_per_ele` é `PetscInt` (4 bytes na compilação testada), mas é transmitido como `MPI_LONG` (8 bytes). O AddressSanitizer acusou a leitura além da variável (`global-buffer-overflow`). Nos processos que recebem, quatro bytes da memória vizinha podem ser sobrescritos. Não foi visto efeito nos resultados.

**Sugestão:**

```cpp
MPI_Bcast(&particles_per_ele,1,MPIU_INT,0,PETSC_COMM_WORLD);
```

Uma verificação automática das chamadas de `MPI_Bcast` encontrou só este caso com tamanho errado. Os `PetscBool` enviados como `MPI_C_BOOL` estão corretos no PETSc 3.24.5, em que `PetscBool` é `bool`. A correção não foi executada.

**Julia:** não se aplica ainda, pois não há leitor de `param.txt` nem `Bcast`.

# 4. Contorno periódico aborta com um processo em x

**Onde:** `src/DM1_2d.cpp:1001`, função `mean_value_periodic_boundary_2d`.

`DMLocalToGlobal` com `INSERT_VALUES` não é aceito pelo PETSc quando há periodicidade e só um processo nessa direção. Com 1 processo, o programa para com *"Available only for boundary none or with parallelism in x direction"*. Com 2 e 4 processos, a mesma entrada (sem arquivo externo de velocidade) rodou. A montagem em `AssembleF_Thermal_2d` usa o mesmo padrão.

**Sugestão:** escrever no vetor global em vez de usar `DMLocalToGlobal(..., INSERT_VALUES, ...)`. Esboço, não executado:

```cpp
DMDAVecGetArray(da,local_F,&ff);   // ghosts, só leitura
DMDAVecGetArray(da,F,&gg);         // vetor global
for (k=sz; k<sz+mmz; k++) for (i=sx; i<sx+mmx; i++) {
    if (i==Nx-1) gg[k][i] = (ff[k][Nx-1] + ff[k][Nx])/esc;
    if (i==0)    gg[k][i] = (ff[k][-1]   + ff[k][0])/esc;
}
```

Alternativa mínima: recusar, com mensagem clara, `periodic_boundary` com um processo em x.

**Julia:** resolvido nas partes portadas (partículas, térmica e Stokes). O elemento da borda soma direto no nó periódico com `ADD_VALUES`. Testado com 1, 2 e 4 processos.

# 5. Com RK4, a superfície não é advectada

**Onde:** `src/main.cpp`, laço de tempo. `sp_move_surface_swarm` só é chamada no ramo Euler.

**Sugestão** (não executada)**:** chamar também no ramo RK4, com o passo completo:

```cpp
if (RK4==1){
    VecCopy(Veloc_fut,Veloc_weight);
    ierr = moveSwarm(dimensions, dt_calor_sec);
    if (dimensions == 2 && sp_surface_tracking) {
        ierr = sp_move_surface_swarm(dimensions, dt_calor_sec); CHKERRQ(ierr);
    }
}
```

Convém conferir como `sp_move_surface_swarm` integra no tempo. **Julia:** a superfície ainda não foi portada.

# 6. `local_divV` é perdido na inicialização do magmatismo

**Onde:** `src/DM1_2d.cpp:799-817` e `src/DMT.cpp:217`.

`create_thermal` cria o vetor com `DMCreateLocalVector`. Depois, `Thermal_init_2d` sobrescreve o ponteiro com `DMGetLocalVector` (o vetor original vaza) e o devolve com `DMRestoreLocalVector`, que o deixa `NULL`. Na busca feita, `local_divV` não é usado em outro lugar; por ora o efeito é só o vazamento, mas um uso futuro receberia `NULL`.

**Sugestão** (não executada)**:** remover o par `DMGetLocalVector`/`DMRestoreLocalVector` e usar o vetor já criado. **Julia:** o magmatismo ainda não foi portado.

# 7. O peso do FSSA não corresponde a `theta_FSSA`

**Onde:** `src/DM1_v.cpp:284`, função `AssembleA_Veloc_2d`.

```cpp
traction_bottom = theta_FSSA*dt_calor_sec*rho_mean_bottom*gravity*dx_const/2.0;
```

Há três fatores envolvidos:

1. `rho_mean_bottom` vem de `dRho`, que já é a força de volume ($-g\rho$, em `calc_drho`). Multiplicar por `gravity` de novo põe um $g$ a mais.
2. `Ke_veloc`, `VCe` e `VfMe` omitem o jacobiano 1/4 e são 4 vezes maiores. O termo do FSSA não tem esse fator, e fica 4 vezes menor em relação a eles.
3. O termo é somado aos nós de baixo e aos de cima do elemento, então a mesma interface é contada duas vezes.

O peso efetivo é $\theta\,g/2$: **5 vezes** o nominal com $g = 10$ e metade com $g = 1$ (adimensional).

**Medição:** superfície perturbada sob uma camada leve e pouco viscosa. A resposta da velocidade define o peso efetivo, por $w_0/w - 1 = \theta\,dt/\tau$. O cálculo foi feito em Julia, com as duas fórmulas, e deu o mesmo com 1, 2 e 4 processos:

| Fórmula | Peso nominal | Peso efetivo medido |
|:------------------|:-:|:-:|
| Como no C++ | 0,5 | 2,481 |
| Consistente (sugestão) | 0,5 | 0,496 |

A fórmula "como no C++" é a do executável: num solve do van Keken gravado pelo C++, as velocidades coincidem até 7e-11 (diferença relativa).

**Sugestão** (mantendo as matrizes como estão):

```cpp
// dRho já contém -g*rho; 1/4 por nó, vezes o fator 4 das matrizes
traction_bottom = theta_FSSA*dt_calor_sec*rho_mean_bottom*dx_const;
```

Executada numa cópia do C++, junto com a correção 8: as velocidades coincidem com as do Julia até 5e-13. **Ela muda resultados com FSSA**, inclusive, provavelmente, os arquivos de referência dos testes, que não foram rodados. Vale confirmar com os autores se o peso atual foi intencional. O 3D usa outra expressão e não foi analisado.

**Julia:** resolvido. `fssa_legacy = true` repete o C++ para comparação.

# 8. O FSSA usa o `dRho` da montagem anterior

**Onde:** `src/DMv.cpp`, função `build_veloc` (`calc_drho()` na linha 541).

`AssembleA_Veloc_2d`, que usa `dRho` no FSSA, é chamada antes de `calc_drho()`. A matriz usa a densidade de uma montagem anterior; a força (`AssembleF_Veloc_2d`) usa a atual.

**Medição:** mesmo solve do van Keken. Com o `dRho` atual no FSSA, o Julia difere do C++ em 4,1e-4; com o `dRho` gravado antes de `AssembleA`, em 7e-11.

**Sugestão:** chamar `calc_drho()` entre `Swarm2Mesh_2d()` e `AssembleA_Veloc_2d`. Executada numa cópia do C++ junto com a correção 7.

**Julia:** resolvido; matriz e força usam o mesmo `dRho`.

# 9. O Uzawa para no limite de iterações sem avisar

**Onde:** `src/DMv.cpp`, função `solve_veloc` (`maxk = 400`).

O laço termina sem mensagem quando não atinge `denok_min`, e a simulação segue com a velocidade não convergida. Foi visto no rift reamostrado, com tolerância apertada e reologia linear: 399 iterações, `denok` = 1,9e-4 e execução terminada normalmente. O caso é extremo e possivelmente mal posto (o balanço de fluxo nas bordas não foi verificado); o ponto é a falta de aviso.

**Sugestão** (não executada):

```cpp
if (denok>denok_min) PetscPrintf(PETSC_COMM_WORLD,
    "WARNING: Uzawa did not converge, denok = %lg\n",denok);
```

**Julia:** emite um aviso.

# Observações que não são bugs do Mandyoc

- **Jacobiano do elemento térmico:** `montaKeThermal_*` omite o fator 1/4. Como ele aparece em todos os termos, a solução não muda.
- **PETSc.jl 0.5.6:** a estrutura `MatStencil` está na ordem `i,j,k,c`, e o header do C usa `k,j,i,c`. Por isso o Julia usa `MatSetValuesLocal`. Visto por leitura.

# Térmica 2D: comparação Julia × C++

Um passo térmico real do rift (201 × 51 nós, com advecção, calor radiogênico e adiabático) foi gravado pelo C++ e repetido em Julia:

| Caso | 1 processo | 2 processos | 4 processos | Mudança de T no passo |
|:----------|:-:|:-:|:-:|:-:|
| Fechado | 1,4e-9 °C | 2,5e-12 °C | 2,8e-9 °C | até 6,6 °C |
| Periódico | 8,4e-10 °C | 2,3e-12 °C | 7,0e-10 °C | até 7,2 °C |

Os valores são a diferença máxima entre Julia e C++. O C++ rodou com 2 processos.

# Stokes 2D: comparação Julia × C++

Um solve gravado pelo C++ (van Keken, 81 × 81 nós, tolerância 1e-32) foi repetido em Julia com `fssa_legacy = true` e o `dRho` defasado do item 8. Diferença relativa máxima da velocidade:

| Caso | 1 processo | 2 processos | 4 processos |
|:----------------------|:-:|:-:|:-:|
| Sem FSSA | 1,2e-13 | 1,1e-13 | 1,1e-13 |
| Com FSSA | 6,5e-11 | 1,5e-11 | 3,7e-11 |
| Periódico, com FSSA | 5,2e-13 | 5,3e-13 | 5,8e-13 |

Partículas e Stokes acoplados (van Keken caso 1a, 80 × 80 elementos, 36 partículas por elemento, sem FSSA):

| Código | Pico de vrms | Tempo do pico | Taxa de crescimento |
|:---------|:-:|:-:|:-:|
| Julia (RK4) | 0,003103 | 210,6 | 0,01191 |
| C++ (Euler) | 0,003111 | 213,8 | 0,01164 |

A diferença cai com o passo de tempo: com o passo 4 vezes menor, as taxas são 0,01200 (Julia) e 0,01191 (C++). Os valores não foram conferidos com a tabela do artigo original.
