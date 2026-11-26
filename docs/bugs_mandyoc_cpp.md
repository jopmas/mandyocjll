---
title: "Bugs encontrados no Mandyoc (C++)"
subtitle: "Levantados ao portar o código para Julia (mandyocjll)"
date: "25 de novembro de 2026"
lang: pt-BR
geometry: margin=2.2cm
fontsize: 10pt
colorlinks: true
---

Código analisado: `mandyoc_jpm`, commit `64dc7d1`. Foi compilado com PETSc 3.24.5 no macOS (arm64), em cópias temporárias. **Nenhum repositório do Mandyoc foi alterado.** As linhas citadas são desse commit. As medições dos itens 10 a 14 e da reologia usaram o `ggciag/mandyoc`, commit `956af67` (o `main` já incorporado ao fork).

**Confirmado** significa que o problema foi reproduzido executando o C++. **Por leitura** significa que ele só foi visto no código.

| Nº | Problema | Quando ocorre | Evidência | Situação no Julia |
|:-:|:-----------------|:------------|:---------|:----------------|
| 1 | `break` no RK4 para a advecção das outras partículas | RK4, 2D, fluxo pelas bordas | Confirmado | Resolvido e testado |
| 2 | RK4 lê fora dos ghosts | RK4, 2D, em paralelo | Confirmado (SEGV) | Resolvido e testado |
| 3 | `MPI_Bcast` de 8 bytes em variável de 4 | Em paralelo | Confirmado (ASan) | Não se aplica |
| 4 | Contorno periódico aborta com 1 processo em x | Periódico | Confirmado | Resolvido e testado |
| 5 | RK4 não move a superfície | RK4 com rastreamento da superfície | Por leitura | Ainda não implementado |
| 6 | `local_divV` é perdido e vira `NULL` | Magmatismo | Por leitura | Ainda não implementado |
| 7 | Peso do FSSA não é `theta_FSSA` (5 vezes maior com g = 10) | FSSA ativo (padrão), 2D | Confirmado | Resolvido e testado |
| 8 | FSSA usa o `dRho` da montagem anterior | FSSA ativo (padrão), 2D | Confirmado | Resolvido e testado |
| 9 | Uzawa para no limite de iterações sem aviso | Sem convergência | Confirmado | Emite aviso |
| 10 | Laço não linear para na 1ª iteração (`Xi = nan`) | Velocidade inicial nula, sem rampa | Confirmado | Resolvido |
| 11 | `Calc_dt_calor` troca `dx` e `dz` | dx diferente de dz | Confirmado | Resolvido |
| 12 | Reologias 9 e 19 viram viscosidade mínima | `non_linear_method = off` | Confirmado | Interrompe com erro |
| 13 | Índice não inicializado em `Swarm_add_remove_2d` | Elemento sem partículas | Por leitura | Resolvido |
| 14 | `sticky_blanket_air` pula elementos na divisa dos processos | Paralelo | Por leitura | Resolvido |
| 15 | Remoção de partículas depende da ordem no processo | Paralelo | Por leitura | Resolvido e testado |
| 16 | `-seed` não valida a camada | Camada inexistente | Por leitura | Resolvido e testado |
| 17 | Padrão de `viscosity_per_element` é uma opção removida | `param.txt` sem a chave | Por leitura | Interrompe com mensagem clara |

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

**Julia:** não se aplica; cada processo lê o `param.txt`, sem `Bcast`.

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

# 10. O laço não linear para na primeira iteração quando a velocidade inicial é nula

**Onde:** `src/veloc_total.cpp`.

`Xi` compara a velocidade anterior com a nova. Se a anterior é nula, `sig1 = 0` e `Xi = 0/0 = nan`. Como `nan > Xi_min` é falso, o laço termina com uma iteração, como se tivesse convergido.

**Medição:** `Schmalholz_2011` reamostrado (101 × 61 nós), com `initial_dynamic_range = False`: o C++ imprime `Xi = nan 0` e segue. Com a rampa inicial ligada, as chamadas seguintes partem de velocidade não nula e o efeito é menor.

**Sugestão** (não executada):

```cpp
Xi = (sig1>0.0 && sig2>0.0) ? 1.0 - vivi/PetscSqrtReal(sig1*sig2) : 1.0;
```

O laço também para em 700 iterações sem aviso. **Julia:** resolvido; avisa se não convergir.

# 11. `Calc_dt_calor` usa a dimensão errada do elemento

**Onde:** `src/main.cpp:405`.

`VecMax` grava o índice do máximo em `ind_v_min`, e `VecMin` grava o do mínimo em `ind_v_max`. Com isso, `ind_v_mod` aponta para o outro extremo, e a escolha entre `dx` e `dz` pode usar a componente errada.

**Medição:** slab reamostrado, passo 1. A maior velocidade em módulo é vertical (w = -6,5e-11 m/s), e o passo deveria usar dz = 11 km: 1,07e6 anos, limitado a 1,0e6 por `dt_max`. O C++ usou dx = 10 km: 9,73e5 anos. O fator de erro é dx/dz ou dz/dx.

**Sugestão** (não executada):

```cpp
VecMax(Veloc_fut,&ind_v_max,&max_v);
VecMin(Veloc_fut,&ind_v_min,&min_v);
```

**Julia:** resolvido.

# 12. Reologias 9 e 19 com `non_linear_method = off`

**Onde:** `src/calc_visc.cpp`.

O cálculo dessas reologias só acontece com `WITH_NON_LINEAR == 1`. Com `off`, a viscosidade fica em `visc_MIN`, sem aviso.

**Medição:** rift reamostrado (reologia 9) com `non_linear_method = off`: viscosidade mínima e mediana de 1e18 Pa s, o valor de `viscosity_min`. O máximo, 5e19, é o piso de `sticky_blanket_air` (50 vezes o mínimo). A medição do item 9 tinha essa mesma configuração.

**Sugestão:** interromper com mensagem ao ler essa combinação. **Julia:** interrompe com erro.

# 13. Índice não inicializado em `Swarm_add_remove_2d`

**Onde:** `src/DMSwarm_move.cpp`.

Se um elemento fica sem nenhuma partícula (`cont_p == 0`), `p_prox` nunca recebe valor, e `p_prox_total` é usado para ler `layer_array`, `rarray` e `strain_fac`. O resultado é indefinido.

**Sugestão** (não executada): se `cont_p == 0`, procurar a partícula mais próxima nos elementos vizinhos ou interromper com mensagem. **Julia:** resolvido; usa a partícula mais próxima dos elementos vizinhos. Esse caminho não foi exercitado num teste.

# 14. `sticky_blanket_air` depende da divisão entre processos

**Onde:** `src/DMSwarm2mesh.cpp:579`.

O laço vai só até `sx+mmx-2` e `sz+mmz-2` em cada processo, então a última linha de elementos de cada processo é pulada. Em z, se a divisa coincidir com a interface do ar, a camada não é aplicada ali. Em x, a atribuição a `[i+1]` cobre a maior parte. No rift com 2 e 4 processos o resultado não mudou.

**Sugestão** (não executada): percorrer todos os elementos do processo e usar os ghosts de `qq_rho`. **Julia:** resolvido; mesmo resultado com 1 e 4 processos.

# 15. A remoção de partículas depende da ordem no processo

**Onde:** `src/DMSwarm_move.cpp`, função `Swarm_add_remove_2d`.

É removida a partícula com a vizinha mais próxima. As duas partículas de um par mais próximo têm sempre a mesma distância, e o empate é decidido pela ordem no array, que muda com o número de processos. Na adição, `cell_seed` também inclui o `rank`. O resultado com 1 e com N processos não é o mesmo.

**Evidência:** o C++ foi só lido. A cópia fiel em Julia deu partículas diferentes com 1, 2 e 4 processos; com desempate pela posição e semente sem o `rank`, o resultado é idêntico, bit a bit.

**Sugestão** (não executada em C++): desempatar pela coordenada e tirar o `rank` de `cell_seed` (passo e índice do elemento já o tornam único). **Julia:** resolvido; `cpp_order = true` repete o C++.

# 16. A opção `-seed` não valida a camada

**Onde:** `src/reader.cpp:687` e `src/main.cpp`.

`weakening_seed[seed_layer[k]]` é escrito sem conferir se a camada existe; uma camada fora do intervalo escreve fora do array. O teste em `main.cpp` só compara a quantidade de valores com `n_interfaces + 1`, depois da escrita, e imprime `Error` sem interromper.

**Sugestão** (não executada): conferir `0 <= seed_layer[k] <= n_interfaces` antes da escrita e interromper. **Julia:** resolvido; interrompe com mensagem.

# 17. O padrão de `viscosity_per_element` é uma opção que não existe mais

**Onde:** `src/header.h:37` e `src/DMSwarm2mesh.cpp:274`.

O padrão é `visc_const_per_element = 0` (`variable`), mas `Swarm2Mesh_2d` imprime que `viscosity_per_element = variable` "is no longer available" e chama `exit(-2)`. Um `param.txt` sem a chave para depois de ler as entradas e montar as partículas, com uma mensagem que não diz qual valor usar. O `exit` dentro do MPI também não passa pelo `PetscFinalize`.

**Sugestão** (não executada): mudar o padrão para `constant` (`visc_const_per_element = 1`) ou recusar `variable` já no `reader.cpp`, dizendo que o valor aceito é `constant`. **Julia:** exige `constant` ao ler o `param.txt` e interrompe com essa mensagem.

# 18. A superfície inicial fica em degraus fora do primeiro processo em x

**Onde:** `src/sp.cpp:158`, função `sp_create_surface_swarm_2d`.

O teste `p < milocal - 1` compara o índice global `p` (de `si` a `si+milocal-1`) com o tamanho local `milocal`. No processo 0 funciona; nos outros processos em x, quase todos os pontos caem no `else` e recebem a profundidade do nó à esquerda, sem interpolar dentro do elemento. A superfície inicial fica em degraus.

**Evidência:** só lido.

**Sugestão** (não executada): usar `if (p < Nx - 1)`. **Julia:** resolvido; a superfície é montada com o índice global.

# 19. Pontos da superfície entre os intervalos z de dois processos são perdidos

**Onde:** `src/sp.cpp:148-165`.

Cada processo guarda os pontos com `z_bottom <= pz <= z_top`, sendo `z_top` o último nó local. Entre `z_top` de um processo e `z_bottom` do processo de cima há um elemento inteiro (`dz`) que nenhum processo cobre, e os pontos com a interface nessa faixa não são criados. Com mais de um processo em z, a superfície pode ter buracos, que a interpolação depois preenche com uma reta.

**Evidência:** só lido.

**Sugestão** (não executada): usar `pz < z_bottom + mjlocal * dz_const` (ou `pz <= z_top` só no último processo em z). **Julia:** resolvido; a superfície é um vetor replicado em todos os processos e nenhum ponto se perde.

# 20. Com mais de um processo em z, a reamostragem duplica pontos

**Onde:** `src/sp.cpp:413`, função `sp_surface_swarm_interpolation`.

Cada processo com pontos locais (`nlocal > 0`) reamostra a sua faixa `si..si+milocal` em x. Se dois processos da mesma coluna têm pontos (a superfície atravessa a divisa em z), os dois geram os mesmos `px`, e a superfície remontada tem pontos repetidos e um tamanho que muda com a partição.

**Evidência:** só lido.

**Sugestão** (não executada): reamostrar só no processo 0 a partir de `seq_array` e redistribuir, ou escolher um único dono por coluna. **Julia:** resolvido; a reamostragem é a mesma em todos os processos.

# 21. Partículas não são atualizadas nos processos sem pontos da superfície

**Onde:** `src/sp.cpp:1110`, função `sp_update_surface_swarm_particles_properties`.

O laço é `for (p = 0; seq_surface_size && p < nlocal; p++)`, com `seq_surface_size` igual ao tamanho **local** do swarm da superfície. Um processo sem pontos da superfície (por exemplo, abaixo dela em z, mas com partículas de ar ou rocha perto da divisa depois da sedimentação ou erosão) não troca ar por sedimento nem rocha por ar. O índice `i = px/sp_dx` também lê `seq_array[2*i+3]` sem conferir o fim do vetor.

**Evidência:** só lido.

**Sugestão** (não executada): usar o tamanho do vetor sequencial completo e limitar `i` a `global_surface_size - 2`. **Julia:** resolvido; todos os processos atualizam as suas partículas com a superfície completa. Com 1 e 4 processos, os resultados de Julia são iguais (superfície 4e-3 m, densidade 4e-6).

# 22. O dPhi do calor latente é só o último incremento das iterações

**Onde:** `src/DMSwarm_move.cpp:448-470` (`moveSwarm`), chamado em cada iteração não linear (`src/DMv.cpp:527`) e em cada um dos 4 subpassos.

A fusão de cada partícula é recalculada em toda chamada de `moveSwarm`, sempre com a mesma temperatura da malha: `dphi = (T - Ts)/(L + facX)`, `Phi += dphi`, e o solidus sobe com a depleção. As chamadas funcionam como iterações de Newton: `Phi` converge para o equilíbrio sem calor latente (`T = Ts`), e `dPhi` guarda só o último incremento, que tende a zero. O calor latente da equação térmica (`-T*dS*dPhi/dt`) usa esse `dPhi`, então depende do número de iterações não lineares e some quando elas convergem bem. A fusão não consome a energia correspondente.

**Evidência:** no rift de teste (121 × 41, `magmatism_Ts0 = 970`, passo 0, Julia, que repete o C++):

| `Xi_min` | Iterações | Soma de Phi nos nós | dPhi máximo |
|:-:|:-:|:-:|:-:|
| 1e-3 | 2 | 2,262 | 1,9e-4 |
| 1e-6 | 10 | 2,296 | 2,9e-6 |
| 1e-9 | 34 | 2,298 | 1,6e-7 |

**Sugestão** (não executada): calcular a fusão uma vez por passo, depois da equação do calor, com `dPhi` igual à fusão do passo, e aplicar o calor latente no passo seguinte (ou baixar a temperatura da partícula em `L*dphi`). **Julia:** repete o C++.

# 23. O calor latente ignora `magmatism_dS`

**Onde:** `src/DM1_2d.cpp:498`; em 3D, também `src/DMSwarm_move.cpp:689-692`.

Na equação do calor, `float dS = 400.0;` é uma variável local que esconde o parâmetro `dS` lido de `magmatism_dS`. Em 3D, `Ts0`, `dTsdP`, `dTsdX` e `dS` também são constantes locais, e os quatro parâmetros `magmatism_*` não têm efeito.

**Sugestão** (não executada): apagar as declarações locais. **Julia:** usa `magmatism_dS` na fusão e no calor latente. Com o padrão (400) o resultado é o mesmo do C++.

# 24. A temperatura da fusão está em °C, e a do calor latente em K

**Onde:** `src/DMSwarm_move.cpp:455` e `src/DM1_2d.cpp:491`.

A fusão usa `L = tp*dS/c_heat_capacity`, com `tp` em °C. O calor latente usa `T_mean + 273`. A energia retirada não é a que a fusão supôs: a razão é `(T + 273)/T`, cerca de 1,2 a 1.300 °C.

**Sugestão** (não executada): usar a temperatura absoluta nos dois (e 273,15). **Julia:** repete o C++.

# 25. Sem fusão na última chamada, o magma vai para x = 0

**Onde:** `src/DMSwarm_move.cpp:251-252` e `464-467`, e `src/magmatic.cpp`.

`x_magma_high` e `z_magma_high` voltam a `(0, -depth)` em cada chamada de `moveSwarm`. A extração, nos passos de saída, usa a posição da última chamada (o último subpasso). Se nenhuma partícula fundiu nessa chamada, o que é comum pelo item 22, o volume é extrudido em x = 0, na borda esquerda. Num empate entre processos, vale o x do processo de maior número.

**Evidência:** no rift de teste com `magmatism_Ts0 = 970`, os 2,26e8 m² do primeiro passo foram extrudidos em x = 0 nos dois códigos.

**Sugestão** (não executada): guardar a posição da partícula mais rasa com fusão desde a última extração (ou no passo inteiro), e não extrudir se não houver. **Julia:** repete o C++, mas desempata pelo menor x (profundidades a menos de 1 cm são empate), e o resultado não depende do número de processos.

# 26. `calc_magmatic_extraction` não retorna valor

**Onde:** `src/magmatic.cpp:33-89`.

A função é `PetscErrorCode`, mas não tem `return`. O compilador avisa (`-Wreturn-type`). Com o clang do macOS, o fim da função vira uma instrução de trap, e o programa para com `SIGTRAP` na primeira extração. Com o gcc (Linux, Hypatia), o comportamento é indefinido e por acaso funciona; o valor também não é conferido em `main.cpp`.

**Evidência:** no macOS, todo cenário com extração parou no passo 1 com `Caught signal number 5 TRAP`. Com `PetscFunctionReturn(0);` no fim, numa cópia de teste, rodou normalmente.

**Sugestão** (não executada): `PetscFunctionReturn(0);` no fim da função e `CHKERRQ(ierr)` na chamada. **Julia:** não se aplica.

# 27. Com `nondimensionalization = True`, o tempo mistura unidades

**Onde:** `src/main.cpp:420-434` (`Calc_dt_calor`), `src/main.cpp:262` (laço do tempo), `src/main.cpp:477` e `530` (`rescaleVeloc`, `multi_veloc_change`), `src/sp.cpp:1162`, `1175` e `1188`, e `src/reader.cpp:574-592` e `783-789`.

No modo adimensional, as velocidades e os comprimentos são adimensionais, então `dt_calor = 0.1*(dh/v)/seg_per_ano` é o passo adimensional dividido por `seg_per_ano`, e o fim de `Calc_dt_calor` o multiplica de novo. O tempo `tempo` acumula assim tempo adimensional (a escala `time0 = h0/veloc0` é cerca de 14 Myr no exemplo `continental_rift_non_dimensional`), mas é comparado com valores em anos ou Myr:

- `dt_max` (anos) é comparado com o passo antes da multiplicação, um número da ordem de 1e-10, e nunca limita o passo;
- `time_max` (anos) é comparado com o tempo adimensional: `time_max = 100e6` corresponde a mais de 1e15 Myr e na prática não para a simulação;
- as mudanças de `scale_bcv.txt`, `multi_veloc.txt`, `sedimentation_rate.txt`, `sediment_layers.txt` e `base_level.txt` comparam o tempo adimensional com `1e6*t` (t em Myr) e acontecem no momento errado.

Além disso, a linha `k` de `interfaces.txt` (difusividade de cada camada) não é dividida por `kappa0_scaled`, mas o termo de difusão é multiplicado por `diffusion_scaled = kappa0*time0/h0^2`; com a linha `k`, a difusividade das camadas fica multiplicada por `kappa0`. O código marca isso como pendente (`//!!! check in the case of variable kappa`). As saídas (`time_*.txt` e os campos) ficam adimensionais.

Por leitura, outras partes também misturam unidades no modo adimensional: os processos superficiais (`sp_d_c` e as taxas em m²/ano com comprimentos e passo adimensionais) e o magmatismo (`magmatism_dTdP`, em K/Pa, multiplicado pela pressão adimensional).

Duas tolerâncias dependem das unidades e fazem a solução mudar com a opção: `denok` (item 28) e o `g` a mais do FSSA (item 7: o peso efetivo é $\theta g/2$, 2,5 $\theta$ em SI e 0,5 $\theta$ com $g = 1$).

**Evidência:** o exemplo `continental_rift_non_dimensional/model1` em 161 × 31, 4 passos, 1 processo, `solver = direct`, `dt_max` grande (o C++ não o aplica no modo adimensional) e sem deformação inicial aleatória. Os campos adimensionais foram convertidos com as escalas de `dimensional_convertion.txt`:

| Comparação (erro relativo da velocidade) | Passo 0 | Passo 1 | Passo 4 |
|:-----------------------------------------------|:-:|:-:|:-:|
| C++ adimensional × C++ dimensional, mesmo `denok = 1e-5` | 1,0 | 1,0 | 1,0 |
| o mesmo, com `denok` equivalente (7,9e-14 em SI) | 1,9e-6 | 1,3e-1 | 9,9e-2 |

Com o `denok` equivalente, os dois modos coincidem no passo 0; a diferença que aparece no passo 1 é a do FSSA. O passo de tempo é o mesmo (84,6 mil anos) nos dois modos.

**Sugestão** (não executada): manter o tempo em segundos dimensionais no laço (`tempo*time0_scaled`) e converter só os termos das equações, ou adimensionalizar também `dt_max`, `time_max` e os tempos dos arquivos; dividir `conductivity[i]` por `kappa0_scaled`; converter `denok` e corrigir o FSSA (item 7). **Julia:** com `nondimensionalization = True`, as equações são resolvidas em SI, que é a mesma solução, e `denok` e o `g` do FSSA de `fssa_legacy` são convertidos. Assim `dt_max`, `time_max`, os arquivos com tempos e a linha `k` funcionam como no modo dimensional. As saídas ficam em SI, e as escalas são escritas em `dimensional_convertion.txt`. Com `fssa_legacy = true`, o Julia repete o C++ adimensional (velocidade 7e-5 e temperatura 2e-5 no passo 4, em SI).

# 28. `denok` é um limite absoluto e depende das unidades

**Onde:** `src/DMv.cpp:608-614`.

As iterações de Uzawa param quando `denok`, o quadrado da norma de $G^T v$ (o divergente integrado, em m²/s em SI), fica abaixo de `denok_min`. O limite é absoluto: em SI esse resíduo é da ordem de 1e-7 a 1e-11 nos modelos de litosfera, então o padrão `denok = 1e-4` (e o 1e-5 do exemplo `continental_rift_non_dimensional`) para o Uzawa já na iteração 0 e deixa a velocidade compressível. Com `nondimensionalization = True` o mesmo valor é 1/(veloc0·h0)² = 1,3e8 vezes mais exigente.

**Evidência:** no mesmo caso do item 27, passo 0:

| Rodada | `denok` | Iterações de Uzawa | $\lVert\nabla\cdot v\rVert/\lVert\nabla v\rVert$ |
|:------------------|:-:|:-:|:-:|
| C++ dimensional | 1e-5 | 9 | 0,12 |
| C++ dimensional | 7,9e-14 | 158 | 2,6e-3 |
| C++ adimensional | 1e-5 | 158 | 2,6e-3 |

Os cenários do tapioca usam `denok = 1e-14` e não são afetados.

**Sugestão** (não executada): usar um critério relativo, por exemplo `denok` dividido pelo resíduo inicial ou por $\lVert G\rVert\,\lVert v\rVert$, ou ao menos mudar o padrão para 1e-14. **Julia:** repete o C++, converte `denok` no modo adimensional e avisa quando `denok` > 1e-10.

# 29. O FSSA não tem efeito em 3D

**Onde:** `src/DM1_v.cpp:748`, `778` e `794`, função `AssembleA_Veloc_3d`.

O termo do FSSA é somado em `Ke_veloc_final` (linha 778), mas a linha 794 faz `Ke_veloc_final[i]=Ke_veloc[i]*volume` para todos os termos e o apaga antes de a matriz ser montada. Em 3D, `free_surface_stab` não muda nada. Além disso, `theta_FSSA = 0.5` é uma variável local (linha 748) que esconde o parâmetro, e o termo tem o mesmo `g` a mais do item 7.

**Evidência:** só lido; o modelo 3D do Julia sem o termo (`fssa_legacy = true`) coincide com o C++ (item "3D: comparação" abaixo).

**Sugestão** (não executada): mover a soma do termo para depois da multiplicação pelo volume, usar o parâmetro `theta_FSSA` e a fórmula do item 7 (`theta*dt*rho_mean_bottom*dx*dy/4`, com o fator das matrizes). **Julia:** o FSSA 3D usa a forma consistente ($\theta\,dt\,\Delta(dRho)\,dx\,dy/8$ por nó); com `fssa_legacy = true` não há termo, como no C++.

# 30. Em 3D, `front` e `back` são faces opostas para a temperatura e para a velocidade

**Onde:** `src/DMT.cpp:311-313` e `src/DMv.cpp:429-435`.

Na temperatura, `front_temperature` vale para y = 0 e `back_temperature` para y = Ly. Na velocidade, `back_*_velocity` vale para y = 0 e `front_*_velocity` para y = Ly. Um `param.txt` que fixa "a frente" aplica as condições em faces diferentes.

**Sugestão** (não executada): usar a mesma convenção nos dois e documentá-la. **Julia:** repete o C++ (para os mesmos `param.txt`) e documenta.

# 31. A pressão do topo é somada em precisão simples

**Onde:** `src/calc_drho.cpp:239`, `257` (2D) e `436`, `455` (3D), funções `shift_pressure_2d/3d`.

A média da pressão no topo, que é subtraída de toda a pressão nos nós, é acumulada em `float` e somada entre processos com `MPI_FLOAT`. A pressão fica com erro relativo da ordem de 1e-7, e o resultado muda com o número de processos.

**Evidência:** em todas as comparações deste documento, a pressão do Julia (que soma em `double`) difere da do C++ em 5e-7 a 6e-7, também quando velocidade e temperatura coincidem em 1e-6.

**Sugestão** (não executada): `PetscReal` e `MPIU_REAL`. **Julia:** soma em precisão dupla.

# Observações a confirmar

- **Nível de base médio:** `sp_evaluate_adjusted_mean_elevation_with_sea_level` (`src/sp.cpp:324`) faz a média só do primeiro quarto dos pontos (`seq_surface_size / 2 / 4`), isto é, da borda esquerda do modelo. Vale confirmar se é intencional. O Julia repete.
- **Altura da superfície nas partículas:** `sp_update_surface_swarm_particles_properties` usa `rx = 0.5` fixo, a média dos dois pontos vizinhos, e não a posição da partícula. O Julia repete.
- **`aggradation_rate`:** é lido, mas não é usado.
- **Volume extraído:** é o aumento da soma de `Phi` nos nós desde a última extração. `Phi` nunca é zerado, então o magma extrudido continua na fonte, e oscilações da média nos nós (advecção, partículas novas) entram no volume. Vale também para a fusão da crosta. O Julia repete.
- **Viscosidade com fusão:** o fator `X_OH*X_m` (Lu & Huismans, 2021, com as constantes fixas no código) é aplicado depois dos limites de viscosidade, então a viscosidade pode passar de `viscosity_max` (até 5 vezes) ou ficar abaixo de `viscosity_min`. A densidade não muda com a fusão nem com a depleção. O Julia repete.
- **Extração sem a superfície:** `magmatism_extraction = on` sem `sp_surface_tracking` não faz nada e não avisa. O Julia avisa.
- **Superfície acima do topo do modelo:** os pontos da superfície que saem da malha (por exemplo, numa extrusão grande) são descartados pela migração do swarm, e a reamostragem preenche o buraco. O Julia mantém esses pontos e avisa.
- **`basal_heat`:** a condutividade é a do manto (`kappa*RHOM*c`), não a da camada da base, e a temperatura fixa da base (`bot_temperature`) é substituída depois de cada solução. O Julia repete.
- **`multi_velocity`:** o novo `input_velocity_<n>.txt` substitui também os valores de contorno, então um fator de `scale_bcv.txt` aplicado antes se perde. Sem `velocity_from_ascii`, a velocidade atual vira o contorno. O Julia repete.
- **3D sem difusividade por camada:** em 3D, a térmica usa a difusividade global (`kappa`, `src/thermal_Ke.cpp:321` e `483`), e `Swarm2Mesh_3d` não calcula `geoq_kappa` (marcado no código como pendente), então a linha `k` de `interfaces.txt` e `high_kappa_in_asthenosphere` não valem. Também faltam em 3D o aquecimento por cisalhamento (marcado como pendente), `basal_heat`, a superfície, os processos superficiais, a extração de magma e a fronteira periódica. O Julia 3D tem a difusividade por camada, `high_kappa_in_asthenosphere`, o aquecimento por cisalhamento e `basal_heat`, como no 2D; superfície, magmatismo e fronteira periódica interrompem a rodada em 3D.
- **`sticky_blanket_air` em 3D:** o laço (`src/DMSwarm2mesh.cpp:1114`) vai até `sx+mmx-2` em cada direção, como no item 14, e escreve a viscosidade mínima também nos nós do último plano em x e y, que não guardam elemento (sem efeito na solução). O Julia usa os elementos do processo.
- **Ponto exatamente em `px` na reamostragem:** é ignorado pelas desigualdades estritas (`cx_aux < px`, `cx_aux > px`). O Julia usa esse ponto; a diferença só aparece nesse caso.
- **Jacobiano do elemento térmico:** `montaKeThermal_*` omite o fator 1/4. Como ele aparece em todos os termos, a solução não muda.
- **Reologia 7:** a viscosidade é `visco_ref*A*exp(-(QE+VE*P)/(R*T))`, que cresce com a temperatura se `QE` é positivo. Vale confirmar se é intencional. O Julia repete a fórmula.
- **Reologia 1:** retorna antes dos limites de viscosidade e da plasticidade. Vale confirmar. O Julia repete.
- **PETSc.jl 0.5.6 (não é do Mandyoc):** a estrutura `MatStencil` está na ordem `i,j,k,c`, e o header do C usa `k,j,i,c`. Por isso o Julia usa `MatSetValuesLocal`. Visto por leitura.

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

# Reologia e modelo acoplado: comparação Julia × C++

A função `calc_visco_ponto` foi comparada em 52.000 pontos aleatórios, com as 13 reologias e as opções de plasticidade, pressão e limites: diferença relativa máxima de 1,1e-14.

O modelo completo (partículas, térmica, Stokes não linear e advecção) foi comparado em dois casos dos testes, reamostrados, sem FSSA e com o `dt` do C++. Erro relativo (norma L2):

| Caso | Passo | Velocidade | log10 da viscosidade | Temperatura |
|:----------------|:-:|:-:|:-:|:-:|
| Rift viscoplástico (201 × 51) | 0 | 1,7e-8 | 4,5e-9 | 0 |
| Rift viscoplástico (201 × 51) | 3 | 2,8e-8 | 1,3e-8 | 8,1e-12 |
| Slab, lei de potência (101 × 61) | 0 | 9,8e-7 | 8,2e-8 | - |
| Slab, lei de potência (101 × 61) | 4 | 4,3e-7 | 6,1e-9 | - |

O rift deu o mesmo com 1 e 4 processos. A temperatura do slab é nula. Nestas linhas do slab, o C++ rodou sem `Swarm_add_remove_2d`.

Com a adição e remoção de partículas nos dois códigos (1 processo, mesma sequência de `rand_r`):

| Caso | Passo | Velocidade | Deformação | Densidade |
|:----------------|:-:|:-:|:-:|:-:|
| Slab (cerca de 100 trocas por passo) | 6 | 3,1e-7 | 6,9e-6 | 8,4e-9 |
| Rift, dt de 50 mil anos | 30 | 8,0e-7 | 1,4e-6 | 6,4e-9 |

No rift, 172 partículas entraram e 103 saíram no último passo.

# Processos superficiais: comparação Julia × C++

`continental_rift_with_surface_processes` do tapioca-svg, 161 × 31 nós, 40 pontos da superfície por elemento, 3 passos de 10 mil anos, 1 processo, `fssa_legacy = true`. A diferença da superfície está no limite da precisão dos `surface_*.txt` do C++ (6 algarismos, 0,1 m em 40 km):

| `sp_mode` | Variação da superfície (C++) | Diferença máxima da superfície | Densidade | Deformação |
|:------------------------------|:-:|:-:|:-:|:-:|
| `diffusion` | 0,2 m | 0,05 m | 3,5e-7 | 1,1e-6 |
| `sedimentation_only` | 1.510 m | 0,03 m | 3,5e-7 | 4,7e-7 |
| `sedimentation_rate_limited` | 518 m | 0,05 m | 3,5e-7 | 1,1e-6 |
| `theunissen_sedimentation` (nível de base variável) | 29 m | 0,12 m | 3,5e-7 | 1,1e-6 |

Densidade e deformação: erro relativo (norma L2) no último passo.

# Magmatismo: comparação Julia × C++

`rifting_sedimentation` do tapioca-svg, 121 × 41 nós, 5 passos, 1 processo, `fssa_legacy = true`, com o solidus baixado para fundir já no início. O C++ com extração rodou numa cópia com o `return` do item 26. Erro relativo (norma L2) no último passo:

| Caso | Phi | X_depletion | Temperatura | Superfície |
|:----------------------------------------|:-:|:-:|:-:|:-:|
| `magmatism_Ts0 = 930`, sem extração | 3,2e-6 | 1,0e-6 | 2,1e-6 | 0,55 m |
| `magmatism_Ts0 = 970`, com extração (2,26e8 m² em x = 0) | 1,2e-5 | 7,7e-7 | 2,1e-6 | 3,5 m |

O volume extraído e a posição batem com o C++ (2,2632e8 m²; com `Ts0 = 930`, 2,6976e9 m² em x = 599,17 km). O `dPhi` difere em cerca de 3e-4, porque é a diferença pequena `T - Ts` em precisão simples. A velocidade e a deformação diferem em 1e-3 também sem magmatismo, porque uma fileira de partículas cai exatamente sobre a interface da crosta em -62,5 km. Com 1 e 4 processos, os resultados de Julia são iguais (Phi 7e-7, superfície 4e-3 m).

# 3D: comparação Julia × C++

`continental_rift_non_dimensional/model2` em 41 × 5 × 11 nós, 15 partículas por elemento, unidades SI, `denok = 1e-14`, `solver = direct`, sem perturbação das partículas e sem deformação aleatória, `fssa_legacy = true`. Erro relativo (norma L2):

| Rodada | Passo | Velocidade | Temperatura | Pressão | Deformação |
|:-----------------------------------------|:-:|:-:|:-:|:-:|:-:|
| 3 passos de 10 mil anos, 1 processo | 3 | 1,2e-6 | 2,2e-6 | 5,7e-7 | 1,2e-6 |
| o mesmo, C++ com 2 processos | 3 | 3,8e-6 | 2,4e-6 | 5,7e-7 | 1,2e-6 |
| 30 passos de 85 mil anos (adição e remoção de partículas) | 30 | 1,0e-3 | 2,2e-6 | | 2,3e-2 |

A viscosidade difere em 5e-4 só nos nós do último plano em x e y (observação sobre `sticky_blanket_air`). Na rodada longa, as diferenças aparecem com a adição e remoção de partículas, cujos sorteios são diferentes nos dois códigos (o C++ usa uma sequência por processo). O Julia com 1 e 4 processos dá o mesmo resultado (temperatura 7e-6, pressão 7e-8, velocidade idêntica fora do ruído de `vy`).

