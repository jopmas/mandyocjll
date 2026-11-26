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

# Observações a confirmar

- **Nível de base médio:** `sp_evaluate_adjusted_mean_elevation_with_sea_level` (`src/sp.cpp:324`) faz a média só do primeiro quarto dos pontos (`seq_surface_size / 2 / 4`), isto é, da borda esquerda do modelo. Vale confirmar se é intencional. O Julia repete.
- **Altura da superfície nas partículas:** `sp_update_surface_swarm_particles_properties` usa `rx = 0.5` fixo, a média dos dois pontos vizinhos, e não a posição da partícula. O Julia repete.
- **`aggradation_rate`:** é lido, mas não é usado.
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

