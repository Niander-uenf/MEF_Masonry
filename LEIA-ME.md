# MEF_Masonry

*English: see [README.md](README.md).*

Elementos finitos térmicos e termomecânicos para alvenaria de blocos em
temperatura, sem ANSYS. Só precisa de Python com NumPy, SciPy e Matplotlib; os
resultados 3D saem em VTK para visualizar no ParaView (gratuito), e o programa
tem interface no navegador com pós-processador 3D próprio.

O pacote nasceu para o ensaio da UFV com o bloco BCEPQ (prisma de 5 fiadas de
bloco e meio, juntas secas, forno a 5 °C/min, 30 min de patamar e carga ao fim
do patamar). Hoje trata de alvenarias de vários blocos e, na versão 1.0, cobre:

| | O que entra |
|---|---|
| Blocos | Maciços; furos redondos verticais (BCEPQ); células retangulares verticais (concreto NBR 6136, cerâmico estrutural NBR 15270, inclusive de paredes vazadas), com ou sem graute; furos horizontais (vedação) |
| Juntas | Seca (contato, folgas e atrito) ou argamassa (camada aderida com dano próprio) |
| Corpos de prova | Bloco isolado, prisma, mini parede e parede em escala real |
| Análises | Térmica 3D; termoelástica 3D; não linear 3D até a ruptura; seção 2D (horizontal, ou vertical nos furos horizontais); parede com segunda ordem; fluência transitória no não linear 3D e na parede; sensibilidade global da parede (Morris) |
| Ensaios | Ambiente, a quente, residual, transiente (carga constante e aquecimento até a ruptura) e incêndio-padrão ISO 834 |
| Biblioteca | 11 blocos, 3 materiais, 4 argamassas e 2 grautes, cada um com a fonte |
| Verificação | 8 roteiros, 240 casos na rodada rápida (o 01 e o 06 rodam com os três elementos), num só comando: `exemplos/00_verificar_tudo.py` |

As seções abaixo descrevem cada etapa, da mais recente (versão 1.0, com o
novo nome) para a mais antiga (modo Parede, 0.5).

## Versão 1.0: o programa passa a se chamar MEF_Masonry

O programa nasceu para o BCEPQ, mas hoje calcula blocos maciços, cerâmicos
estruturais e de vedação, blocos de concreto com e sem graute, com junta seca
ou de argamassa. Por isso, na versão 1.0, ele passa a se chamar MEF_Masonry.
O BCEPQ continua como um dos blocos da biblioteca e como projeto padrão.

O que mudou de nome:

| Antes (até 0.10) | Agora (1.0) |
|---|---|
| Pasta `mef_bcepq` | Pasta `MEF_Masonry` |
| Pacote Python `mef_bcepq` | Pacote `mef_masonry` |
| `Iniciar_MEF_BCEPQ.bat` | `Iniciar_MEF_Masonry.bat` |
| Executável `dist/MEF_BCEPQ/MEF_BCEPQ.exe` | `dist/MEF_Masonry/MEF_Masonry.exe` |

O que continua funcionando:

- **Projetos salvos** (`.json`) das versões anteriores abrem normalmente; os
  campos novos entram com os valores padrão.
- **O último projeto aberto no navegador** volta sozinho na primeira vez que
  a 1.0 abre.
- **Roteiros próprios** que fazem `import mef_bcepq` continuam rodando: o
  nome antigo é um apelido de `mef_masonry`, que avisa a troca e será retirado
  numa versão futura. Troque `mef_bcepq` por `mef_masonry` nos seus roteiros.
- **Histórico de análises:** fica na pasta `resultados_app` da instalação
  antiga. Para vê-lo na página Análises anteriores, copie essa pasta para
  dentro de `MEF_Masonry`.

Os cálculos não mudaram: as 240 verificações dão os mesmos resultados da
versão 0.10.

## Instalação

```
python -m pip install -r requirements.txt
```

Python 3.9 ou mais recente. Opcional: `pip install pypardiso` (Windows/Linux)
acelera o solver mecânico 3D e reduz a memória.

## Programa com interface

```
python app.py          (ou dois cliques em Iniciar_MEF_Masonry.bat no Windows)
```

No Windows, o `Iniciar_MEF_Masonry.bat` procura o Python (`py` ou `python`),
instala o NumPy, o SciPy e o Matplotlib na primeira vez e abre o programa;
se faltar algo, a janela fica aberta com a explicação.

Abre no navegador uma interface local (nada vai para a internet) com:

- **Geometria paramétrica:** dimensões do bloco, diâmetro e número de furos
  por bloco, blocos por fiada (múltiplos de 0,5), número de fiadas,
  amarração (contrafiada ou a prumo), junta e capeamento. A malha pode ser
  pré-visualizada.
- **Propriedades** da cerâmica, com tabelas opcionais em função da
  temperatura (k, fator de E, deformação térmica, ft e fc).
- **Modelos e elementos:**
  - análise térmica 3D, termoelástica linear 3D, não linear 3D até a
    ruptura (dano de Mazars com fechamento de fissuras; juntas com folgas e
    atrito), seção 2D ou parede (faixa vertical com segunda ordem);
  - elementos 3D W6 (cunha linear, padrão do artigo), H8 (hexaedro linear
    sobre os quadriláteros da malha em "O") e W15 (cunha quadrática); T3 na
    seção 2D. Veja "Elementos" abaixo.
- **Ensaios:** compressão a 30 °C, a quente (rampa, patamar e carga),
  residual (aquece, resfria e carrega) e incêndio-padrão ISO 834 numa face,
  sob carga constante. No ISO 834 não linear, a ruptura segue os limites de
  deformação axial da EN 1363-1 (C = h/100 e dC/dt = 3h/1000 por minuto);
  esse caso é experimental, porque a carga entra como pressão no topo, sem
  placa rígida.
- **Execução** com registro ao vivo e cancelamento.
- **Resultados:** indicadores, gráficos, tabelas, CSV e VTK.
- **Pós-processador 3D** (WebGL, funciona sem internet):
  - mapas de cor de temperatura, deslocamentos, tensões principais e
    vertical, dano de tração, de compressão e efetivo, e folgas das juntas;
  - estados do aquecimento e da carga;
  - forma deformada com escala;
  - cortes em x, y e z;
  - espelhamento da simetria;
  - materiais ocultáveis, valor no ponto (duplo clique) e imagem PNG.
- **Projetos** salvos e abertos em .json, além de projetos prontos (os
  ensaios da UFV).

As análises ficam em `resultados_app/`, uma pasta por execução, com
`projeto.json`, `resultados.json`, figuras, CSV, VTK e `visual.json`.

**Executável (.exe):** no computador de destino, rode `python empacotar.py`
(ou dois cliques em `Empacotar_Windows.bat`). O PyInstaller gera
`dist/MEF_Masonry/MEF_Masonry.exe`, que pode ser copiado para outras máquinas
sem Python. O empacotamento não foi testado no Windows a partir deste
ambiente; se faltar algum módulo, acrescente um `--hidden-import` em
`empacotar.py`.

**Pelo código:** `mef_masonry.projeto.executar(projeto, pasta)` roda o mesmo
projeto sem interface (veja `projeto_padrao()` e `ESQUEMA`).

## Sensibilidade global da parede (versão 0.10)

Antes de calibrar o modelo e de planejar ensaios, vale saber quais entradas
mais mudam a resposta. A versão 0.10 traz uma análise de sensibilidade global
da parede no ISO 834 pelo **método de Morris** (efeitos elementares): o
programa roda a parede muitas vezes, variando as entradas incertas dentro de
faixas, e ordena quais mais pesam nos tempos de resistência (R) e de
isolamento (I) e na capacidade a frio.

**Como usar:** em Modelos e elementos, escolha a análise "Sensibilidade da
parede no ISO 834 (Morris)". Aparece a página 6, Sensibilidade, com os
parâmetros e as faixas. O resto do projeto (geometria, materiais, carga,
apoios) são os valores nominais em torno dos quais tudo varia. O projeto
pronto "Parede BCEPQ 2,94 × 3,00 m — sensibilidade no ISO 834" já vem
montado.

**O método:**

- Cada parâmetro tem uma faixa, mapeada em [0, 1] e dividida numa grade de
  p níveis (padrão 4).
- Uma trajetória parte de um ponto sorteado da grade e muda um parâmetro de
  cada vez, em ordem aleatória, de Δ = p / [2 (p − 1)] (2/3 da faixa com
  4 níveis). Com k parâmetros, cada trajetória custa k + 1 rodadas.
- O efeito elementar do parâmetro i é EE = [y(x + Δ eᵢ) − y(x)] / Δ, na
  escala da faixa inteira. Com r trajetórias (padrão 10), cada parâmetro tem
  r efeitos.
- μ* = média de |EE| (Campolongo et al., 2007) é quanto a saída muda, em média,
  quando o parâmetro percorre a faixa: é a medida que ordena os parâmetros.
  μ, a média com sinal, diz para que lado a saída anda. σ, o desvio-padrão,
  mede a não linearidade e as interações.
- O intervalo de 95% de μ* sai de um bootstrap nas trajetórias.
- A razão σ/μ* classifica o efeito (Garcia Sanchez et al., 2014): quase
  linear (< 0,1), monotônico (< 0,5), quase monotônico (< 1) ou não linear e
  com interações (≥ 1).

| Parâmetro | Faixa padrão | Por quê |
|---|---|---|
| Condutividade térmica do bloco | ±20% | Incerteza da medida e variação entre lotes; vale também para a tabela k(T), se houver |
| Calor específico do bloco | ±15% | Inclui a incerteza do pico de umidade; vale também para a tabela c(T) |
| Emissividade das paredes dos furos | 0,7 a 0,95 | Cerâmica: cerca de 0,85 a 0,95; o limite baixo cobre superfícies mais lisas |
| Convecção na face exposta | 15 a 35 W/m²·K | A EN 1991-1-2 usa 25 W/m²·K na face exposta ao incêndio-padrão |
| Emissividade da superfície | 0,6 a 0,95 | A EN 1991-1-2 recomenda 0,8 para a superfície do elemento |
| Dilatação térmica do bloco | 4 a 8 ×10⁻⁶/°C | Cerâmica vermelha: cerca de 4 a 8 ×10⁻⁶/°C; substitui o α e escala a tabela de deformação térmica |
| Fluência transitória do bloco (k_tr) | 0 a 2,35 | De sem fluência (0) ao limite superior do concreto (2,35; Anderberg e Thelandersson, 1976) |
| Resistência da alvenaria (fator na lei) | ±15% | Dispersão da resistência dos prismas; multiplica as tensões da lei da parede |
| Rigidez da alvenaria (fator na lei) | ±25% | Dispersão do módulo dos prismas; multiplica o módulo da lei da parede |
| Imperfeição inicial (H dividida por) | H/1000 a H/250 | Flecha inicial H/n no meio da altura; a amplitude varia linearmente entre os dois extremos |
| Excentricidade da carga (topo e base) | -10 a 10 mm | Aplicada igual no topo e na base; positiva para o lado exposto ao fogo |

**Custo:** r (k + 1) + 1 rodadas, com a rodada nominal. Com os 11 parâmetros e
r = 10, são 121 rodadas. Na malha rápida da triagem (seção com 10 mm e 16
elementos de viga), cada rodada leva de 10 a 30 s num núcleo (a fluência
transitória é a mais cara); a malha recomendada é cerca de 6 vezes mais
lenta. As rodadas correm em paralelo (núcleos menos um,
por padrão).

**Saída:** no resumo, os valores nominais e o parâmetro de maior efeito em
cada tempo; uma tabela de efeitos por saída (μ*, intervalo, μ, σ, σ/μ* e a
leitura); as figuras de μ* e do plano μ*–σ; e os arquivos
`sensibilidade_R.csv`, `sensibilidade_I.csv`, `sensibilidade_Pmax.csv` e
`sensibilidade_rodadas.csv` (o plano completo, com os valores físicos e as
saídas de cada rodada).

**Parede BCEPQ** (projeto pronto: 2,94 × 3,00 m, 150 kN, ISO 834 por 300 min,
11 parâmetros, r = 10, 121 rodadas na malha rápida, 11,5 min em 2 núcleos;
as saídas estão em `resultados/sensibilidade_parede_bcepq/`).
Nominais: R = 211 min, I = 112 min e capacidade a frio de 554 kN.

| Tempo de resistência R | μ* (min) | IC 95% | μ (min) | σ/μ* |
|---|---|---|---|---|
| Dilatação térmica α (4 a 8 ×10⁻⁶/°C) | 253 | 181 a 328 | −253 | 0,50 |
| Excentricidade da carga (−10 a +10 mm) | 157 | 74 a 238 | +157 | 0,91 |
| Imperfeição inicial (H/1000 a H/250) | 55 | 21 a 102 | −55 | 1,28 |
| k_tr do bloco (0 a 2,35) | 24 | 4 a 51 | +24 | 1,72 |
| Rigidez da alvenaria (±25%) | 20 | 4 a 40 | +20 | 1,59 |
| Demais (k, ε da face, c, ε dos furos, resistência, h) | 0,3 a 10 | | | |

| Tempo de isolamento I | μ* (min) | IC 95% | μ (min) | σ/μ* |
|---|---|---|---|---|
| Condutividade k (±20%) | 42 | 39 a 47 | −42 | 0,15 |
| Calor específico c (±15%) | 27 | 26 a 29 | +27 | 0,12 |
| Emissividade das faces (0,6 a 0,95) | 8,0 | 6,0 a 10,3 | +8,0 | 0,46 |
| Emissividade dos furos (0,7 a 0,95) | 7,6 | 6,8 a 8,4 | −7,6 | 0,18 |
| Convecção na face exposta (15 a 35 W/m²·K) | 0,8 | 0,6 a 1,1 | −0,8 | 0,52 |

Na capacidade a frio, pesam a resistência da alvenaria e a excentricidade
(148 kN cada), a imperfeição (87 kN) e a rigidez (25 kN); as entradas
térmicas dão efeito zero, como devem, e as mecânicas dão zero no isolamento.

O que isso diz para o programa experimental:

- **α do bloco é a primeira medida a fazer** (dilatometria até 1000 °C). O
  α = 6 ×10⁻⁶/°C do BCEPQ é provisório e cai bem num limiar: nenhuma rodada
  com α ≥ 6,7 passou de 300 min (medianas de 90 e 37 min), e cerca de 70% das
  rodadas com α ≤ 5,3 passaram. O mecanismo é o arqueamento térmico para o
  lado do fogo somado à segunda ordem.
- **A excentricidade e a imperfeição do ensaio precisam ser medidas.** Carga
  deslocada para a face exposta compensa o arqueamento (mediana de R acima de
  300 min com +10 mm, contra 64 min com −10 mm). As 12 rodadas que romperam
  antes de 30 min tinham α ≥ 6,7 e carga deslocada para a face fria.
- **k_tr e a rigidez vêm em segundo plano**, com efeito não linear e com
  interações (σ/μ* > 1): o ensaio transiente vale a pena, mas depois de α.
- **A resistência da alvenaria quase não muda R** (2 min), embora mande na
  capacidade a frio: a parede rompe por estabilidade, não por esmagamento.
- **No isolamento, k(T) e c(T) mandam**; as emissividades pesam pouco e a
  convecção, nada. Na calibração com os termopares, o foco é k e c (com o
  pico de umidade); as condições de contorno pesam menos.

Em 59 das 121 rodadas a parede não rompeu até 300 min, e o R entrou como 300
min (limite inferior): os efeitos no R, sobretudo o de α, são maiores do que
os da tabela. Os números valem para as faixas acima, para os apoios
articulados do projeto e para a lei do prisma BCEPQ.

**Verificação:** `exemplos/15_verificacao_sensibilidade.py`, 31 casos (32 com
`--completo`, mais 4 resultados informativos), todos dentro da tolerância:

- plano de Morris (um parâmetro por passo, passo Δ, pontos na grade);
- função linear (μ = aᵢ, σ = 0, exatos), função aditiva quadrática (μ*
  igual à conta feita à parte pelos pares de níveis; com r = 2000, a 0,4%
  do valor esperado), interação x₁x₂ (σ > 0) e função g de Sobol (ordem
  conhecida dos parâmetros);
- rodada que falha fica de fora só nos dois efeitos vizinhos; a mesma
  semente repete o plano; o intervalo de 95% contém μ*;
- faixas físicas, tabelas k(T) e ε_th(T) escaladas junto, projeto original
  intocado, fatores da lei (resistência de pico e módulo exatos);
- parede: a rodada do estudo dá os mesmos R, I e Pmax da análise Parede; em
  paralelo e em sequência, resultados idênticos (`--completo`).

Com a malha recomendada (5 mm e 32 elementos), o R nominal do BCEPQ cai de
210,9 para 205,1 min e o I passa de 112,0 para 112,5 min; cada rodada leva
cerca de 6,6 vezes mais. Nas 3 primeiras trajetórias do estudo (36 rodadas)
rodadas também com a malha recomendada, o R mudou 0,9 min em média (6,6 min
no máximo) e o I, no máximo 0,5 min; os cinco primeiros parâmetros saíram na
mesma ordem nas duas malhas, com μ* a menos de 5%.

## Fluência transitória e ensaio transiente (versão 0.9)

Num ensaio a quente (aquecer e depois carregar), o material já chega quente e
sem carga. Num incêndio, a parede está carregada enquanto esquenta, e o
material encurta além da dilatação livre e da deformação elástica: é a
**fluência transitória**, ou deformação térmica induzida pela carga. Ela é
irreversível e só aparece no primeiro aquecimento sob compressão.

**O modelo** (página Modelos e elementos, Fluência transitória):

- Anderberg e Thelandersson (1976): ε_tr = φ(T) σ / fc20, com
  φ(T) = k_tr ε_th(T). σ é a tensão de compressão e fc20 a resistência à
  compressão do material a 20 °C. No lugar de k_tr, o bloco aceita φ(T)
  tabelado, tirado do ensaio transiente.
- Incremental: Δε_tr = Δφ σ / fc20 só quando φ passa do maior valor já
  atingido. Não volta no resfriamento nem soma ao reaquecer abaixo da
  temperatura máxima.
- No 3D, Δε_tr = Δφ / fc20 · [(1 + ν) σ⁻ − ν tr(σ⁻) I], com σ⁻ a parte
  compressiva da tensão: no uniaxial, a fluência encurta na direção da carga
  e expande na transversal com o próprio Poisson do material. Cada intervalo
  entre campos térmicos é percorrido em subpassos com preditor e corretor
  (regra do trapézio), e o incremento de ε_tr por subpasso fica abaixo de
  2 × 10⁻⁴.
- Na parede, bloco e junta de leito estão em série na altura da fiada:
  c(T) = (1 − f) φ_b / fc_b + f φ_a / fc_a, com f a fração da junta na altura
  da fiada. O graute, contínuo, tem o seu próprio coeficiente.

**Valores padrão** (a fluência vem ligada):

| Material | k_tr | Por quê |
|---|---|---|
| Argamassa de cimento | 2,0 | Prakash, Azenha, Pereira e Lourenço (2020), na faixa 1,8 a 2,5 do concreto |
| Bloco cerâmico (BCEPQ e outros) | 0 | A cerâmica queimada não tem pasta de cimento e não há dado publicado; Prakash et al. (2020) também desprezam a fluência na unidade cerâmica |
| Bloco de concreto e graute com as curvas da EN 1992-1-2 | 0 | As curvas da norma já incluem a fluência (Gernay e Franssen, 2012); somar outra seria contar duas vezes |

Com esses padrões, o BCEPQ de junta seca não muda: os resultados das
versões anteriores continuam valendo. A fluência aparece nas alvenarias com
argamassa e onde você informar k_tr ou φ(T).

**Ensaio transiente** (Ensaio = transiente, só no não linear 3D): a carga
constante é aplicada a frio e o forno aquece pela rampa e pelo patamar até a
ruptura. A ruptura é a perda de equilíbrio ou os limites de deformação axial
da EN 1363-1. A saída é a curva deformação × temperatura, com a dilatação
livre e a parcela de fluência, e as temperaturas na ruptura. É o ensaio que
mede φ(T): com a deformação total ε, a livre ε_th e a elástica σ/E(T),

    φ(T) = (ε − ε_th − σ/E(T)) fc20 / σ

(encurtamentos positivos). Os resultados de dois projetos prontos novos:

| Projeto pronto | Sem fluência | Com fluência |
|---|---|---|
| Tijolo maciço com argamassa 1:1:6, 46 kN (40%) até 800 °C | 69% do tijolo com dano > 0,5 no fim | 17% (a fluência da argamassa alivia as tensões térmicas entre tijolo e junta) |
| BCEPQ, 49 kN (50%) até 800 °C, k_tr hipotético de 2 no bloco | 54% do bloco com dano > 0,5 no fim | 44%, com 4,0‰ de fluência média; sem ruptura nos dois (o BCEPQ não perde resistência a quente) |

**Verificação:** `exemplos/14_verificacao_fluencia.py`, 19 casos (20 com
`--completo`), todos dentro da tolerância:

- deformação sob tensão constante e aquecimento uniforme, axial e lateral,
  com os três elementos (exata);
- aquecimento restringido contra a solução analítica da relaxação,
  σ = fc/k [exp(−E k ε_th / fc) − 1] (0,6% com o passo padrão e 0,05% com
  10⁻⁵);
- irreversibilidade no resfriamento e no reaquecimento, ausência de fluência
  sob tração e φ(T) tabelado (exatos);
- parede: encurtamento sob carga constante (exato) e coeficiente em série de
  bloco e junta;
- parede BCEPQ: com k_tr = 0, ligar a fluência não muda nada.

**Sensibilidade da parede BCEPQ** (projeto pronto do ISO 834 sob 150 kN, 32
elementos, 240 min), com k_tr hipotéticos no bloco:

| k_tr do bloco | 0 | 0,5 | 1,0 | 1,5 | 2,0 |
|---|---|---|---|---|---|
| Tempo de resistência (min) | 205,1 | 179,4 | 167,6 | 176,2 | 200,9 |

O efeito não é monótono. A fluência encurta as fibras comprimidas que
esquentam: tira rigidez da seção e, ao mesmo tempo, reduz a flecha térmica
para o lado do fogo. Com k_tr perto de 1, o primeiro efeito domina e o tempo
cai 18%; com k_tr = 2, os dois quase se compensam. A flecha na ruptura fica
entre 43 e 46 mm em todos os casos, e o modo de ruptura é o mesmo. Como a
cerâmica não tem k_tr publicado, o valor do BCEPQ precisa sair de um ensaio
transiente (prisma carregado e aquecido), e o programa já reproduz esse
ensaio para comparar.

## Furos horizontais e consolidação (versão 0.8)

**Blocos de furos horizontais** (vedação), com Furação = Furos horizontais:

- **Geometria:** furos na espessura e na altura do bloco (o de 8 furos é
  2 × 4), espessura das paredes externas, dos septos (entre os furos, nas duas
  direções) e das paredes de leito (em cima e embaixo). Os furos atravessam o
  bloco ao longo do comprimento; as juntas verticais de argamassa os fecham.
- **3D:** a planta segue as bordas das colunas de furos, e cada fiada se divide
  em subcamadas nas bordas das fileiras de furos; os furos abrem nas fileiras,
  fora das juntas verticais. A radiação nos furos usa a largura do furo na
  direção do fluxo de calor.
- **Seção 2D vertical:** para estes blocos, a seção 2D é a vertical de uma
  fiada (espessura × altura, com meia junta em cima e embaixo). Como os furos
  atravessam o bloco, essa seção vale ao longo de todo o comprimento da parede.
  As bordas de cima e de baixo são cortes periódicos (sem troca de calor), e o
  estado plano generalizado fica ao longo dos furos.
- **Parede:** a térmica sai da mesma seção vertical, e as fibras são as paredes
  que descem contínuas pela fiada (faces e septos verticais, com as juntas
  sobre elas). Os septos horizontais entram só na térmica, porque não carregam
  compressão vertical.
- **Biblioteca:** blocos de 8 furos 90 × 190 × 190 mm (SINAPI 7271) e de 6
  furos 90 × 140 × 240 mm (paredes medidas por Guedes et al., 2019), com
  argamassa 1:2:9. E, fc e ft são estimativas a partir da resistência mínima de
  1,5 MPa da NBR 15270 (veja a nota de cada item).
- **Projetos prontos:** isolamento no ISO 834 da seção vertical do bloco de 8
  furos (24 min sem revestimento, nas propriedades da biblioteca) e prisma não
  linear a 20 °C.
- **Verificação:** `exemplos/13_verificacao_furos_horizontais.py`, 12 casos (14
  com `--completo`), todos dentro da tolerância:
  - volume dos furos no prisma 3D, área da seção horizontal mínima e áreas da
    seção vertical (exatos);
  - térmica: a seção vertical 2D reproduz o 3D extrudado ao longo dos furos
    (diferença de 0,3 °C aos 60 min de ISO 834);
  - fibras da parede: área por metro igual à largura menos os furos, e centro
    no plano médio;
  - prisma não linear de 8 furos: o pico é o das paredes contínuas
    (fc × área delas, dentro de 5%).

**Consolidação:**

- `exemplos/00_verificar_tudo.py` roda todas as verificações (01, 06 e 10 a 13,
  com os três elementos onde se aplica) e grava
  `resultados/resumo_verificacoes.csv`. Sem opção, leva cerca de 3 min; com
  `--completo`, inclui os casos longos (convergência da parede e prismas de
  Camacho et al.) e leva cerca de 30 min.
- A seção 2D, as figuras, o VTK e o pós-processador mostram a seção vertical em
  pé (espessura na horizontal, altura na vertical).

Os resultados do BCEPQ não mudaram: o prisma a 30 °C continua com 98,68 kN e
1212 MPa.

## Células retangulares e graute (versão 0.7)

Além dos furos redondos, o bloco pode ter **células retangulares**, como os
blocos de concreto (NBR 6136) e os cerâmicos estruturais (NBR 15270):

- **Geometria:** número de células por bloco, espessura das paredes externas
  (faces), dos septos e das paredes de topo. Com "células pequenas por
  parede externa" > 0, cada face ganha uma fileira de células pequenas, como
  nos blocos cerâmicos de paredes vazadas.
- **Malha:** a grade estruturada segue todas as bordas das células de todas as
  fiadas. Com contrafiada, as células mudam de lugar de uma fiada para a
  outra; cada elemento é classificado fiada a fiada (bloco, furo ou graute).
  As juntas de leito ficam com argamassa onde há bloco dos dois lados e com
  graute onde a coluna grauteada passa.
- **Graute:** em todos os furos ou células, ou alternados (um sim, um não),
  em toda a altura. Tem dano próprio, como a argamassa, e propriedades
  próprias. A quente, segue por padrão o concreto silicoso da EN 1992-1-2:
  fc(θ) e εc1(θ) da Tabela 3.1 (o fator de E é fc(θ)/εc1(θ)), ft(θ) e a
  deformação térmica. Vale também para os furos redondos.
- **Onde entra:** no 3D (térmico, linear e não linear), na seção 2D (o graute
  conduz calor no lugar do ar do furo) e no modo Parede (as fibras de graute
  somam à seção, com a sua dilatação térmica).
- **Biblioteca:**
  - blocos de concreto 14 × 19 × 39 (Mohamad et al., 2006), 14 × 19 × 29
    (Camacho et al., 2015), M-20 19 × 19 × 39 e vedação classe C 9 × 19 × 39
    (NBR 6136);
  - cerâmicos estruturais 14 × 19 × 29 de paredes maciças e de paredes
    vazadas e 14 × 19 × 39 de 3 células (NBR 15270);
  - grautes G1 (17,7 MPa) e G2 (26,7 MPa) de Camacho et al. (2015);
  - o material "concreto" agora traz as tabelas da EN 1992-1-2 (E, fc, ft e
    deformação térmica com a temperatura).

  Nos itens sem ensaio, E, fc e ft são estimativas, e a nota de cada item diz
  de onde vieram.
- **Verificação:** `exemplos/12_verificacao_celulas_graute.py`, 15 casos
  (16 com `--completo`), todos dentro da tolerância:
  - área do bloco fiada a fiada (com contrafiada e com paredes vazadas, exata);
  - volume de graute no prisma (exato);
  - rigidez em paralelo bloco + graute (exata);
  - propriedades do graute nos solvers, com as tabelas da EN 1992-1-2;
  - graute na térmica 2D e nas fibras da parede;
  - consistência da biblioteca.

  Com `--completo`, roda também o prisma de Camacho et al. (bloco B1, 3
  fiadas, cerca de 10 minutos) e compara com os ensaios, sem calibrar nada:

  | Prisma (área bruta) | Modelo | Ensaio | Diferença |
  |---|---|---|---|
  | Oco | 4,27 MPa | 5,63 MPa | −24% |
  | Grauteado com G1 | 13,86 MPa | 11,29 MPa | +23% |

  No prisma oco, a ruptura do modelo vem do esmagamento da argamassa (6,1
  MPa). O modelo de dano não tem dilatância, e por isso confina pouco a
  argamassa. No grauteado, E, ft e ν do graute são estimativas, e o modelo não
  inclui a retração do graute nem a sua incompatibilidade com o bloco. As
  duas coisas explicam boa parte das diferenças. Para um estudo, calibre E
  da junta e do graute com um prisma do seu lote.

  Os resultados do BCEPQ não mudaram: o prisma a 30 °C continua com
  98,68 kN e 1212 MPa, e as verificações 01, 06, 10 e 11 passam.

## Juntas de argamassa e biblioteca (versão 0.6)

A partir da versão 0.6, a junta pode ser **seca** (como no BCEPQ: camada fina
com contato, folgas e atrito) ou de **argamassa**.

**Junta de argamassa:**

- É uma camada contínua de espessura real, com dano próprio (o mesmo modelo
  de Mazars do bloco, com os parâmetros da argamassa) e perfeitamente aderida
  ao bloco.
- O assentamento pode ser em toda a face ou só nas paredes externas (furos e
  septos sem argamassa). As juntas verticais podem ser preenchidas ou vazias.
- Comprimento e altura do bloco passam a ser medidas de eixo a eixo (bloco +
  junta).
- Quando a argamassa expande lateralmente mais que o bloco (ν/E maior), ela
  traciona o bloco junto à junta (mecanismo de Hilsdorf), e o prisma rompe
  entre a resistência da argamassa e a do bloco.

**Outras novidades:**

- **Calor específico e massa específica em função da temperatura** (tabelas),
  pelo método da capacidade aparente. Servem, por exemplo, para o pico de
  umidade do concreto e do solo-cimento.
- **Biblioteca** (`mef_masonry/biblioteca.py`), com a fonte de cada item:
  - blocos: BCEPQ e tijolo maciço 100 × 60 × 220 mm (os de células
    retangulares entraram na versão 0.7; os de furos horizontais aparecem
    desabilitados até a fase 3);
  - materiais do bloco: BCEPQ, cerâmica e concreto pela NBR 15220-2, e
    concreto com k(T) e c(T) da EN 1992-1-2;
  - argamassas: os traços 1:0,25:3, 1:0,5:4,5, 1:1:6 e 1:2:9 (Mohamad,
    Lourenço e Roman, 2006).

  Na interface, escolher o item preenche os campos e mostra a fonte e o que
  falta conferir. Ao trocar um valor à mão, o item passa a "Personalizado".
- **Verificação:** `exemplos/11_verificacao_materiais_argamassa.py`, 12 casos,
  todos dentro da tolerância:
  - c(T) contra a solução exata de um corpo aquecido (erro de 0,06%);
  - coluna tijolo + argamassa em série (exata);
  - sinal das tensões laterais de Hilsdorf;
  - ruptura pela argamassa sem confinamento (0,4%) e ganho com confinamento;
  - propriedades da argamassa nos solvers;
  - consistência da biblioteca.

  Os resultados do BCEPQ não mudaram: o prisma a 30 °C continua com
  98,68 kN, e as verificações 01, 06 e 10 passam.

## Modo Parede (versão 0.5)

Para paredes em escala real, o 3D sólido fica pesado demais (uma parede de
1,26 × 2,58 m tem cerca de 300 mil nós) e não tem os efeitos que governam o
ensaio de parede em incêndio: arqueamento térmico, segunda ordem e
excentricidade. O modo Parede (`mef_masonry/parede.py`) resolve isso com um
modelo de faixa vertical de 1 m de comprimento:

- **Elementos de viga 2D corrotacionais** (Crisfield 1991; Battini 2002),
  com grandes deslocamentos e rotações: segunda ordem (P-delta), flambagem e
  arqueamento térmico. O padrão é de 32 elementos ao longo da altura.
- **Seção de fibras:** os triângulos de cerâmica da malha 2D de um bloco de
  comprimento, com os furos. Cada fibra tem a sua temperatura, tirada da
  térmica 2D com ISO 834 numa face (o mesmo cálculo da seção 2D).
- **Material homogeneizado:** curva tensão (área líquida) × deformação do
  prisma, calculada pelo modelo não linear 3D calibrado a 30, 600, 800 e
  1.000 °C (`mef_masonry/dados_parede.py`) e interpolada na temperatura. A
  alvenaria seca não resiste à tração, porque as juntas abrem. A descarga é
  secante, e o ramo pós-pico é regularizado pelo comprimento do elemento.
  Opcionalmente, a curva de 30 °C sai do modelo do prisma com as entradas
  atuais.
- **Deformação térmica livre** da cerâmica (α ou tabela).
- **Entradas:** carga total, apoios articulados ou engastados no topo e na
  base, excentricidades (positivas para o lado exposto) e imperfeição
  inicial senoidal (H/450 por padrão).
- **Saídas:**
  - capacidade a frio e fator de redução por esbeltez;
  - tempo de resistência (perda de equilíbrio) e critério de contração
    axial da EN 1363-1;
  - tempo de isolamento (140/180 °C);
  - flecha e deslocamento do topo ao longo do tempo;
  - deformadas;
  - perfis de temperatura;
  - curva tempo de resistência × nível de carga;
  - no pós-processador, a seção com a temperatura e a tensão vertical nas
    fibras na meia altura.
- **Tempo:** a capacidade leva cerca de 5 s; o ISO 834 de 4 h com 5 níveis
  de carga, cerca de 5 min.

Verificação (`exemplos/10_verificacao_parede.py`, 14 casos, todos dentro da
tolerância):

- amplificação da imperfeição de Euler (0,1–0,3%);
- arqueamento térmico livre (0,15%);
- fórmula da secante (0,03%);
- área e inércia das fibras (0,04%);
- capacidade da parede curta igual à resistência do prisma a 20 e a 600 °C
  (0,04%);
- dilatação livre (0,15%);
- convergência com 16, 32 e 64 elementos: 32 elementos ficam a 0,8% do tempo
  de resistência com 64.

Limitações do modo Parede:

- **Bordas verticais:** são livres, sem restrição à flexão.
- **Lascamento:** não está incluído. A fluência transitória entra desde a
  versão 0.9, mas só com k_tr > 0 ou φ(T) informado (a cerâmica queimada não
  tem valor publicado).
- **Estanqueidade (E):** não é verificada.
- **Temperatura ao longo da altura:** é a mesma em toda a parede.
- **α:** é provisório, e a flecha térmica e o tempo de resistência dependem
  diretamente dele.
- **Apoios engastados:** a dilatação restringida faz a face exposta esmagar
  junto ao apoio, e o modelo tende a perder o equilíbrio cedo (resultado
  conservador).

## Ordem de uso

| Passo | Script | O que faz | Tempo e memória* |
| --- | --- | --- | --- |
| 1 | `exemplos/01_verificacao.py` | Faz 36 comparações com soluções analíticas e de referência; todas devem dar OK | 11 s |
| 2 | `exemplos/02_calibrar_juntas.py` | Ajusta a rigidez das juntas secas ao módulo medido do prisma e compara com a fórmula em série e a NBR 16868 | 45 s, 1 GB |
| 3 | `exemplos/03_ensaio_ufv_3d.py` | Simula os ensaios de 600, 800 e 1000 °C no prisma 3D | ~1 min por temperatura, 1,6 GB |
| 4 | `exemplos/04_parametrico_secao.py` | Diâmetro dos furos × condutividade: temperatura na carga, tração térmica e isolamento ISO 834 | 1 min |
| 5 | `exemplos/05_otimizacao_secao.py` | Modelo de otimização (evolução diferencial) e fronteira de Pareto massa × isolamento | 2 min |
| 6 | `exemplos/06_verificacao_nao_linear.py` | Verifica o modelo não linear em 11 casos com solução conhecida | 2 min |
| 7 | `exemplos/08_calibrar_imperfeicao.py` | Calibra a folga das juntas (pelo E do prisma) e a fc efetiva (pela carga a 30 °C) | 15 min |
| 8 | `exemplos/07_ensaio_ufv_nao_linear.py` | Leva o prisma à ruptura a 30, 600, 800 e 1.000 °C (modelo não linear) | ~5 min a 30 °C; 10–15 min aquecido |
| 9 | `exemplos/09_estudo_convergencia.py` | Estudos de malha, incremento, tolerâncias e Δt (Tabela S2 do artigo) | 5–60 min por caso |
| 10 | `exemplos/10_verificacao_parede.py` | Verifica o modo Parede em 14 casos (`--rapido` pula a convergência) | 20 s; 4 min completo |
| 11 | `exemplos/11_verificacao_materiais_argamassa.py` | Verifica c(T), a junta de argamassa e a biblioteca em 12 casos | 1 min |
| 12 | `exemplos/12_verificacao_celulas_graute.py` | Verifica células retangulares e graute em 15 casos (`--completo` roda também os prismas de Camacho et al.) | 5 s; 10 min completo |
| 13 | `exemplos/13_verificacao_furos_horizontais.py` | Verifica os furos horizontais em 12 casos (`--completo` roda também o prisma de 8 furos) | 10 s; 1 min completo |
| 14 | `exemplos/14_verificacao_fluencia.py` | Verifica a fluência transitória em 19 casos (`--completo` roda também a sensibilidade da parede BCEPQ e o ensaio transiente do tijolo maciço) | 2 min; 10 min completo |
| — | `exemplos/00_verificar_tudo.py` | Roda todas as verificações e grava o resumo | 3 min; 30 min completo |

\* Medidos num computador com 2 núcleos. Malha padrão do 3D linear: h = 8 mm,
dz = 12 mm (com h = 6 mm e dz = 10 mm, pede cerca de 5 GB de memória). O não
linear usa h = 12 mm e dz = 20 mm, cerca de 3 vezes mais rápido; a calibração
vale para essa malha.

Os resultados vão para `resultados/` (CSV, PNG e `vtk/`). Rode os scripts a
partir da pasta do pacote, por exemplo `python exemplos/03_ensaio_ufv_3d.py --temperaturas 600`.
Cada script aceita `--help`.

## O modelo

**Geometria.** A seção de bloco e meio (420 × 140 mm, três furos de 84 mm)
tem malha estruturada fora dos furos e malha em "O" em volta de cada furo, que
segue a circunferência. O prisma 3D é a extrusão da seção: placa de aço,
capeamento, cinco fiadas de 60 mm com juntas secas de 2 mm, capeamento e placa.
As juntas verticais entre bloco e meio-bloco alternam entre x = 280 e 140 mm.
O plano de simetria y = 70 mm reduz o modelo à metade.

**Materiais** (`mef_masonry/parametros.py`, com as fontes):

| Propriedade | Valor | Fonte |
| --- | --- | --- |
| k da cerâmica | 0,52 W/m·K | Cerâmica 2026, n = 20 |
| ρ, c | 1.210 kg/m³; 1.285 J/kg·K | BPQ; c derivado da efusividade (898) e de k |
| E do bloco (área bruta) | 2.801 MPa | BPQ, n = 6 |
| E do material (área líquida) | 3.907 MPa | 2.801 / 0,717 |
| ν | 0,18 | BPQ |
| α | 6×10⁻⁶ /°C | **provisório**, sem medição |
| E da junta seca (2 mm) | 84,4 MPa | calibrado (passo 2) para E_prisma = 1.100 MPa |

Propriedades dependentes da temperatura entram por tabelas (`tabela_k`,
`tabela_E_fator`, `tabela_eps_termica`), por exemplo a curva de dilatometria
com o salto do quartzo perto de 573 °C.

**Térmica.** Condução transiente com Euler implícito. As faces expostas trocam
calor por convecção (h = 10 W/m²·K) e radiação (ε = 0,8) com o ar do forno.
O ar dos furos conduz por uma condutividade radiativa equivalente,
k_rad = 4σT³Dε/(2−ε).

**Mecânica.** Termoelasticidade linear com deformação térmica. Os pratos da
prensa são apoios rígidos sem atrito: base com u_z = 0 e, na carga, topo com
deslocamento uniforme ajustado para a reação igualar a carga do ensaio. O
caso térmico (topo livre) e o de carga se somam. Saídas: tensões principais,
dilatação e encurtamento entre placas, deformação no terço central da face
externa (como um extensômetro ou DIC) e tração máxima na cerâmica.

**Elementos.** Todos isoparamétricos, com jacobiano completo
(`mef_masonry/elementos3d.py`):

| Elemento | Nós | Integração | Na mesma malha | Uso |
| --- | --- | --- | --- | --- |
| W6, cunha linear | 6 | 3 × 2 pontos | referência | padrão do artigo |
| H8, hexaedro linear | 8 | 2 × 2 × 2 | mesmos nós, metade dos elementos | conferir a influência do elemento |
| W15, cunha quadrática | 15 | 6 × 3 pontos | cerca de 4,7 vezes mais nós | flexão e gradientes fortes, com malha grossa |

O dano e o estado das juntas valem para o elemento inteiro e são avaliados
no centro dele; a rigidez é integrada com a regra completa. Na W15, a
largura da banda de fissuração usa a cunha linear dos cantos, e a térmica
usa matrizes consistentes de capacidade e de troca de calor (concentrar nos
nós das faces quadráticas desvia o fluxo de entrada). As cargas de pressão
usam o vetor consistente. A W15 pesa: num prisma de 12 × 20 mm, 71.674 nós
contra 15.080 da W6, e o modelo linear a 600 °C leva cerca de 7 min e 4,6 GB
de memória sem o pypardiso. Para a não linear, use a malha mais grossa
(16 × 30 mm, 35.869 nós) e, se puder, instale o pypardiso; o programa avisa
quando a malha da W15 passa de 40 mil nós na não linear.

Escolha do elemento, ensaio não linear a 30 °C (Tabela S2c do artigo; 2 núcleos):

| Elemento | Malha (mm) | Nós | Pico (kN) | E de corda (MPa) | Tempo |
| --- | --- | --- | --- | --- | --- |
| W6 | 12 × 20 | 15.080 | 98,68 | 1.204 | 2,6 min |
| H8 | 12 × 20 | 15.080 | 97,74 | 1.212 | 2,9 min |
| W6 | 16 × 30 | 7.638 | 97,38 | 1.333 | 0,8 min |
| W15 | 16 × 30 | 35.869 | 96,34 | 1.326 | 26 min (1,7 GB) |

O elemento muda o pico em até 1,1%; a malha pesa mais no módulo de corda.

**Seção 2D.** Mesma térmica, com tensões em estado plano generalizado. Serve
para estudos rápidos: forno com todas as faces expostas, ou parede com ISO 834
numa face e ambiente na outra. O isolamento usa o critério de 140/180 °C na
face não exposta.

## Modelo não linear

Para chegar à carga de ruptura, a cerâmica e as juntas deixam de ser
elásticas (`mef_masonry/naolinear.py`).

**Cerâmica: dano com fechamento de fissuras** (Mazars, 1984; Mazars e
Pijaudier-Cabot, 1989). σ = (1 − d) D (ε − ε_t), com
1 − d = (1 − r·d_t)(1 − d_c):

- d_t (tração) começa pelo critério de Rankine, quando a tensão principal
  efetiva chega a ft;
- d_c (compressão) começa pelo critério de von Mises das tensões de
  compressão, quando chega a fc;
- r é o peso da tração no estado atual: 1 em tração pura, 0 em compressão
  pura. Uma fissura térmica comprimida pela prensa fecha e volta a transmitir
  compressão.

Os dois danos amolecem de forma exponencial. A energia dissipada é
regularizada pela banda de fissuração (Bažant e Oh, 1983), com a largura do
elemento medida na direção normal à fissura (Oliver, 1989), para o resultado
não depender da malha.

**Juntas secas e capeamento: contato com folgas e atrito.** São camadas finas
(Desai et al., 1984) que não transmitem tração e transmitem cisalhamento até o
atrito de Coulomb sem coesão, μ|σ_n|. As faces de assentamento dos blocos não
são planas. Sem carga, só uma parte da área nominal se toca: 23% em blocos sem
tratamento (Ngapeya e Waldmann, 2020) e cerca de 15% subindo a 95% perto da
ruptura (Zahra et al., 2016). As folgas concentram a carga e governam a
fissuração e a capacidade do prisma seco (Ngapeya, Waldmann e Scholzen, 2018).
No modelo, cada junta de leito tem um campo aleatório de folgas iniciais
(gaussiano, correlação de 40 mm, 20% da área tocando sem carga). Onde há
contato, a junta é tão rígida quanto a cerâmica, de modo que a flexibilidade
da junta vem só das folgas que fecham. A folga entra como deformação inicial
da camada fechada: σ_n = E (ε_n + g/t).

**Carregamento:** o prisma é aquecido campo a campo, com os campos térmicos
reais da rampa e do patamar, sob a carga de assentamento de 2 kN. Depois, o
prato de topo desce em deslocamento controlado até depois do pico.

| Parâmetro | Valor | Referência |
| --- | --- | --- |
| Resistência à tração ft | 0,66 MPa | Gooch et al. (2023): tração direta de tijolo cerâmico prensado, área líquida |
| Resistência à compressão fc (literatura) | 5,76 MPa | BPQ: fb médio de 4,13 MPa em área bruta, dividido por 0,717 |
| Energia de fratura em tração G_f | 0,018 mm × ft ≈ 12 N/m | índice de ductilidade de unidades cerâmicas (van der Pluijm, 1992; Lourenço, 1996; via Dauda et al., 2020). É coerente com os 0,05–0,085 N/mm medidos para ft ≈ 3 MPa (Almeida et al., 2002) |
| Energia de fratura em compressão G_c | 0,33 mm × fc | índice de ductilidade em compressão de unidades cerâmicas (Angelillo et al., 2014; via Dauda et al., 2020) |
| Atrito em junta seca μ | 0,60, coesão nula | 0,61–0,64 em juntas secas serradas (Lourenço et al., 2004) |
| Área em contato sem carga | 20% | 23% (Ngapeya e Waldmann, 2020); ~15% (Zahra et al., 2016) |
| Desvio-padrão das folgas | **0,028 mm** (calibrado) | E de corda do prisma (5–30% da carga máxima) = 1.100 MPa na média de 3 sementes (BPQ: 1.001–1.216 MPa). O script 08 dá 0,029 mm; com 0,028 mm, o E médio é cerca de 1.130 MPa |
| fc efetiva da cerâmica | **2,48 MPa** (calibrada) | carga máxima a 30 °C (98,7 kN); só com `--calibrado` |

**Calibração** (`08_calibrar_imperfeicao.py`), em duas etapas independentes:
(1) a folga sai só do trecho inicial da curva, pelo módulo do prisma, na
média de três sementes do campo aleatório (a rigidez varia cerca de 20% entre
sementes); (2) com essa folga, a fc efetiva sai da carga máxima a 30 °C.
O script deu 0,029 mm e 2,48 MPa (valores padrão: 0,028 mm e 2,48 MPa).

Com as resistências da literatura, o modelo dá 230 kN a 30 °C, mais que o
dobro do ensaio. Para reproduzir o ensaio, a fc efetiva precisa ser 0,43 vez
a do bloco. Explicações prováveis:

- o ensaio da unidade baixa e larga (60 × 140 mm), capeada e confinada pelo
  atrito dos pratos, superestima a resistência do material;
- as asperezas menores que a malha não estão no modelo;
- a carga a 30 °C ainda precisa ser confirmada.

Os ensaios aquecidos são previsões, sem novo ajuste.

**Resultados** (`resultados/nl_resumo_calibrado.csv` e `nl_ensaios_calibrado.png`;
folga de 0,028 mm, fc efetiva de 2,48 MPa; malha h = 12 mm, dz = 20 mm):

| Ensaio | Carga máxima do modelo | Ensaio (UFV) | Diferença | Cerâmica com d > 0,5 após o aquecimento | E de corda 5–30% |
| --- | --- | --- | --- | --- | --- |
| 30 °C | 98,6 kN | 98,7 kN | calibrado | 0,06% | 1.210 MPa |
| 600 °C | 97,1 kN | 97,9 kN | −1% | 13% | 882 MPa |
| 800 °C | 98,1 kN | 98,9 kN | −1% | 17% | 964 MPa |
| 1.000 °C | 99,0 kN | 114,4 kN | −13% | 18% | 1.072 MPa |

Com a fc da literatura (5,76 MPa), o mesmo modelo dá 230 kN a 30 °C.

No aquecimento, 13–18% da cerâmica fissura por tração. Ainda assim, a carga
máxima quase não muda até 800 °C, como nos ensaios, porque as fissuras
térmicas fecham sob a compressão da prensa. Sem o fechamento (dano
isotrópico, versão anterior, com fc de 2,58 MPa), o prisma a 600 °C rompia
com 80 kN.

A 1.000 °C o ensaio ganha resistência, e o modelo não acompanha porque as
propriedades são constantes com a temperatura. Uma explicação provável é a
continuação da sinterização perto da temperatura de queima.

A rigidez cai depois do aquecimento, o que é uma previsão a conferir com o DIC
ou os extensômetros dos ensaios aquecidos.

Com outras sementes do campo de folgas (2 e 3), a carga máxima a 30 °C foi de
99,5 e 98,8 kN (variação < 1%), e o E de corda, de 964 e 1.162 MPa. A
rigidez depende mais da distribuição das folgas que a resistência; a
dispersão é parecida com a dos prismas do BPQ.

Tempo: 3–4 min a 30 °C e 11–14 min nos ensaios aquecidos (2 núcleos).
Incrementos aceitos sem convergência completa (passo mínimo): 0, 3, 1 e 8,
respectivamente; o script informa essa contagem.

**Dependência da temperatura.** As tabelas `tabela_ft_fator`, `tabela_fc_fator`
(DanoCeramica), `tabela_E_fator` e `tabela_eps_termica` (Ceramica) recebem os
valores. O padrão é constante, e α = 6×10⁻⁶ /°C é provisório. Para tijolos
cerâmicos vazados, a literatura indica que o módulo sobe até cerca de 1,75 vez
e cai acima de 800 °C, e que a dilatação tem pico perto de 600 °C (Nguyen e
Meftah, 2014, conforme Prieler et al., 2025). Os números devem vir do original
ou do Anexo D da EN 1996-1-2 (acesso pago); não foram inventados aqui.

**Verificação** (`exemplos/06_verificacao_nao_linear.py`, 11 casos):

- picos de tração e de compressão, com erro < 0,4%;
- energia de fratura com dois tamanhos de elemento: erro de 5–6% e diferença
  entre as malhas < 2%;
- fissura de tração comprimida recupera a rigidez intacta;
- junta seca: sem tração, rigidez em compressão exata, cisalhamento
  saturando em μN, e folga que só transmite força depois de fechar.

**Limitações:**

- dano escalar: a fissura não tem orientação própria, e o fechamento depende
  só do sinal das tensões;
- atrito tratado como rigidez secante, sem histórico de deslizamento;
- dano local com banda de fissuração (sensibilidade moderada à malha);
- folgas aleatórias: o E de corda varia cerca de 20% entre sementes; a
  carga máxima, menos de 1%;
- fluência transitória pelo modelo de Anderberg e Thelandersson, com k_tr
  (ou φ(T)) por material e Poisson transitório igual ao elástico;
- fc efetiva calibrada, e não medida.

### Referências do modelo não linear

- Al-Fakih, A.; Mohammed, B.S.; Liew, M.S. (2018). Behavior of the dry bed joint in the mortarless interlocking masonry system: an overview. *Civil Engineering Research Journal*, 4(3), 555639.
- Almeida, J.C.; Lourenço, P.B.; Barros, J.A. (2002). Characterization of brick and brick–mortar interface under uniaxial tension. Seminário, Belo Horizonte.
- Bažant, Z.P.; Oh, B.H. (1983). Crack band theory for fracture of concrete. *Materials and Structures*, 16, 155–177.
- Dauda, J.A.; Iuorio, O.; Lourenço, P.B. (2020). Numerical analysis and experimental characterisation of brick masonry. *Int. J. Masonry Research and Innovation*, 5(3), 321–347.
- Desai, C.S.; Zaman, M.M.; Lightner, J.G.; Siriwardane, H.J. (1984). Thin-layer element for interfaces and joints. *Int. J. Numer. Anal. Methods Geomech.*, 8, 19–43.
- Gooch, L.J.; Masia, M.J.; Stewart, M.G.; Lam, C.Y. (2023). Statistical assessment of tensile and shear properties of unreinforced clay brick masonry. *Construction and Building Materials*, 386, 131578.
- Jirásek, M.; Bažant, Z.P. (2002). *Inelastic Analysis of Structures*. Wiley.
- Lourenço, P.B.; Ramos, L.F.; Vasconcelos, G. (2004). On the cyclic behaviour of stone dry masonry joints. 13th Int. Brick and Block Masonry Conference, Amsterdã.
- Mazars, J. (1984). *Application de la mécanique de l'endommagement au comportement non linéaire et à la rupture du béton de structure*. Tese, Université Paris VI.
- Mazars, J.; Pijaudier-Cabot, G. (1989). Continuum damage theory — application to concrete. *J. Eng. Mech.*, 115(2), 345–365.
- Ngapeya, G.G.C.; Waldmann, D.; Scholzen, F. (2018). Impact of the height imperfections of masonry blocks on the load bearing capacity of dry-stack masonry walls. *Construction and Building Materials*, 165, 898–913.
- Ngapeya, G.G.C.; Waldmann, D. (2020). Overcome of bed-joint imperfections and improvement of actual contact in dry-stacked masonry. *Construction and Building Materials*, 233, 117173.
- Nguyen, T.-D.; Meftah, F. (2014). Behavior of hollow clay brick masonry walls during fire. Part 2: 3D finite element modeling and spalling assessment. *Fire Safety Journal*, 66, 35–45.
- Oliver, J. (1989). A consistent characteristic length for smeared cracking models. *Int. J. Numer. Methods Eng.*, 28, 461–474.
- Prieler, R. et al. (2025). Thermo-mechanical analysis of masonry brick walls with embedded test specimen under fire exposure. *Fire Technology*, 61, 875–926.
- Zahra, T. et al. (2016), conforme Al-Fakih et al. (2018): área de contato de juntas secas de ~15% a ~95% com a carga.

## Verificação

Para conferir tudo de uma vez (por exemplo, depois de instalar ou atualizar):

```
python exemplos/00_verificar_tudo.py            (cerca de 7 min)
python exemplos/00_verificar_tudo.py --completo (cerca de 45 min)
```

| Roteiro | O que verifica | Casos |
|---|---|---|
| 01 (W6, H8, W15) | Térmica e termoelástica linear contra soluções exatas | 36 por elemento |
| 06 (W6, H8, W15) | Dano e juntas com atrito | 11 por elemento |
| 10 | Modo Parede | 10 (14 completo) |
| 11 | c(T), junta de argamassa e biblioteca | 12 |
| 12 | Células retangulares e graute | 15 (16 completo, mais 2 comparações com ensaios) |
| 13 | Furos horizontais | 12 (14 completo) |
| 14 | Fluência transitória | 19 (20 completo, mais 3 resultados informativos) |
| 15 | Sensibilidade (Morris) e rodadas da parede | 31 (32 completo, mais 4 resultados informativos) |

`01_verificacao.py` grava `resultados/verificacao.csv`; com
`--elemento H8` ou `--elemento W15`, roda os mesmos casos com esse elemento
(`verificacao_H8.csv`, `verificacao_W15.csv`). Os três passam nos 36 casos.
`06_verificacao_nao_linear.py` aceita a mesma opção, e os três passam nos 11
casos do modelo não linear. Casos de `01`:

- condução transiente numa placa (série de Fourier), em 2D e 3D, com dois passos de tempo;
- compressão uniaxial, com prato rígido e com pressão;
- dilatação livre uniforme;
- coluna em camadas, cerâmica + juntas (associação em série);
- placa com perfil parabólico de temperatura, em 3D e 2D;
- seção com furos, comparando a térmica 2D com o 3D extrudado.

Erros: menores que 0,5% nas tensões e de 1 °C na condução com dt = 15 s.
No último caso, H8 e W15 têm tolerância de 1 °C (e não 0,5 °C), porque a
seção 3D não é a mesma malha T3 do modelo 2D.

## Dados experimentais

Coloque em `dados/` arquivos `termopares_600.csv` (e 800, 1000) com as colunas
`tempo_min,T_C`. O passo 3 sobrepõe essas medições às curvas do modelo.

## Limitações

- Modelo linear (passos 1 a 5): material elástico, sem fissuração nem dano, e
  juntas secas como camadas finas equivalentes. Para a carga de ruptura, use o
  modelo não linear (passos 6 a 8).
- Furos tratados como cavidades fechadas, sem convecção vertical.
- α, E(T) e resistência à tração ainda precisam ser medidos.
- A rigidez calibrada das juntas depende da malha; recalibre se mudar h ou dz.
- Junta de argamassa perfeitamente aderida ao bloco: sem descolamento nem
  deslizamento na interface. Isso basta para a compressão de prismas, mas não
  para cisalhamento ou flexão fora do plano governados pela aderência. O α da
  argamassa é provisório.
- Na seção 2D e no modo Parede, a seção é só a do bloco (as juntas verticais de
  argamassa não entram na térmica 2D). Com células retangulares e contrafiada,
  a seção 2D é a da primeira fiada.
- Graute perfeitamente aderido ao bloco, sem retração. Nos blocos de concreto e
  cerâmicos estruturais da biblioteca, E, fc e ft do bloco e do graute são
  estimativas (veja a nota de cada item). As paredes das células são
  prismáticas (sem a conicidade dos blocos reais); a biblioteca usa espessuras
  equivalentes que dão a área líquida medida.
- O dano em compressão não tem dilatância: prismas cuja ruptura é governada
  pelo esmagamento da argamassa confinada tendem a ficar abaixo do ensaio.
- Fluência transitória: o Poisson transitório é o elástico (não há dado para
  outro valor); k_tr da cerâmica queimada é desconhecido (padrão 0); a
  fluência básica (de longa duração) não entra.
- Sensibilidade (Morris): é uma triagem, que ordena as entradas mas não
  reparte a variância da saída (para isso, índices de Sobol, com muito mais
  rodadas). Os efeitos dependem das faixas escolhidas; o tempo de resistência
  sem ruína até o fim da exposição entra como a duração (limite inferior).
  Só a parede no ISO 834 tem esse estudo.
- Furos horizontais: furos abertos nas pontas do prisma, sem convecção ao longo
  deles. Na parede, a seção vertical é a de uma fiada (a mesma em toda a
  altura), e as fibras somam as paredes contínuas; a lei do material continua
  sendo a do prisma (calibre pelo prisma do seu bloco). Nos blocos de vedação
  da biblioteca, E, fc e ft são estimativas.

## Para o artigo experimental × numérico × teórico

| Comparação | Experimental | Numérico | Teórico |
| --- | --- | --- | --- |
| Campo térmico | termopares | passo 3 | condução em placa (verificação) |
| Rigidez a 30 °C | E do prisma (BPQ, DIC) | passo 2 | série bloco + junta; NBR 16868 (600 fpk) |
| Tensões térmicas | fissuração observada, ultrassom | passos 3 e 4 | placa com perfil parabólico |
| Geometria | — | passos 4 e 5 | — |
| Carga de ruptura × temperatura | carga máxima dos ensaios | passo 8 | — |
