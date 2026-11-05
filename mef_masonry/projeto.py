"""Projeto do programa: entradas em JSON -> modelo -> análise -> resultados.

Um projeto é um dicionário com cinco grupos (geometria, ceramica, modelo,
ensaio, malha). As unidades das entradas são as de laboratório (mm, MPa,
°C, min); a conversão para SI fica aqui. ESQUEMA descreve cada campo
(rótulo, unidade, ajuda, opções) e serve tanto à interface quanto à
validação.

    from mef_masonry.projeto import projeto_padrao, executar
    p = projeto_padrao()
    p["ensaio"]["tipo"] = "quente"; p["ensaio"]["T_alvo"] = 600
    res = executar(p, "resultados/app/teste")
"""
import copy
import json
import math
import os
import time
from dataclasses import replace
from typing import Callable, Dict, List

import numpy as np

from . import malha as ml, termico as tm, mecanico as mc, curvas, estudos
from .malha import AR, CERAMICA, JUNTA, CAPEAMENTO, ACO, ARGAMASSA, GRAUTE
from .parametros import Geometria, DanoCeramica, JuntaAtrito, CALIBRACAO_30C
from .termico import Materiais, Contorno
from .saida import salvar_vtk, salvar_csv
from . import visual
from . import biblioteca as _bib

VERSAO = "1.0"

# ----------------------------------------------------------------------------
# Esquema das entradas
# ----------------------------------------------------------------------------
def _n(id, rotulo, padrao, unidade="", ajuda="", minimo=None, maximo=None, passo=None, **kw):
    d = {"id": id, "rotulo": rotulo, "tipo": "numero", "padrao": padrao, "unidade": unidade, "ajuda": ajuda}
    if minimo is not None: d["min"] = minimo
    if maximo is not None: d["max"] = maximo
    if passo is not None: d["passo"] = passo
    d.update(kw)
    return d


def _o(id, rotulo, padrao, opcoes, ajuda="", **kw):
    return {"id": id, "rotulo": rotulo, "tipo": "opcao", "padrao": padrao, "opcoes": opcoes, "ajuda": ajuda, **kw}


def _b(id, rotulo, padrao, ajuda="", **kw):
    return {"id": id, "rotulo": rotulo, "tipo": "logico", "padrao": padrao, "ajuda": ajuda, **kw}


def _t(id, rotulo, colunas, ajuda="", **kw):
    return {"id": id, "rotulo": rotulo, "tipo": "tabela", "padrao": [], "colunas": colunas, "ajuda": ajuda, **kw}


NL = {"modelo.analise": ["nao_linear"]}
LIN = {"modelo.analise": ["linear"]}
AQ = {"ensaio.tipo": ["quente", "residual", "transiente"]}
NLP = {"modelo.analise": ["nao_linear", "parede", "sensibilidade"]}
FLU = {**NLP, "modelo.fluencia": [True]}
PAR = {"modelo.analise": ["parede", "sensibilidade"]}
PAREDES = ("parede", "sensibilidade")
P3D = {"modelo.analise": ["termica", "linear", "nao_linear"]}
ARG = {"geometria.tipo_junta": ["argamassa"]}
SECA = {"geometria.tipo_junta": ["seca"]}
CIRC = {"geometria.forma_furo": ["circular"]}
RET = {"geometria.forma_furo": ["retangular"]}
HOR = {"geometria.forma_furo": ["horizontal"]}
VERT = {"geometria.forma_furo": ["circular", "retangular"]}
RET_HOR = {"geometria.forma_furo": ["retangular", "horizontal"]}
GRA = {"geometria.grauteamento": ["todos", "alternados"]}

_AJUDA_SENS = {
    "k": "Incerteza da medida e variação entre lotes; vale também para a tabela k(T), se houver",
    "c": "Inclui a incerteza do pico de umidade; vale também para a tabela c(T)",
    "eps_furos": "Cerâmica: cerca de 0,85 a 0,95; o limite baixo cobre superfícies mais lisas",
    "h": "A EN 1991-1-2 usa 25 W/m²·K na face exposta ao incêndio-padrão",
    "eps": "A EN 1991-1-2 recomenda 0,8 para a superfície do elemento",
    "alfa": "Cerâmica vermelha: cerca de 4 a 8 ×10⁻⁶/°C; substitui o α e escala a tabela de deformação térmica",
    "ktr": "De sem fluência (0) ao limite superior do concreto (2,35; Anderberg e Thelandersson, 1976)",
    "res": "Dispersão da resistência dos prismas; multiplica as tensões da lei da parede",
    "rig": "Dispersão do módulo dos prismas; multiplica o módulo da lei da parede",
    "imp": "Flecha inicial H/n no meio da altura; a amplitude varia linearmente entre os dois extremos",
    "exc": "Aplicada igual no topo e na base; positiva para o lado exposto ao fogo",
}


def _campos_sensibilidade():
    from .sensibilidade import PARAMETROS
    cs = [
        _n("trajetorias", "Trajetórias r", 10, "", "10 a 20 é o usual; cada trajetória custa (parâmetros + 1) "
           "rodadas da parede", minimo=2, maximo=100, passo=1),
        _n("niveis", "Níveis da grade p", 4, "", "Par; com 4 níveis, cada passo percorre 2/3 da faixa",
           minimo=4, maximo=10, passo=1),
        _n("semente", "Semente do sorteio das trajetórias", 1, "", "A mesma semente repete o mesmo plano",
           minimo=0, passo=1),
        _n("processos", "Processos em paralelo", 0, "", "0 = automático (núcleos do computador menos um)",
           minimo=0, maximo=64, passo=1),
        _b("malha_rapida", "Malha rápida na triagem (seção com 10 mm e 16 elementos de viga)", True,
           "Cerca de 6 vezes mais rápido que a malha recomendada (5 mm e 32 elementos). Na parede BCEPQ, os "
           "tempos mudaram menos de 7 min e a ordem dos parâmetros foi a mesma nas duas malhas"),
    ]
    for d in PARAMETROS:
        ch = d["chave"]
        se = {"sensibilidade.var_" + ch: [True]}
        cs.append(_b("var_" + ch, f"Variar: {d['nome'][0].lower() + d['nome'][1:]}", True, _AJUDA_SENS.get(ch, ""),
                     quebra=True))
        if d["tipo"] == "rel":
            cs.append(_n("faixa_" + ch, f"{d['nome']}: faixa", d["padrao"][0], "± %", minimo=0, maximo=90,
                         mostrar_se=se))
        else:
            lim = ({"minimo": 0, "maximo": 1} if ch in ("eps_furos", "eps") else
                   {"minimo": 0} if ch in ("h", "alfa", "ktr", "imp") else {})
            cs.append(_n("min_" + ch, f"{d['nome']}: de", d["padrao"][0], d["unidade"], mostrar_se=se, **lim))
            cs.append(_n("max_" + ch, f"{d['nome']}: até", d["padrao"][1], d["unidade"], mostrar_se=se, **lim))
    return cs


ESQUEMA = [
    {"id": "geometria", "titulo": "Geometria", "descricao":
        "Corpo de prova de blocos maciços, com furos circulares verticais (alinhados no meio da largura), com "
        "células retangulares verticais (blocos de concreto e cerâmicos estruturais, com ou sem graute) ou com furos "
        "horizontais (blocos de vedação): bloco isolado (1 bloco, 1 fiada), prisma (padrão) ou mini parede (mais "
        "blocos por fiada e mais fiadas), com junta seca ou de argamassa.", "campos": [
        _o("bloco_biblioteca", "Bloco da biblioteca", "bcepq", _bib.opcoes(_bib.BLOCOS),
           "Preenche a geometria, a junta e os materiais; depois troque o que tiver medido"),
        _n("comprimento_bloco", "Comprimento do bloco", 280, "mm", "De eixo a eixo: com argamassa, bloco + junta",
           minimo=40),
        _n("largura", "Largura do bloco", 140, "mm", minimo=40),
        _n("altura_fiada", "Altura da fiada", 60, "mm", "De eixo a eixo: com argamassa, bloco + junta", minimo=20),
        _o("forma_furo", "Furação", "circular",
           [["circular", "Furos redondos verticais (ou bloco maciço)"],
            ["retangular", "Células retangulares verticais (concreto, cerâmico estrutural)"],
            ["horizontal", "Furos horizontais (blocos de vedação)"]],
           "Os furos horizontais atravessam o bloco ao longo do comprimento"),
        _n("diametro_furo", "Diâmetro dos furos", 84, "mm", "0 = bloco maciço", minimo=0, mostrar_se=CIRC),
        _n("furos_por_bloco", "Furos ou células por bloco", 2, "", "Par quando há meio-bloco ou contrafiada; "
           "0 = bloco maciço", minimo=0, passo=1, mostrar_se=VERT),
        _n("furos_largura", "Furos na espessura do bloco", 2, "", "Colunas de furos horizontais", minimo=1, passo=1,
           mostrar_se=HOR),
        _n("furos_altura", "Furos na altura do bloco", 4, "", "Fileiras de furos horizontais (8 furos = 2 × 4)",
           minimo=1, passo=1, mostrar_se=HOR),
        _n("parede_externa", "Espessura das paredes externas (faces)", 25, "mm",
           "Paredes longitudinais do bloco", minimo=3, mostrar_se=RET_HOR),
        _n("septo", "Espessura dos septos", 25, "mm", "Paredes internas entre as células ou entre os furos",
           minimo=3, mostrar_se=RET_HOR),
        _n("parede_leito", "Espessura das paredes de leito", 8, "mm",
           "Paredes de cima e de baixo, onde a argamassa assenta", minimo=3, mostrar_se=HOR),
        _n("parede_topo", "Espessura das paredes de topo", 25, "mm", "Paredes transversais nas pontas do bloco",
           minimo=3, mostrar_se=RET),
        _n("celulas_parede", "Células pequenas por bloco em cada parede externa", 0, "",
           "Blocos cerâmicos de paredes vazadas; 0 = paredes maciças", minimo=0, passo=1, mostrar_se=RET),
        _n("parede_minima", "Parede em volta das células pequenas", 8, "mm", minimo=3,
           mostrar_se={"geometria.forma_furo": ["retangular"]}),
        _o("grauteamento", "Graute", "nenhum",
           [["nenhum", "Sem graute"], ["todos", "Todos os furos ou células"],
            ["alternados", "Furos ou células alternados (um sim, um não)"]],
           "O graute preenche os furos ou as células principais em toda a altura", mostrar_se=VERT),
        _n("blocos_por_fiada", "Blocos por fiada", 1.5, "", "Múltiplo de 0,5 (1,5 = bloco e meio)", minimo=0.5, passo=0.5),
        _n("n_fiadas", "Número de fiadas", 5, "", minimo=1, passo=1),
        _o("amarracao", "Amarração", "contrafiada",
           [["contrafiada", "Contrafiada (juntas verticais desencontradas)"], ["prumo", "Junta a prumo"]]),
        _o("tipo_junta", "Tipo de junta", "seca",
           [["seca", "Seca (camada fina com contato, folgas e atrito)"],
            ["argamassa", "Argamassa (camada contínua de espessura real, aderida ao bloco)"]]),
        _n("espessura_junta", "Espessura da junta", 2, "mm",
           "Seca: camada fina (2 mm). Argamassa: espessura real (cerca de 10 mm)", minimo=0.5),
        _o("assentamento", "Assentamento da argamassa", "total",
           [["total", "Em toda a face do bloco"], ["lateral", "Só nas paredes externas (furos e septos sem argamassa)"]],
           mostrar_se=ARG),
        _b("juntas_verticais_preenchidas", "Juntas verticais preenchidas com argamassa", True, mostrar_se=ARG),
        _b("capeamento", "Capeamento no topo e na base", True),
        _n("espessura_placa", "Espessura das placas de aço", 25, "mm", "Entram só na análise térmica", minimo=0),
    ]},
    {"id": "ceramica", "titulo": "Materiais", "descricao":
        "Propriedades do bloco e da argamassa. As tabelas em função da temperatura são opcionais.", "campos": [
        _o("material_biblioteca", "Material do bloco (biblioteca)", "bcepq", _bib.opcoes(_bib.MATERIAIS),
           "Preenche as propriedades com fonte; as demais ficam como estão"),
        _n("k", "Condutividade térmica", 0.52, "W/m·K", minimo=0.01),
        _n("rho", "Massa específica", 1210, "kg/m³", minimo=100),
        _n("c", "Calor específico", 1285, "J/kg·K", minimo=100),
        _n("E_bloco_bruto", "Módulo do bloco (área bruta)", 2801, "MPa",
           "O módulo do material sai dividindo pela razão área líquida/bruta da geometria", minimo=1),
        _n("nu", "Coeficiente de Poisson", 0.18, "", minimo=0, maximo=0.49),
        _n("alfa", "Coeficiente de dilatação térmica", 6.0, "×10⁻⁶/°C", "Provisório no BCEPQ", minimo=0),
        _t("tabela_k", "k em função de T", ["T (°C)", "k (W/m·K)"]),
        _t("tabela_E_fator", "Fator de E em função de T", ["T (°C)", "fator"]),
        _t("tabela_eps_termica", "Deformação térmica livre (substitui α)", ["T (°C)", "deformação (‰)"]),
        _t("tabela_c", "Calor específico em função de T", ["T (°C)", "c (J/kg·K)"],
           "Use para o pico de umidade do concreto e do solo-cimento"),
        _t("tabela_rho", "Massa específica em função de T", ["T (°C)", "ρ (kg/m³)"]),
        _o("argamassa_biblioteca", "Argamassa (biblioteca)", "1:1:6", _bib.opcoes(_bib.ARGAMASSAS),
           "Mohamad, Lourenço e Roman (2006); térmicas da NBR 15220-2", mostrar_se=ARG),
        _n("arg_fc", "Resistência à compressão da argamassa", 4.2, "MPa", minimo=0.1, mostrar_se=ARG),
        _n("arg_ft", "Resistência à tração da argamassa", 0.78, "MPa", "Estimativa: 0,3 fc^(2/3)", minimo=0.01,
           mostrar_se=ARG),
        _n("arg_E", "Módulo da argamassa", 4033, "MPa", "Tangente inicial", minimo=1, mostrar_se=ARG),
        _n("arg_nu", "Poisson da argamassa", 0.20, "", minimo=0, maximo=0.49, mostrar_se=ARG),
        _n("arg_k", "Condutividade térmica da argamassa", 1.15, "W/m·K", minimo=0.01, mostrar_se=ARG),
        _n("arg_rho", "Massa específica da argamassa", 2000, "kg/m³", minimo=100, mostrar_se=ARG),
        _n("arg_c", "Calor específico da argamassa", 1000, "J/kg·K", minimo=100, mostrar_se=ARG),
        _n("arg_alfa", "Dilatação térmica da argamassa", 10.0, "×10⁻⁶/°C", "Provisório", minimo=0, mostrar_se=ARG),
        _o("graute_biblioteca", "Graute (biblioteca)", "G1", _bib.opcoes(_bib.GRAUTES),
           "Camacho et al. (2015); E pela EN 1992-1-1; térmicas da NBR 15220-2", mostrar_se=GRA),
        _n("gra_fc", "Resistência à compressão do graute", 17.7, "MPa", minimo=0.5, mostrar_se=GRA),
        _n("gra_ft", "Resistência à tração do graute", 2.04, "MPa", "Estimativa: 0,3 fc^(2/3)", minimo=0.05,
           mostrar_se=GRA),
        _n("gra_E", "Módulo do graute", 26100, "MPa", "Estimativa: 22 (fc/10)^0,3 GPa (EN 1992-1-1)", minimo=1,
           mostrar_se=GRA),
        _n("gra_nu", "Poisson do graute", 0.20, "", minimo=0, maximo=0.49, mostrar_se=GRA),
        _n("gra_k", "Condutividade térmica do graute", 1.75, "W/m·K", minimo=0.01, mostrar_se=GRA),
        _n("gra_rho", "Massa específica do graute", 2400, "kg/m³", minimo=100, mostrar_se=GRA),
        _n("gra_c", "Calor específico do graute", 1000, "J/kg·K", minimo=100, mostrar_se=GRA),
        _b("gra_en1992", "Graute a quente como o concreto silicoso da EN 1992-1-2", True,
           "Deformação térmica e redução de E, fc e ft com a temperatura (Tabela 3.1)", mostrar_se=GRA),
        _n("gra_alfa", "Dilatação térmica do graute", 10.0, "×10⁻⁶/°C", "Provisório", minimo=0,
           mostrar_se={**GRA, "ceramica.gra_en1992": [False]}),
    ]},
    {"id": "modelo", "titulo": "Modelos e elementos", "descricao":
        "Tipo de análise, elemento finito e modelos constitutivos.", "campos": [
        _o("analise", "Análise", "nao_linear", [
            ["termica", "Térmica 3D (campos de temperatura)"],
            ["linear", "Termoelástica linear 3D (tensões e deformações)"],
            ["nao_linear", "Não linear 3D até a ruptura (dano e juntas com atrito)"],
            ["secao_2d", "Seção 2D (térmica e tensões em estado plano generalizado)"],
            ["parede", "Parede: faixa vertical com segunda ordem (seção 2D + vigas com fibras)"],
            ["sensibilidade", "Sensibilidade da parede no ISO 834 (Morris: quais entradas mais pesam)"]],
           "Parede: a altura vem do número de fiadas e o comprimento, dos blocos por fiada; a seção "
           "de um bloco (com os furos) é aquecida em 2D e cada triângulo de cerâmica vira uma fibra. "
           "Sensibilidade: roda a parede muitas vezes, variando as entradas incertas nas faixas da página "
           "Sensibilidade, e ordena quais mais mudam os tempos de resistência e de isolamento"),
        _o("elemento", "Elemento finito", "W6", [
            ["W6", "Cunha linear de 6 nós (W6)"],
            ["H8", "Hexaedro linear de 8 nós (H8)"],
            ["W15", "Cunha quadrática de 15 nós (W15)"],
            ["T3", "Triângulo linear de 3 nós (seção 2D)"]],
           "W6: padrão do artigo. H8: hexaedros sobre os quadriláteros da malha em O, metade dos "
           "elementos. W15: cunha quadrática, cerca de 4,7 vezes mais nós que a W6 na mesma malha; "
           "use malhas mais grossas (por exemplo, 16 × 30 mm). Na seção 2D o elemento é sempre o T3",
           mostrar_se={"modelo.analise": ["termica", "linear", "nao_linear", "secao_2d"]}),
        _o("lei_parede", "Material da parede", "bcepq", [
            ["bcepq", "Curvas do prisma BCEPQ a 30, 600, 800 e 1000 °C (modelo não linear calibrado)"],
            ["prisma", "Calibrar agora pelo modelo não linear do prisma a 30 °C (alguns minutos a mais)"]],
           "Lei tensão × deformação homogeneizada, sem tração (as juntas abrem). Na calibração pelo prisma, "
           "a curva de 30 °C sai das entradas atuais e as demais temperaturas seguem as razões do BCEPQ",
           mostrar_se=PAR),
        _n("fator_resistencia_parede", "Fator na resistência da alvenaria (lei da parede)", 1.0, "",
           "Multiplica as tensões da lei em todas as temperaturas; 1 = lei sem ajuste", minimo=0.1, mostrar_se=PAR),
        _n("fator_rigidez_parede", "Fator na rigidez da alvenaria (lei da parede)", 1.0, "",
           "Multiplica o módulo da lei (as deformações de pico mudam na razão inversa); 1 = lei sem ajuste",
           minimo=0.1, mostrar_se=PAR),
        _n("E_junta", "E da camada de junta (modelo linear)", 84.4, "MPa",
           "Calibrado para E_prisma = 1.100 MPa na malha de 8 × 12 mm", minimo=0.1, mostrar_se={"modelo.analise": ["linear"], "geometria.tipo_junta": ["seca"]}),
        _b("calibrar_junta", "Calibrar E da junta pelo módulo do prisma", False, mostrar_se={"modelo.analise": ["linear"], "geometria.tipo_junta": ["seca"]}),
        _n("E_prisma_alvo", "Módulo do prisma para a calibração", 1100, "MPa", minimo=10, mostrar_se={"modelo.analise": ["linear"], "geometria.tipo_junta": ["seca"]}),
        _o("modelo_ceramica", "Modelo de dano do bloco, da argamassa e do graute", "dano_mazars", [
            ["dano_mazars", "Dano de tração e compressão com fechamento de fissuras (Mazars)"],
            ["dano_isotropico", "Dano sem fechamento de fissuras"]], mostrar_se=NL),
        _n("ft", "Resistência à tração do bloco ft", 0.66, "MPa", "Gooch et al. (2023)", minimo=0.01, mostrar_se=NL),
        _n("fc", "Resistência à compressão do bloco fc", 2.48, "MPa",
           "2,48 = efetiva calibrada a 30 °C; 5,76 = derivada do bloco (BPQ)", minimo=0.1, mostrar_se=NLP),
        _n("du_t", "Índice de ductilidade em tração G_f/ft", 0.018, "mm", minimo=0.001, mostrar_se=NL),
        _n("du_c", "Índice de ductilidade em compressão G_c/fc", 0.33, "mm", minimo=0.01, mostrar_se=NL),
        _t("tabela_ft_fator", "Fator de ft em função de T", ["T (°C)", "fator"], mostrar_se=NL),
        _t("tabela_fc_fator", "Fator de fc em função de T", ["T (°C)", "fator"], mostrar_se=NL),
        _n("mu", "Atrito da junta seca μ", 0.60, "", "Coesão nula", minimo=0, mostrar_se={"modelo.analise": ["nao_linear"], "geometria.tipo_junta": ["seca"]}),
        _n("desvio_folga", "Desvio-padrão das folgas das juntas", 0.028, "mm", "0 = contato perfeito", minimo=0,
           mostrar_se={"modelo.analise": ["nao_linear"], "geometria.tipo_junta": ["seca"]}),
        _n("contato_inicial", "Área em contato sem carga", 20, "%", minimo=1, maximo=100, mostrar_se={"modelo.analise": ["nao_linear"], "geometria.tipo_junta": ["seca"]}),
        _n("comprimento_corr", "Comprimento de correlação das folgas", 40, "mm", minimo=5, mostrar_se={"modelo.analise": ["nao_linear"], "geometria.tipo_junta": ["seca"]}),
        _n("E_contato", "E da junta seca e do capeamento onde há contato", 3907, "MPa", minimo=1, mostrar_se=NL),
        _n("semente", "Semente do campo de folgas", 1, "", minimo=0, passo=1, mostrar_se={"modelo.analise": ["nao_linear"], "geometria.tipo_junta": ["seca"]}),
        _b("fluencia", "Fluência transitória (deformação térmica induzida pela carga)", True,
           "No primeiro aquecimento sob compressão, ε_tr = φ(T) σ / fc20, com φ(T) = k_tr ε_th(T) "
           "(Anderberg e Thelandersson, 1976) ou uma tabela medida no ensaio transiente", mostrar_se=NLP),
        _n("k_tr_bloco", "Coeficiente k_tr do bloco", 0.0, "",
           "Cerâmica queimada: sem dados publicados (0). Concreto com as curvas da EN 1992-1-2: 0, porque a "
           "fluência já está nelas. Concreto sem essas curvas: 1,8 a 2,35", minimo=0, mostrar_se=FLU),
        _t("tabela_fluencia_bloco", "φ(T) do bloco, do ensaio transiente (substitui k_tr)", ["T (°C)", "φ (‰)"],
           "φ = (deformação do ensaio transiente − livre − elástica) × fc20 / σ", mostrar_se=FLU),
        _n("k_tr_argamassa", "Coeficiente k_tr da argamassa", 2.0, "",
           "Prakash et al. (2020), argamassa de cimento", minimo=0, mostrar_se={**FLU, **ARG}),
        _n("k_tr_graute", "Coeficiente k_tr do graute", 0.0, "",
           "0 com as curvas da EN 1992-1-2 (fluência já implícita)", minimo=0, mostrar_se={**FLU, **GRA}),
    ]},
    {"id": "ensaio", "titulo": "Ensaio e carga", "descricao":
        "Programa térmico do forno e forma de carregamento.", "campos": [
        _o("tipo", "Tipo de ensaio", "quente", [
            ["ambiente", "Compressão a temperatura ambiente"],
            ["quente", "Compressão a quente (rampa, patamar e carga no forno)"],
            ["residual", "Compressão residual (aquece, resfria e depois carrega)"],
            ["iso834", "Incêndio-padrão ISO 834 numa face, sob carga constante"],
            ["transiente", "Ensaio transiente: carga constante e aquecimento no forno até a ruptura"]]),
        _n("T_inicial", "Temperatura inicial", 30, "°C"),
        _n("taxa", "Taxa de aquecimento", 5, "°C/min", minimo=0.1, mostrar_se=AQ),
        _n("T_alvo", "Temperatura-alvo do forno", 600, "°C", mostrar_se=AQ),
        _n("patamar", "Patamar na temperatura-alvo", 30, "min", minimo=0, mostrar_se=AQ),
        _n("taxa_resfriamento", "Taxa de resfriamento", 2, "°C/min", minimo=0.1, mostrar_se={"ensaio.tipo": ["residual"]}),
        _n("patamar_final", "Espera após o resfriamento", 60, "min", minimo=0, mostrar_se={"ensaio.tipo": ["residual"]}),
        _n("duracao_iso", "Duração da exposição ISO 834", 120, "min", minimo=5, mostrar_se={"ensaio.tipo": ["iso834"]}),
        _n("h_conveccao", "Coeficiente de convecção do forno", 10, "W/m²·K", "ISO 834: 25 na face exposta", minimo=0),
        _n("emissividade", "Emissividade da superfície", 0.8, "", minimo=0, maximo=1),
        _n("emissividade_furos", "Emissividade das paredes dos furos", 0.9, "", minimo=0, maximo=1),
        _n("carga", "Carga", 97.9, "kN",
           "Linear: carga aplicada no fim do programa térmico. ISO 834 e ensaio transiente: carga constante. "
           "Não linear (demais ensaios): só a referência para o módulo de corda; a ruptura sai do cálculo. "
           "Parede: carga total no comprimento da parede (ISO 834); a 30 °C, só a referência", minimo=0),
        _n("carga_assentamento", "Carga de assentamento no aquecimento", 2, "kN",
           "No ensaio transiente e no ISO 834 do 3D, vale a carga constante", minimo=0, mostrar_se=P3D),
        _o("apoio_base", "Apoio na base", "articulado",
           [["articulado", "Articulado (gira livremente)"], ["engastado", "Engastado (rotação impedida)"]],
           mostrar_se=PAR),
        _o("apoio_topo", "Apoio no topo", "articulado",
           [["articulado", "Articulado (gira livremente)"], ["engastado", "Engastado (rotação impedida)"]],
           "O topo sempre desloca na vertical sob a carga", mostrar_se=PAR),
        _n("excentricidade_topo", "Excentricidade da carga no topo", 0, "mm",
           "Positiva para o lado exposto ao fogo", mostrar_se=PAR),
        _n("excentricidade_base", "Excentricidade da reação na base", 0, "mm",
           "Positiva para o lado exposto ao fogo", mostrar_se=PAR),
        _n("imperfeicao", "Imperfeição inicial: altura dividida por", 450, "",
           "Flecha inicial senoidal no meio da altura, para o lado exposto (H/450 é um valor usual); "
           "0 = parede perfeita", minimo=0, mostrar_se=PAR),
        _n("niveis_carga", "Níveis de carga na curva tempo × carga", 5, "",
           "Repete o ISO 834 com 10% a 90% da capacidade a frio; 0 = só a carga informada", minimo=0, maximo=12,
           passo=1, mostrar_se={"modelo.analise": ["parede"], "ensaio.tipo": ["iso834"]}),
    ]},
    {"id": "malha", "titulo": "Malha e solução", "descricao":
        "Discretização e controle da solução numérica.", "campos": [
        _n("h", "Tamanho do elemento em planta", 12, "mm", "8 mm no térmico e linear do artigo; 12 mm no não linear; "
           "5 mm na seção 2D e na parede",
           minimo=2),
        _n("dz", "Altura máxima das camadas de elementos", 20, "mm", minimo=2, mostrar_se=P3D),
        _n("n_elementos", "Elementos de viga ao longo da altura", 32, "",
           "Número par; com 32, o tempo de resistência fica a menos de 1% do obtido com 64", minimo=4, passo=1,
           mostrar_se=PAR),
        _n("dt", "Passo de tempo da térmica", 30, "s", minimo=1),
        _n("delta_max", "Deslocamento máximo do prato", 0.8, "mm", minimo=0.05, mostrar_se=NL),
        _n("incremento", "Incremento de deslocamento", 0.02, "mm", minimo=0.001, mostrar_se=NL),
        _n("tol_u", "Tolerância em deslocamento", 1e-4, "", minimo=1e-8, mostrar_se=NL),
        _n("tol_d", "Tolerância em dano", 5e-3, "", minimo=1e-6, mostrar_se=NL),
        _n("n_campos_rampa", "Campos térmicos na rampa", 6, "", minimo=1, passo=1, mostrar_se=NL),
        _n("n_campos_patamar", "Campos térmicos no patamar", 2, "", minimo=1, passo=1, mostrar_se=NL),
        _n("intervalo_iso", "Intervalo entre campos no ISO 834", 5, "min", minimo=1, mostrar_se=NL),
    ]},
    {"id": "sensibilidade", "titulo": "Sensibilidade", "descricao":
        "Análise de sensibilidade global da parede no ISO 834 pelo método de Morris (efeitos elementares). "
        "Cada parâmetro marcado varia na faixa informada, em torno dos valores das outras páginas; o programa "
        "roda a parede r × (parâmetros + 1) vezes e ordena os parâmetros pelo efeito médio μ* nos tempos de "
        "resistência e de isolamento. As faixas padrão são pontos de partida: troque-as pela incerteza que "
        "você conhece em cada entrada.", "campos": _campos_sensibilidade()},
]


def projeto_padrao() -> Dict:
    p = {"nome": "Prisma BCEPQ", "versao": VERSAO}
    for g in ESQUEMA:
        p[g["id"]] = {c["id"]: copy.deepcopy(c["padrao"]) for c in g["campos"]}
    return p


def _pred_bloco(chave, extra):
    """Projeto pronto a partir de um bloco da biblioteca e das mudanças extras."""
    q = copy.deepcopy(_bib.BLOCOS[chave]["mudancas"])
    q.setdefault("geometria", {})["bloco_biblioteca"] = chave
    for g, d in extra.items():
        q.setdefault(g, {}).update(d)
    return q


PREDEFINIDOS = {
    "BCEPQ — 30 °C, não linear": {"ensaio": {"tipo": "ambiente", "carga": 98.7}},
    "BCEPQ — 600 °C a quente, não linear": {"ensaio": {"tipo": "quente", "T_alvo": 600, "carga": 97.9}},
    "BCEPQ — 800 °C a quente, não linear": {"ensaio": {"tipo": "quente", "T_alvo": 800, "carga": 98.9}},
    "BCEPQ — 1000 °C a quente, não linear": {"ensaio": {"tipo": "quente", "T_alvo": 1000, "carga": 114.4}},
    "BCEPQ — 600 °C, térmica e tensões (linear, malha do artigo)": {
        "modelo": {"analise": "linear"}, "ensaio": {"tipo": "quente", "T_alvo": 600, "carga": 97.9},
        "malha": {"h": 8, "dz": 12}},
    "BCEPQ — residual após 600 °C, não linear": {"ensaio": {"tipo": "residual", "T_alvo": 600}},
    "Parede BCEPQ — ISO 834, seção 2D (isolamento)": {
        "modelo": {"analise": "secao_2d", "elemento": "T3"},
        "ensaio": {"tipo": "iso834", "T_inicial": 20, "h_conveccao": 25, "duracao_iso": 240}, "malha": {"h": 5}},
    "Tijolo maciço com argamassa 1:1:6 — prisma a 30 °C, não linear": _pred_bloco("macico_10x6x22", {
        "modelo": {"analise": "nao_linear"}, "ensaio": {"tipo": "ambiente", "carga": 50},
        "malha": {"h": 8, "dz": 10, "delta_max": 1.2, "incremento": 0.02}}),
    "Tijolo maciço com argamassa 1:1:6 — 600 °C, térmica e tensões (linear)": _pred_bloco("macico_10x6x22", {
        "modelo": {"analise": "linear"}, "ensaio": {"tipo": "quente", "T_alvo": 600, "carga": 50},
        "malha": {"h": 8, "dz": 10}}),
    "Bloco de concreto 14 × 19 × 29 (Camacho et al., 2015) — prisma oco a 20 °C, não linear": _pred_bloco(
        "concreto_m15_29", {"geometria": {"grauteamento": "nenhum"}, "modelo": {"analise": "nao_linear"},
                            "ensaio": {"tipo": "ambiente", "T_inicial": 20, "carga": 180},
                            "malha": {"h": 12, "dz": 20, "delta_max": 1.5, "incremento": 0.03}}),
    "Bloco de concreto 14 × 19 × 29 (Camacho et al., 2015) — prisma grauteado (G1) a 20 °C, não linear": _pred_bloco(
        "concreto_m15_29", {"geometria": {"grauteamento": "todos"}, "modelo": {"analise": "nao_linear"},
                            "ensaio": {"tipo": "ambiente", "T_inicial": 20, "carga": 470},
                            "malha": {"h": 12, "dz": 20, "delta_max": 2.0, "incremento": 0.03}}),
    "Bloco de concreto 14 × 19 × 39 grauteado — ISO 834, seção 2D (isolamento)": _pred_bloco(
        "concreto_m15_39", {"geometria": {"grauteamento": "alternados", "blocos_por_fiada": 1},
                            "ceramica": {**_bib.GRAUTES["G1"]["mudancas"]["ceramica"], "graute_biblioteca": "G1"},
                            "modelo": {"analise": "secao_2d", "elemento": "T3"},
                            "ensaio": {"tipo": "iso834", "T_inicial": 20, "h_conveccao": 25, "duracao_iso": 240},
                            "malha": {"h": 5}}),
    "Bloco de vedação de 8 furos — ISO 834, seção 2D vertical (isolamento)": _pred_bloco(
        "vedacao_8_furos", {"modelo": {"analise": "secao_2d", "elemento": "T3"},
                            "ensaio": {"tipo": "iso834", "T_inicial": 20, "h_conveccao": 25, "duracao_iso": 240},
                            "malha": {"h": 4}}),
    "Bloco de vedação de 8 furos — prisma a 20 °C, não linear": _pred_bloco(
        "vedacao_8_furos", {"modelo": {"analise": "nao_linear"},
                            "ensaio": {"tipo": "ambiente", "T_inicial": 20, "carga": 30},
                            "malha": {"h": 12, "dz": 20, "delta_max": 1.2, "incremento": 0.03}}),
    "BCEPQ — ensaio transiente sob 49 kN (50% da carga a 30 °C), não linear": {
        "ensaio": {"tipo": "transiente", "carga": 49, "T_alvo": 800, "taxa": 5, "patamar": 30}},
    "Tijolo maciço com argamassa 1:1:6 — ensaio transiente sob 46 kN (40%), não linear": _pred_bloco(
        "macico_10x6x22", {"modelo": {"analise": "nao_linear"},
                           "ensaio": {"tipo": "transiente", "carga": 46, "T_alvo": 800, "taxa": 5, "patamar": 30,
                                      "T_inicial": 20}}),
    "Parede BCEPQ 2,94 × 3,00 m — ISO 834 sob carga": {
        "geometria": {"blocos_por_fiada": 10.5, "n_fiadas": 50},
        "modelo": {"analise": "parede"},
        "ensaio": {"tipo": "iso834", "T_inicial": 20, "h_conveccao": 25, "duracao_iso": 240, "carga": 150},
        "malha": {"h": 5}},
    "Parede BCEPQ 2,94 × 3,00 m — capacidade a 20 °C": {
        "geometria": {"blocos_por_fiada": 10.5, "n_fiadas": 50},
        "modelo": {"analise": "parede"},
        "ensaio": {"tipo": "ambiente", "T_inicial": 20, "carga": 150}, "malha": {"h": 5}},
    "Parede BCEPQ 2,94 × 3,00 m — sensibilidade no ISO 834 (Morris, 11 parâmetros)": {
        "geometria": {"blocos_por_fiada": 10.5, "n_fiadas": 50},
        "modelo": {"analise": "sensibilidade"},
        "ensaio": {"tipo": "iso834", "T_inicial": 20, "h_conveccao": 25, "duracao_iso": 300, "carga": 150},
        "malha": {"h": 5}},
}


def aplicar(p: Dict, mudancas: Dict) -> Dict:
    q = copy.deepcopy(p)
    for g, d in mudancas.items():
        if isinstance(d, dict):
            q.setdefault(g, {}).update(d)
        else:
            q[g] = d
    return q


def completar(p: Dict) -> Dict:
    """Preenche campos ausentes com os padrões (projetos de versões anteriores)."""
    return aplicar(projeto_padrao(), {k: v for k, v in p.items()})


# ----------------------------------------------------------------------------
# Conversão para os objetos do modelo
# ----------------------------------------------------------------------------
def _tab(t, fy=1.0):
    if not t:
        return None
    a = [(float(r[0]), float(r[1]) * fy) for r in t if r is not None and len(r) >= 2
         and r[0] not in ("", None) and r[1] not in ("", None)]
    a.sort()
    return a or None


def geometria(p) -> Geometria:
    g = p["geometria"]
    return Geometria(comprimento_bloco=g["comprimento_bloco"] / 1e3, largura=g["largura"] / 1e3,
                     altura_fiada=g["altura_fiada"] / 1e3, diametro_furo=g["diametro_furo"] / 1e3,
                     n_fiadas=int(g["n_fiadas"]), espessura_placa=g["espessura_placa"] / 1e3,
                     furos_por_bloco=int(g["furos_por_bloco"]), blocos_por_fiada=float(g["blocos_por_fiada"]),
                     amarracao=g["amarracao"], tipo_junta=g.get("tipo_junta", "seca"),
                     assentamento=g.get("assentamento", "total"),
                     juntas_verticais_preenchidas=bool(g.get("juntas_verticais_preenchidas", True)),
                     espessura_junta=g.get("espessura_junta", 2) / 1e3,
                     forma_furo=g.get("forma_furo", "circular"), parede_externa=g.get("parede_externa", 25) / 1e3,
                     septo=g.get("septo", 25) / 1e3, parede_topo=g.get("parede_topo", 25) / 1e3,
                     celulas_parede=int(g.get("celulas_parede", 0)), parede_minima=g.get("parede_minima", 8) / 1e3,
                     grauteamento=("nenhum" if g.get("forma_furo") == "horizontal"
                                   else g.get("grauteamento", "nenhum")),
                     furos_largura=int(g.get("furos_largura", 2)), furos_altura=int(g.get("furos_altura", 4)),
                     parede_leito=g.get("parede_leito", 8) / 1e3)


def razao_area_liquida(geo: Geometria) -> float:
    """Área do bloco (sem furos) / área bruta, na fiada 0: converte E_bloco_bruto em E do material."""
    return 1.0 - geo.area_furos(0) / (geo.comprimento * geo.largura)


def razao_area_solida(geo: Geometria) -> float:
    """Área resistente (bloco + graute) / área bruta, na fiada 0."""
    return 1.0 - geo.area_furos(0, graute=True) / (geo.comprimento * geo.largura)


def materiais(p, geo: Geometria) -> Materiais:
    c, mo, en = p["ceramica"], p["modelo"], p["ensaio"]
    M = Materiais()
    M.ceramica.k, M.ceramica.rho, M.ceramica.c = c["k"], c["rho"], c["c"]
    M.ceramica.E_bloco_bruto = c["E_bloco_bruto"] * 1e6
    M.ceramica.razao_area_liquida = razao_area_liquida(geo)
    M.ceramica.nu, M.ceramica.alfa = c["nu"], c["alfa"] * 1e-6
    M.ceramica.T_ref = en["T_inicial"]
    M.ceramica.tabela_k = _tab(c.get("tabela_k"))
    M.ceramica.tabela_E_fator = _tab(c.get("tabela_E_fator"))
    eps = _tab(c.get("tabela_eps_termica"), 1e-3)
    M.ceramica.tabela_eps_termica = eps
    M.ceramica.tabela_c = _tab(c.get("tabela_c"))
    M.ceramica.tabela_rho = _tab(c.get("tabela_rho"))
    a = M.argamassa
    a.k, a.rho, a.c = c["arg_k"], c["arg_rho"], c["arg_c"]
    a.E_bloco_bruto, a.razao_area_liquida = c["arg_E"] * 1e6, 1.0
    a.nu, a.alfa, a.T_ref = c["arg_nu"], c["arg_alfa"] * 1e-6, en["T_inicial"]
    gr = M.graute
    gr.k, gr.rho, gr.c = c.get("gra_k", 1.75), c.get("gra_rho", 2400), c.get("gra_c", 1000)
    gr.E_bloco_bruto, gr.razao_area_liquida = c.get("gra_E", 26100) * 1e6, 1.0
    gr.nu, gr.alfa, gr.T_ref = c.get("gra_nu", 0.2), c.get("gra_alfa", 10) * 1e-6, en["T_inicial"]
    if c.get("gra_en1992", True):
        from .parametros import EN1992_E_FATOR, EN1992_EPS_TERMICA
        gr.tabela_E_fator, gr.tabela_eps_termica = list(EN1992_E_FATOR), list(EN1992_EPS_TERMICA)
    t = p["geometria"]["espessura_junta"] / 1e3
    for j in (M.junta, M.capeamento):
        j.espessura, j.E, j.nu = t, mo["E_junta"] * 1e6, c["nu"]
    M.ar.emissividade_paredes = en["emissividade_furos"]
    M.diametro_furo = max(geo.dimensao_furo, 1e-6)
    return M


def dano_arg(p) -> DanoCeramica:
    """Dano da argamassa: fc e ft informados; índices de ductilidade de Dauda et al. (2020)."""
    from .parametros import dano_argamassa
    c, mo = p["ceramica"], p["modelo"]
    d = dano_argamassa(fc=c["arg_fc"] * 1e6, ft=c["arg_ft"] * 1e6)
    d.fechamento_fissuras = mo["modelo_ceramica"] == "dano_mazars"
    return d


def dano_gra(p) -> DanoCeramica:
    """Dano do graute: fc e ft informados; mesmos índices de ductilidade da argamassa."""
    from .parametros import dano_graute
    c, mo = p["ceramica"], p["modelo"]
    d = dano_graute(fc=c.get("gra_fc", 17.7) * 1e6, ft=c.get("gra_ft", 2.04) * 1e6)
    d.fechamento_fissuras = mo["modelo_ceramica"] == "dano_mazars"
    if c.get("gra_en1992", True):
        from .parametros import EN1992_FC_FATOR, EN1992_FT_FATOR
        d.tabela_fc_fator, d.tabela_ft_fator = list(EN1992_FC_FATOR), list(EN1992_FT_FATOR)
    return d


def dano(p) -> DanoCeramica:
    mo = p["modelo"]
    return DanoCeramica(ft=mo["ft"] * 1e6, fc=mo["fc"] * 1e6, du_t=mo["du_t"] * 1e-3, du_c=mo["du_c"] * 1e-3,
                        fechamento_fissuras=(mo["modelo_ceramica"] == "dano_mazars"),
                        tabela_ft_fator=_tab(mo.get("tabela_ft_fator")), tabela_fc_fator=_tab(mo.get("tabela_fc_fator")))


def junta_atrito(p) -> JuntaAtrito:
    mo, en = p["modelo"], p["ensaio"]
    return JuntaAtrito(mu=mo["mu"], desvio_folga=mo["desvio_folga"] * 1e-3, contato_inicial=mo["contato_inicial"] / 100,
                       comprimento_corr=mo["comprimento_corr"] * 1e-3, E_contato=mo["E_contato"] * 1e6,
                       semente=int(mo["semente"]), pre_carga_N=en["carga_assentamento"] * 1e3)


def validar(p) -> List[str]:
    erros = []
    try:
        geo = geometria(p)
        erros += geo.validar()
    except Exception as e:  # noqa
        return [f"Geometria inválida: {e}"]
    mo, en, ma = p["modelo"], p["ensaio"], p["malha"]
    if mo["analise"] == "secao_2d":
        if en["tipo"] == "ambiente":
            erros.append("A seção 2D só calcula efeitos térmicos; escolha um ensaio com aquecimento.")
    elif mo["analise"] in PAREDES:
        if mo["analise"] == "sensibilidade":
            erros += _validar_sensibilidade(p)
        elif en["tipo"] not in ("ambiente", "iso834"):
            erros.append("Na parede, escolha a compressão a temperatura ambiente (capacidade) ou o ISO 834 sob carga.")
        if en["tipo"] == "iso834" and en["carga"] <= 0:
            erros.append("No ISO 834 da parede, informe a carga constante (total no comprimento da parede).")
        t2 = p["geometria"]["largura"] / 2
        if max(abs(en["excentricidade_topo"]), abs(en["excentricidade_base"])) >= t2:
            erros.append("A excentricidade deve ser menor que metade da espessura da parede.")
        if ma["h"] > 10:
            erros.append("Na parede, use elementos de no máximo 10 mm na seção (recomendado: 5 mm).")
        if ma["n_elementos"] < 4:
            erros.append("Use pelo menos 4 elementos de viga ao longo da altura.")
    else:
        if mo["elemento"] == "T3":
            erros.append("O triângulo T3 é só para a seção 2D; no 3D use W6, H8 ou W15.")
    if en["tipo"] in ("quente", "residual", "transiente") and en["T_alvo"] <= en["T_inicial"]:
        erros.append("A temperatura-alvo deve ser maior que a inicial.")
    if en["tipo"] == "transiente":
        if mo["analise"] not in ("nao_linear", "termica"):
            erros.append("O ensaio transiente é do modelo não linear 3D (ou só da térmica).")
        elif mo["analise"] == "nao_linear" and en["carga"] <= 0:
            erros.append("No ensaio transiente, informe a carga constante.")
    tab = mo.get("tabela_fluencia_bloco") or []
    if any(len(r) != 2 or r[1] < 0 for r in tab):
        erros.append("A tabela φ(T) do bloco precisa de pares (T, φ) com φ ≥ 0.")
    if mo["analise"] == "linear" and en["tipo"] != "iso834" and en["carga"] <= 0:
        erros.append("A análise linear precisa de uma carga positiva.")
    if mo["analise"] == "nao_linear" and en["tipo"] == "iso834" and en["carga"] <= 0:
        erros.append("No ISO 834 com análise não linear, informe a carga constante.")
    g = p["geometria"]
    if g.get("tipo_junta") == "argamassa" and g["espessura_junta"] > 0.3 * g["altura_fiada"]:
        erros.append("A junta de argamassa deve ter no máximo 30% da altura da fiada (de eixo a eixo).")
    if g.get("bloco_biblioteca") in _bib.BLOCOS and _bib.BLOCOS[g["bloco_biblioteca"]].get("fase", 0) > _bib.FASE_ATUAL:
        erros.append("O bloco escolhido na biblioteca ainda não roda nesta versão; escolha outro ou 'Personalizado'.")
    if mo["analise"] not in PAREDES and ma["h"] > p["geometria"]["largura"] / 3:
        erros.append("Tamanho de elemento grande demais para a largura do bloco.")
    if mo["analise"] == "nao_linear" and ma["incremento"] >= ma["delta_max"]:
        erros.append("O incremento deve ser menor que o deslocamento máximo.")
    return erros


# ----------------------------------------------------------------------------
# Programa térmico
# ----------------------------------------------------------------------------
def programa_termico(p):
    """Devolve (f(t), t_final em s, marcos {nome: t}) do gás do forno."""
    en = p["ensaio"]
    Ti = en["T_inicial"]
    if en["tipo"] in ("quente", "transiente"):
        f, tf, tr = curvas.rampa_patamar(Ti, en["taxa"], en["T_alvo"], en["patamar"] * 60)
        return f, tf, {"fim da rampa": tr, "carga": tf}
    if en["tipo"] == "residual":
        f, tf, (t1, t2, t3) = curvas.rampa_patamar_resfriamento(Ti, en["taxa"], en["T_alvo"], en["patamar"] * 60,
                                                               en["taxa_resfriamento"], en["patamar_final"] * 60)
        return f, tf, {"fim da rampa": t1, "fim do patamar": t2, "fim do resfriamento": t3, "carga": tf}
    if en["tipo"] == "iso834":
        return curvas.iso834(Ti), en["duracao_iso"] * 60, {"fim": en["duracao_iso"] * 60}
    return curvas.constante(Ti), 0.0, {"carga": 0.0}


def _contornos_3d(p, m):
    en = p["ensaio"]
    f, tf, _ = programa_termico(p)
    if en["tipo"] == "iso834":
        return [Contorno(tm.areas_laterais_3d(m, faces=("y0",)), f, en["h_conveccao"], en["emissividade"], "exposta"),
                Contorno(tm.areas_laterais_3d(m, faces=("y1",)), curvas.constante(en["T_inicial"]), 4.0,
                         en["emissividade"], "não exposta")]
    return [Contorno(tm.areas_laterais_3d(m), f, en["h_conveccao"], en["emissividade"], "laterais")]


def _sondas_3d(m, geo):
    H = geo.altura_prisma
    zm = m.zs[np.argmin(np.abs(m.zs - H / 2))]
    W, Wm = m.m2d.largura, m.m2d.largura_modelo
    cf = geo.centros_furos
    x_septo = 0.5 * (cf[0] + cf[1]) if len(cf) >= 2 else (cf[0] if cf else geo.comprimento / 2)
    if geo.retangular:
        cel = sorted(c for c in geo.celulas(0, m.m2d.largura_faixa) if c[4])
        if len(cel) >= 2:
            x_septo = 0.5 * (cel[0][1] + cel[1][0])
    y_meio = min(W / 2, Wm)
    if len(cf) == 1:            # um furo: núcleo da parede entre o furo e a face
        y_meio = 0.5 * (W / 2 - geo.diametro_furo / 2)
    s = {"superfície (y = 0, meia altura)": ml.no_mais_proximo(m.xyz, (geo.comprimento / 3, 0.0, zm)),
         "núcleo do septo (meia altura)": ml.no_mais_proximo(m.xyz, (x_septo, y_meio, zm))}
    if not m.m2d.simetria:
        s["face oposta (y = W, meia altura)"] = ml.no_mais_proximo(m.xyz, (geo.comprimento / 3, W, zm))
    cer = np.unique(m.cun[(m.mat == CERAMICA) & (np.abs(m.xyz[m.cun[:, 0], 2] - H / 2) < H / 6)])
    s["média da cerâmica no terço central"] = cer
    return s


# ----------------------------------------------------------------------------
# Execução
# ----------------------------------------------------------------------------
class Cancelado(Exception):
    pass


class _Saida:
    """Resultados acumulados de uma execução."""

    def __init__(self, pasta):
        self.pasta = pasta
        self.resumo, self.tabelas, self.figuras, self.arquivos, self.avisos = [], [], [], [], []
        self.visual = None

    def pos(self, m, estados, descricao="", material=None):
        visual.exportar(os.path.join(self.pasta, "visual.json"), m, estados, descricao, material)
        self.visual = "visual.json"

    def item(self, rotulo, valor, unidade=""):
        if isinstance(valor, (float, np.floating)):
            valor = float(f"{valor:.4g}")
        self.resumo.append([rotulo, valor, unidade])

    def tabela(self, titulo, colunas, linhas, arquivo=None):
        self.tabelas.append({"titulo": titulo, "colunas": colunas, "linhas": linhas})
        if arquivo:
            salvar_csv(os.path.join(self.pasta, arquivo), linhas, colunas)
            self.arquivos.append({"nome": arquivo, "descricao": titulo})

    def figura(self, fig, arquivo, titulo):
        fig.savefig(os.path.join(self.pasta, arquivo), dpi=130, facecolor=fig.get_facecolor())
        self.figuras.append({"arquivo": arquivo, "titulo": titulo})
        import matplotlib.pyplot as plt
        plt.close(fig)

    def vtk(self, arquivo, m, **kw):
        salvar_vtk(os.path.join(self.pasta, arquivo), m, **kw)
        self.arquivos.append({"nome": arquivo, "descricao": "campos 3D para o ParaView"})

    def dict(self):
        return {"resumo": self.resumo, "tabelas": self.tabelas, "figuras": self.figuras,
                "arquivos": self.arquivos, "avisos": self.avisos, "visual": self.visual}


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e4e3df", "axes.edgecolor": "#52514e"})
    return plt


CORES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6"]


ELEMENTOS_3D = {"W6": "cunhas W6", "H8": "hexaedros H8", "W15": "cunhas W15"}


def malha_3d(p, geo=None):
    """Malha 3D do prisma com o elemento escolhido."""
    geo = geometria(p) if geo is None else geo
    ma = p["malha"]
    el = p["modelo"]["elemento"] if p["modelo"]["elemento"] in ELEMENTOS_3D else "W6"
    return ml.malha_prisma(geo, h=ma["h"] / 1e3, dz=ma["dz"] / 1e3, t_junta=p["geometria"]["espessura_junta"] / 1e3,
                           simetria=(p["ensaio"]["tipo"] != "iso834"), capeamento=p["geometria"]["capeamento"],
                           elemento=el)


def figura_malha(p, pasta=None, arquivo="malha.png"):
    """Figura da seção e do corte vertical da malha; devolve (caminho ou bytes, resumo)."""
    _plt()
    from matplotlib.figure import Figure          # sem pyplot: pode rodar em paralelo a uma análise
    from matplotlib.collections import PolyCollection
    geo = geometria(p)
    ma, en, mo = p["malha"], p["ensaio"], p["modelo"]
    if mo["analise"] in PAREDES:
        return _figura_parede(p, geo, pasta, arquivo)
    if mo["analise"] == "secao_2d":
        m2 = ml.malha_secao(geo, h=ma["h"] / 1e3, simetria=False, faixas=())
        info = ml.resumo(m2)
    else:
        m = malha_3d(p, geo)
        m2 = m.m2d
        info = ml.resumo(m)
    if m2.vertical:                          # seção vertical (furos horizontais), espessura × altura
        fig = Figure(figsize=(10, 4.2))
        a0 = fig.subplots(1, 1)
        mat2 = tm.material_secao(m2, geo)
        cores2 = {AR: "#cde2fb", CERAMICA: "#d9c7b0", ARGAMASSA: "#7f8fa6", JUNTA: "#eb6834"}
        pc = PolyCollection(m2.xy[m2.tri][:, :, ::-1] * 1e3, facecolors=[cores2[int(k)] for k in mat2],
                            edgecolors="#52514e", linewidths=0.25)
        a0.add_collection(pc); a0.autoscale(); a0.set_aspect("equal"); a0.grid(False)
        a0.set_title("Seção vertical de uma fiada (vale ao longo de todo o comprimento)\n"
                     "furos em azul, meias juntas em cinza-azulado", fontsize=9)
        a0.set_xlabel("y, espessura (mm; y = 0 é a face exposta)"); a0.set_ylabel("altura na fiada (mm)")
        fig.tight_layout()
        if pasta:
            caminho = os.path.join(pasta, arquivo)
            fig.savefig(caminho, dpi=130)
            return caminho, info
        import io
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110)
        return buf.getvalue(), info
    fig = Figure(figsize=(10, 3.6))
    ax = fig.subplots(1, 1 if mo["analise"] == "secao_2d" else 2,
                      gridspec_kw=None if mo["analise"] == "secao_2d" else {"width_ratios": [2.2, 1]})
    ax = np.atleast_1d(ax)
    pol, marca = m2.tri, m2.marca
    quad = mo["analise"] != "secao_2d" and m.tipo == "H8"
    if quad:
        pol, marca = m2.quad, m2.marca_quad
    _, gr0 = ml.estados_fiada(m2, geo, 0, quad=quad)
    cor = np.where(gr0, "#5b5f66", np.where(marca == 0, "#cde2fb", np.where(marca > 0, "#eb6834", "#d9c7b0")))
    pc = PolyCollection(m2.xy[pol] * 1e3, facecolors=cor, edgecolors="#52514e", linewidths=0.25)
    ax[0].add_collection(pc); ax[0].autoscale(); ax[0].set_aspect("equal")
    ax[0].set_title("Planta (faixas de junta vertical em laranja; os furos horizontais aparecem no corte)"
                    if geo.horizontal else
                    "Seção na fiada 1 (furos em azul, graute em cinza, faixas de junta vertical em laranja)"
                    if gr0.any() else "Seção na fiada 1 (furos em azul, faixas de junta vertical em laranja)",
                    fontsize=9, loc="left")
    ax[0].set_xlabel("x (mm)"); ax[0].set_ylabel("y (mm)"); ax[0].grid(False)
    if mo["analise"] != "secao_2d":
        cores = {AR: "#cde2fb", CERAMICA: "#d9c7b0", JUNTA: "#eb6834", CAPEAMENTO: "#eda100", ACO: "#8a8984",
                 ARGAMASSA: "#7f8fa6", GRAUTE: "#5b5f66"}
        # corte vertical: uma coluna de elementos em y ~ 0
        quads = []
        cs = []
        meio = geo.retangular or geo.grauteamento != "nenhum"     # corte pelos furos
        # furos horizontais: face x = 0 (espessura × altura), onde aparecem as fileiras de furos
        eixo, coords = (0, [1, 2]) if geo.horizontal else (1, [0, 2])
        v_corte = 0.0 if geo.horizontal else (geo.largura / 2 if meio else 0.0)
        for f in m.elem.faces_vis:
            if len(f) != 4:
                continue
            fc = m.cun[:, f]
            sel = np.isclose(m.xyz[fc, eixo], v_corte).all(1)
            if geo.horizontal:                        # só a primeira fiada, para os furos aparecerem
                sel &= (m.xyz[fc, 2] <= geo.altura_fiada + 1e-9).all(1)
            for e in np.nonzero(sel)[0]:
                quads.append(m.xyz[fc[e]][:, coords])
                cs.append(cores[int(m.mat[e])])
        pc2 = PolyCollection(np.array(quads) * 1e3, facecolors=cs, edgecolors="#52514e", linewidths=0.2)
        ax[1].add_collection(pc2); ax[1].autoscale(); ax[1].set_aspect("equal")
        ax[1].set_title("Face x = 0, fiada 1" if geo.horizontal else "Corte em y = W/2" if meio
                        else "Face y = 0 (juntas e placas)", fontsize=9, loc="left")
        ax[1].set_xlabel("y (mm)" if geo.horizontal else "x (mm)"); ax[1].set_ylabel("z (mm)"); ax[1].grid(False)
    fig.tight_layout()
    if pasta:
        caminho = os.path.join(pasta, arquivo)
        fig.savefig(caminho, dpi=130)
        return caminho, info
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    return buf.getvalue(), info


def executar(p: Dict, pasta: str, log: Callable[[str], None] = print) -> Dict:
    """Roda o projeto e grava tudo em `pasta`. Devolve o dicionário de resultados."""
    p = completar(p)
    os.makedirs(pasta, exist_ok=True)
    with open(os.path.join(pasta, "projeto.json"), "w", encoding="utf-8") as f:
        json.dump(p, f, ensure_ascii=False, indent=2)
    erros = validar(p)
    if erros:
        raise ValueError(" ".join(erros))
    t0 = time.time()
    S = _Saida(pasta)
    geo = geometria(p)
    M = materiais(p, geo)
    an = p["modelo"]["analise"]
    log(f"mef_masonry {VERSAO}: análise '{an}', ensaio '{p['ensaio']['tipo']}'")
    S.item("Razão área líquida/bruta", razao_area_liquida(geo))
    if geo.grauteamento != "nenhum":
        S.item("Razão área resistente (bloco + graute)/bruta", razao_area_solida(geo))
    if an not in PAREDES:
        S.item("E do material (área líquida)", M.ceramica.E_material / 1e6, "MPa")
    if an == "secao_2d":
        _executar_2d(p, geo, M, S, log)
    elif an == "parede":
        figura_malha(p, pasta)
        S.figuras.append({"arquivo": "malha.png", "titulo": "Seção de fibras e esquema da parede"})
        _executar_parede(p, geo, M, S, log)
    elif an == "sensibilidade":
        figura_malha(p, pasta)
        S.figuras.append({"arquivo": "malha.png", "titulo": "Seção de fibras e esquema da parede (valores nominais)"})
        _executar_sensibilidade(p, geo, S, log)
    else:
        ma = p["malha"]
        m = malha_3d(p, geo)
        log("Malha: " + ml.resumo(m))
        S.item("Nós da malha", m.n_nos)
        S.item(f"Elementos ({ELEMENTOS_3D[m.tipo]})", len(m.cun))
        if m.tipo == "W15" and m.n_nos > (40000 if an == "nao_linear" else 80000):
            try:
                import pypardiso  # noqa: F401
                solver = ""
            except ImportError:
                solver = " e instalar o pypardiso"
            S.avisos.append(f"Malha grande para a W15 ({m.n_nos} nós): o cálculo pode levar horas e muita memória. "
                            f"Considere a malha de 16 × 30 mm{solver}.")
            log("Aviso: " + S.avisos[-1])
        figura_malha(p, pasta)
        S.figuras.append({"arquivo": "malha.png", "titulo": "Malha"})
        if an == "termica":
            _executar_termica(p, geo, M, m, S, log)
        elif an == "linear":
            _executar_linear(p, geo, M, m, S, log)
        else:
            _executar_nao_linear(p, geo, M, m, S, log)
    S.item("Tempo de cálculo", (time.time() - t0) / 60, "min")
    res = S.dict()
    res["projeto"] = p
    with open(os.path.join(pasta, "resultados.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2, default=float)
    log("Concluído.")
    return res


def _termica_3d(p, geo, M, m, S, log, instantes=()):
    f, tf, marcos = programa_termico(p)
    if tf <= 0:
        return None, {}, marcos
    sondas = _sondas_3d(m, geo)
    log(f"Térmica: {tf / 60:.0f} min de ensaio, Δt = {p['malha']['dt']:.0f} s")
    rt = tm.resolver_3d(m, M, _contornos_3d(p, m), tf, dt=p["malha"]["dt"], T_inicial=p["ensaio"]["T_inicial"],
                        sondas=sondas, registrar_cada=60.0, instantes_campo=sorted(set(instantes)))
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 3.8))
    t = rt.t / 60
    ax.plot(t, rt.T_gas, color="#9a9893", ls=(0, (4, 2)), lw=1.5, label="gás do forno")
    for k, (nome, v) in enumerate(rt.sondas.items()):
        ax.plot(t, v, color=CORES[k % len(CORES)], lw=1.8, label=nome)
    ax.set_xlabel("tempo (min)"); ax.set_ylabel("temperatura (°C)"); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    S.figura(fig, "temperaturas.png", "Temperaturas ao longo do ensaio")
    nomes = list(rt.sondas)
    S.tabela("Temperaturas ao longo do ensaio", ["tempo (min)", "gás (°C)"] + [f"{n} (°C)" for n in nomes],
             [[round(a, 2), round(b, 1)] + [round(rt.sondas[n][i], 1) for n in nomes]
              for i, (a, b) in enumerate(zip(t, rt.T_gas))][:: max(1, len(t) // 60)],
             arquivo="temperaturas.csv")
    for n in nomes:
        S.item(f"T final — {n}", rt.sondas[n][-1], "°C")
    return rt, sondas, marcos


def _executar_termica(p, geo, M, m, S, log):
    _, tf, _ = programa_termico(p)
    dt = p["malha"]["dt"]
    inst = [round(tf * k / 12 / dt) * dt for k in range(1, 12)] if tf > 0 else []
    rt, _, _ = _termica_3d(p, geo, M, m, S, log, instantes=inst)
    if rt is None:
        S.avisos.append("Ensaio a temperatura ambiente: não há campo térmico a calcular.")
        return
    S.vtk("temperatura_final.vtk", m, dados_nos={"T_C": rt.T_final}, dados_elem={"material": m.mat.astype(float)})
    est = [{"rotulo": f"t = {t / 60:.0f} min", "campos": {"T": rt.campos[t]}} for t in sorted(rt.campos)]
    est.append({"rotulo": f"t = {tf / 60:.0f} min (fim)", "campos": {"T": rt.T_final}})
    S.pos(m, est, "Térmica 3D")


def _executar_linear(p, geo, M, m, S, log):
    en, mo = p["ensaio"], p["modelo"]
    if mo["calibrar_junta"]:
        log(f"Calibrando E da junta para E_prisma = {mo['E_prisma_alvo']:.0f} MPa")
        Ej, _ = estudos.calibrar_junta(m, M, mo["E_prisma_alvo"], geo, verboso=False)
        S.item("E da junta calibrado", Ej / 1e6, "MPa")
    if en["tipo"] == "iso834":
        n = max(2, int(en["duracao_iso"] // 30))
        inst = [en["duracao_iso"] * 60 * (k + 1) / n for k in range(n)]
        rt, _, _ = _termica_3d(p, geo, M, m, S, log, instantes=inst)
        linhas, estados_iso = [], []
        for t in sorted(rt.campos):
            R = mc.resolver_3d(m, M, T_nos=rt.campos[t], carga_N=en["carga"] * 1e3)
            Rt = mc.combinar(R["termico"], R["carga"]) if en["carga"] > 0 else R["termico"]
            s, p99, _ = mc.tracao_maxima(m, Rt)
            d = mc.encurtamento_entre_placas(m, Rt, geo.altura_prisma) * 1e3
            linhas.append([round(t / 60, 1), round(s / 1e6, 3), round(p99 / 1e6, 3), round(d, 3)])
            estados_iso.append({"rotulo": f"t = {t / 60:.0f} min", "u": Rt.u,
                                "campos": {"T": rt.campos[t], **visual.campos_u(Rt.u), "s1": Rt.princ[:, 2] / 1e6,
                                           "s3": Rt.princ[:, 0] / 1e6, "szz": Rt.sig[:, 2] / 1e6}})
            log(f"  t = {t / 60:.0f} min: tração máx {s / 1e6:.2f} MPa")
        S.tabela("Tensões e deformação ao longo do ISO 834 (termoelástico)",
                 ["tempo (min)", "tração máxima (MPa)", "tração p99 (MPa)", "variação de altura (mm)"], linhas,
                 arquivo="iso834_linear.csv")
        S.vtk("campo_final.vtk", m, dados_nos={"T_C": rt.T_final, "u_mm": Rt.u * 1e3},
              dados_elem={"material": m.mat.astype(float), "s1_MPa": Rt.princ[:, 2] / 1e6})
        S.pos(m, estados_iso, "Termoelástico, ISO 834")
        return
    rt, sondas, marcos = _termica_3d(p, geo, M, m, S, log)
    T = rt.T_final if rt is not None else None
    log("Termoelástico: caso térmico e de carga")
    R = mc.resolver_3d(m, M, T_nos=T, carga_N=en["carga"] * 1e3, incluir_termico=T is not None)
    Rc = R["carga"]
    Rt = R.get("termico")
    Rtot = mc.combinar(Rt, Rc) if Rt is not None else Rc
    H = geo.altura_prisma
    A = geo.comprimento * geo.largura
    eps, _, _ = mc.deformacao_extensometro(m, Rc, H / 3, 2 * H / 3)
    S.item("Módulo do corpo de prova, terço central", en["carga"] * 1e3 / A / -eps / 1e6, "MPa")
    S.item("Encurtamento pela carga entre placas", -mc.encurtamento_entre_placas(m, Rc, H) * 1e3, "mm")
    linhas = []
    for nome, Rr in (("térmico", Rt), ("carga", Rc), ("total", Rtot)):
        if Rr is None:
            continue
        s, p99, pos = mc.tracao_maxima(m, Rr)
        x, y, z = (float(v) for v in np.round(pos * 1e3, 1))
        linhas.append([nome, round(s / 1e6, 3), round(p99 / 1e6, 3), x, y, z])
    S.tabela("Tração principal máxima na cerâmica",
             ["caso", "máxima (MPa)", "percentil 99 (MPa)", "x (mm)", "y (mm)", "z (mm)"],
             linhas, arquivo="tracao_maxima.csv")
    if Rt is not None:
        S.item("Dilatação térmica entre placas", mc.encurtamento_entre_placas(m, Rt, H) * 1e3, "mm")
    dn = {"u_total_mm": Rtot.u * 1e3}
    if T is not None:
        dn["T_C"] = T
    S.vtk("tensoes.vtk", m, dados_nos=dn, dados_elem={"material": m.mat.astype(float),
                                                      "s1_total_MPa": Rtot.princ[:, 2] / 1e6,
                                                      "s3_total_MPa": Rtot.princ[:, 0] / 1e6})
    est = []
    for nome, Rr in (("só térmico", Rt), ("só carga", Rc), ("térmico + carga", Rtot)):
        if Rr is None:
            continue
        c = {**visual.campos_u(Rr.u), "s1": Rr.princ[:, 2] / 1e6, "s3": Rr.princ[:, 0] / 1e6, "szz": Rr.sig[:, 2] / 1e6}
        if T is not None:
            c["T"] = T
        est.append({"rotulo": nome, "campos": c, "u": Rr.u})
    S.pos(m, est, "Termoelástico 3D")


def _campos_nl(p, geo, M, m, S, log):
    """Campos de temperatura para o aquecimento do modelo não linear."""
    en, ma = p["ensaio"], p["malha"]
    f, tf, marcos = programa_termico(p)
    if tf <= 0:
        return None
    if en["tipo"] in ("quente", "transiente"):
        tr = marcos["fim da rampa"]
        inst = list(np.linspace(0, tr, int(ma["n_campos_rampa"]) + 1)[1:]) + \
            list(np.linspace(tr, tf, int(ma["n_campos_patamar"]) + 1)[1:])
    else:  # residual
        t1, t2, t3 = marcos["fim da rampa"], marcos["fim do patamar"], marcos["fim do resfriamento"]
        nr, npat = int(ma["n_campos_rampa"]), int(ma["n_campos_patamar"])
        inst = list(np.linspace(0, t1, nr + 1)[1:]) + list(np.linspace(t1, t2, npat + 1)[1:]) + \
            list(np.linspace(t2, t3, nr + 1)[1:]) + list(np.linspace(t3, tf, npat + 1)[1:])
    dt = ma["dt"]
    inst = sorted({round(t / dt) * dt for t in inst})
    rt, _, _ = _termica_3d(p, geo, M, m, S, log, instantes=inst)
    return rt.campos


def _coef_fluencia(p):
    """(ligada, k_bloco, k_argamassa, k_graute, tabela_bloco [(T, φ)]) da fluência transitória."""
    mo, g = p["modelo"], p["geometria"]
    tab = [(float(T), float(v) * 1e-3) for T, v in (mo.get("tabela_fluencia_bloco") or [])]
    kb = float(mo.get("k_tr_bloco", 0.0))
    km = float(mo.get("k_tr_argamassa", 2.0)) if g.get("tipo_junta") == "argamassa" else 0.0
    kg = float(mo.get("k_tr_graute", 0.0)) if g.get("grauteamento", "nenhum") != "nenhum" else 0.0
    ligada = bool(mo.get("fluencia", True)) and (kb > 0 or km > 0 or kg > 0 or bool(tab))
    return ligada, kb, km, kg, tab


def _ligar_fluencia(mod, p, S=None):
    ligada, kb, km, kg, tab = _coef_fluencia(p)
    if ligada:
        mod.definir_fluencia(kb, km, kg, tab or None)
        if S is not None:
            S.item("Fluência transitória: k_tr do bloco, da argamassa e do graute",
                   f"{'tabela' if tab else f'{kb:g}'}; {km:g}; {kg:g}")
    elif S is not None:
        S.item("Fluência transitória", "desligada" if not p["modelo"].get("fluencia", True) else "nula (k_tr = 0)")
    return ligada


def _executar_nao_linear(p, geo, M, m, S, log):
    from .naolinear import ModeloNaoLinear
    en, ma = p["ensaio"], p["malha"]
    Pd, Pj = dano(p), junta_atrito(p)
    H = geo.altura_prisma
    if en["tipo"] == "iso834":
        return _nao_linear_iso(p, geo, M, m, S, log, Pd, Pj)
    if en["tipo"] == "transiente":
        return _nao_linear_transiente(p, geo, M, m, S, log, Pd, Pj)
    campos = _campos_nl(p, geo, M, m, S, log)
    n_passos = max(2, int(round(ma["delta_max"] / ma["incremento"])))
    mod = ModeloNaoLinear(m, M, Pd, Pj, verboso=True, dano_argamassa=dano_arg(p), dano_graute=dano_gra(p))
    mod.tol_u, mod.tol_d = ma["tol_u"], ma["tol_d"]
    mod.instantaneos = []
    mod.definir_extensometro(H / 3, 2 * H / 3)
    _ligar_fluencia(mod, p, S)
    if campos is None:
        campos = {0.0: np.full(m.n_nos, M.ceramica.T_ref)}
    log(f"Aquecimento: {len(campos)} campo(s) sob {en['carga_assentamento']:.1f} kN")
    mod.aquecer(campos)
    vc = mod.vol[mod.ic]
    frac = float(vc[mod.d_aquecimento > 0.5].sum() / vc.sum())
    log(f"Carregamento até {ma['delta_max']} mm em incrementos de {ma['incremento']} mm")
    d, P = mod.carregar(ma["delta_max"] / 1e3, n_passos)
    i = int(np.argmax(P))
    hist = [q for q in mod.historico if q.fase == "carga"]
    S.item("Carga máxima", P[i] / 1e3, "kN")
    S.item("Resistência à compressão (área bruta)", P[i] / (geo.comprimento * geo.largura) / 1e6, "MPa")
    S.item("Deslocamento do prato no pico", d[i] * 1e3, "mm")
    S.item("Cerâmica com dano > 0,5 após o aquecimento", 100 * frac, "%")
    if en["carga"] > 0:
        S.item("Módulo de corda 5–30% da carga de referência", estudos.modulo_corda(mod, en["carga"] * 1e3, geo), "MPa")
    nc = sum(1 for q in hist if not q.convergiu)
    S.item("Incrementos aceitos sem convergência completa", nc)
    if nc:
        S.avisos.append(f"{nc} incremento(s) aceitos no passo mínimo sem convergência completa (ver log).")
    if P[-1] > 0.95 * P[i] and i == len(P) - 1:
        S.avisos.append("A carga ainda subia no fim: aumente o deslocamento máximo para captar o pico.")
    S.tabela("Curva carga × deslocamento",
             ["deslocamento (mm)", "carga (kN)", "deformação no terço central (‰)", "contato nas juntas (%)",
              "cerâmica com dano > 0,5 (%)"],
             [[round(q.rotulo * 1e3, 4), round(q.carga_N / 1e3, 3), round(-q.eps_ext * 1e3, 4),
               round(100 * q.contato, 1), round(100 * q.vol_danificado, 2)] for q in hist],
             arquivo="carga_deslocamento.csv")
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(d * 1e3, P / 1e3, color=CORES[0], lw=2)
    ax.plot(d[i] * 1e3, P[i] / 1e3, "o", color=CORES[0], ms=7, mec="white", mew=1.5)
    ax.annotate(f"{P[i] / 1e3:.1f} kN", (d[i] * 1e3, P[i] / 1e3), xytext=(6, 6), textcoords="offset points", fontsize=8)
    ax.set_xlabel("deslocamento do prato (mm)"); ax.set_ylabel("carga (kN)"); ax.set_xlim(0, None); ax.set_ylim(0, None)
    fig.tight_layout()
    S.figura(fig, "carga_deslocamento.png", "Curva carga × deslocamento")
    folga = np.full(len(m.cun), np.nan); folga[mod.id_glob[mod.ij]] = mod.folga * 1e3
    d_aq = np.full(len(m.cun), np.nan); d_aq[mod.id_glob[mod.ic]] = mod.d_aquecimento
    S.vtk("dano.vtk", m, dados_nos={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
          dados_elem={"material": m.mat.astype(float), "dano_efetivo": mod.dano_global("d"),
                      "dano_tracao": mod.dano_global("dt"), "dano_compressao": mod.dano_global("dc"),
                      "dano_apos_aquecimento": d_aq, "folga_inicial_mm": folga})
    est = visual.estados_nao_linear(mod)
    if est:
        est[0]["campos"]["folga"] = folga
    S.pos(m, est, "Não linear 3D")


def _nao_linear_iso(p, geo, M, m, S, log, Pd, Pj):
    """ISO 834 numa face sob carga constante; ruptura pelos limites de
    deformação axial (EN 1363-1: C = h/100 e dC/dt = 3h/1000 por minuto,
    h em mm) ou pela perda de equilíbrio do modelo."""
    from .naolinear import ModeloNaoLinear
    en, ma = p["ensaio"], p["malha"]
    H = geo.altura_prisma
    dtc = ma["intervalo_iso"] * 60
    inst = list(np.arange(dtc, en["duracao_iso"] * 60 + 1e-6, dtc))
    rt, _, _ = _termica_3d(p, geo, M, m, S, log, instantes=inst)
    mod = ModeloNaoLinear(m, M, Pd, Pj, verboso=True, dano_argamassa=dano_arg(p), dano_graute=dano_gra(p))
    mod.tol_u, mod.tol_d = ma["tol_u"], ma["tol_d"]
    mod.instantaneos = []
    _ligar_fluencia(mod, p, S)
    P = en["carga"] * 1e3
    log(f"Aplicando a carga constante de {en['carga']:.1f} kN")
    mod.aquecer({0.0: np.full(m.n_nos, M.ceramica.T_ref)}, pre_carga_N=P)
    u0 = float(mod.u[3 * mod.topo + 2].mean())
    hmm = H * 1e3
    C_lim, v_lim = hmm / 100, 3 * hmm / 1000
    linhas, t_rup, ant = [[0.0, 0.0, 0.0, 0.0]], None, (0.0, 0.0)
    for t in sorted(rt.campos):
        mod.aquecer({t: rt.campos[t]}, pre_carga_N=P)
        C = -(float(mod.u[3 * mod.topo + 2].mean()) - u0) * 1e3          # contração positiva, mm
        v = (C - ant[1]) / ((t - ant[0]) / 60)
        frac = float(mod.vol[mod.ic][np.maximum(mod.dt, mod.dc) > 0.5].sum() / mod.vol[mod.ic].sum())
        linhas.append([round(t / 60, 1), round(C, 4), round(v, 4), round(100 * frac, 2)])
        log(f"  t = {t / 60:.0f} min: contração {C:.3f} mm, taxa {v:.3f} mm/min")
        if t_rup is None and C >= C_lim and v >= v_lim:
            t_rup = t / 60
            break
        if not mod.historico[-1].convergiu and abs(C) > 10 * C_lim:
            t_rup = t / 60
            S.avisos.append("O modelo perdeu o equilíbrio (deslocamentos sem limite).")
            break
        ant = (t, C)
    S.item("Limite de contração axial C = h/100", C_lim, "mm")
    S.item("Limite de taxa de contração 3h/1000", v_lim, "mm/min")
    S.item("Tempo até a perda da capacidade portante", t_rup if t_rup else f"> {en['duracao_iso']:.0f}", "min")
    if mod.fluencia is not None:
        S.item("Deformação transitória média no fim (vertical)", -mod.deformacao_transitoria_media * 1e3, "‰")
    S.tabela("Contração axial ao longo do ISO 834 (positiva = encurtamento)",
             ["tempo (min)", "contração (mm)", "taxa (mm/min)", "cerâmica com dano > 0,5 (%)"], linhas,
             arquivo="iso834_contracao.csv")
    plt = _plt()
    a = np.array(linhas, float)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(a[:, 0], a[:, 1], color=CORES[0], lw=2, label="contração do modelo")
    ax.axhline(C_lim, color="#9a9893", ls=(0, (4, 2)), lw=1, label="limite h/100")
    ax.set_xlabel("tempo de exposição (min)"); ax.set_ylabel("contração axial (mm)"); ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    S.figura(fig, "iso834_contracao.png", "Contração axial sob carga constante")
    S.avisos.append("Critério de ruptura pelos limites de deformação da EN 1363-1 para elementos verticais; "
                    "sem tabelas de E(T), ft(T) e fc(T), o modelo tende a não romper.")
    S.vtk("dano.vtk", m, dados_nos={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
          dados_elem={"material": m.mat.astype(float), "dano_tracao": mod.dano_global("dt"),
                      "dano_compressao": mod.dano_global("dc")})
    S.pos(m, visual.estados_nao_linear(mod), "Não linear 3D, ISO 834")


def _nao_linear_transiente(p, geo, M, m, S, log, Pd, Pj):
    """Ensaio transiente: carga constante aplicada a frio e aquecimento no forno
    (rampa e patamar) até a ruptura. A ruptura é a perda de equilíbrio do modelo
    ou os limites de deformação axial da EN 1363-1 (C = h/100 e dC/dt = 3h/1000
    por minuto). A curva deformação × temperatura é a saída clássica do ensaio:
    a diferença para a dilatação livre é a deformação induzida pela carga."""
    from .naolinear import ModeloNaoLinear
    en, ma = p["ensaio"], p["malha"]
    H = geo.altura_prisma
    A = geo.comprimento * geo.largura
    campos = _campos_nl(p, geo, M, m, S, log)
    f_gas, tf, _ = programa_termico(p)
    mod = ModeloNaoLinear(m, M, Pd, Pj, verboso=True, dano_argamassa=dano_arg(p), dano_graute=dano_gra(p))
    mod.tol_u, mod.tol_d = ma["tol_u"], ma["tol_d"]
    mod.instantaneos = []
    mod.definir_extensometro(H / 3, 2 * H / 3)
    _ligar_fluencia(mod, p, S)
    P = en["carga"] * 1e3
    S.item("Carga constante", en["carga"], "kN")
    S.item("Tensão na área bruta", P / A / 1e6, "MPa")
    log(f"Ensaio transiente: {en['carga']:.1f} kN a frio e aquecimento a {en['taxa']:g} °C/min até "
        f"{en['T_alvo']:.0f} °C")
    mod.aquecer({0.0: np.full(m.n_nos, M.ceramica.T_ref)}, pre_carga_N=P)
    u0 = float(mod.u[3 * mod.topo + 2].mean())
    e0 = mod.deformacao_extensometro()
    vc = mod.vol[mod.ic]
    bl = ~(mod.arg | mod.gra)

    def T_bloco():
        Tc = mod.T[mod.cun[mod.ic]].mean(1)
        return float((Tc * vc)[bl].sum() / vc[bl].sum())

    hmm = H * 1e3
    C_lim, v_lim = hmm / 100, 3 * hmm / 1000
    Ti = en["T_inicial"]
    linhas = [[0.0, round(Ti, 1), round(T_bloco(), 1), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    rup, ant = None, (0.0, 0.0)
    for t in sorted(campos):
        mod.aquecer({t: campos[t]}, pre_carga_N=P)
        conv = mod.historico[-1].convergiu
        Tb = T_bloco()
        C = -(float(mod.u[3 * mod.topo + 2].mean()) - u0) * 1e3            # contração positiva, mm
        v = (C - ant[1]) / ((t - ant[0]) / 60)
        eps = -(mod.deformacao_extensometro() - e0) * 1e3                   # encurtamento positivo, ‰
        livre = -float(M.ceramica.eps_termica(np.array([Tb]))[0]) * 1e3
        tr = -mod.deformacao_transitoria_media * 1e3
        frac = float(vc[np.maximum(mod.dt, mod.dc) > 0.5].sum() / vc.sum())
        linhas.append([round(t / 60, 2), round(float(f_gas(t)), 1), round(Tb, 1), round(eps, 4), round(livre, 4),
                       round(eps - livre, 4), round(tr, 4), round(C, 4), round(100 * frac, 2)])
        log(f"  t = {t / 60:.0f} min: bloco {Tb:.0f} °C, deformação {eps:+.3f} ‰ (livre {livre:+.3f} ‰), "
            f"contração {C:.3f} mm")
        if (C >= C_lim and v >= v_lim) or (not conv and abs(C) > C_lim):
            rup = (t, Tb, float(f_gas(t)))
            break
        ant = (t, C)
    if rup:
        S.item("Tempo até a ruptura", rup[0] / 60, "min")
        S.item("Temperatura do forno na ruptura", rup[2], "°C")
        S.item("Temperatura média do bloco na ruptura", rup[1], "°C")
    else:
        S.item("Tempo até a ruptura", f"> {tf / 60:.0f} (sem ruptura)", "min")
    if mod.fluencia is not None:
        S.item("Deformação transitória média no fim (vertical)", linhas[-1][6], "‰")
    S.tabela("Ensaio transiente (deformações de encurtamento positivas, no terço central)",
             ["tempo (min)", "forno (°C)", "bloco, média (°C)", "deformação total (‰)", "livre (‰)",
              "induzida pela carga (‰)", "fluência transitória (‰)", "contração do topo (mm)",
              "bloco com dano > 0,5 (%)"], linhas, arquivo="ensaio_transiente.csv")
    plt = _plt()
    a = np.array(linhas, float)
    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.plot(a[:, 2], a[:, 3], color=CORES[0], lw=2, label=f"total, sob {en['carga']:g} kN")
    ax.plot(a[:, 2], a[:, 4], color="#9a9893", lw=1.5, ls=(0, (4, 2)), label="dilatação livre")
    if mod.fluencia is not None:
        ax.plot(a[:, 2], a[:, 6], color=CORES[1], lw=1.5, label="fluência transitória")
    ax.axhline(0, color="#d6d4cf", lw=0.8)
    ax.set_xlabel("temperatura média do bloco (°C)"); ax.set_ylabel("encurtamento (‰)")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    S.figura(fig, "ensaio_transiente.png", "Deformação × temperatura no ensaio transiente")
    S.avisos.append("Deformação no terço central da face y = 0; encurtamento positivo e dilatação negativa.")
    S.vtk("dano.vtk", m, dados_nos={"T_C": mod.T, "u_mm": mod.u.reshape(-1, 3) * 1e3},
          dados_elem={"material": m.mat.astype(float), "dano_tracao": mod.dano_global("dt"),
                      "dano_compressao": mod.dano_global("dc")})
    S.pos(m, visual.estados_nao_linear(mod), "Não linear 3D, ensaio transiente")


def _executar_2d(p, geo, M, S, log):
    from .malha import FURO
    en, ma = p["ensaio"], p["malha"]
    m2 = ml.malha_secao(geo, h=ma["h"] / 1e3, simetria=False, faixas=())
    mat2 = tm.material_secao(m2, geo)
    log("Malha 2D: " + ml.resumo(m2))
    S.item("Nós da malha", m2.n_nos)
    S.item("Elementos (triângulos T3)", len(m2.tri))
    f, tf, _ = programa_termico(p)
    face = None
    if en["tipo"] == "iso834":
        exp = tm.comprimentos_contorno_2d(m2, faces=("y0",))
        nexp = tm.comprimentos_contorno_2d(m2, faces=("y1",))
        ctr = [Contorno(exp, f, en["h_conveccao"], en["emissividade"], "exposta"),
               Contorno(nexp, curvas.constante(en["T_inicial"]), 4.0, en["emissividade"], "não exposta")]
        face = np.nonzero(nexp > 0)[0]
        w = nexp[face] / nexp[face].sum()
        Ti = en["T_inicial"]
        sondas = {"face não exposta (média)": face}

        def parar(t, T):
            return (T[face] @ w - Ti >= 140.0) or (T[face].max() - Ti >= 180.0)
    else:
        ctr = [Contorno(tm.comprimentos_contorno_2d(m2), f, en["h_conveccao"], en["emissividade"])]
        cen = m2.xy[m2.tri].mean(1)
        sol = np.unique(m2.tri[mat2 == CERAMICA])
        sondas = {"mínima da cerâmica": sol[np.argsort(np.linalg.norm(m2.xy[sol] - [geo.comprimento / 2, geo.largura / 2], axis=1))[:1]]}
        parar = None
    log(f"Térmica 2D: {tf / 60:.0f} min")
    r = tm.resolver_2d(m2, M, ctr, tf, dt=ma["dt"], T_inicial=en["T_inicial"], sondas=sondas, registrar_cada=60.0, mat=mat2,
                       parar_se=parar)
    if en["tipo"] == "iso834":
        S.item("Tempo de isolamento (140/180 °C na face não exposta)",
               r.parou_em / 60 if r.parou_em else f"> {tf / 60:.0f}", "min")
    s = mc.resolver_secao_gps(m2, M, r.T_final, mat=mat2)
    sol = mat2 != AR
    S.item("Tração máxima (estado plano generalizado)", np.nanmax(s["s1"]) / 1e6, "MPa")
    S.item("Tração p99", np.nanpercentile(s["s1"][sol], 99) / 1e6, "MPa")
    from .saida import figura_secao
    Te = r.T_final[m2.tri].mean(1)
    figura_secao(m2, Te, "Temperatura no fim (°C)", os.path.join(S.pasta, "secao_T.png"), "°C")
    S.figuras.append({"arquivo": "secao_T.png", "titulo": "Temperatura na seção no fim"})
    figura_secao(m2, s["s1"] / 1e6, "Tensão principal máxima (MPa)", os.path.join(S.pasta, "secao_s1.png"), "MPa",
                 mascara=sol, cmap="RdBu_r")
    S.figuras.append({"arquivo": "secao_s1.png", "titulo": "Tensão principal máxima na seção"})
    nomes = list(r.sondas)
    S.tabela("Temperaturas ao longo do ensaio", ["tempo (min)", "gás (°C)"] + [f"{n} (°C)" for n in nomes],
             [[round(a / 60, 2), round(b, 1)] + [round(float(np.mean(r.sondas[n][i])), 1) for n in nomes]
              for i, (a, b) in enumerate(zip(r.t, r.T_gas))], arquivo="temperaturas_2d.csv")
    S.vtk("secao.vtk", m2, dados_nos={"T_C": r.T_final}, dados_elem={"s1_MPa": s["s1"] / 1e6})
    S.pos(m2, [{"rotulo": f"fim, t = {r.t[-1] / 60:.0f} min", "campos": {"T": r.T_final, "s1": s["s1"] / 1e6}}],
          "Seção 2D", material=mat2)


# ----------------------------------------------------------------------------
# Parede (faixa vertical com segunda ordem)
# ----------------------------------------------------------------------------
def _secao_parede(p, geo):
    """Malha 2D de um bloco de comprimento (período da seção) e as fibras."""
    from .parede import SecaoFibras
    geo1 = replace(geo, blocos_por_fiada=1.0)
    m2 = ml.malha_secao(geo1, h=p["malha"]["h"] / 1e3, simetria=False, faixas=())
    return geo1, m2, SecaoFibras.da_malha(m2, geo1)


def _n_elementos(p):
    n = max(4, int(p["malha"]["n_elementos"]))
    return n + (n % 2)


def _figura_parede(p, geo, pasta=None, arquivo="malha.png"):
    _plt()
    from matplotlib.figure import Figure
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Rectangle, FancyArrowPatch
    en = p["ensaio"]
    geo1, m2, sec = _secao_parede(p, geo)
    H, t, L = geo.altura_prisma, geo.largura, geo.comprimento
    n = _n_elementos(p)
    info = (f"Parede {L * 1e3:.0f} × {H * 1e3:.0f} × {t * 1e3:.0f} mm (H/t = {H / t:.1f}); "
            + ("seção vertical de uma fiada" if m2.vertical else "seção de um bloco")
            + f" com {len(sec.y)} fibras; {n} elementos de viga")
    fig = Figure(figsize=(10, 3.8))
    ax = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1.3, 1] if m2.vertical else [2.4, 1]})
    if m2.vertical:
        # espessura na horizontal (face exposta à esquerda), altura da fiada na vertical
        fib = np.zeros(len(m2.tri), bool); fib[sec.tri] = True
        mat2 = tm.material_secao(m2, geo1)
        cor = np.where(mat2 == AR, "#cde2fb", np.where(~fib, "#efe6da",
                       np.where(mat2 == CERAMICA, "#d9c7b0", "#7f8fa6")))
        pc = PolyCollection(m2.xy[m2.tri][:, :, ::-1] * 1e3, facecolors=cor, edgecolors="#52514e", linewidths=0.25)
        ax[0].add_collection(pc); ax[0].autoscale(); ax[0].set_aspect("equal"); ax[0].grid(False)
        ax[0].axvline(0, color="#eb6834", lw=3)
        ax[0].set_title("Seção vertical de uma fiada (fibras: paredes contínuas)", fontsize=9, loc="left")
        ax[0].set_xlabel("y, espessura (mm; y = 0 é a face " + ("exposta)" if en["tipo"] == "iso834" else "y = 0)"))
        ax[0].set_ylabel("altura na fiada (mm)")
    else:
        _, gr0 = ml.estados_fiada(m2, geo1, 0)
        cor = np.where(gr0, "#5b5f66", np.where(m2.marca == 0, "#cde2fb", "#d9c7b0"))
        pc = PolyCollection(m2.xy[m2.tri] * 1e3, facecolors=cor, edgecolors="#52514e", linewidths=0.25)
        ax[0].add_collection(pc); ax[0].autoscale(); ax[0].set_aspect("equal"); ax[0].grid(False)
        ax[0].axhline(0, color="#eb6834", lw=3)
        ax[0].text(m2.comprimento * 1e3 / 2, -6, "face exposta (ISO 834)" if en["tipo"] == "iso834" else "face y = 0",
                   color="#eb6834", ha="center", va="top", fontsize=8)
        ax[0].set_ylim(-18, t * 1e3 + 4)
        ax[0].set_title("Seção de fibras: um bloco de comprimento (bloco em bege, furos em azul"
                        + (", graute em cinza)" if gr0.any() else ")"), fontsize=9, loc="left")
        ax[0].set_xlabel("x (mm)"); ax[0].set_ylabel("y (mm)")
    # elevação esquemática: y horizontal (exposto à esquerda), altura vertical
    a = ax[1]; a.grid(False)
    imp = H / en["imperfeicao"] if en["imperfeicao"] > 0 else 0.0
    xs = np.linspace(0, H, n + 1)
    esc = 0.25 * t / max(imp, 1e-9) if imp > 0 else 0
    ys = -imp * np.sin(np.pi * xs / H) * esc
    a.add_patch(Rectangle((-t / 2 * 1e3, 0), t * 1e3, H * 1e3, facecolor="#efe6da", edgecolor="#9a9893", lw=0.8))
    a.plot(ys * 1e3, xs * 1e3, "-o", color="#2a78d6", ms=2.5, lw=1.2)
    if en["tipo"] == "iso834":
        a.axvline(-t / 2 * 1e3, color="#eb6834", lw=3)
    for z, apo, ex in ((0.0, en["apoio_base"], en["excentricidade_base"]), (H, en["apoio_topo"], en["excentricidade_topo"])):
        yy = -ex
        if apo == "engastado":
            a.plot([-t / 2 * 1e3 * 1.6, t / 2 * 1e3 * 1.6], [z * 1e3] * 2, color="#52514e", lw=3)
        else:
            a.plot(yy, z * 1e3, marker="^" if z == 0 else "v", color="#52514e", ms=9)
    yt = -en["excentricidade_topo"]
    a.add_patch(FancyArrowPatch((yt, H * 1e3 * 1.16), (yt, H * 1e3 * 1.02), arrowstyle="-|>", mutation_scale=12,
                                color="#eb6834" if en["tipo"] == "iso834" else "#52514e", lw=1.5))
    a.text(yt, H * 1e3 * 1.18, "P", ha="center", va="bottom", fontsize=9)
    a.set_xlim(-t * 1e3 * 1.4, t * 1e3 * 1.4); a.set_ylim(-0.06 * H * 1e3, 1.26 * H * 1e3)
    a.set_title("Faixa vertical (eixo e imperfeição ampliados)", fontsize=9, loc="left")
    a.set_xlabel("y (mm)"); a.set_ylabel("altura (mm)")
    fig.tight_layout()
    if pasta:
        caminho = os.path.join(pasta, arquivo)
        fig.savefig(caminho, dpi=130)
        return caminho, info
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    return buf.getvalue(), info


def _curva_prisma_30(p, log):
    """Roda o modelo não linear do prisma (5 fiadas de bloco e meio) a 30 °C com
    as entradas atuais; devolve (delta mm, carga kN), H do prisma e área líquida."""
    from .naolinear import ModeloNaoLinear
    q = completar(copy.deepcopy(p))
    q["geometria"].update({"blocos_por_fiada": 1.5, "n_fiadas": 5, "amarracao": "contrafiada"})
    q["modelo"].update({"analise": "nao_linear", "elemento": "W6"})
    q["ensaio"].update({"tipo": "ambiente"})
    q["malha"].update({"h": 12, "dz": 20, "delta_max": 0.8, "incremento": 0.02})
    geo = geometria(q)
    M = materiais(q, geo)
    m = malha_3d(q, geo)
    log("Calibração pelo prisma a 30 °C: " + ml.resumo(m))
    mod = ModeloNaoLinear(m, M, dano(q), junta_atrito(q), verboso=False, dano_argamassa=dano_arg(q), dano_graute=dano_gra(q))
    mod.aquecer({0.0: np.full(m.n_nos, M.ceramica.T_ref)})
    if geo.tipo_junta == "argamassa":       # prisma argamassado: pico mais tardio
        d, P = mod.carregar(0.006 * geo.altura_prisma, 60)
    else:
        d, P = mod.carregar(0.8e-3, 40)
    i = int(np.argmax(P))
    log(f"  prisma: pico {P[i] / 1e3:.1f} kN com {d[i] * 1e3:.3f} mm")
    A_liq = geo.comprimento * geo.largura * razao_area_solida(geo)
    return np.c_[d * 1e3, P / 1e3], geo.altura_prisma, A_liq, P[i] / (geo.comprimento * geo.largura)


def _lei_parede(p, geo, log, S):
    from .parede import LeiParede
    from .dados_parede import CURVAS_PRISMA_BCEPQ, H_REF, A_BRUTA_REF
    ref = LeiParede.de_prisma(CURVAS_PRISMA_BCEPQ, H_REF, A_BRUTA_REF * razao_area_liquida(Geometria()))
    if p["modelo"]["lei_parede"] != "prisma":
        return ref
    curva, H_p, A_liq, f_bruta = _curva_prisma_30(p, log)
    nova = LeiParede.de_prisma({30.0: curva}, H_p, A_liq)
    k_s = nova.s_pico[0] / ref.s_pico[0]
    k_e = nova.e_pico[0] / ref.e_pico[0]
    S.item("Resistência do prisma calibrado a 30 °C (área bruta)", f_bruta / 1e6, "MPa")
    log(f"  razões em relação ao BCEPQ: tensão {k_s:.3f}, deformação {k_e:.3f}")
    return ref.escalada(k_s, k_e).com_curva(30.0, nova.e[0], nova.s[0])


def _fluencia_parede(p, geo, M, sec, S=None):
    """c(T) por fibra (1/Pa) da fluência transitória na parede, ou None.

    Bloco e junta de leito estão em série na altura da fiada: c = (1 - f) φ_b/fc_b
    + f φ_a/fc_a, com f = espessura da junta / altura da fiada (só com argamassa).
    O graute é contínuo: c = φ_g/fc_g."""
    ligada, kb, km, kg, tab = _coef_fluencia(p)
    if S is not None:
        S.item("Fluência transitória: k_tr do bloco, da argamassa e do graute",
               f"{'tabela' if tab else f'{kb:g}'}; {km:g}; {kg:g}" if ligada else
               ("desligada" if not p["modelo"].get("fluencia", True) else "nula (k_tr = 0)"))
    if not ligada:
        return None
    c, mo = p["ceramica"], p["modelo"]
    fc_b = mo["fc"] * 1e6
    fc_a = c.get("arg_fc", 4.2) * 1e6
    fc_g = c.get("gra_fc", 17.7) * 1e6
    f_a = geo.espessura_junta / geo.altura_fiada if geo.tipo_junta == "argamassa" else 0.0
    gr = np.asarray(sec.graute, bool)

    def phi_bloco(T):
        if tab:
            Ts, vs = zip(*sorted(tab))
            return np.interp(T, Ts, vs)
        return kb * np.maximum(M.ceramica.eps_termica(T), 0.0)

    def coef(T):
        T = np.asarray(T, float)
        cb = (1 - f_a) * phi_bloco(T) / fc_b + f_a * km * np.maximum(M.argamassa.eps_termica(T), 0.0) / fc_a
        cg = kg * np.maximum(M.graute.eps_termica(T), 0.0) / fc_g
        return np.where(gr, cg, cb)
    return coef


def _executar_parede(p, geo, M, S, log, lei=None, rapido=False):
    """Parede no ISO 834 ou a frio. `lei`: lei da parede já pronta (sem os fatores
    do projeto); `rapido`: sem figuras nem níveis de carga, devolve R, I e Pmax."""
    from .parede import ModeloParede
    from .saida import figura_secao
    en, ma = p["ensaio"], p["malha"]
    geo1, m2, sec = _secao_parede(p, geo)
    H, t, L = geo.altura_prisma, geo.largura, geo.comprimento
    Ti = en["T_inicial"]
    n = _n_elementos(p)
    log(f"Parede {L * 1e3:.0f} × {H * 1e3:.0f} × {t * 1e3:.0f} mm; seção com {len(sec.y)} fibras; "
        f"{n} elementos de viga")
    S.item("Altura da parede", H * 1e3, "mm")
    S.item("Esbeltez H/t", H / t)
    S.item("Fibras na seção", len(sec.y))
    if lei is None:
        lei = _lei_parede(p, geo, log, S)
    fr = p["modelo"].get("fator_resistencia_parede", 1.0)
    fk = p["modelo"].get("fator_rigidez_parede", 1.0)
    if fr != 1.0 or fk != 1.0:
        lei = lei.escalada(fr, fr / fk)
        S.item("Fatores na resistência e na rigidez da lei", f"{fr:.3g}; {fk:.3g}")
    if not rapido and p["modelo"]["lei_parede"] != "prisma" and (geo.tipo_junta != "seca"
                                                  or p["ceramica"].get("material_biblioteca", "bcepq") != "bcepq"):
        S.avisos.append("A lei da parede usada é a do prisma BCEPQ; para outra alvenaria, escolha 'Calibrar agora "
                        "pelo modelo não linear do prisma' em Modelos e elementos.")
    razao = sec.area / t
    S.item(f"Resistência da alvenaria a {Ti:.0f} °C (área bruta)", lei.resistencia(Ti) * razao / 1e6, "MPa")
    imp = H / en["imperfeicao"] if en["imperfeicao"] > 0 else 0.0
    e_t, e_b = en["excentricidade_topo"] / 1e3, en["excentricidade_base"] / 1e3
    gr_fib, arg_fib = sec.graute, sec.argamassa

    def eth(T):
        """Deformação térmica livre por fibra (bloco, graute ou junta)."""
        return np.where(gr_fib, M.graute.eps_termica(T),
                        np.where(arg_fib, M.argamassa.eps_termica(T), M.ceramica.eps_termica(T)))

    flu = _fluencia_parede(p, geo, M, sec, S)

    def novo():
        return ModeloParede(H, sec, lei, eth, n_el=n, apoio_base=en["apoio_base"], apoio_topo=en["apoio_topo"],
                            imperfeicao=imp, T_inicial=Ti, fluencia=flu)

    if imp == 0 and e_t == 0 and e_b == 0:
        S.avisos.append("Parede perfeita e carga centrada: sem imperfeição, a flambagem não é detectada a frio.")
    # ---------------------------------------------------------------- capacidade a frio
    log(f"Capacidade a {Ti:.0f} °C (carga crescente até perder o equilíbrio)")
    m0 = novo()
    curva = [(0.0, 0.0, 0.0)]
    Pmax = m0.capacidade(e_t, e_b, registrar=lambda P: curva.append((P, m0.flecha_meio, m0.deslocamento_topo)))
    N_sec = lei.resistencia(Ti) * sec.area
    S.item(f"Capacidade da parede a {Ti:.0f} °C", Pmax * L / 1e3, "kN")
    S.item(f"Capacidade por metro a {Ti:.0f} °C", Pmax / 1e3, "kN/m")
    S.item("Fator de redução (capacidade / resistência da seção)", Pmax / N_sec)
    log(f"  capacidade {Pmax * L / 1e3:.1f} kN ({Pmax / 1e3:.1f} kN/m); fator de redução {Pmax / N_sec:.3f}")
    plt = _plt()

    def campo_tensao(mod):
        a = np.full(len(m2.tri), np.nan)
        sg = mod.fibras_meia_altura()
        if sg is not None:
            a[sec.tri] = sg / 1e6
        return a

    if en["tipo"] != "iso834":
        c = np.array(curva)
        S.tabela("Carga × flecha a frio", ["carga (kN)", "carga (kN/m)", "flecha no meio (mm)", "deslocamento do topo (mm)"],
                 [[round(a * L / 1e3, 2), round(a / 1e3, 3), round(b * 1e3, 3), round(d * 1e3, 4)] for a, b, d in c],
                 arquivo="parede_capacidade.csv")
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.8))
        ax[0].plot(-c[:, 1] * 1e3, c[:, 0] * L / 1e3, color=CORES[0], lw=2)
        ax[0].set_xlabel("flecha no meio da altura, para o lado y < 0 (mm)"); ax[0].set_ylabel("carga (kN)")
        ax[0].axhline(Pmax * L / 1e3, color="#9a9893", ls=(0, (4, 2)), lw=1)
        ax[0].annotate(f"{Pmax * L / 1e3:.0f} kN", (0, Pmax * L / 1e3), xytext=(4, 4), textcoords="offset points", fontsize=8)
        ax[1].plot(m0.flecha * 1e3, m0.x0 * 1e3, color=CORES[0], lw=2)
        ax[1].axvline(0, color="#9a9893", lw=0.8)
        ax[1].set_xlabel("deslocamento horizontal no pico (mm)"); ax[1].set_ylabel("altura (mm)")
        fig.tight_layout()
        S.figura(fig, "parede_capacidade.png", "Carga × flecha e deformada no pico")
        S.pos(m2, [{"rotulo": "pico de carga (meia altura)", "campos": {"szz": campo_tensao(m0)}}],
              "Parede: tensão vertical nas fibras na meia altura", material=tm.material_secao(m2, geo1))
        S.avisos.append("Modelo de faixa vertical: as bordas verticais da parede não restringem a flexão. "
                        "Material homogeneizado pelo prisma, sem tração (as juntas abrem).")
        return
    # ---------------------------------------------------------------- térmica 2D
    f, tf, _ = programa_termico(p)
    exp = tm.comprimentos_contorno_2d(m2, faces=("y0",))
    nexp = tm.comprimentos_contorno_2d(m2, faces=("y1",))
    ctr = [Contorno(exp, f, en["h_conveccao"], en["emissividade"], "exposta"),
           Contorno(nexp, curvas.constante(Ti), 4.0, en["emissividade"], "não exposta")]
    face = np.nonzero(nexp > 0)[0]
    w = nexp[face] / nexp[face].sum()
    face_exp = np.nonzero(exp > 0)[0]
    we = exp[face_exp] / exp[face_exp].sum()
    t_iso = []

    def isolamento(tt, T):
        if not t_iso and (T[face] @ w - Ti >= 140.0 or T[face].max() - Ti >= 180.0):
            t_iso.append(tt)
        return False

    inst = list(np.arange(60.0, tf + 1e-6, 60.0))
    log(f"Térmica 2D da seção: {tf / 60:.0f} min de ISO 834, Δt = {ma['dt']:.0f} s")
    r = tm.resolver_2d(m2, M, ctr, tf, dt=ma["dt"], T_inicial=Ti, registrar_cada=60.0, instantes_campo=inst,
                       mat=tm.material_secao(m2, geo1),
                       parar_se=isolamento)
    tempos = [0.0] + sorted(r.campos)
    Tnos = {0.0: np.full(m2.n_nos, Ti), **r.campos}
    campos = [(tt, sec.temperaturas(m2, Tnos[tt])) for tt in tempos]
    S.item("Tempo de isolamento (140/180 °C na face não exposta)", t_iso[0] / 60 if t_iso else f"> {tf / 60:.0f}", "min")
    # ---------------------------------------------------------------- incêndio sob carga
    P = en["carga"] * 1e3 / L
    S.item("Carga aplicada", en["carga"], "kN")
    S.item("Carga aplicada por metro", P / 1e3, "kN/m")
    S.item(f"Carga / capacidade a {Ti:.0f} °C", 100 * P / Pmax, "%")
    passo_foto = max(15.0, round(tf / 60 / 8 / 15) * 15) * 60
    fotos = {}
    hist = []
    mod = novo()

    def registrar(tt):
        hist.append((tt, mod.flecha_meio, mod.flecha_max, mod.deslocamento_topo))
        k = round(tt / 60)
        if abs(tt - 60 * k) < 1e-6 and (tt % passo_foto < 1e-6) and tt > 0:
            fotos[tt] = (mod.flecha.copy(), campo_tensao(mod))
        if abs(tt - 60 * k) < 1e-6 and k % 10 == 0 and k > 0:
            log(f"  t = {k} min: flecha {mod.flecha_max * 1e3:.1f} mm, topo {mod.deslocamento_topo * 1e3:+.2f} mm")

    log(f"ISO 834 sob {en['carga']:.1f} kN ({P / 1e3:.1f} kN/m, {100 * P / Pmax:.0f}% da capacidade a frio)")
    if P >= Pmax or not mod.carregar(P, e_t, e_b):
        t_col, t_ult = 0.0, 0.0
        S.avisos.append("A carga informada é maior que a capacidade da parede a frio.")
    else:
        hist.append((0.0, mod.flecha_meio, mod.flecha_max, mod.deslocamento_topo))
        fotos[0.0] = (mod.flecha.copy(), campo_tensao(mod))
        t_col, t_ult = mod.incendio(P, e_t, e_b, campos, registrar=registrar)
        fotos[t_ult] = (mod.flecha.copy(), campo_tensao(mod))
    R = t_col / 60 if t_col is not None else None
    S.item("Tempo de resistência (perda de equilíbrio)", R if R is not None else f"> {tf / 60:.0f}", "min")
    # critério de deformação axial da EN 1363-1 (h em mm)
    hmm = H * 1e3
    C_lim, v_lim = hmm / 100, 3 * hmm / 1000
    t_C = None
    if hist:
        hm = [(a, d) for a, b, c, d in hist if abs(a / 60 - round(a / 60)) < 1e-6]
        u0 = hm[0][1]
        for (ta, ua), (tb, ub) in zip(hm[:-1], hm[1:]):
            C = -(ub - u0) * 1e3
            v = -(ub - ua) * 1e3 / ((tb - ta) / 60)
            if C >= C_lim and v >= v_lim:
                t_C = tb / 60
                break
    S.item("Critério de contração axial (EN 1363-1)", t_C if t_C else "não atingido", "min" if t_C else "")
    if rapido:
        return {"R": R if R is not None else tf / 60, "R_censurado": R is None,
                "I": t_iso[0] / 60 if t_iso else tf / 60, "I_censurado": not t_iso,
                "Pmax": Pmax * L / 1e3, "flecha": hist[-1][2] * 1e3 if hist else float("nan"),
                "carga_capacidade": 100 * P / Pmax, "tf": tf / 60}
    if hist:
        S.item("Flecha máxima no fim", hist[-1][2] * 1e3, "mm")
        xc, rc = mod.compressao_critica()
        S.item("Altura da fibra mais comprimida no fim", xc * 1e3, "mm")
    if R is not None:
        log(f"  perda de equilíbrio aos {R:.1f} min")
    a = np.array(hist) if hist else np.zeros((1, 4))
    S.tabela("Histórico do ISO 834 sob carga",
             ["tempo (min)", "flecha no meio (mm)", "flecha máxima (mm)", "deslocamento do topo (mm)"],
             [[round(x / 60, 3), round(b * 1e3, 3), round(c * 1e3, 3), round(d * 1e3, 4)] for x, b, c, d in a
              if abs(x / 60 - round(x / 60)) < 1e-6 or x == a[-1, 0]],
             arquivo="parede_historico.csv")
    # ---------------------------------------------------------------- níveis de carga
    niveis = []
    k = int(en["niveis_carga"])
    if k > 0:
        log(f"Curva tempo de resistência × carga: {k} níveis")
        for frac in np.linspace(0.1, 0.9, k):
            mk = novo()
            Pk = frac * Pmax
            if not mk.carregar(Pk, e_t, e_b):
                niveis.append((frac, Pk, 0.0)); continue
            tc, _ = mk.incendio(Pk, e_t, e_b, campos)
            Rk = tc / 60 if tc is not None else np.inf
            niveis.append((frac, Pk, Rk))
            log(f"  {100 * frac:.0f}% ({Pk * L / 1e3:.0f} kN): " + (f"{Rk:.1f} min" if np.isfinite(Rk) else f"> {tf / 60:.0f} min"))
        S.tabela("Tempo de resistência × carga", ["carga / capacidade a frio (%)", "carga (kN)", "carga (kN/m)",
                                                  "tempo de resistência (min)"],
                 [[round(100 * fr, 1), round(Pk * L / 1e3, 1), round(Pk / 1e3, 2),
                   round(Rk, 1) if np.isfinite(Rk) else f"> {tf / 60:.0f}"] for fr, Pk, Rk in niveis],
                 arquivo="parede_tempo_carga.csv")
    # ---------------------------------------------------------------- figuras
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    ax[0].plot(a[:, 0] / 60, -a[:, 2] * 1e3, color=CORES[0], lw=2)
    ax[0].set_xlabel("tempo de exposição (min)"); ax[0].set_ylabel("flecha máxima para o fogo (mm)")
    ax[1].plot(a[:, 0] / 60, a[:, 3] * 1e3, color=CORES[1], lw=2)
    ax[1].axhline(0, color="#9a9893", lw=0.8)
    ax[1].set_xlabel("tempo de exposição (min)"); ax[1].set_ylabel("deslocamento vertical do topo (mm)")
    for x in ax:
        if R is not None:
            x.axvline(R, color="#9a9893", ls=(0, (4, 2)), lw=1)
    if R is not None:
        ax[0].annotate(f"perda de equilíbrio: {R:.0f} min", (R, 0), xytext=(-4, 6), textcoords="offset points",
                       ha="right", fontsize=8)
    fig.tight_layout()
    S.figura(fig, "parede_flecha.png", "Flecha e deslocamento do topo ao longo do ISO 834")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.2))
    tf_ord = sorted(fotos)
    cmap = plt.get_cmap("inferno")
    for i, tt in enumerate(tf_ord):
        cor = cmap(0.1 + 0.75 * i / max(1, len(tf_ord) - 1))
        ax[0].plot(fotos[tt][0] * 1e3, mod.x0 * 1e3, color=cor, lw=1.6, label=f"{tt / 60:.0f} min")
    ax[0].axvline(0, color="#9a9893", lw=0.8)
    ax[0].set_xlabel("deslocamento horizontal (mm; < 0 = para o fogo)"); ax[0].set_ylabel("altura (mm)")
    if ax[0].get_legend_handles_labels()[0]:
        ax[0].legend(fontsize=7, frameon=False)
    # temperatura média através da espessura
    yb = np.linspace(-t / 2, t / 2, 29)
    ib = np.clip(np.digitize(sec.y, yb) - 1, 0, len(yb) - 2)
    yc = 0.5 * (yb[1:] + yb[:-1])
    Ab = np.bincount(ib, sec.A, len(yc))
    for i, tt in enumerate([x for x in tf_ord if x > 0 and x in Tnos] or [tempos[-1]]):
        Tf = sec.temperaturas(m2, Tnos[tt])
        Tm = np.bincount(ib, sec.A * Tf, len(yc)) / np.maximum(Ab, 1e-30)
        ok = Ab > 0
        cor = cmap(0.1 + 0.75 * i / max(1, len(tf_ord) - 1))
        ax[1].plot((yc[ok] + t / 2) * 1e3, Tm[ok], color=cor, lw=1.6, label=f"{tt / 60:.0f} min")
    ax[1].set_xlabel("distância à face exposta (mm)"); ax[1].set_ylabel("temperatura média da cerâmica (°C)")
    ax[1].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    S.figura(fig, "parede_deformada.png", "Deformada e temperaturas através da espessura")
    if niveis:
        fr = np.array([x[0] for x in niveis]); Rs = np.array([x[2] for x in niveis])
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        ok = np.isfinite(Rs)
        ax.plot(Rs[ok], 100 * fr[ok], "o-", color=CORES[0], lw=2, mec="white", mew=1.2)
        if (~ok).any():
            ax.plot(np.full((~ok).sum(), tf / 60), 100 * fr[~ok], ">", color=CORES[0], ms=8)
        if R is not None:
            ax.plot([R], [100 * P / Pmax], "s", color=CORES[1], ms=7, mec="white", label="carga informada")
            ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel("tempo de resistência (min)"); ax.set_ylabel("carga / capacidade a frio (%)")
        ax.set_xlim(0, tf / 60 * 1.05); ax.set_ylim(0, 100)
        fig.tight_layout()
        S.figura(fig, "parede_tempo_carga.png", "Tempo de resistência × nível de carga")
    Te = Tnos[tempos[-1]][m2.tri].mean(1)
    figura_secao(m2, Te, f"Temperatura aos {tempos[-1] / 60:.0f} min (°C)", os.path.join(S.pasta, "secao_T.png"), "°C")
    S.figuras.append({"arquivo": "secao_T.png", "titulo": "Temperatura na seção no fim do ISO 834"})
    est = []
    for tt in tf_ord:
        tn = min(Tnos, key=lambda x: abs(x - tt))
        est.append({"rotulo": f"t = {tt / 60:.0f} min", "campos": {"T": Tnos[tn], "szz": fotos[tt][1]}})
    S.pos(m2, est, "Parede: temperatura e tensão vertical nas fibras na meia altura", material=tm.material_secao(m2, geo1))
    S.avisos.append("Modelo de faixa vertical: as bordas verticais da parede não restringem a flexão. Material "
                    "homogeneizado pelas curvas do prisma aquecido a 5 °C/min (ensaio em regime permanente); "
                    + ("com fluência transitória pelos k_tr informados; " if flu is not None else
                       "sem fluência transitória (k_tr = 0 ou desligada); ")
                    + "não inclui lascamento nem a estanqueidade (E).")
    if M.ceramica.tabela_eps_termica is None:
        S.avisos.append(f"Deformação térmica por α = {M.ceramica.alfa * 1e6:.1f}×10⁻⁶/°C (provisório): a flecha "
                        "térmica e o tempo de resistência dependem diretamente desse valor.")
    if "engastado" in (en["apoio_base"], en["apoio_topo"]):
        S.avisos.append("Com apoio engastado, a dilatação da face exposta fica restringida e a cerâmica dessa face "
                        "esmaga junto ao apoio"
                        + ("; sem fluência transitória, o resultado tende a ser conservador." if flu is None else
                           ". A fluência transitória alivia parte dessa compressão."))


# ----------------------------------------------------------------------------
# Sensibilidade da parede no ISO 834 (Morris)
# ----------------------------------------------------------------------------
class _SaidaMuda(_Saida):
    """Saída que só guarda o resumo (rodadas do estudo de sensibilidade)."""

    def __init__(self):
        super().__init__(None)

    def pos(self, *a, **k):
        pass

    def tabela(self, *a, **k):
        pass

    def figura(self, fig, *a, **k):
        import matplotlib.pyplot as plt
        plt.close(fig)

    def vtk(self, *a, **k):
        pass


def avaliar_parede(q: Dict, lei=None) -> Dict:
    """Uma rodada da parede no ISO 834 sem figuras: R e I (min) e Pmax (kN)."""
    q = completar(q)
    geo = geometria(q)
    M = materiais(q, geo)
    return _executar_parede(q, geo, M, _SaidaMuda(), lambda s: None, lei=lei, rapido=True)


def _validar_sensibilidade(p) -> List[str]:
    from . import sensibilidade as sb
    en, se = p["ensaio"], p["sensibilidade"]
    erros = []
    if en["tipo"] != "iso834":
        erros.append("A sensibilidade é da parede no ISO 834: escolha esse ensaio em Ensaio e carga.")
    if en["carga"] <= 0:
        erros.append("Informe a carga constante da parede (total no comprimento).")
    pars = sb.ativos(se)
    if len(pars) < 2:
        erros.append("Marque pelo menos dois parâmetros na página Sensibilidade.")
    t2 = p["geometria"]["largura"] / 2
    for d in pars:
        lo, hi, ch = d["lo"], d["hi"], d["chave"]
        if d["tipo"] == "rel":
            if not 0 < hi < 90:
                erros.append(f"{d['nome']}: a faixa deve ficar entre 0 e 90%.")
            continue
        if lo == hi:
            erros.append(f"{d['nome']}: os dois extremos da faixa são iguais.")
        if ch == "imp" and min(lo, hi) <= 0:
            erros.append("Imperfeição: use H dividida por valores positivos (por exemplo, 1000 e 250).")
        if ch == "exc" and max(abs(lo), abs(hi)) >= t2:
            erros.append("Excentricidade: a faixa deve ficar dentro de metade da espessura da parede.")
        if ch in ("eps", "eps_furos") and not (0 <= min(lo, hi) and max(lo, hi) <= 1):
            erros.append(f"{d['nome']}: use valores entre 0 e 1.")
        if ch in ("h", "alfa", "ktr") and min(lo, hi) < 0:
            erros.append(f"{d['nome']}: use valores positivos.")
    if int(se["trajetorias"]) < 2:
        erros.append("Use pelo menos 2 trajetórias.")
    if int(se["niveis"]) % 2:
        erros.append("O número de níveis da grade deve ser par (4, 6 ou 8).")
    return erros


def _base_sensibilidade(p):
    """Projeto nominal das rodadas: parede no ISO 834, sem níveis de carga, malha rápida."""
    b = completar(copy.deepcopy(p))
    b["modelo"]["analise"] = "parede"
    b["ensaio"].update({"tipo": "iso834", "niveis_carga": 0})
    if b["sensibilidade"].get("malha_rapida", True):
        b["malha"]["h"] = max(b["malha"]["h"], 10)
        b["malha"]["n_elementos"] = min(int(b["malha"]["n_elementos"]), 16)
    return b


def _faixa_texto(d):
    if d["tipo"] == "rel":
        return f"±{d['hi']:g}%"
    if d["chave"] == "imp":
        return f"H/{d['lo']:g} a H/{d['hi']:g}"
    u = f" {d['unidade']}" if d["unidade"] else ""
    return f"{d['lo']:g} a {d['hi']:g}{u}".replace(".", ",")


def _comportamento(ms, sg, mmax):
    """Leitura do par (μ*, σ): Garcia Sanchez et al. (2014)."""
    if not np.isfinite(ms) or ms < 0.05 * mmax:
        return "efeito pequeno"
    r = sg / ms if ms > 0 else np.inf
    return ("quase linear" if r < 0.1 else "monotônico" if r < 0.5 else "quase monotônico" if r < 1.0
            else "não linear ou com interações")


def _executar_sensibilidade(p, geo, S, log):
    from . import sensibilidade as sb
    se = p["sensibilidade"]
    base = _base_sensibilidade(p)
    pars = sb.ativos(se)
    k, r, nv = len(pars), int(se["trajetorias"]), int(se["niveis"])
    X = sb.trajetorias(k, r, nv, int(se["semente"]))
    projetos = [base] + [sb.aplicar(base, pars, x) for x in X.reshape(-1, k)]
    procs = int(se["processos"]) or sb.processos_padrao()
    ma = base["malha"]
    log(f"Sensibilidade (Morris): {k} parâmetros, {r} trajetórias, grade de {nv} níveis; "
        f"{len(projetos)} rodadas da parede (h = {ma['h']:g} mm, {_n_elementos(base)} elementos de viga), "
        f"{procs} processo(s)")
    S.item("Parâmetros variados", k)
    S.item("Rodadas da parede", len(projetos))
    geo_b = geometria(base)
    lei = _lei_parede(base, geo_b, log, S)
    t0 = time.time()
    res = sb.rodar(projetos, lei, procs, log)
    S.item("Tempo médio por rodada (com os processos em paralelo)", (time.time() - t0) / len(projetos), "s")
    falhas = [(i, x["erro"]) for i, x in enumerate(res) if "erro" in x]
    nom = res[0]
    if "erro" in nom:
        raise RuntimeError("A rodada nominal da parede falhou: " + nom["erro"])
    corpo = res[1:]

    def col(ch):
        return np.array([x.get(ch, np.nan) for x in corpo], float).reshape(r, k + 1)

    Y = {"R": col("R"), "I": col("I"), "Pmax": col("Pmax")}
    cens = {"R": int(sum(x.get("R_censurado", False) for x in corpo)),
            "I": int(sum(x.get("I_censurado", False) for x in corpo))}
    ef = {s: sb.efeitos(X, Y[s]) for s in Y}
    tf = nom["tf"]
    # ------------------------------------------------------------ resumo
    S.item("Tempo de resistência nominal (malha da triagem)", nom["R"] if not nom["R_censurado"] else f"> {tf:.0f}",
           "min")
    S.item("Tempo de isolamento nominal", nom["I"] if not nom["I_censurado"] else f"> {tf:.0f}", "min")
    S.item("Capacidade a frio nominal", nom["Pmax"], "kN")
    for s, rot in (("R", "tempo de resistência"), ("I", "isolamento")):
        e = ef[s]
        if np.all(~np.isfinite(e["mu_estrela"])):
            continue
        i = int(np.nanargmax(e["mu_estrela"]))
        S.item(f"Maior efeito no {rot}", f"{pars[i]['nome']} (μ* = {e['mu_estrela'][i]:.0f} min)")
    vR = Y["R"][np.isfinite(Y["R"])]
    if vR.size:
        S.item("Tempo de resistência nas rodadas (mínimo a máximo)", f"{vR.min():.0f} a {vR.max():.0f}", "min")
    # ------------------------------------------------------------ tabelas e arquivos
    nomes = [d["nome"] for d in pars]
    for s, tit, u in sb.SAIDAS:
        e = ef[s]
        mmax = np.nanmax(e["mu_estrela"]) if np.any(np.isfinite(e["mu_estrela"])) else 0.0
        ordem = np.argsort(-np.nan_to_num(e["mu_estrela"], nan=-1))
        linhas = [[nomes[i], _faixa_texto(pars[i]), round(float(e["mu_estrela"][i]), 2),
                   f"{e['ic_lo'][i]:.1f} a {e['ic_hi'][i]:.1f}".replace(".", ","), round(float(e["mu"][i]), 2),
                   round(float(e["sigma"][i]), 2),
                   round(float(e["sigma"][i] / e["mu_estrela"][i]), 2) if e["mu_estrela"][i] > 1e-9 else "—",
                   _comportamento(e["mu_estrela"][i], e["sigma"][i], mmax)] for i in ordem]
        S.tabela(f"Efeitos de Morris {'na' if tit.startswith('capacidade') else 'no'} {tit} ({u})",
                 ["parâmetro", "faixa", f"μ* ({u})", "IC 95% de μ*", f"μ ({u})", f"σ ({u})", "σ/μ*", "comportamento"],
                 linhas, arquivo=f"sensibilidade_{s}.csv")
    cab = ["rodada", "trajetória", "passo"] + [f"{d['curto']} ({d['unidade']})" if d["unidade"] else d["curto"]
                                                for d in pars] + ["R (min)", "I (min)", "Pmax (kN)",
                                                                  "flecha no fim (mm)", "observação"]
    linhas = []
    for j, (x, rr) in enumerate(zip([None] + list(X.reshape(-1, k)), res)):
        tj, ps = ("nominal", "") if j == 0 else ((j - 1) // (k + 1) + 1, (j - 1) % (k + 1))
        vals = (["nominal"] * k if x is None else
                [round(sb.valor(d, xi), 4) if d["tipo"] == "abs" else f"{sb.valor(d, xi):+.1f}%"
                 for d, xi in zip(pars, x)])
        obs = rr.get("erro", "") or "; ".join(n for n, c in (("sem ruína até o fim", rr.get("R_censurado")),
                                                               ("sem perda de isolamento", rr.get("I_censurado"))) if c)
        linhas.append([j, tj, ps] + vals + [round(rr.get("R", np.nan), 2), round(rr.get("I", np.nan), 2),
                                             round(rr.get("Pmax", np.nan), 1), round(rr.get("flecha", np.nan), 1),
                                             obs])
    S.tabela("Rodadas do estudo (plano de Morris)", cab, linhas, arquivo="sensibilidade_rodadas.csv")
    # ------------------------------------------------------------ figuras
    _figuras_sensibilidade(S, pars, ef)
    # ------------------------------------------------------------ avisos
    S.avisos.append("Morris é um método de triagem: ordena as entradas pelo efeito médio μ*, a mudança média da "
                    "saída quando o parâmetro percorre a faixa inteira (μ com sinal: + quando a saída cresce do "
                    "primeiro para o segundo extremo). σ alto indica efeito não linear ou interação com outras "
                    "entradas. Os resultados valem para as faixas informadas.")
    if base["sensibilidade"].get("malha_rapida", True):
        S.avisos.append("Triagem com a malha rápida: os tempos diferem alguns minutos dos da malha recomendada "
                        "(5 mm e 32 elementos). Na parede BCEPQ, a ordem dos parâmetros foi a mesma nas duas "
                        "malhas; para outra parede, confirme repetindo o estudo sem a malha rápida (cerca de 6 "
                        "vezes mais demorado).")
    if cens["R"]:
        S.avisos.append(f"Em {cens['R']} rodada(s) a parede não perdeu o equilíbrio até {tf:.0f} min; o tempo de "
                        f"resistência entrou como {tf:.0f} min (limite inferior), o que reduz os efeitos no R. "
                        "Aumente a duração do ISO 834 se quiser medi-los por inteiro.")
    if cens["I"]:
        S.avisos.append(f"Em {cens['I']} rodada(s) a face não exposta não atingiu 140/180 °C até {tf:.0f} min; o "
                        f"isolamento entrou como {tf:.0f} min.")
    if falhas:
        S.avisos.append(f"{len(falhas)} rodada(s) falharam e ficaram de fora dos efeitos (primeira: rodada "
                        f"{falhas[0][0]}, {falhas[0][1]}).")
    if any(d["chave"] == "ktr" for d in pars) and not p["modelo"].get("fluencia", True):
        S.avisos.append("A fluência transitória estava desligada no projeto; nas rodadas ela foi ligada para "
                        "variar o k_tr do bloco.")


def _figuras_sensibilidade(S, pars, ef):
    plt = _plt()
    from . import sensibilidade as sb
    tinta, tinta2, grade = "#2b2a28", "#6b6a66", "#e4e3df"
    curtos = [d["curto"] for d in pars]
    nomes = [d["nome"] for d in pars]
    oR = np.argsort(np.nan_to_num(ef["R"]["mu_estrela"], nan=-1))       # ordem do R, do menor para o maior
    # ---------------------------------------------- μ* com IC (pequenos múltiplos, mesma ordem)
    fig, ax = plt.subplots(1, 3, figsize=(11, 0.34 * len(pars) + 1.6), sharey=True,
                           gridspec_kw={"width_ratios": [1.15, 1.15, 0.9]})
    y = np.arange(len(pars))
    for a, (s, tit, u) in zip(ax, sb.SAIDAS):
        e = ef[s]
        ms = np.nan_to_num(e["mu_estrela"][oR])
        a.barh(y, ms, height=0.56, color=CORES[0], zorder=2)
        lo, hi = np.nan_to_num(e["ic_lo"][oR]), np.nan_to_num(e["ic_hi"][oR])
        a.hlines(y, lo, hi, color=tinta, lw=1.0, zorder=3)
        top = np.argsort(-ms)[:3]
        xmax = max(hi.max(), ms.max(), 1e-9)
        for i in top:
            if ms[i] > 0.02 * xmax:
                a.text(hi[i] + 0.02 * xmax, y[i], f"{ms[i]:.0f}" if xmax >= 10 else f"{ms[i]:.1f}",
                       va="center", fontsize=8, color=tinta)
        a.set_xlim(0, xmax * 1.18)
        a.set_title(f"{tit} ({u})", fontsize=9, loc="left", color=tinta)
        a.grid(axis="y", visible=False); a.grid(axis="x", color=grade)
        a.tick_params(axis="y", length=0)
        a.set_xlabel(f"μ* ({u})", color=tinta2)
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([nomes[i] for i in oR], fontsize=8)
    fig.tight_layout()
    S.figura(fig, "sensibilidade_mu.png", "Efeito médio μ* de cada parâmetro, com intervalo de 95% (bootstrap); "
             "mesma ordem nos três painéis, pelo tempo de resistência")
    # ---------------------------------------------- plano μ*–σ
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for a, s, tit in ((ax[0], "R", "tempo de resistência"), (ax[1], "I", "tempo de isolamento")):
        e = ef[s]
        ms, sg = np.nan_to_num(e["mu_estrela"]), np.nan_to_num(e["sigma"])
        lim = max(ms.max(), sg.max(), 1e-9) * 1.15
        xs = np.array([0, lim])
        for rz, rot in ((1.0, "σ/μ* = 1"), (0.5, "0,5"), (0.1, "0,1")):
            a.plot(xs, rz * xs, color="#b9b7b1", lw=0.9, ls=(0, (4, 3)), zorder=1)
            xe = lim if rz <= 1 else lim / rz
            a.text(xe * 0.98, rz * xe * 0.98, rot, fontsize=7, color=tinta2, ha="right", va="bottom")
        a.scatter(ms, sg, s=46, color=CORES[0], edgecolor="white", linewidth=1.5, zorder=3)
        for i in range(len(pars)):
            if ms[i] >= 0.08 * ms.max():
                a.annotate(curtos[i], (ms[i], sg[i]), xytext=(5, 3), textcoords="offset points", fontsize=8,
                           color=tinta)
        a.set_xlim(0, lim); a.set_ylim(0, lim)
        a.set_xlabel("μ* (min)", color=tinta2); a.set_ylabel("σ (min)", color=tinta2)
        a.set_title(tit, fontsize=9, loc="left", color=tinta)
    fig.tight_layout()
    S.figura(fig, "sensibilidade_mu_sigma.png", "Plano μ*–σ: à direita, efeito maior; acima da linha σ/μ* = 1, "
             "efeito não linear ou com interações (parâmetros pequenos sem rótulo)")
