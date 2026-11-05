"""Biblioteca de blocos, materiais e argamassas.

Cada item traz as mudanças que aplica ao projeto (grupo -> campos, nas
unidades da interface), a fonte e uma nota do que falta conferir. Os itens com
"fase" > FASE_ATUAL ainda não rodam (geometria em desenvolvimento) e aparecem
desabilitados na interface. Fontes e discussão: documento "Biblioteca de
blocos e tijolos do MEF".

Convenção de dimensões: comprimento_bloco e altura_fiada são de eixo a eixo
(bloco + junta); com junta seca, a junta é uma camada fina dentro do bloco.
"""
from typing import Dict

from .parametros import EN1992_FC_FATOR, EN1992_FT_FATOR, EN1992_E_FATOR, EN1992_EPS_TERMICA

FASE_ATUAL = 3


def _ft_argamassa(fc):
    """ft = 0,3 fc^(2/3) em MPa (Eurocode 2; usada por Dauda, Iuorio e Lourenço, 2020)."""
    return round(0.3 * fc ** (2 / 3), 2)


# ---------------------------------------------------------------------------
# Argamassas (Mohamad, Lourenço e Roman, 2006: fa e Ea tangente inicial, bloco M1)
# ---------------------------------------------------------------------------
_TERM_ARG = {"arg_k": 1.15, "arg_rho": 2000, "arg_c": 1000}     # NBR 15220-2, via LabEEE/UFSC

ARGAMASSAS: Dict[str, Dict] = {
    "1:0,25:3": {"nome": "1 : 0,25 : 3 (cimento : cal : areia), fa = 19,9 MPa",
                 "mudancas": {"ceramica": {"arg_fc": 19.9, "arg_E": 11230, **_TERM_ARG}}},
    "1:0,5:4,5": {"nome": "1 : 0,5 : 4,5, fa = 8,6 MPa",
                  "mudancas": {"ceramica": {"arg_fc": 8.63, "arg_E": 6409, **_TERM_ARG}}},
    "1:1:6": {"nome": "1 : 1 : 6, fa = 4,2 MPa",
              "mudancas": {"ceramica": {"arg_fc": 4.20, "arg_E": 4033, **_TERM_ARG}}},
    "1:2:9": {"nome": "1 : 2 : 9, fa = 2,3 MPa",
              "mudancas": {"ceramica": {"arg_fc": 2.29, "arg_E": 2042, **_TERM_ARG}}},
    "personalizada": {"nome": "Personalizada (valores informados abaixo)", "mudancas": {}},
}
for _k, _v in ARGAMASSAS.items():
    c = _v["mudancas"].get("ceramica")
    if c:
        c["arg_ft"] = _ft_argamassa(c["arg_fc"])
        c["arg_nu"] = 0.20
        _v["fonte"] = "Mohamad, Lourenço e Roman (2006); térmicas da NBR 15220-2; ft = 0,3 fc^(2/3)"
        _v["nota"] = ("fa e Ea medidos em cilindros 5 × 10 cm (bloco M1). ν = 0,20 é nominal e α da argamassa "
                      "é provisório.")
    else:
        _v["fonte"], _v["nota"] = "", "Informe as propriedades medidas da argamassa."


# ---------------------------------------------------------------------------
# Grautes (Camacho et al., 2015: traços em massa, abatimento 200 ± 30 mm)
# ---------------------------------------------------------------------------
def _E_concreto(fc):
    """Ecm = 22 (fcm/10)^0,3 GPa (EN 1992-1-1), em MPa."""
    return round(22000 * (fc / 10) ** 0.3, -2)


_TERM_GRA = {"gra_k": 1.75, "gra_rho": 2400, "gra_c": 1000}       # concreto, NBR 15220-2
GRAUTES: Dict[str, Dict] = {
    "G1": {"nome": "G1: 1 : 1,83 : 2,17 (cimento : areia : brita, massa), a/c 0,67, fgm = 17,7 MPa",
           "mudancas": {"ceramica": {"gra_fc": 17.73, **_TERM_GRA}}},
    "G2": {"nome": "G2: 1 : 1,20 : 1,80, a/c 0,55, fgm = 26,7 MPa",
           "mudancas": {"ceramica": {"gra_fc": 26.66, **_TERM_GRA}}},
    "personalizado": {"nome": "Personalizado (valores informados abaixo)", "mudancas": {}},
}
for _k, _v in GRAUTES.items():
    c = _v["mudancas"].get("ceramica")
    if c:
        c["gra_ft"] = _ft_argamassa(c["gra_fc"])
        c["gra_E"] = _E_concreto(c["gra_fc"])
        c["gra_nu"] = 0.20
        _v["fonte"] = "Camacho et al. (2015); E pela EN 1992-1-1; ft = 0,3 fc^(2/3); térmicas da NBR 15220-2"
        _v["nota"] = ("Só fgm foi medido; E, ft e ν são estimativas. A quente, o graute segue o concreto silicoso "
                      "da EN 1992-1-2 (opção abaixo).")
    else:
        _v["fonte"], _v["nota"] = "", "Informe as propriedades medidas do graute."


# ---------------------------------------------------------------------------
# Materiais do bloco
# ---------------------------------------------------------------------------
def _concreto_k_superior():
    """Condutividade do concreto, limite superior da EN 1992-1-2 (W/m.K)."""
    return [[T, round(2 - 0.2451 * (T / 100) + 0.0107 * (T / 100) ** 2, 3)] for T in range(20, 1201, 100)]


def _concreto_c_seco():
    """Calor específico do concreto seco, EN 1992-1-2 (J/kg.K), sem o pico de umidade."""
    return [[20, 900], [100, 900], [200, 1000], [400, 1100], [1200, 1100]]


def _tab(t, f=1.0):
    """Tabela [(T, v)] -> [[T, v*f]] (formato da interface)."""
    return [[T, round(v * f, 6)] for T, v in t]


MATERIAIS: Dict[str, Dict] = {
    "bcepq": {"nome": "Cerâmica BCEPQ (medida)",
              "mudancas": {"ceramica": {"k": 0.52, "rho": 1210, "c": 1285, "E_bloco_bruto": 2801, "nu": 0.18,
                                        "alfa": 6.0, "tabela_k": [], "tabela_c": [], "tabela_rho": [],
                                        "tabela_E_fator": [], "tabela_eps_termica": []},
                           "modelo": {"tabela_ft_fator": [], "tabela_fc_fator": [], "k_tr_bloco": 0.0,
                                      "tabela_fluencia_bloco": []}},
              "fonte": "Manuscrito BCEPQ e Cerâmica 2026",
              "nota": "α ainda provisório (sem dilatometria)."},
    "ceramica_nbr15220": {"nome": "Cerâmica, NBR 15220-2 (1600 kg/m³)",
                          "mudancas": {"ceramica": {"k": 0.90, "rho": 1600, "c": 920, "tabela_k": [],
                                                    "tabela_c": [], "tabela_rho": [], "tabela_E_fator": [],
                                                    "tabela_eps_termica": []},
                                       "modelo": {"tabela_ft_fator": [], "tabela_fc_fator": [], "k_tr_bloco": 0.0,
                                      "tabela_fluencia_bloco": []}},
                          "fonte": "NBR 15220-2, via catálogo do LabEEE/UFSC",
                          "nota": "Só as propriedades térmicas vêm da norma: informe E, ν e α do bloco ou calibre."},
    "concreto_nbr15220": {"nome": "Concreto (bloco), NBR 15220-2 + EN 1992-1-2",
                          "mudancas": {"ceramica": {"k": 1.75, "rho": 2400, "c": 1000, "nu": 0.20, "alfa": 10.0,
                                                    "tabela_k": _concreto_k_superior(),
                                                    "tabela_c": _concreto_c_seco(), "tabela_rho": [],
                                                    "tabela_E_fator": _tab(EN1992_E_FATOR),
                                                    "tabela_eps_termica": _tab(EN1992_EPS_TERMICA, 1e3)},
                                       "modelo": {"tabela_fc_fator": _tab(EN1992_FC_FATOR),
                                                  "tabela_ft_fator": _tab(EN1992_FT_FATOR),
                                                  "k_tr_bloco": 0.0, "tabela_fluencia_bloco": []}},
                          "fonte": ("NBR 15220-2 (20 °C); EN 1992-1-2, agregado silicoso: k (limite superior), "
                                    "c seco, fc(θ) e εc1(θ) da Tabela 3.1, ft(θ) e deformação térmica"),
                          "nota": ("Curvas do concreto estrutural, aplicadas ao bloco como aproximação. O pico de calor "
                                   "específico da umidade (100 a 200 °C) não está incluído. O fator de E é fc(θ)/εc1(θ) "
                                   "em relação a 20 °C. A deformação térmica não foi conferida numa cópia aberta da "
                                   "norma. Como as curvas da EN já incluem a fluência transitória, k_tr do bloco fica "
                                   "0. Informe E, fc e ft do bloco.")},
    "personalizado": {"nome": "Personalizado (valores informados abaixo)", "mudancas": {},
                      "fonte": "", "nota": "Informe as propriedades medidas do bloco."},
}


# ---------------------------------------------------------------------------
# Blocos e tijolos
# ---------------------------------------------------------------------------
def _geo(L, H, C, junta, **kw):
    """Geometria de eixo a eixo a partir das dimensões reais (mm) e da junta."""
    seca = junta == 0
    g = {"largura": L, "altura_fiada": H + junta, "comprimento_bloco": C + junta,
         "tipo_junta": "seca" if seca else "argamassa", "espessura_junta": 2 if seca else junta}
    g.update(kw)
    return g


_CONCRETO = {**MATERIAIS["concreto_nbr15220"]["mudancas"]["ceramica"]}
_CERAMICA_EST = {**MATERIAIS["ceramica_nbr15220"]["mudancas"]["ceramica"], "nu": 0.15, "alfa": 6.0}
# índices de ductilidade: concreto pelo Model Code 2010 (G_f = 73 fcm^0,18 N/m; G_f/ft ≈ 0,1 mm) e G_c/fc
# como o da argamassa; cerâmica com os valores das unidades cerâmicas (Dauda et al., 2020)
_MODELO_CONCRETO = {"du_t": 0.10, "du_c": 1.6, **MATERIAIS["concreto_nbr15220"]["mudancas"]["modelo"]}
_MODELO_CERAMICA = {"du_t": 0.018, "du_c": 0.33, "tabela_ft_fator": [], "tabela_fc_fator": [],
                    "k_tr_bloco": 0.0, "tabela_fluencia_bloco": []}

BLOCOS: Dict[str, Dict] = {
    "bcepq": {
        "nome": "BCEPQ (bloco de encaixe, UFV): 140 × 60 × 280 mm, 2 furos de 84 mm, junta seca",
        "fase": 0,
        "mudancas": {"geometria": {"comprimento_bloco": 280, "largura": 140, "altura_fiada": 60, "diametro_furo": 84,
                                   "furos_por_bloco": 2, "blocos_por_fiada": 1.5, "n_fiadas": 5,
                                   "amarracao": "contrafiada", "tipo_junta": "seca", "espessura_junta": 2,
                                   "capeamento": True, "forma_furo": "circular", "grauteamento": "nenhum"},
                     "ceramica": {**MATERIAIS["bcepq"]["mudancas"]["ceramica"], "material_biblioteca": "bcepq"},
                     "modelo": {**_MODELO_CERAMICA, "ft": 0.66, "fc": 2.48}},
        "fonte": "Ensaios da UFV e manuscrito BCEPQ",
        "nota": "fc = 2,48 MPa é a efetiva calibrada no prisma a 30 °C."},
    "macico_10x6x22": {
        "nome": "Tijolo maciço cerâmico 100 × 60 × 220 mm, argamassa 1:1:6 de 10 mm",
        "fase": 1,
        "mudancas": {"geometria": _geo(100, 60, 220, 10, diametro_furo=0, furos_por_bloco=0, blocos_por_fiada=1,
                                       n_fiadas=5, amarracao="prumo", assentamento="total", forma_furo="circular",
                                       grauteamento="nenhum",
                                       juntas_verticais_preenchidas=True, capeamento=True),
                     "ceramica": {**MATERIAIS["ceramica_nbr15220"]["mudancas"]["ceramica"],
                                  **ARGAMASSAS["1:1:6"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "ceramica_nbr15220", "argamassa_biblioteca": "1:1:6"},
                     "modelo": {**_MODELO_CERAMICA, "ft": 0.66, "fc": 5.76}},
        "fonte": "Dimensões e térmicas: catálogo do LabEEE/UFSC (NBR 15220-2); argamassa: Mohamad et al. (2006)",
        "nota": ("E, ν, α, fc e ft do tijolo ficaram com os valores da cerâmica BCEPQ (fc do material = 5,76 MPa, "
                 "ft = 0,66 MPa de Gooch et al., 2023): troque pelos do seu lote.")},
    "solo_cimento_modular": {
        "nome": "Tijolo modular de solo-cimento 125 × 62,5 × 250 mm, 2 furos — em desenvolvimento (medir furos)",
        "fase": 9, "mudancas": {}, "fonte": "Floresta/Sahara; Garcia Jr. (2025)",
        "nota": "O diâmetro dos furos e as propriedades do solo-cimento não têm fonte aberta."},
    "concreto_m15_39": {
        "nome": "Bloco de concreto 140 × 190 × 390 mm, 2 células, argamassa 1:0,5:4,5 (Mohamad et al., 2006)",
        "fase": 2,
        "mudancas": {"geometria": _geo(140, 190, 390, 10, forma_furo="retangular", furos_por_bloco=2,
                                       parede_externa=32, septo=30, parede_topo=30, celulas_parede=0,
                                       diametro_furo=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="nenhum"),
                     "ceramica": {**_CONCRETO, "E_bloco_bruto": 9160,
                                  **ARGAMASSAS["1:0,5:4,5"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "concreto_nbr15220", "argamassa_biblioteca": "1:0,5:4,5"},
                     "modelo": {**_MODELO_CONCRETO, "ft": 1.14, "fc": 18.5}},
        "fonte": "Mohamad, Lourenço e Roman (2006); NBR 6136; térmicas da NBR 15220-2 e EN 1992-1-2",
        "nota": ("fb = 18,5 MPa na área líquida (razão 0,58). As paredes têm 25 mm no topo e são cônicas: o modelo usa "
                 "paredes prismáticas equivalentes (32 mm nas faces, 30 mm nos septos e topos), que dão a mesma área "
                 "líquida. ft = 0,1137 fb^0,789 (mesma fonte). E do bloco não medido: E líquido de 15,8 GPa sai do "
                 "módulo do prisma (15,15 GPa) pela associação em série. Confira o número de fiadas e o "
                 "assentamento do prisma na fonte.")},
    "concreto_m15_29": {
        "nome": "Bloco de concreto 140 × 190 × 290 mm (B1), prisma de 3 blocos, graute G1 (Camacho et al., 2015)",
        "fase": 2,
        "mudancas": {"geometria": _geo(140, 190, 290, 10, forma_furo="retangular", furos_por_bloco=2,
                                       parede_externa=28, septo=28, parede_topo=28, celulas_parede=0,
                                       diametro_furo=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="todos"),
                     "ceramica": {**_CONCRETO, "E_bloco_bruto": 5300,
                                  **ARGAMASSAS["1:0,5:4,5"]["mudancas"]["ceramica"],
                                  "arg_fc": 6.10, "arg_ft": _ft_argamassa(6.10), "arg_E": 5050,
                                  **GRAUTES["G1"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "concreto_nbr15220", "argamassa_biblioteca": "personalizada",
                                  "graute_biblioteca": "G1"},
                     "modelo": {**_MODELO_CONCRETO, "ft": 0.98, "fc": 15.39}},
        "fonte": "Camacho, Logullo, Parsekian e Soudais (2015); NBR 6136",
        "nota": ("Bloco B1: fb = 8,64 MPa bruto e 15,39 MPa líquido (razão 0,57; paredes equivalentes de 28 mm). "
                 "Prismas medidos (área bruta): oco 5,63 MPa; com graute G1 11,29 MPa; com G2 10,50 MPa. Argamassa "
                 "1:0,5:4,5 com 6,10 MPa (E = 5050 MPa, interpolado entre as argamassas de Mohamad et al., 2006). E do bloco estimado pela "
                 "deformação de ruptura (3,3‰): E líquido ≈ 2 fb/εu = 9,3 GPa. ft = 0,1137 fb^0,789.")},
    "concreto_m20_39": {
        "nome": "Bloco de concreto M-20 190 × 190 × 390 mm, classe A (NBR 6136), 2 células",
        "fase": 2,
        "mudancas": {"geometria": _geo(190, 190, 390, 10, forma_furo="retangular", furos_por_bloco=2,
                                       parede_externa=32, septo=25, parede_topo=25, celulas_parede=0,
                                       diametro_furo=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="nenhum"),
                     "ceramica": {**_CONCRETO, "E_bloco_bruto": 7340,
                                  **ARGAMASSAS["1:0,5:4,5"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "concreto_nbr15220", "argamassa_biblioteca": "1:0,5:4,5"},
                     "modelo": {**_MODELO_CONCRETO, "ft": 1.07, "fc": 17.2}},
        "fonte": "NBR 6136 (paredes mínimas da classe A); demais valores estimados",
        "nota": ("Paredes nas espessuras mínimas (32 mm nas faces, 25 mm nos septos e topos; área líquida 0,46). "
                 "fc = fbk mínimo da classe A (8 MPa bruto) na área líquida; ft = 0,1137 fb^0,789; E líquido de "
                 "15,8 GPa como no M-15 de Mohamad et al. Troque pelos valores do seu lote.")},
    "concreto_c_9_39": {
        "nome": "Bloco de concreto de vedação classe C 90 × 190 × 390 mm, 2 células (NBR 6136)",
        "fase": 2,
        "mudancas": {"geometria": _geo(90, 190, 390, 10, forma_furo="retangular", furos_por_bloco=2,
                                       parede_externa=18, septo=18, parede_topo=18, celulas_parede=0,
                                       diametro_furo=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="nenhum"),
                     "ceramica": {**_CONCRETO, "E_bloco_bruto": 3000,
                                  **ARGAMASSAS["1:1:6"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "concreto_nbr15220", "argamassa_biblioteca": "1:1:6"},
                     "modelo": {**_MODELO_CONCRETO, "ft": 0.48, "fc": 6.2}},
        "fonte": "NBR 6136 (paredes de 18 mm, classe C); demais valores estimados",
        "nota": ("fc = fbk mínimo da classe C (3 MPa bruto) na área líquida (razão 0,48); ft = 0,1137 fb^0,789; "
                 "E líquido ≈ 1000 fb (ordem de grandeza). Troque pelos valores do seu lote.")},
    "est_ceramico_14x19x29_pm": {
        "nome": "Estrutural cerâmico 140 × 190 × 290 mm, paredes maciças, 2 células (NBR 15270)",
        "fase": 2,
        "mudancas": {"geometria": _geo(140, 190, 290, 10, forma_furo="retangular", furos_por_bloco=2,
                                       parede_externa=20, septo=8, parede_topo=20, celulas_parede=0,
                                       diametro_furo=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="nenhum"),
                     "ceramica": {**_CERAMICA_EST, "E_bloco_bruto": 6000,
                                  **ARGAMASSAS["1:1:6"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "ceramica_nbr15220", "argamassa_biblioteca": "1:1:6"},
                     "modelo": {**_MODELO_CERAMICA, "ft": 1.49, "fc": 14.9}},
        "fonte": "NBR 15270-2 (paredes externas ≥ 20 mm, septos ≥ 8 mm); térmicas da NBR 15220-2",
        "nota": ("Geometria nas espessuras mínimas (área líquida 0,40); blocos reais têm septos mais grossos: meça o "
                 "seu. fc de um bloco de 6 MPa bruto na área líquida; ft = 0,10 fc e E líquido ≈ 1000 fb são "
                 "estimativas de ordem de grandeza.")},
    "est_ceramico_14x19x29_pv": {
        "nome": "Estrutural cerâmico 140 × 190 × 290 mm, paredes vazadas (NBR 15270)",
        "fase": 2,
        "mudancas": {"geometria": _geo(140, 190, 290, 10, forma_furo="retangular", furos_por_bloco=2,
                                       parede_externa=30, septo=8, parede_topo=8, celulas_parede=6,
                                       parede_minima=8, diametro_furo=0, blocos_por_fiada=1, n_fiadas=3,
                                       amarracao="prumo", assentamento="total", juntas_verticais_preenchidas=True,
                                       capeamento=True, grauteamento="nenhum"),
                     "ceramica": {**_CERAMICA_EST, "E_bloco_bruto": 4700,
                                  **ARGAMASSAS["1:1:6"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "ceramica_nbr15220", "argamassa_biblioteca": "1:1:6"},
                     "modelo": {**_MODELO_CERAMICA, "ft": 1.91, "fc": 19.1}},
        "fonte": "NBR 15270-2 (paredes ≥ 8 mm, septos ≥ 7 mm); térmicas da NBR 15220-2",
        "nota": ("Ilustrativo: paredes externas de 30 mm com uma fileira de 6 células de 14 mm, paredes internas de "
                 "8 mm (área líquida 0,31). fc de um bloco de 6 MPa bruto na área líquida; ft = 0,10 fc e E "
                 "líquido ≈ 1000 fb são estimativas. Meça o seu bloco.")},
    "est_ceramico_14x19x39": {
        "nome": "Estrutural cerâmico 140 × 190 × 390 mm, paredes maciças, 3 células (arranjo típico)",
        "fase": 2,
        "mudancas": {"geometria": _geo(140, 190, 390, 10, forma_furo="retangular", furos_por_bloco=3,
                                       parede_externa=20, septo=8, parede_topo=20, celulas_parede=0,
                                       diametro_furo=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="nenhum"),
                     "ceramica": {**_CERAMICA_EST, "E_bloco_bruto": 6000,
                                  **ARGAMASSAS["1:1:6"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "ceramica_nbr15220", "argamassa_biblioteca": "1:1:6"},
                     "modelo": {**_MODELO_CERAMICA, "ft": 1.49, "fc": 14.9}},
        "fonte": "NBR 15270-2; arranjo de 3 células a conferir no catálogo do fabricante",
        "nota": ("Espessuras mínimas da norma; com 3 células as juntas de contrafiada não alinham as células "
                 "(o modelo aceita, sem graute). Propriedades estimadas como no bloco de 29 cm.")},
    "vedacao_8_furos": {
        "nome": "Vedação cerâmico 90 × 190 × 190 mm, 8 furos horizontais (2 × 4), argamassa 1:2:9",
        "fase": 3,
        "mudancas": {"geometria": _geo(90, 190, 190, 10, forma_furo="horizontal", furos_largura=2, furos_altura=4,
                                       parede_externa=8, septo=7, parede_leito=8, diametro_furo=0,
                                       furos_por_bloco=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="nenhum"),
                     "ceramica": {**_CERAMICA_EST, "E_bloco_bruto": 1500,
                                  **ARGAMASSAS["1:2:9"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "ceramica_nbr15220", "argamassa_biblioteca": "1:2:9"},
                     "modelo": {**_MODELO_CERAMICA, "ft": 0.59, "fc": 5.87}},
        "fonte": "SINAPI 7271 (dimensões e furos); NBR 15270-1 (paredes e resistência mínimas); NBR 15220-2",
        "nota": ("Paredes externas e de leito de 8 mm e septos de 7 mm (acima dos mínimos de 7 e 6 mm): furos de "
                 "33,5 × 38,3 mm. As paredes que descem contínuas (faces e septo vertical) somam 26% da largura; "
                 "fc vem da resistência mínima de 1,5 MPa na área bruta, levada a essa área. ft = 0,10 fc e E "
                 "líquido ≈ 1000 fb são estimativas. Meça o seu bloco.")},
    "vedacao_6_furos": {
        "nome": "Vedação cerâmico 90 × 140 × 240 mm, 6 furos horizontais (2 × 3), argamassa 1:2:9 (Guedes et al., 2019)",
        "fase": 3,
        "mudancas": {"geometria": _geo(90, 140, 240, 10, forma_furo="horizontal", furos_largura=2, furos_altura=3,
                                       parede_externa=7.8, septo=7.1, parede_leito=7.8, diametro_furo=0,
                                       furos_por_bloco=0, blocos_por_fiada=1, n_fiadas=3, amarracao="prumo",
                                       assentamento="total", juntas_verticais_preenchidas=True, capeamento=True,
                                       grauteamento="nenhum"),
                     "ceramica": {**_CERAMICA_EST, "E_bloco_bruto": 1500,
                                  **ARGAMASSAS["1:2:9"]["mudancas"]["ceramica"],
                                  "material_biblioteca": "ceramica_nbr15220", "argamassa_biblioteca": "1:2:9"},
                     "modelo": {**_MODELO_CERAMICA, "ft": 0.60, "fc": 5.95}},
        "fonte": "Guedes et al. (2019): paredes e septos medidos em três olarias; NBR 15270-1; NBR 15220-2",
        "nota": ("Paredes externas de 7,8 mm e septos de 7,1 mm (médias das faixas medidas, 7,4 a 8,3 e 6,3 a "
                 "8,0 mm); paredes de leito iguais às externas. fc vem da resistência mínima de 1,5 MPa bruta, "
                 "levada à área das paredes contínuas (25%); ft = 0,10 fc e E líquido ≈ 1000 fb são estimativas.")},
    "personalizado": {"nome": "Personalizado (geometria informada abaixo)", "fase": 0, "mudancas": {},
                      "fonte": "", "nota": "Informe a geometria do bloco."},
}


def opcoes(tabela: Dict[str, Dict]):
    return [[k, v["nome"]] for k, v in tabela.items()]


def para_interface():
    """Mudanças e notas de cada item, por campo seletor (enviado à interface)."""
    def pack(tab):
        return {k: {"mudancas": v["mudancas"], "nota": v.get("nota", ""), "fonte": v.get("fonte", "")}
                for k, v in tab.items()}
    return {"geometria.bloco_biblioteca": pack(BLOCOS),
            "ceramica.material_biblioteca": pack(MATERIAIS),
            "ceramica.argamassa_biblioteca": pack(ARGAMASSAS),
            "ceramica.graute_biblioteca": pack(GRAUTES)}


def aplicar(p: Dict, seletor: str, chave: str) -> Dict:
    """Aplica um item ao projeto p (cópia), como a interface faz."""
    import copy
    tab = {"bloco": BLOCOS, "material": MATERIAIS, "argamassa": ARGAMASSAS, "graute": GRAUTES}[seletor]
    q = copy.deepcopy(p)
    for g, d in tab[chave]["mudancas"].items():
        q.setdefault(g, {}).update(copy.deepcopy(d))
    campo = {"bloco": ("geometria", "bloco_biblioteca"), "material": ("ceramica", "material_biblioteca"),
             "argamassa": ("ceramica", "argamassa_biblioteca"), "graute": ("ceramica", "graute_biblioteca")}[seletor]
    q[campo[0]][campo[1]] = chave
    return q
