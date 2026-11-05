"""Parâmetros do modelo, com a origem de cada valor.

Unidades internas: SI (m, s, K ou °C, Pa, W). As entradas de geometria
também estão em metros.
"""
from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np

Tabela = Optional[Sequence[Tuple[float, float]]]  # [(T em °C, valor), ...]


def interpolar(T, tabela: Tabela, padrao):
    """Interpola linearmente uma tabela (T, valor); fora dela usa os extremos.
    Sem tabela, devolve o valor padrão (constante)."""
    if tabela is None:
        return np.full(np.shape(T), padrao, dtype=float) if np.ndim(T) else float(padrao)
    t = np.asarray(tabela, dtype=float)
    return np.interp(T, t[:, 0], t[:, 1])


@dataclass
class Ceramica:
    """Cerâmica BCEPQ (bloco extrudado, prensado e queimado)."""
    k: float = 0.52            # W/m.K   Cerâmica 2026 (média de 20 medições: 0,5188)
    rho: float = 1210.0        # kg/m³   manuscrito BPQ (massa específica aparente)
    c: float = 1285.0          # J/kg.K  derivado: efusividade 898 W.s^0,5/m².K e k (Cerâmica 2026)
    E_bloco_bruto: float = 2801e6  # Pa  módulo do bloco (BPQ, n = 6), referido à área bruta
    razao_area_liquida: float = 0.717  # área líquida / bruta da seção com furos de 84 mm
    nu: float = 0.18           # BPQ (0,182)
    alfa: float = 6e-6         # 1/K     PROVISÓRIO, sem medição (o artigo usa 6e-6)
    T_ref: float = 30.0        # °C      temperatura sem deformação térmica
    # Dependência da temperatura (opcional). Ex.: tabela_k = [(20, 0.52), (1000, 0.60)]
    tabela_k: Tabela = None
    tabela_E_fator: Tabela = None   # fator multiplicador de E: [(20, 1.0), (800, 0.6)]
    tabela_eps_termica: Tabela = None  # deformação térmica livre total: [(30, 0), (573, 0.0035), ...]
    tabela_c: Tabela = None         # calor específico (J/kg.K) com a temperatura (pico da umidade, por exemplo)
    tabela_rho: Tabela = None       # massa específica (kg/m³) com a temperatura (perda de água)

    @property
    def E_material(self):
        """E do material cerâmico (área líquida). Se o E do bloco foi medido com
        tensão na área bruta, o material é mais rígido: E_mat = E_bruto / (A_liq/A_bruta)."""
        return self.E_bloco_bruto / self.razao_area_liquida

    def k_T(self, T):
        return interpolar(T, self.tabela_k, self.k)

    def c_T(self, T):
        return interpolar(T, self.tabela_c, self.c)

    def rho_T(self, T):
        return interpolar(T, self.tabela_rho, self.rho)

    @property
    def capacidade_variavel(self):
        return self.tabela_c is not None or self.tabela_rho is not None

    def E_T(self, T):
        return self.E_material * interpolar(T, self.tabela_E_fator, 1.0)

    def eps_termica(self, T):
        if self.tabela_eps_termica is not None:
            return interpolar(T, self.tabela_eps_termica, 0.0)
        return self.alfa * (np.asarray(T, float) - self.T_ref)


@dataclass
class JuntaSeca:
    """Junta seca entre fiadas (camada fina equivalente).
    E_junta deve ser calibrado para o prisma reproduzir o módulo medido
    (exemplo 02). O valor padrão sai da calibração com E_prisma = 1.100 MPa."""
    espessura: float = 0.002   # m
    E: float = 84.4e6          # Pa: calibrado (exemplo 02) para E_prisma = 1.100 MPa, malha h = 8 mm, dz = 12 mm
    nu: float = 0.18
    fator_k: float = 1.0       # condutividade da junta = fator_k * k da cerâmica


@dataclass
class Argamassa(Ceramica):
    """Argamassa de assentamento (camada contínua de espessura real, aderida ao
    bloco). Térmicas: NBR 15220-2 (via catálogo do LabEEE/UFSC). Mecânicas do
    traço 1:1:6 (cimento : cal : areia, em volume): Mohamad, Lourenço e Roman
    (2006), fa = 4,2 MPa e Ea = 4033 MPa (módulo tangente inicial). E_material
    = E_bloco_bruto (razao_area_liquida = 1)."""
    k: float = 1.15
    rho: float = 2000.0
    c: float = 1000.0
    E_bloco_bruto: float = 4033e6
    razao_area_liquida: float = 1.0
    nu: float = 0.20            # nominal
    alfa: float = 10e-6         # PROVISÓRIO: ordem de grandeza de argamassas de cimento


# ---------------------------------------------------------------------------
# Concreto de agregado silicoso a quente (EN 1992-1-2). Tabela 3.1: fc(θ)/fck e
# εc1(θ); 3.2.2.2: kc,t(θ) = 1 até 100 °C, caindo linearmente a 0 em 600 °C;
# 3.3.1: deformação térmica -1,8e-4 + 9e-6 θ + 2,3e-11 θ³ até 700 °C e 14e-3
# acima. O fator de E é a razão dos módulos secantes no pico, fc(θ)/εc1(θ),
# com piso de 0,001 (rigidez residual).
# ---------------------------------------------------------------------------
EN1992_FC_FATOR = [(20, 1.0), (100, 1.0), (200, 0.95), (300, 0.85), (400, 0.75), (500, 0.60), (600, 0.45),
                   (700, 0.30), (800, 0.15), (900, 0.08), (1000, 0.04), (1100, 0.01), (1200, 0.0)]
_EN1992_EC1 = [0.0025, 0.0040, 0.0055, 0.0070, 0.0100, 0.0150, 0.0250, 0.0250, 0.0250, 0.0250, 0.0250, 0.0250, 0.0250]
EN1992_FT_FATOR = [(20, 1.0), (100, 1.0), (600, 0.0), (1200, 0.0)]
EN1992_E_FATOR = [(T, round(max(f * 0.0025 / e, 0.001), 4)) for (T, f), e in zip(EN1992_FC_FATOR, _EN1992_EC1)]


def _eps_silicoso(T):
    return -1.8e-4 + 9e-6 * T + 2.3e-11 * T ** 3 if T <= 700 else 14e-3


EN1992_EPS_TERMICA = [(T, round(_eps_silicoso(T), 7)) for T in list(range(20, 701, 20)) + [1200]]
EN1992_EPS_TERMICA = sorted(dict(EN1992_EPS_TERMICA).items())


@dataclass
class Graute(Ceramica):
    """Graute dos furos ou células (concreto fino fluido). Padrão: graute G1 de
    Camacho et al. (2015), fgk = 17,7 MPa; E pela EN 1992-1-2/EN 1992-1-1,
    Ecm = 22 (fcm/10)^0,3 GPa (estimativa, sem medida); térmicas de concreto
    normal da NBR 15220-2. E_material = E_bloco_bruto (razao_area_liquida = 1)."""
    k: float = 1.75
    rho: float = 2400.0
    c: float = 1000.0
    E_bloco_bruto: float = 26.1e9
    razao_area_liquida: float = 1.0
    nu: float = 0.20
    alfa: float = 10e-6         # PROVISÓRIO: concreto silicoso, ordem de grandeza


@dataclass
class Aco:
    k: float = 45.0
    rho: float = 7850.0
    c: float = 500.0
    E: float = 200e9
    nu: float = 0.30
    alfa: float = 12e-6


@dataclass
class ArFuros:
    """Ar nos furos: condução + radiação entre as paredes do furo,
    representadas por uma condutividade equivalente k_rad = 4 σ T³ D ε/(2-ε)."""
    emissividade_paredes: float = 0.9
    rho_c: float = 1.0e3       # J/m³.K (desprezível)


@dataclass
class Geometria:
    """Geometria paramétrica do prisma de blocos com furos circulares verticais,
    alinhados no meio da largura. O padrão é o prisma BCEPQ da UFV: cinco
    fiadas de bloco e meio (420 x 140 mm), dois furos de 84 mm por bloco,
    juntas verticais desencontradas (fiada par: bloco 0-280 e meio-bloco
    280-420 mm; ímpar: meio-bloco 0-140 e bloco 140-420)."""
    comprimento_bloco: float = 0.280
    largura: float = 0.140
    altura_fiada: float = 0.060
    diametro_furo: float = 0.084       # 0 = bloco maciço
    n_fiadas: int = 5
    espessura_placa: float = 0.025
    furos_por_bloco: int = 2
    blocos_por_fiada: float = 1.5      # comprimento do prisma em blocos (múltiplo de 0,5)
    amarracao: str = "contrafiada"     # "contrafiada" (juntas desencontradas) ou "prumo"
    # Posições x das juntas verticais em cada fiada; None = calculadas pela amarração.
    # Aceita um número por fiada (formato antigo) ou uma tupla por fiada.
    juntas_verticais: Optional[Tuple] = None
    # Juntas: "seca" (camada fina com contato, folgas e atrito) ou "argamassa"
    # (camada contínua de espessura real). Com argamassa, comprimento_bloco e
    # altura_fiada são medidas de eixo a eixo (bloco + junta).
    tipo_junta: str = "seca"
    assentamento: str = "total"         # argamassa em toda a face ("total") ou só nas paredes externas ("lateral")
    juntas_verticais_preenchidas: bool = True
    espessura_junta: float = 0.002      # m (define a extensão real do bloco nas células retangulares)
    # Furos: "circular" (furos redondos verticais alinhados, diametro_furo),
    # "retangular" (células verticais: furos_por_bloco células principais ao
    # longo do comprimento, paredes externas, septos e paredes de topo; com
    # celulas_parede > 0, cada parede externa tem uma fileira de células
    # pequenas, como nos blocos cerâmicos de paredes vazadas) ou "horizontal"
    # (grade de furos_largura × furos_altura furos retangulares que atravessam
    # o bloco ao longo do comprimento, como nos blocos de vedação: paredes
    # externas, septos entre os furos e paredes de leito em cima e embaixo).
    forma_furo: str = "circular"
    parede_externa: float = 0.025       # m, paredes longitudinais (faces do bloco)
    septo: float = 0.025                # m, paredes transversais internas
    parede_topo: float = 0.025          # m, paredes transversais das pontas do bloco
    celulas_parede: int = 0             # células pequenas por bloco em cada parede externa
    parede_minima: float = 0.008        # m, paredes em volta das células pequenas
    grauteamento: str = "nenhum"        # "nenhum", "todos" ou "alternados" (furos ou células principais)
    furos_largura: int = 2              # furos horizontais na espessura do bloco
    furos_altura: int = 4               # furos horizontais na altura do bloco
    parede_leito: float = 0.008         # m, paredes de cima e de baixo dos blocos de furos horizontais

    @property
    def retangular(self):
        return self.forma_furo == "retangular" and self.furos_por_bloco > 0

    @property
    def horizontal(self):
        return self.forma_furo == "horizontal" and self.furos_largura > 0 and self.furos_altura > 0

    @property
    def tem_furos(self):
        """Furos verticais (redondos ou células): os que recebem graute."""
        return self.retangular or (self.forma_furo == "circular" and self.diametro_furo > 0
                                   and self.furos_por_bloco > 0)

    @property
    def espessura_casca(self):
        """Espessura das paredes externas (para o assentamento lateral)."""
        if self.retangular or self.horizontal:
            return self.parede_externa
        if self.diametro_furo <= 0:
            return self.largura / 2
        return (self.largura - self.diametro_furo) / 2

    @property
    def dimensao_furo(self):
        """Dimensão do furo na direção da espessura (radiação no ar dos furos)."""
        if self.retangular:
            return self.largura - 2 * self.parede_externa
        if self.horizontal:
            return self.largura_furo_h
        return self.diametro_furo

    # ---------------------------------------------------------------- furos horizontais
    @property
    def largura_furo_h(self):
        """Largura (na espessura do bloco) de cada furo horizontal."""
        n = self.furos_largura
        return (self.largura - 2 * self.parede_externa - (n - 1) * self.septo) / n

    def altura_furo_h(self, t_junta=None):
        """Altura de cada furo horizontal; o bloco tem altura_fiada - junta."""
        tj = self.espessura_junta if t_junta is None else t_junta
        n = self.furos_altura
        return (self.altura_fiada - tj - 2 * self.parede_leito - (n - 1) * self.septo) / n

    def furos_h_y(self):
        """Faixas (y0, y1) das colunas de furos horizontais."""
        if not self.horizontal:
            return []
        w, te, tw = self.largura_furo_h, self.parede_externa, self.septo
        return [(te + k * (w + tw), te + k * (w + tw) + w) for k in range(self.furos_largura)]

    def furos_h_z(self, zc, t_junta=None):
        """Faixas (z0, z1) das fileiras de furos horizontais de um bloco centrado em zc."""
        if not self.horizontal:
            return []
        tj = self.espessura_junta if t_junta is None else t_junta
        hz, tb, tw = self.altura_furo_h(tj), self.parede_leito, self.septo
        z0 = zc - (self.altura_fiada - tj) / 2 + tb
        return [(z0 + k * (hz + tw), z0 + k * (hz + tw) + hz) for k in range(self.furos_altura)]

    def blocos_da_fiada(self, i):
        """Intervalos (a, b) dos blocos e meios-blocos da fiada i, de eixo a eixo."""
        xs = [0.0] + list(self.juntas_por_fiada[i]) + [self.comprimento]
        return list(zip(xs[:-1], xs[1:]))

    def celulas(self, i, t_junta=None):
        """Células retangulares da fiada i: [(x0, x1, y0, y1, principal)], em m."""
        if not self.retangular:
            return []
        tj = self.espessura_junta if t_junta is None else t_junta
        W, Lt = self.largura, self.comprimento
        te, tw, tt = self.parede_externa, self.septo, self.parede_topo
        res = []
        for a, b in self.blocos_da_fiada(i):
            x0 = a + (tj / 2 if a > 1e-9 else 0.0)
            x1 = b - (tj / 2 if b < Lt - 1e-9 else 0.0)
            n = max(1, int(round(self.furos_por_bloco * (b - a) / self.comprimento_bloco)))
            Lc = (x1 - x0 - 2 * tt - (n - 1) * tw) / n
            for k in range(n):
                cx0 = x0 + tt + k * (Lc + tw)
                res.append((cx0, cx0 + Lc, te, W - te, True))
            if self.celulas_parede > 0:
                ms = max(1, int(round(self.celulas_parede * (b - a) / self.comprimento_bloco)))
                tm = self.parede_minima
                Ls = (x1 - x0 - (ms + 1) * tm) / ms
                wy = te - 2 * tm
                for k in range(ms):
                    sx0 = x0 + tm + k * (Ls + tm)
                    res.append((sx0, sx0 + Ls, tm, tm + wy, False))
                    res.append((sx0, sx0 + Ls, W - tm - wy, W - tm, False))
        return res

    def colunas_grauteadas(self, t_junta=None):
        """Faixas x (x0, x1) dos furos ou células principais grauteados, na fiada 0."""
        if self.grauteamento == "nenhum" or not self.tem_furos:
            return []
        if self.retangular:
            cols = sorted((c[0], c[1]) for c in self.celulas(0, t_junta) if c[4])
        else:
            R = self.diametro_furo / 2
            cols = [(x - R, x + R) for x in self.centros_furos]
        if self.grauteamento == "alternados":
            cols = cols[::2]
        return cols

    def area_furos(self, i=0, graute=False, t_junta=None):
        """Área dos furos (m²) na fiada i: graute=False conta todos os furos (área do
        bloco = bruta - furos); graute=True conta só os furos sem graute."""
        if self.horizontal:          # seção horizontal mínima: passa pelas fileiras de furos
            return self.furos_largura * self.largura_furo_h * self.comprimento
        if not self.tem_furos:
            return 0.0
        cols = self.colunas_grauteadas(t_junta) if graute else []

        def grauteado(x0, x1):
            return any(min(x1, b) - max(x0, a) > 0.5 * (x1 - x0) for a, b in cols)
        if self.retangular:
            return sum((c[1] - c[0]) * (c[3] - c[2]) for c in self.celulas(i, t_junta)
                       if not (c[4] and grauteado(c[0], c[1])))
        R = self.diametro_furo / 2
        return sum(np.pi * R * R for x in self.centros_furos if not grauteado(x - R, x + R))

    @property
    def comprimento(self):
        return self.blocos_por_fiada * self.comprimento_bloco

    @property
    def espacamento_furos(self):
        return self.comprimento_bloco / self.furos_por_bloco

    @property
    def centros_furos(self):
        if self.forma_furo != "circular" or self.diametro_furo <= 0 or self.furos_por_bloco <= 0:
            return ()
        s = self.espacamento_furos
        n = int(round(self.comprimento / s))
        return tuple(s / 2 + k * s for k in range(n))

    @property
    def juntas_por_fiada(self) -> Tuple[Tuple[float, ...], ...]:
        if self.juntas_verticais is not None:
            return tuple(tuple(j) if isinstance(j, (tuple, list)) else (float(j),) for j in self.juntas_verticais)
        Lb, Lt = self.comprimento_bloco, self.comprimento
        res = []
        for i in range(self.n_fiadas):
            off = Lb / 2 if (self.amarracao == "contrafiada" and i % 2 == 1) else 0.0
            xs = [round(off + k * Lb, 9) for k in range(int(Lt / Lb) + 2)]
            res.append(tuple(x for x in xs if 1e-6 < x < Lt - 1e-6))
        return tuple(res)

    @property
    def posicoes_juntas(self) -> Tuple[float, ...]:
        return tuple(sorted({x for f in self.juntas_por_fiada for x in f}))

    @property
    def altura_prisma(self):
        return self.n_fiadas * self.altura_fiada

    def validar(self):
        """Lista de problemas da geometria (vazia se estiver tudo certo)."""
        p = []
        if self.comprimento_bloco <= 0 or self.largura <= 0 or self.altura_fiada <= 0 or self.n_fiadas < 1:
            p.append("Dimensões e número de fiadas devem ser positivos.")
        if abs(self.blocos_por_fiada * 2 - round(self.blocos_por_fiada * 2)) > 1e-9 or self.blocos_por_fiada < 0.5:
            p.append("Blocos por fiada deve ser múltiplo de 0,5.")
        if self.amarracao not in ("contrafiada", "prumo"):
            p.append("Amarração deve ser 'contrafiada' ou 'prumo'.")
        if self.tipo_junta not in ("seca", "argamassa"):
            p.append("Tipo de junta deve ser 'seca' ou 'argamassa'.")
        if self.assentamento not in ("total", "lateral"):
            p.append("Assentamento deve ser 'total' ou 'lateral'.")
        if self.grauteamento not in ("nenhum", "todos", "alternados"):
            p.append("Grauteamento deve ser 'nenhum', 'todos' ou 'alternados'.")
        if self.grauteamento != "nenhum" and not self.tem_furos:
            p.append("Só furos verticais recebem graute (bloco maciço ou de furos horizontais).")
        if self.forma_furo not in ("circular", "retangular", "horizontal"):
            p.append("Forma dos furos deve ser 'circular', 'retangular' ou 'horizontal'.")
        if self.forma_furo == "horizontal":
            if self.furos_largura < 1 or self.furos_altura < 1:
                p.append("Com furos horizontais, use pelo menos 1 furo na espessura e 1 na altura.")
            elif self.largura_furo_h < 0.003:
                p.append("Furos horizontais estreitos demais (menos de 3 mm): reveja paredes externas e septos.")
            elif self.altura_furo_h() < 0.003:
                p.append("Furos horizontais baixos demais (menos de 3 mm): reveja paredes de leito, septos e junta.")
        elif self.retangular:
            n = self.furos_por_bloco
            meio = abs(self.blocos_por_fiada % 1) > 1e-9
            if (self.amarracao == "contrafiada" or meio) and n % 2 and self.grauteamento != "nenhum":
                p.append("Com graute e meio-bloco ou contrafiada, o número de células por bloco deve ser par "
                         "(as células precisam se alinhar de uma fiada para a outra).")
            if 2 * self.parede_externa >= 0.9 * self.largura:
                p.append("As paredes externas ocupam quase toda a largura do bloco.")
            if self.celulas_parede > 0 and self.parede_externa <= 2 * self.parede_minima + 1e-4:
                p.append("Com paredes vazadas, a parede externa deve ser maior que duas vezes a parede mínima.")
            for i in range(self.n_fiadas):
                if any(c[1] - c[0] < 0.004 for c in self.celulas(i)):
                    p.append("Células curtas demais (menos de 4 mm): reveja septos, paredes de topo e células.")
                    break
        elif self.diametro_furo > 0 and self.forma_furo == "circular":
            n = self.furos_por_bloco
            meio = abs(self.blocos_por_fiada % 1) > 1e-9
            if n < 1:
                p.append("Com furos, use pelo menos 1 furo por bloco.")
            elif (self.amarracao == "contrafiada" or meio) and n % 2:
                p.append("Com meio-bloco ou contrafiada, o número de furos por bloco deve ser par "
                         "(os furos precisam ficar alinhados entre fiadas).")
            R, W = self.diametro_furo / 2, self.largura
            if R >= 0.9 * W / 2:
                p.append("O furo é largo demais para a largura do bloco (parede menor que 10% da meia-largura).")
            elif n >= 1 and R >= 0.9 * self.espacamento_furos / 2:
                p.append("Furos próximos demais entre si ou das juntas; reduza o diâmetro ou o número de furos.")
        return p


@dataclass
class Ensaio:
    """Protocolo da UFV (dissertação do Jonathan): forno a 5 °C/min até a
    temperatura-alvo, 30 min de patamar e compressão ainda dentro do forno."""
    T_inicial: float = 30.0
    taxa_C_min: float = 5.0
    T_alvo: float = 600.0
    patamar_s: float = 1800.0
    h_conveccao: float = 10.0      # W/m².K (forno elétrico, ar parado) — testar 10 a 25
    emissividade: float = 0.8      # superfície cerâmica — testar 0,8 a 0,9
    carga_kgf: float = 9980.9      # carga máxima do ensaio

    @property
    def carga_N(self):
        return self.carga_kgf * 9.80665


@dataclass
class DanoCeramica:
    """Modelo de dano isotrópico da cerâmica (análise não linear).

    Critério de Rankine em tração e de von Mises em compressão (tensões
    efetivas), amolecimento exponencial e regularização pela banda de
    fissuração (Bažant e Oh, 1983). As energias de fratura saem dos índices de
    ductilidade d_u = G/f usados para unidades cerâmicas em micromodelos."""
    ft: float = 0.66e6      # Pa: tração direta de tijolo cerâmico prensado, área líquida (Gooch et al., 2023)
    fc: float = 5.76e6      # Pa: fb médio 4,13 MPa em área bruta (BPQ) / 0,717 -> área líquida
    du_t: float = 0.018e-3  # m: G_f/ft da unidade cerâmica (van der Pluijm 1992; Lourenço 1996; via Dauda et al., 2020)
    du_c: float = 0.33e-3   # m: G_c/fc da unidade cerâmica (Angelillo et al. 2014; via Dauda et al., 2020)
    d_max: float = 0.9999   # dano máximo (rigidez residual mínima, para estabilidade numérica)
    fechamento_fissuras: bool = True  # fissura de tração comprimida volta a transmitir compressão (Mazars, 1984)
    tabela_ft_fator: Tabela = None  # fatores de ft com a temperatura: [(T, fator), ...]
    tabela_fc_fator: Tabela = None

    @property
    def G_f(self):
        return self.du_t * self.ft   # N/m

    @property
    def G_c(self):
        return self.du_c * self.fc

    def ft_T(self, T):
        return self.ft * interpolar(T, self.tabela_ft_fator, 1.0)

    def fc_T(self, T):
        return self.fc * interpolar(T, self.tabela_fc_fator, 1.0)


def dano_argamassa(fc: float = 4.2e6, ft: Optional[float] = None) -> "DanoCeramica":
    """Parâmetros de dano da argamassa. ft = 0,3 fc^(2/3) (MPa; fórmula do
    Eurocode 2 usada por Dauda, Iuorio e Lourenço, 2020). Índices de
    ductilidade da argamassa: G_f/ft = 0,065 mm e G_c/fc = 1,6 mm (mesma fonte,
    a partir de Lourenço, 1996)."""
    if ft is None:
        ft = 0.3 * (fc / 1e6) ** (2 / 3) * 1e6
    return DanoCeramica(ft=ft, fc=fc, du_t=0.065e-3, du_c=1.6e-3)


def dano_graute(fc: float = 17.7e6, ft: Optional[float] = None) -> "DanoCeramica":
    """Parâmetros de dano do graute: ft = 0,3 fc^(2/3) (MPa) e os mesmos índices
    de ductilidade da argamassa (G_f/ft = 0,065 mm, G_c/fc = 1,6 mm). Para
    fc = 17,7 MPa, G_f = 0,13 N/mm, perto de G_f = 73 fcm^0,18 N/m do Model
    Code 2010 (0,12 N/mm)."""
    if ft is None:
        ft = 0.3 * (fc / 1e6) ** (2 / 3) * 1e6
    return DanoCeramica(ft=ft, fc=fc, du_t=0.065e-3, du_c=1.6e-3)


@dataclass
class JuntaAtrito:
    """Junta seca como camada fina não linear: não transmite tração e
    transmite cisalhamento até o limite de atrito de Coulomb sem coesão."""
    mu: float = 0.60             # atrito em juntas secas serradas: 0,61-0,64 (Lourenço, Ramos e Vasconcelos, 2004)
    fator_aberta: float = 1e-4   # rigidez normal residual da junta aberta (fração de E_junta)
    fator_plano: float = 1e-3    # rigidez no plano da camada (só para estabilidade)
    fator_G_min: float = 1e-3    # rigidez ao cisalhamento mínima em deslizamento
    eps_tol: float = 1e-6        # faixa de histerese do contato (deformação normal) para evitar oscilação
    pre_carga_N: float = 2e3     # carga de assentamento durante o aquecimento (prisma inteiro)
    # Irregularidade das faces de assentamento (juntas de leito). Sem carga, só
    # parte da área nominal toca: 23% em blocos sem tratamento (Ngapeya e
    # Waldmann, 2020) e ~15% subindo a ~95% perto da ruptura (Zahra et al.,
    # 2016). O resto tem folgas que fecham com a carga e concentram tensões
    # (Ngapeya, Waldmann e Scholzen, 2018). As folgas são um campo gaussiano
    # com desvio-padrão desvio_folga e comprimento de correlação comprimento_corr.
    desvio_folga: float = 0.028e-3  # m; calibrado pelo módulo do prisma (E de corda 5-30% = 1.100 MPa na
                                    # média de 3 sementes; script 08 dá 0,029 mm; malha h = 12, dz = 20 mm)
                                    # 0 = contato perfeito
    E_contato: Optional[float] = 3907e6  # Pa: E das juntas e do capeamento onde há contato; igual ao da
                                         # cerâmica (a flexibilidade da junta vem das folgas). None = Materiais.junta.E
    comprimento_corr: float = 0.04  # m (hipótese: cerca de metade do diâmetro dos furos)
    contato_inicial: float = 0.20   # fração da área em contato sem carga
    semente: int = 1                # semente do campo aleatório


# Resistência à compressão efetiva da cerâmica que reproduz a carga máxima
# do prisma a 30 °C com as folgas calibradas (script 08, malha h = 12 mm,
# dz = 20 mm). É 0,43 vez o valor derivado do ensaio do bloco. Explicações
# prováveis: o ensaio da unidade baixa e larga (60 x 140 mm), capeada e
# confinada pelo atrito dos pratos, superestima a resistência do material; as
# asperezas menores que a malha não estão no modelo; a carga a 30 °C ainda
# precisa ser confirmada. Use só como calibração (opção --calibrado do script 07).
CALIBRACAO_30C = {"fc": 2.48e6}

# Cargas máximas dos prismas da UFV (kgf); a de 30 °C deve ser confirmada.
CARGAS_UFV_KGF = {30: 10066.6, 600: 9980.9, 800: 10087.3, 1000: 11662.1}

# Módulo secante de prismas secos de 5 fiadas (BPQ, n = 4 por grupo), MPa
E_PRISMA_SECO_MPA = {"5ID (blocos inteiros)": 1216.0, "5HD (meio-bloco alternado)": 1001.0}

SIGMA = 5.670374419e-8  # Stefan-Boltzmann
