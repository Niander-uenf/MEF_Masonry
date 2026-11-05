"""Malhas: seção horizontal 2D do prisma (triângulos e quadriláteros) e prisma
3D por extrusão: cunhas de 6 nós (W6), hexaedros de 8 nós (H8) ou cunhas
quadráticas de 15 nós (W15).

A seção é constante ao longo da altura (furos alinhados entre fiadas), por isso
o prisma 3D é obtido extrudando a malha 2D. A malha 2D é estruturada e os nós
na fronteira dos furos são projetados sobre a circunferência, o que dá furos
lisos e malha conforme sem gerador externo.
"""
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .parametros import Geometria

# Materiais
AR, CERAMICA, JUNTA, CAPEAMENTO, ACO, ARGAMASSA, GRAUTE = 0, 1, 2, 3, 4, 5, 6
NOMES = {AR: "ar (furos)", CERAMICA: "bloco", JUNTA: "junta seca",
         CAPEAMENTO: "capeamento", ACO: "aço (placas)", ARGAMASSA: "argamassa", GRAUTE: "graute"}

# Marcadores dos triângulos da seção
SOLIDO, FURO = -1, 0  # faixas de junta vertical: 1, 2, ...


@dataclass
class Malha2D:
    xy: np.ndarray          # (n, 2) coordenadas, m
    tri: np.ndarray         # (ne, 3) conectividade, sentido anti-horário
    marca: np.ndarray       # (ne,) SOLIDO, FURO ou índice da faixa (1, 2, ...)
    faixas: Tuple[float, ...]
    comprimento: float
    largura_modelo: float   # largura modelada (metade se houver simetria)
    largura: float
    simetria: bool
    quad: Optional[np.ndarray] = None        # (nq, 4) quadriláteros de origem (anti-horários)
    marca_quad: Optional[np.ndarray] = None  # (nq,)
    largura_faixa: float = 0.002             # espessura das juntas verticais (e recuo das células)
    retangular: bool = False                 # células retangulares (furos podem mudar de fiada para fiada)
    vertical: bool = False                   # seção vertical de uma fiada (furos horizontais): a primeira
                                             # coordenada é a altura na fiada, a segunda a espessura

    @property
    def n_nos(self):
        return len(self.xy)

    def areas(self):
        p = self.xy[self.tri]
        return 0.5 * ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                      - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))

    def arestas_contorno(self):
        """Arestas que pertencem a um único triângulo: (a, b)."""
        e = np.vstack([self.tri[:, [0, 1]], self.tri[:, [1, 2]], self.tri[:, [2, 0]]])
        es = np.sort(e, axis=1)
        _, inv, cnt = np.unique(es, axis=0, return_inverse=True, return_counts=True)
        return e[cnt[inv.ravel()] == 1]


def _coordenadas(a, b, h, obrigatorias=()):
    n = max(1, int(np.ceil((b - a) / h - 1e-9)))
    c = np.linspace(a, b, n + 1)
    obr = np.asarray(sorted(obrigatorias), float)
    if obr.size:
        manter = np.ones(c.size, bool)
        for o in obr:
            manter &= np.abs(c - o) > 0.35 * h
        manter[0] = manter[-1] = True
        c = np.unique(np.round(np.concatenate([c[manter], obr]), 12))
    return c


def _quads_em_triangulos(q):
    """q: (nq, 4) nós de quadriláteros; divide alternando a diagonal."""
    alt = (np.arange(len(q)) % 2 == 0)[:, None]
    t1 = np.where(alt, q[:, [0, 1, 2]], q[:, [0, 1, 3]])
    t2 = np.where(alt, q[:, [0, 2, 3]], q[:, [1, 2, 3]])
    return np.vstack([t1, t2])


def _grade(pontos, h):
    """Coordenadas que passam por todos os pontos dados, com passo <= h entre eles."""
    p = np.unique(np.round(np.asarray(pontos, float), 9))
    out = [p[0]]
    for a, b in zip(p[:-1], p[1:]):
        n = max(1, int(np.ceil((b - a) / h - 1e-9)))
        out.extend(a + (b - a) * np.arange(1, n + 1) / n)
    return np.asarray(out)


def malha_secao_retangular(geo: Geometria, h: float = 0.005, simetria: bool = True,
                           faixas: Optional[Sequence[float]] = None,
                           t_junta: Optional[float] = None) -> Malha2D:
    """Seção (planta) de blocos com células retangulares (concreto, cerâmico
    estrutural) ou com furos horizontais (vedação: a grade segue as bordas das
    colunas de furos, e a extrusão abre os furos nas fileiras).

    Grade estruturada alinhada a todas as bordas de células de todas as fiadas
    (com contrafiada as células mudam de posição de uma fiada para a outra) e
    às faixas das juntas verticais. A marca FURO segue a fiada 0 (seção 2D);
    a extrusão classifica cada elemento fiada a fiada (estados_fiada).
    """
    L, W = geo.comprimento, geo.largura
    tj = geo.espessura_junta if t_junta is None else t_junta
    faixas = tuple(geo.posicoes_juntas if faixas is None else faixas)
    px = [0.0, L] + [p for x in faixas for p in (x - tj / 2, x + tj / 2)]
    py = [0.0, W / 2, W] + [v for f in geo.furos_h_y() for v in f]
    for i in range(geo.n_fiadas):
        for c in geo.celulas(i, tj):
            px += [c[0], c[1]]
            py += [c[2], c[3]]
    xs, ys = _grade(px, h), _grade(py, h)
    if simetria:
        ys = ys[ys <= W / 2 + 1e-12]
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nx, ny = X.shape
    xy = np.c_[X.ravel(), Y.ravel()]
    idx = np.arange(nx * ny).reshape(nx, ny)
    i, j = (a.ravel() for a in np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), indexing="ij"))
    n00, n10, n11, n01 = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
    alt = ((i + j) % 2 == 0)[:, None]
    tri = np.vstack([np.where(alt, np.c_[n00, n10, n11], np.c_[n00, n10, n01]),
                     np.where(alt, np.c_[n00, n11, n01], np.c_[n10, n11, n01])])
    quad = np.c_[n00, n10, n11, n01]

    def marcar(pol):
        cen = xy[pol].mean(1)
        mk = np.full(len(pol), SOLIDO)
        for k, xf in enumerate(faixas, start=1):
            mk[np.abs(cen[:, 0] - xf) < tj / 2] = k
        for c in geo.celulas(0, tj):
            mk[(cen[:, 0] > c[0]) & (cen[:, 0] < c[1]) & (cen[:, 1] > c[2]) & (cen[:, 1] < c[3])] = FURO
        return mk
    return Malha2D(xy=xy, tri=tri, marca=marcar(tri), faixas=faixas, comprimento=L,
                   largura_modelo=W / 2 if simetria else W, largura=W, simetria=simetria,
                   quad=quad, marca_quad=marcar(quad), largura_faixa=tj, retangular=True)


def malha_secao_vertical(geo: Geometria, h: float = 0.005, t_junta: Optional[float] = None) -> Malha2D:
    """Seção vertical de uma fiada de blocos de furos horizontais (espessura ×
    altura, de eixo a eixo), para a térmica 2D e o modo Parede. A primeira
    coordenada é a altura na fiada (0 a altura_fiada) e a segunda, a espessura
    (y = 0 é a face exposta). Meia junta embaixo e meia em cima (marcas 1 e 2);
    furos com a marca FURO. Como os furos atravessam o bloco, a seção vale ao
    longo de todo o comprimento da parede."""
    W, hf = geo.largura, geo.altura_fiada
    tj = geo.espessura_junta if t_junta is None else t_junta
    linhas_z = [0.0, tj / 2, hf - tj / 2, hf] + [v for f in geo.furos_h_z(hf / 2, tj) for v in f]
    ys = _grade([0.0, W / 2, W] + [v for f in geo.furos_h_y() for v in f], h)
    zs = _grade(linhas_z, h)
    Z, Y = np.meshgrid(zs, ys, indexing="ij")
    nz, ny = Z.shape
    xy = np.c_[Z.ravel(), Y.ravel()]
    idx = np.arange(nz * ny).reshape(nz, ny)
    i, j = (a.ravel() for a in np.meshgrid(np.arange(nz - 1), np.arange(ny - 1), indexing="ij"))
    n00, n10, n11, n01 = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
    alt = ((i + j) % 2 == 0)[:, None]
    tri = np.vstack([np.where(alt, np.c_[n00, n10, n11], np.c_[n00, n10, n01]),
                     np.where(alt, np.c_[n00, n11, n01], np.c_[n10, n11, n01])])
    quad = np.c_[n00, n10, n11, n01]

    def marcar(pol):
        cen = xy[pol].mean(1)
        mk = np.full(len(pol), SOLIDO)
        mk[cen[:, 0] < tj / 2] = 1
        mk[cen[:, 0] > hf - tj / 2] = 2
        for z0, z1 in geo.furos_h_z(hf / 2, tj):
            for y0, y1 in geo.furos_h_y():
                mk[(cen[:, 0] > z0) & (cen[:, 0] < z1) & (cen[:, 1] > y0) & (cen[:, 1] < y1)] = FURO
        return mk
    return Malha2D(xy=xy, tri=tri, marca=marcar(tri), faixas=(0.0, hf), comprimento=hf, largura_modelo=W,
                   largura=W, simetria=False, quad=quad, marca_quad=marcar(quad), largura_faixa=tj, vertical=True)


def estados_fiada(m2d: Malha2D, geo: Geometria, i: int, quad: bool = False):
    """Furo e graute por elemento 2D (triângulo, ou quadrilátero se quad) na fiada i.

    Furos circulares: alinhados em todas as fiadas (marca FURO); grauteados os
    furos das colunas escolhidas. Células retangulares: classificadas pelas
    células da fiada i; uma célula principal é grauteada se cobre mais da
    metade de uma coluna grauteada (definida na fiada 0)."""
    pol = m2d.quad if quad else m2d.tri
    marca = m2d.marca_quad if quad else m2d.marca
    cen = m2d.xy[pol].mean(1)
    cols = geo.colunas_grauteadas(m2d.largura_faixa)
    if m2d.retangular and geo.retangular:
        furo = np.zeros(len(pol), bool)
        graute = np.zeros(len(pol), bool)
        for c in geo.celulas(i, m2d.largura_faixa):
            dentro = (cen[:, 0] > c[0]) & (cen[:, 0] < c[1]) & (cen[:, 1] > c[2]) & (cen[:, 1] < c[3])
            furo |= dentro
            if c[4] and any(min(c[1], b) - max(c[0], a) > 0.5 * (c[1] - c[0]) for a, b in cols):
                graute |= dentro
        return furo, graute
    furo = marca == FURO
    graute = np.zeros(len(pol), bool)
    for a, b in cols:
        graute |= furo & (cen[:, 0] > a - 1e-9) & (cen[:, 0] < b + 1e-9)
    return furo, graute


def malha_secao(geo: Geometria = Geometria(), h: float = 0.005, simetria: bool = True,
                faixas: Optional[Sequence[float]] = (0.14, 0.28), largura_faixa: float = 0.002,
                furos: Optional[Sequence[Tuple[float, float]]] = None,
                raio: Optional[float] = None) -> Malha2D:
    """Malha da seção horizontal do prisma de bloco e meio (420 x 140 mm).

    Fora dos furos a malha é retangular estruturada; em volta de cada furo há
    uma caixa quadrada com malha em "O" (anéis concêntricos), que segue a
    circunferência exatamente. O interior do furo (ar) também é malhado,
    para a análise térmica.

    h        : tamanho de elemento (m)
    simetria : modela só y <= W/2 (plano de simetria pelos centros dos furos)
    faixas   : posições x das juntas verticais (faixas finas onde a junta pode existir)
    furos    : centros (x, y) dos furos, todos com o mesmo y; padrão = BCEPQ
    raio     : raio dos furos; padrão = diametro_furo/2 (0 = seção maciça)
    """
    if geo.horizontal and furos is None and raio is None:
        return malha_secao_vertical(geo, h)
    if geo.retangular and furos is None and raio is None:
        return malha_secao_retangular(geo, h, simetria, faixas or (),
                                      largura_faixa if faixas else geo.espessura_junta)
    L, W = geo.comprimento, geo.largura
    R = geo.diametro_furo / 2 if raio is None else raio
    if furos is None:
        furos = [(x, W / 2) for x in geo.centros_furos]
    furos = sorted(furos) if R > 0 else []
    faixas = tuple(faixas or ())
    bordas = [p for x in faixas for p in (x - largura_faixa / 2, x + largura_faixa / 2)]

    if furos:
        cy = furos[0][1]
        if any(abs(f[1] - cy) > 1e-12 for f in furos):
            raise ValueError("Os furos devem estar alinhados (mesmo y).")
        xs_f = [f[0] for f in furos]
        folga = [min(cy, W - cy), xs_f[0], L - xs_f[-1]]
        if len(xs_f) > 1:
            folga.append(float(np.min(np.diff(xs_f))) / 2)
        a = R + (min(folga) - R) / 2               # meia-largura da caixa em O
        n = max(8, 2 * int(round(a / h)))          # segmentos por lado (par)
        lado = lambda c: c - a + 2 * a * np.arange(n + 1) / n
    # coordenadas x e y da grade retangular
    partes, ini = [], 0.0
    for cx, _ in furos:
        partes += [_coordenadas(ini, cx - a, h, [b for b in bordas if ini < b < cx - a]), lado(cx)]
        ini = cx + a
    partes.append(_coordenadas(ini, L, h, [b for b in bordas if ini < b < L]))
    xs = np.unique(np.round(np.concatenate(partes), 12))
    if furos:
        ys = np.unique(np.round(np.concatenate(
            [_coordenadas(0.0, cy - a, h), lado(cy), _coordenadas(cy + a, W, h)]), 12))
    else:
        ys = _coordenadas(0.0, W, h, [W / 2])

    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nx, ny = X.shape
    pts = [np.c_[X.ravel(), Y.ravel()]]
    nn = nx * ny
    idx = np.arange(nn).reshape(nx, ny)
    i, j = (v.ravel() for v in np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), indexing="ij"))
    cxr = 0.5 * (xs[i] + xs[i + 1]); cyr = 0.5 * (ys[j] + ys[j + 1])
    fora = np.ones(len(i), bool)
    for cx, _ in furos:
        fora &= ~((np.abs(cxr - cx) < a) & (np.abs(cyr - cy) < a))
    i, j = i[fora], j[fora]
    q = np.c_[idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]]
    tris = [_quads_em_triangulos(q)]
    marcas = [np.full(len(tris[0]), SOLIDO)]
    quads, marcas_q = [q], [np.full(len(q), SOLIDO)]

    for cx, _ in furos:
        ix0 = int(np.argmin(np.abs(xs - (cx - a)))); iy0 = int(np.argmin(np.abs(ys - (cy - a))))
        # contorno da caixa, anti-horário a partir do canto inferior esquerdo
        per = np.concatenate([idx[ix0:ix0 + n, iy0], idx[ix0 + n, iy0:iy0 + n],
                              idx[ix0 + n:ix0:-1, iy0 + n], idx[ix0, iy0 + n:iy0:-1]])
        Q = pts[0][per]
        c = np.array([cx, cy])
        d = Q - c
        C = c + R * d / np.linalg.norm(d, axis=1)[:, None]
        b = 0.45 * R
        I = c + d * (b / a)
        m1 = max(2, int(round(((a - R) + (a * np.sqrt(2) - R)) / 2 / h)))
        m2 = max(2, int(round((R - b * 1.2) / h)))
        # anéis: interior (I) -> circunferência (C) -> caixa (Q)
        aneis = [I + (C - I) * s for s in np.arange(m2) / m2] + \
                [C + (Q - C) * s for s in np.arange(m1) / m1]
        ids_aneis = []
        for anel in aneis:
            ids_aneis.append(nn + np.arange(len(anel))); pts.append(anel); nn += len(anel)
        ids_aneis.append(per)
        for k in range(len(ids_aneis) - 1):
            p0, p1 = ids_aneis[k], ids_aneis[k + 1]
            qq = np.c_[p0, np.roll(p0, -1), np.roll(p1, -1), p1]
            t = _quads_em_triangulos(qq)
            tris.append(t)
            marcas.append(np.full(len(t), FURO if k < m2 else SOLIDO))
            quads.append(qq); marcas_q.append(np.full(len(qq), FURO if k < m2 else SOLIDO))
        # quadrado interno estruturado (ar)
        g = -np.ones((n + 1, n + 1), int)
        anel0 = ids_aneis[0]
        pos = np.concatenate([np.c_[np.arange(n), np.zeros(n, int)],
                              np.c_[np.full(n, n), np.arange(n)],
                              np.c_[np.arange(n, 0, -1), np.full(n, n)],
                              np.c_[np.zeros(n, int), np.arange(n, 0, -1)]])
        g[pos[:, 0], pos[:, 1]] = anel0
        ii, jj = np.nonzero(g < 0)
        novos = c + np.c_[-b + 2 * b * ii / n, -b + 2 * b * jj / n]
        g[ii, jj] = nn + np.arange(len(ii)); pts.append(novos); nn += len(ii)
        a_, b_ = (v.ravel() for v in np.meshgrid(np.arange(n), np.arange(n), indexing="ij"))
        qq = np.c_[g[a_, b_], g[a_ + 1, b_], g[a_ + 1, b_ + 1], g[a_, b_ + 1]]
        t = _quads_em_triangulos(qq)
        tris.append(t); marcas.append(np.full(len(t), FURO))
        quads.append(qq); marcas_q.append(np.full(len(qq), FURO))

    xy = np.vstack(pts)
    tri = np.vstack(tris)
    marca = np.concatenate(marcas)
    quad = np.vstack(quads)
    marca_q = np.concatenate(marcas_q)
    cen = xy[tri].mean(axis=1)
    cenq = xy[quad].mean(axis=1)
    for k, xf in enumerate(faixas, start=1):
        marca[(marca == SOLIDO) & (np.abs(cen[:, 0] - xf) < largura_faixa / 2)] = k
        marca_q[(marca_q == SOLIDO) & (np.abs(cenq[:, 0] - xf) < largura_faixa / 2)] = k
    if simetria:
        manter = cen[:, 1] < W / 2
        tri, marca = tri[manter], marca[manter]
        mq = cenq[:, 1] < W / 2
        quad, marca_q = quad[mq], marca_q[mq]
    # remove nós órfãos e renumera
    usados = np.unique(tri)
    novo = -np.ones(len(xy), int); novo[usados] = np.arange(len(usados))
    xy, tri, quad = xy[usados], novo[tri], novo[quad]
    # orientação anti-horária
    p = xy[tri]
    ar = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    tri[ar < 0] = tri[ar < 0][:, [0, 2, 1]]
    pq = xy[quad]
    arq = 0.5 * np.sum(pq[:, :, 0] * np.roll(pq[:, :, 1], -1, axis=1) - np.roll(pq[:, :, 0], -1, axis=1) * pq[:, :, 1], axis=1)
    quad[arq < 0] = quad[arq < 0][:, ::-1]
    m = Malha2D(xy=xy, tri=tri, marca=marca, faixas=faixas, comprimento=L,
                largura_modelo=W / 2 if simetria else W, largura=W, simetria=simetria,
                quad=quad, marca_quad=marca_q, largura_faixa=largura_faixa)
    if (np.abs(m.areas()) < 1e-14).any():
        raise RuntimeError("Malha 2D com elementos degenerados.")
    return m


def malha_retangulo(Lx: float, Ly: float, h: float, hy: Optional[float] = None) -> Malha2D:
    """Retângulo maciço (usado nas verificações); hy permite malha anisotrópica."""
    xs = _coordenadas(0.0, Lx, h)
    ys = _coordenadas(0.0, Ly, hy or h)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    nx, ny = X.shape
    xy = np.c_[X.ravel(), Y.ravel()]
    idx = np.arange(nx * ny).reshape(nx, ny)
    i, j = (a.ravel() for a in np.meshgrid(np.arange(nx - 1), np.arange(ny - 1), indexing="ij"))
    n00, n10, n11, n01 = idx[i, j], idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]
    alt = ((i + j) % 2 == 0)[:, None]
    t1 = np.where(alt, np.c_[n00, n10, n11], np.c_[n00, n10, n01])
    t2 = np.where(alt, np.c_[n00, n11, n01], np.c_[n10, n11, n01])
    tri = np.vstack([t1, t2])
    quad = np.c_[n00, n10, n11, n01]
    return Malha2D(xy=xy, tri=tri, marca=np.full(len(tri), SOLIDO), faixas=(),
                   comprimento=Lx, largura_modelo=Ly, largura=Ly, simetria=False,
                   quad=quad, marca_quad=np.full(len(quad), SOLIDO))


# ----------------------------------------------------------------------------
# Extrusão 3D
# ----------------------------------------------------------------------------
@dataclass
class Camada:
    z0: float
    z1: float
    tipo: str            # 'placa', 'capeamento', 'junta', 'unidade'
    fiada: int = -1      # 'unidade': a fiada; 'junta': a fiada de baixo; 'capeamento': a fiada vizinha
    furo_h: int = -1     # 'unidade' de bloco com furos horizontais: fileira de furos nesta subcamada (-1 = parede)
    material: int = CERAMICA  # usado quando tipo == 'uniforme'


@dataclass
class Malha3D:
    xyz: np.ndarray        # (n, 3)
    cun: np.ndarray        # (ne, 6): 3 nós de baixo + 3 de cima, mesma ordem
    mat: np.ndarray        # (ne,) material
    camada: np.ndarray     # (ne,) índice da subcamada
    tri_de: np.ndarray     # (ne,) triângulo 2D de origem (quadrilátero, no H8)
    zs: np.ndarray         # cotas dos planos de nós
    m2d: Malha2D
    camadas: List[Camada] = field(default_factory=list)
    tipo: str = "W6"       # elemento: W6, H8 ou W15

    @property
    def n_nos(self):
        return len(self.xyz)

    @property
    def elem(self):
        from .elementos3d import tipo
        return tipo(self.tipo)

    @property
    def n2(self):
        return self.m2d.n_nos


def camadas_prisma(geo: Geometria = Geometria(), t_junta: float = 0.002,
                   capeamento: bool = True, placas: bool = True,
                   t_capeamento: Optional[float] = None) -> List[Camada]:
    """Camadas do prisma: placa inferior, capeamento, fiadas com juntas
    entre elas, capeamento e placa superior. As juntas ficam centradas nas
    interfaces entre fiadas (z = 60, 120, ... mm). O capeamento tem a espessura
    da junta, salvo t_capeamento."""
    hf, n = geo.altura_fiada, geo.n_fiadas
    H = hf * n
    tc = t_junta if t_capeamento is None else t_capeamento
    cam: List[Camada] = []
    if placas:
        cam.append(Camada(-geo.espessura_placa, 0.0, "placa"))
    z = 0.0
    if capeamento:
        cam.append(Camada(0.0, tc, "capeamento", fiada=0))
        z = tc
    for i in range(n):
        topo = (i + 1) * hf - (t_junta / 2 if i < n - 1 else (tc if capeamento else 0.0))
        # blocos de furos horizontais: a unidade se divide nas bordas das fileiras de furos
        cortes = [(z0, z1, k) for k, (z0, z1) in enumerate(geo.furos_h_z((i + 0.5) * hf, t_junta))]
        for z0, z1, k in cortes:
            if z0 > z + 1e-9:
                cam.append(Camada(z, z0, "unidade", fiada=i))
            cam.append(Camada(z0, z1, "unidade", fiada=i, furo_h=k))
            z = z1
        cam.append(Camada(z, topo, "unidade", fiada=i))
        z = topo
        if i < n - 1:
            cam.append(Camada(z, z + t_junta, "junta", fiada=i))
            z += t_junta
    if capeamento:
        cam.append(Camada(z, H, "capeamento", fiada=n - 1))
    if placas:
        cam.append(Camada(H, H + geo.espessura_placa, "placa"))
    return cam


def _arestas_tri(tri):
    """Arestas únicas dos triângulos e, por triângulo, os índices das arestas (0-1, 1-2, 2-0)."""
    e = np.stack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]], axis=1)      # (nt, 3, 2)
    es = np.sort(e, axis=2).reshape(-1, 2)
    unicas, inv = np.unique(es, axis=0, return_inverse=True)
    return unicas, inv.reshape(-1, 3)


def extrudar(m2d: Malha2D, camadas: Sequence[Camada], dz: float = 0.010,
             geo: Geometria = Geometria(), tipo: str = "W6") -> Malha3D:
    """Extruda a seção nas camadas dadas (subdividindo cada uma com passo <= dz).
    tipo: 'W6' (triângulos -> cunhas), 'H8' (quadriláteros -> hexaedros) ou
    'W15' (triângulos -> cunhas quadráticas, com nós no meio das arestas)."""
    zs = [camadas[0].z0]
    info = []  # (camada, índice)
    for ic, c in enumerate(camadas):
        nsub = max(1, int(np.ceil((c.z1 - c.z0) / dz - 1e-9)))
        for s in range(nsub):
            zs.append(c.z0 + (c.z1 - c.z0) * (s + 1) / nsub)
            info.append(ic)
    zs = np.asarray(zs)
    n2 = m2d.n_nos
    nz = len(zs) - 1
    if tipo == "H8":
        if m2d.quad is None:
            raise ValueError("A seção não tem quadriláteros para o H8.")
        base2d, marca2d = m2d.quad, m2d.marca_quad
    else:
        base2d, marca2d = m2d.tri, m2d.marca
    nt = len(base2d)
    k = np.repeat(np.arange(nz), nt)
    t = np.tile(base2d, (nz, 1))
    if tipo == "W15":
        arestas, eid = _arestas_tri(m2d.tri)
        ne2 = len(arestas)
        n2q = n2 + ne2
        xy2 = np.vstack([m2d.xy, 0.5 * (m2d.xy[arestas[:, 0]] + m2d.xy[arestas[:, 1]])])
        zm = 0.5 * (zs[:-1] + zs[1:])
        xyz = np.vstack([np.c_[np.tile(xy2, (nz + 1, 1)), np.repeat(zs, n2q)],
                         np.c_[np.tile(m2d.xy, (nz, 1)), np.repeat(zm, n2)]])
        e = np.tile(eid, (nz, 1))
        kb, kt = (k * n2q)[:, None], ((k + 1) * n2q)[:, None]
        km = ((nz + 1) * n2q + k * n2)[:, None]
        cun = np.c_[t + kb, t + kt, n2 + e + kb, n2 + e + kt, t + km]
    else:
        xyz = np.c_[np.tile(m2d.xy, (nz + 1, 1)), np.repeat(zs, n2)]
        cun = np.c_[t + (k * n2)[:, None], t + ((k + 1) * n2)[:, None]]
    mat = np.empty(nz * nt, dtype=np.int8)
    juntas = geo.juntas_por_fiada
    argamassa = geo.tipo_junta == "argamassa"
    mat_vertical = (ARGAMASSA if geo.juntas_verticais_preenchidas else AR) if argamassa else JUNTA
    # assentamento lateral: argamassa só nas faixas das paredes externas (|y - W/2| >= W/2 - casca)
    yc = m2d.xy[base2d].mean(1)[:, 1]
    leito = np.ones(len(base2d), bool)
    if argamassa and geo.assentamento == "lateral" and geo.tem_furos:
        leito = np.abs(yc - geo.largura / 2) >= geo.largura / 2 - geo.espessura_casca - 1e-9
    # estado de cada elemento 2D em cada fiada: 0 sólido, 1 furo (ar), 2 furo grauteado
    nf = max(1, geo.n_fiadas)
    estados = []
    for i in range(nf):
        f, g = estados_fiada(m2d, geo, i, quad=(tipo == "H8"))
        estados.append(np.where(g, 2, np.where(f, 1, 0)))

    def est(i):
        return estados[min(max(i, 0), nf - 1)]
    for kk in range(nz):
        c = camadas[info[kk]]
        sl = slice(kk * nt, (kk + 1) * nt)
        mk = marca2d
        if c.tipo == "placa":
            mat[sl] = ACO
        elif c.tipo == "uniforme":
            mat[sl] = np.where(mk == FURO, AR, c.material)
        elif c.tipo == "capeamento":
            e = est(c.fiada)
            mat[sl] = np.where(e == 1, AR, CAPEAMENTO)
        elif c.tipo == "junta":
            # junta entre as fiadas i e i+1: graute onde há graute e nenhum furo vazio
            # dos dois lados; argamassa (ou junta seca) onde é sólido nos dois lados
            e0, e1 = est(c.fiada), est(c.fiada + 1 if c.fiada >= 0 else 0)
            vazio = (e0 == 1) | (e1 == 1)
            gr = ~vazio & ((e0 == 2) | (e1 == 2))
            sol = (e0 == 0) & (e1 == 0)
            base = np.where(sol & leito, ARGAMASSA, AR) if argamassa else np.where(sol, JUNTA, AR)
            mat[sl] = np.where(gr, GRAUTE, base)
        else:  # unidade
            e = est(c.fiada)
            m = np.where(e == 2, GRAUTE, np.where(e == 1, AR, CERAMICA))
            faixa = np.zeros(len(base2d), bool)            # juntas verticais desta fiada
            if m2d.faixas and 0 <= c.fiada < len(juntas):
                for xj in juntas[c.fiada]:
                    kf = int(np.argmin([abs(xj - f) for f in m2d.faixas])) + 1
                    faixa |= mk == kf
                m = np.where(faixa & (e == 0), mat_vertical, m)
            if c.furo_h >= 0:                              # fileira de furos horizontais
                col = np.zeros(len(base2d), bool)
                for y0, y1 in geo.furos_h_y():
                    col |= (yc > y0) & (yc < y1)
                m = np.where(col & ~faixa, AR, m)
            mat[sl] = m
    return Malha3D(xyz=xyz, cun=cun, mat=mat, camada=np.repeat(np.asarray(info), nt),
                   tri_de=np.tile(np.arange(nt), nz), zs=zs, m2d=m2d, camadas=list(camadas), tipo=tipo)


def malha_prisma(geo: Geometria = Geometria(), h: float = 0.006, dz: float = 0.010,
                 t_junta: float = 0.002, simetria: bool = True, placas: bool = True,
                 capeamento: bool = True, elemento: str = "W6") -> Malha3D:
    """Atalho: malha 3D do prisma (padrão: 5 fiadas de bloco e meio, juntas secas)."""
    problemas = geo.validar()
    if problemas:
        raise ValueError(" ".join(problemas))
    if geo.retangular or geo.horizontal:
        m2d = malha_secao_retangular(geo, h=h, simetria=simetria, faixas=geo.posicoes_juntas, t_junta=t_junta)
    else:
        m2d = malha_secao(geo, h=h, simetria=simetria, largura_faixa=t_junta, faixas=geo.posicoes_juntas)
    t_cap = t_junta if geo.tipo_junta == "seca" else min(t_junta, 0.003)
    cam = camadas_prisma(geo, t_junta=t_junta, capeamento=capeamento, placas=placas, t_capeamento=t_cap)
    return extrudar(m2d, cam, dz=dz, geo=geo, tipo=elemento)


def no_mais_proximo(xyz: np.ndarray, ponto) -> int:
    return int(np.argmin(np.linalg.norm(xyz - np.asarray(ponto)[None, :], axis=1)))


def resumo(m) -> str:
    if isinstance(m, Malha3D):
        cont = {NOMES[k]: int((m.mat == k).sum()) for k in NOMES if (m.mat == k).any()}
        nome = {"W6": "cunhas W6", "H8": "hexaedros H8", "W15": "cunhas W15"}.get(m.tipo, m.tipo)
        return (f"{m.n_nos} nós, {len(m.cun)} {nome}, {len(m.zs) - 1} camadas; "
                + ", ".join(f"{k}: {v}" for k, v in cont.items()))
    a = m.areas()
    return f"{m.n_nos} nós, {len(m.tri)} triângulos; área mín/med = {a.min() / np.median(a):.2f}"
