"""Verificação da fase 3: blocos com furos horizontais (vedação).

V1  Volume dos furos no prisma 3D (bloco e meio em contrafiada, 3 fiadas,
    juntas de argamassa): furos × comprimento de bloco de cada fiada (as
    juntas verticais fecham os furos), exato.
V2  Seção horizontal mínima: numa fileira de furos, a área sólida é a área
    bruta menos as colunas de furos (fora das juntas verticais), exata.
V3  Seção vertical (2D): área dos furos e das meias juntas, exatas.
V4  Térmica: a seção vertical 2D reproduz o 3D extrudado ao longo dos furos
    (uma fiada aquecida pelo ISO 834 numa face, 60 min; mesmas propriedades
    no bloco e na junta para as duas malhas terem o mesmo problema).
V5  Parede: as fibras ficam nas paredes contínuas (faces e septos verticais),
    com área por metro igual à largura menos as colunas de furos, e centro de
    área no plano médio.
V6  Radiação nos furos: a dimensão usada é a largura do furo na direção do
    fluxo de calor.
V7  Prisma não linear do bloco de 8 furos a 20 °C: roda até o pico, e o pico
    é o das paredes que descem contínuas (faces e septo vertical), fc × área
    delas, com tolerância de 5% (os septos horizontais quase não carregam).
    Só com --completo (cerca de 1 min).
V8  Biblioteca: os blocos de furos horizontais passam na validação; o tijolo
    de solo-cimento (diâmetro dos furos ainda sem fonte) é recusado.

Uso: python exemplos/13_verificacao_furos_horizontais.py [--completo]
Grava resultados/verificacao_furos_horizontais.csv.
"""
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import malha as ml, termico as tm, projeto as pj, biblioteca as bib, curvas
from mef_masonry.malha import AR, ARGAMASSA, CERAMICA
from mef_masonry.parametros import Geometria
from mef_masonry.termico import Materiais, Contorno
from mef_masonry.parede import SecaoFibras
from mef_masonry.saida import salvar_csv

COMPLETO = "--completo" in sys.argv
linhas = []


def registrar(caso, modelo, ref, tol_rel, unidade="", tol_abs=0.0):
    erro = abs(modelo - ref) / max(abs(ref), 1e-30)
    ok = erro <= tol_rel or abs(modelo - ref) <= tol_abs
    tol = f"{100 * tol_rel:.1f}%" if tol_abs <= 0 else f"{tol_abs:g} {unidade}"
    linhas.append([caso, round(modelo, 6), round(ref, 6), unidade, f"{100 * erro:.3f}%", tol,
                   "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: modelo {modelo:.6g} {unidade}, referência {ref:.6g} ({100 * erro:.3f}%)")


def condicao(caso, ok, a, b, unidade, criterio):
    linhas.append([caso, round(a, 6), round(b, 6), unidade, "—", criterio, "ok" if ok else "FALHOU"])
    print(f"{'ok ' if ok else 'XX '} {caso}: {a:.4g} e {b:.4g} {unidade} ({criterio})")


def volumes(m):
    xyz = m.xyz[m.cun]
    p = xyz[:, :3, :2]
    a = 0.5 * np.abs((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
                     - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1]))
    return a * (xyz[:, 3, 2] - xyz[:, 0, 2])


tj = 0.010
g = Geometria(forma_furo="horizontal", largura=0.090, altura_fiada=0.200, comprimento_bloco=0.200,
              espessura_junta=tj, tipo_junta="argamassa", parede_externa=0.008, septo=0.007, parede_leito=0.008,
              furos_largura=2, furos_altura=4, blocos_por_fiada=1.5, amarracao="contrafiada", n_fiadas=3)
assert not g.validar(), g.validar()
wy, hz = g.largura_furo_h, g.altura_furo_h(tj)

# ---------------------------------------------------------------------------- V1 e V2
m = ml.malha_prisma(g, h=0.008, dz=0.010, t_junta=tj, simetria=False, placas=False, capeamento=True)
vol = volumes(m)
ref = sum(g.furos_largura * g.furos_altura * wy * hz * (g.comprimento - tj * len(j)) for j in g.juntas_por_fiada)
registrar("V1 volume dos furos horizontais no prisma 3D", vol[m.mat == AR].sum() * 1e6, ref * 1e6, 1e-9, "cm³")

c_fileira = [k for k, c in enumerate(m.camadas) if c.tipo == "unidade" and c.furo_h == 1 and c.fiada == 1][0]
sel = m.camada == c_fileira
dz = (m.xyz[m.cun[sel][:, 3], 2] - m.xyz[m.cun[sel][:, 0], 2])
nsub = len(np.unique(np.round(m.xyz[m.cun[sel][:, 0], 2], 9)))
a_sol = (vol[sel][m.mat[sel] != AR] / dz[m.mat[sel] != AR]).sum() / nsub
ref2 = g.comprimento * g.largura - g.furos_largura * wy * (g.comprimento - tj * len(g.juntas_por_fiada[1]))
registrar("V2 área sólida numa fileira de furos (fiada 2)", a_sol * 1e4, ref2 * 1e4, 1e-9, "cm²")

# ---------------------------------------------------------------------------- V3
m2 = ml.malha_secao(g, h=0.004, simetria=False, faixas=())
a2 = np.abs(m2.areas())
mat2 = tm.material_secao(m2, g)
registrar("V3 área dos furos na seção vertical", a2[mat2 == AR].sum() * 1e4,
          g.furos_largura * g.furos_altura * wy * hz * 1e4, 1e-9, "cm²")
registrar("V3 área das meias juntas na seção vertical", a2[mat2 == ARGAMASSA].sum() * 1e4, g.largura * tj * 1e4,
          1e-9, "cm²")

# ---------------------------------------------------------------------------- V4
g4 = replace(g, comprimento_bloco=0.040, blocos_por_fiada=1, amarracao="prumo", n_fiadas=1)
M4 = Materiais()
M4.ceramica.k, M4.ceramica.rho, M4.ceramica.c = 0.90, 1600.0, 920.0
M4.argamassa.k, M4.argamassa.rho, M4.argamassa.c = 0.90, 1600.0, 920.0
M4.diametro_furo = g4.dimensao_furo
Ti, tf, h4 = 20.0, 3600.0, 0.004
f_iso = curvas.iso834(Ti)
m2v = ml.malha_secao(g4, h=h4, simetria=False, faixas=())
ctr2 = [Contorno(tm.comprimentos_contorno_2d(m2v, faces=("y0",)), f_iso, 25.0, 0.8),
        Contorno(tm.comprimentos_contorno_2d(m2v, faces=("y1",)), curvas.constante(Ti), 4.0, 0.8)]
face2 = np.nonzero(tm.comprimentos_contorno_2d(m2v, faces=("y1",)) > 0)[0]
r2 = tm.resolver_2d(m2v, M4, ctr2, tf, dt=15.0, T_inicial=Ti, mat=tm.material_secao(m2v, g4))
m3 = ml.malha_prisma(g4, h=h4, dz=h4, t_junta=tj, simetria=False, placas=False, capeamento=False)
ctr3 = [Contorno(tm.areas_laterais_3d(m3, faces=("y0",)), f_iso, 25.0, 0.8),
        Contorno(tm.areas_laterais_3d(m3, faces=("y1",)), curvas.constante(Ti), 4.0, 0.8)]
r3 = tm.resolver_3d(m3, M4, ctr3, tf, dt=15.0, T_inicial=Ti)
face3 = np.nonzero(np.isclose(m3.xyz[:, 1], g4.largura))[0]
registrar("V4 face não exposta aos 60 min: seção vertical 2D × 3D", float(r2.T_final[face2].mean()),
          float(r3.T_final[face3].mean()), 0.0, "°C", tol_abs=1.0)
meio2 = np.nonzero(np.isclose(m2v.xy[:, 1], g4.largura / 2))[0]
meio3 = np.nonzero(np.isclose(m3.xyz[:, 1], g4.largura / 2))[0]
registrar("V4 septo vertical (meio da espessura) aos 60 min: 2D × 3D", float(r2.T_final[meio2].mean()),
          float(r3.T_final[meio3].mean()), 0.0, "°C", tol_abs=1.0)

# ---------------------------------------------------------------------------- V5
sec = SecaoFibras.da_malha(m2, g)
registrar("V5 área das fibras por metro de parede = largura - colunas de furos", sec.area * 1e3,
          (g.largura - g.furos_largura * wy) * 1e3, 1e-9, "mm²/mm")
cy = float((sec.A * sec.y).sum() / sec.area)
condicao("V5 centro de área das fibras no plano médio", abs(cy) < 1e-9, cy * 1e3, 0.0, "mm", "|y| < 1e-6 mm")
cen = m2.xy[m2.tri[sec.tri]].mean(1)
dentro = sum(((cen[:, 1] > y0) & (cen[:, 1] < y1)).sum() for y0, y1 in g.furos_h_y())
condicao("V5 nenhuma fibra nas colunas de furos", dentro == 0, float(dentro), 0.0, "fibras", "zero")

# ---------------------------------------------------------------------------- V6
p6 = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", "vedacao_8_furos"))
geo6 = pj.geometria(p6)
registrar("V6 dimensão do furo na radiação = largura do furo", pj.materiais(p6, geo6).diametro_furo * 1e3,
          geo6.largura_furo_h * 1e3, 1e-12, "mm")

# ---------------------------------------------------------------------------- V7
if COMPLETO:
    from mef_masonry.naolinear import ModeloNaoLinear
    q = pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", "vedacao_8_furos"))
    q["modelo"]["analise"] = "nao_linear"
    q["ensaio"].update({"tipo": "ambiente", "T_inicial": 20})
    q["malha"].update({"h": 12, "dz": 20})
    geo = pj.geometria(q)
    M = pj.materiais(q, geo)
    m7 = pj.malha_3d(q, geo)
    mod = ModeloNaoLinear(m7, M, pj.dano(q), pj.junta_atrito(q), verboso=False, dano_argamassa=pj.dano_arg(q),
                          dano_graute=pj.dano_gra(q))
    mod.aquecer({0.0: np.full(m7.n_nos, M.ceramica.T_ref)})
    d, P = mod.carregar(0.006 * geo.altura_prisma, 60)
    i = int(np.argmax(P))
    A_paredes = geo.comprimento * (geo.largura - geo.furos_largura * geo.largura_furo_h)
    teto = q["modelo"]["fc"] * 1e6 * A_paredes
    condicao("V7 prisma de 8 furos: pico alcançado (a carga cai depois)", i < len(P) - 1, P[i] / 1e3, P[-1] / 1e3,
             "kN", "pico antes do fim")
    registrar("V7 pico do prisma de 8 furos = fc × área das paredes contínuas", P[i] / 1e3, teto / 1e3, 0.05, "kN")

# ---------------------------------------------------------------------------- V8
ok8 = []
for k in ("vedacao_8_furos", "vedacao_6_furos"):
    ok8.append(not pj.validar(pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", k))))
recusa = bool(pj.validar(pj.completar(bib.aplicar(pj.projeto_padrao(), "bloco", "solo_cimento_modular"))))
registrar("V8 blocos de furos horizontais que passam na validação", sum(ok8), 2, 0.0, "itens")
condicao("V8 tijolo de solo-cimento (sem diâmetro dos furos) recusado", recusa, float(recusa), 1.0, "", "recusado")

saida = os.path.join(os.path.dirname(__file__), "..", "resultados", "verificacao_furos_horizontais.csv")
salvar_csv(saida, linhas, ["caso", "modelo", "referência", "unidade", "erro", "tolerância", "situação"])
nf = sum(1 for l in linhas if l[-1] == "FALHOU")
print(f"\n{len(linhas) - nf}/{len(linhas)} casos dentro da tolerância. Gravado em {os.path.abspath(saida)}")
sys.exit(1 if nf else 0)
