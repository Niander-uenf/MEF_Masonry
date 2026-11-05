"""Comparação com o ANSYS: prisma BCEPQ do projeto de pesquisa da UENF (Sá, 2025).

O modelo ANSYS (Mechanical 2024 R1, arquivos de entrada ds.dat de 06/01/2025) é o
prisma de 5 fiadas de bloco e meio, com o encaixe macho-fêmea real, elementos
SOLID187, todas as superfícies em contato coladas (bonded, MPC), base engastada,
temperatura uniforme no corpo (T_ref = 22 °C) e pressão no ressalto da fiada de
topo, fora do encaixe. Material linear: E = 2801 MPa (aplicado ao material),
ν = 0,182, α = 6e-6 /°C. A análise térmica transitória do mesmo projeto não tem
calor específico, então devolve o regime permanente: o prisma inteiro na
temperatura imposta nas faces (o campo importado no ds.dat é 600,0 °C em todos
os 23 070 nós).

O MEF_Masonry modela o mesmo prisma como um corpo só (as faces coladas fazem do
conjunto um sólido contínuo), de 422,03 × 140 × 301,1 mm (da base ao ressalto
carregado), com os três furos de 84 mm alinhados nas posições médias do ANSYS
(nas fiadas do ANSYS os furos se deslocam até 7 mm) e o encaixe omitido. A
variante "com sulco" retira o sulco de 5 mm × 7 mm da face inferior, que no
ANSYS fica fora do apoio, e a carga vai só no ressalto (fora da projeção do
encaixe). Como o problema é linear, cada malha resolve dois casos unitários
(temperatura e carga) e os seis casos do ANSYS saem por superposição.

Saída: resultados/artigo/benchmark_ansys.json e uma tabela no terminal."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from mef_masonry.parametros import Geometria, Ceramica
from mef_masonry.termico import Materiais
from mef_masonry import malha as mh, mecanico as mc

AQUI = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "resultados", "artigo")

# ---------------------------------------------------------------- dados do ANSYS
L, W, H = 0.42203, 0.140, 0.3011          # m: comprimento, largura, base -> ressalto do topo
FUROS_X = (0.07368, 0.21125, 0.34850)     # m: média das posições nas fiadas ímpares e pares
D_FURO = 0.084
T_REF = 22.0
E, NU, ALFA = 2801e6, 0.182, 6e-6
# sulco da face inferior (fiada 1: meio-bloco 0-142,3 mm e bloco 142,03-422,03 mm):
# retângulo externo do sulco em cada unidade (x0, x1, y0, y1) e largura do sulco
SULCO = [(0.0287, 0.1250, 0.0217, 0.1183), (0.14203 + 0.0215, 0.14203 + 0.2575, 0.0217, 0.1183)]
B_SULCO, P_SULCO = 0.00507, 0.007
# encaixe da fiada de topo (projeção do macho, que fica sem carga): meio-bloco à esquerda, bloco à direita
ENCAIXE_TOPO = [(0.0287, 0.1250, 0.0217, 0.1183), (0.14203 + 0.0215, 0.14203 + 0.2575, 0.0217, 0.1183)]
# casos do relatório (T em °C, força em N) e respostas do ANSYS: deslocamento total máximo (mm)
# e tensão de von Mises máxima (MPa)
CASOS = [
    ("600 °C", 600.0, 97806.38, 1.1151, 55.858),
    ("800 °C", 800.0, 98956.41, 1.5566, 73.889),
    ("1000 °C", 1000.0, 114405.20, 1.9745, 92.489),
    ("carga de 600 °C, sem temperatura", T_REF, 97806.38, 0.30700, 34.694),
    ("carga de 800 °C, sem temperatura", T_REF, 98956.41, 0.31061, 35.102),
    ("carga de 1000 °C, sem temperatura", T_REF, 114405.20, 0.35911, 40.582),
]


def _no_retangulo(x, y, r, borda=0.0):
    x0, x1, y0, y1 = r
    return (x > x0 - borda) & (x < x1 + borda) & (y > y0 - borda) & (y < y1 + borda)


def _na_faixa_sulco(x, y, retangulos=SULCO):
    """Faixa de largura B_SULCO junto ao contorno interno de cada retângulo (sulco ou macho)."""
    s = np.zeros(np.shape(x), bool)
    for r in retangulos:
        x0, x1, y0, y1 = r
        dentro = _no_retangulo(x, y, r)
        miolo = (x > x0 + B_SULCO) & (x < x1 - B_SULCO) & (y > y0 + B_SULCO) & (y < y1 - B_SULCO)
        s |= dentro & ~miolo
    return s


def modelo(h, elemento, sulco):
    geo = Geometria(comprimento_bloco=L / 1.5, largura=W, altura_fiada=H / 5, diametro_furo=D_FURO,
                    n_fiadas=5, blocos_por_fiada=1.5, tipo_junta="seca", espessura_junta=0.002)
    m2 = mh.malha_secao(geo, h=h, simetria=True, faixas=None, furos=[(x, W / 2) for x in FUROS_X],
                        raio=D_FURO / 2)
    cam = [mh.Camada(0.0, P_SULCO, "uniforme", material=mh.CERAMICA),
           mh.Camada(P_SULCO, H, "uniforme", material=mh.CERAMICA)]
    m3 = mh.extrudar(m2, cam, dz=h, geo=geo, tipo=elemento)
    if sulco:
        base2d = m2.quad if elemento == "H8" else m2.tri
        c2 = m2.xy[base2d].mean(1)
        faixa = _na_faixa_sulco(c2[:, 0], c2[:, 1])
        zc = m3.xyz[m3.cun[:, :6] if elemento != "H8" else m3.cun].mean(1)[:, 2]
        tira = (zc < P_SULCO) & faixa[m3.tri_de]
        m3.mat[tira & (m3.mat == mh.CERAMICA)] = mh.AR
    return m3


def resolver(h, elemento, sulco):
    m3 = modelo(h, elemento, sulco)
    M = Materiais(ceramica=Ceramica(E_bloco_bruto=E, razao_area_liquida=1.0, nu=NU, alfa=ALFA, T_ref=T_REF))
    T600 = np.full(m3.n_nos, 600.0)
    regiao = None
    if sulco:   # a pressão do ANSYS cobre todo o ressalto, menos a base do macho (faixa de ~5 mm)
        regiao = lambda x, y: ~_na_faixa_sulco(x, y, ENCAIXE_TOPO)
    t0 = time.time()
    R = mc.resolver_3d(m3, M, T_nos=T600, carga_N=CASOS[0][2], placas="aco", base_engastada=True,
                       regiao_carga=regiao)
    dt = time.time() - t0
    uT, uF = R["termico"].u, R["carga"].u          # 600 °C e 97,8 kN
    sT, sF = R["termico"].sig, R["carga"].sig
    ativos = R["termico"].ativos
    res = []
    for nome, T, F, _, _ in CASOS:
        a, b = (T - T_REF) / (600.0 - T_REF), F / CASOS[0][2]
        u = a * uT + b * uF
        s = (a * sT + b * sF)[ativos]
        vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2)
                     + 3 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
        umax = np.linalg.norm(u, axis=1)
        i = int(np.argmax(umax))
        res.append({"caso": nome, "u_max_mm": float(umax[i] * 1e3), "no": m3.xyz[i].tolist(),
                    "von_mises_max_MPa": float(vm.max() / 1e6)})
    nel = int(ativos.sum())
    return {"h_mm": h * 1e3, "elemento": elemento, "sulco": sulco, "elementos": nel, "nos": int(m3.n_nos),
            "tempo_s": round(dt, 1), "casos": res, "ajuste": ajuste(uT, uF)}


def ajuste(uT, uF):
    """Fatores a (parte térmica) e s (parte da carga) que levam o MEF_Masonry aos seis
    valores do ANSYS: u = a uT (T - T_ref)/578 + s uF F/97,8 kN. Mostra em que parte da
    resposta está a diferença (a = s = 1: mesma resposta)."""
    from scipy.optimize import least_squares
    alvo = np.array([c[3] for c in CASOS])
    kT = np.array([(c[1] - T_REF) / (600.0 - T_REF) for c in CASOS])
    kF = np.array([c[2] / CASOS[0][2] for c in CASOS])

    def umax(p):
        a, s = p
        return np.array([np.linalg.norm(a * t * uT + s * f * uF, axis=1).max() * 1e3 for t, f in zip(kT, kF)])
    r = least_squares(lambda p: umax(p) / alvo - 1, [1.0, 1.0])
    fit = umax(r.x)
    return {"a_termico": float(r.x[0]), "s_carga": float(r.x[1]),
            "u_max_ajustado_mm": fit.tolist(), "erro_max_pct": float(np.abs(fit / alvo - 1).max() * 100)}


if __name__ == "__main__":
    rapido = "--rapido" in sys.argv
    # malhas que cabem em ~2 GB com o SuperLU (sem pypardiso)
    malhas = [(0.010, "H8", False), (0.010, "H8", True)] if rapido else \
        [(0.015, "H8", False), (0.010, "H8", False), (0.0075, "H8", False), (0.015, "W15", False),
         (0.010, "H8", True), (0.0075, "H8", True), (0.015, "W15", True), (0.010, "W15", True)]
    saida = {"ansys": [{"caso": c[0], "T": c[1], "F_N": c[2], "u_max_mm": c[3], "von_mises_max_MPa": c[4]}
                       for c in CASOS], "mef_masonry": []}
    for h, el, sulco in malhas:
        r = resolver(h, el, sulco)
        saida["mef_masonry"].append(r)
        print(f"\n{el} h = {h * 1e3:.1f} mm, sulco = {sulco}: {r['elementos']} elementos, {r['tempo_s']} s")
        for c, a in zip(r["casos"], CASOS):
            d = (c["u_max_mm"] - a[3]) / a[3] * 100
            print(f"  {c['caso']:36s} u_max {c['u_max_mm']:.4f} mm (ANSYS {a[3]:.4f}, {d:+.2f}%)"
                  f"  von Mises máx {c['von_mises_max_MPa']:.1f} MPa (ANSYS {a[4]:.1f})", flush=True)
        aj = r["ajuste"]
        print(f"  ajuste: a = {aj['a_termico']:.4f} (térmica), s = {aj['s_carga']:.4f} (carga), "
              f"erro máximo {aj['erro_max_pct']:.2f}%", flush=True)
    os.makedirs(AQUI, exist_ok=True)
    json.dump(saida, open(os.path.join(AQUI, "benchmark_ansys.json"), "w"), indent=1, ensure_ascii=False)
