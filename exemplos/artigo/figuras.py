"""Figures 1-3 of the paper, drawn from the results in this repository.

Fig. 1: workflow of the program (diagram).
Fig. 2: Morris screening of the BCEPQ wall, from resultados/sensibilidade_parede_bcepq/sensibilidade_*.csv.
Fig. 3: design charts, from resultados/artigo/abacos/H*_a*/parede_tempo_carga.csv.

Writes SVG to resultados/artigo/figuras/; if rsvg-convert is on the PATH, also PNG (600 dpi) and PDF.
Usage: python exemplos/artigo/figuras.py
"""
import csv, glob, json, os, re, shutil, subprocess
from xml.sax.saxutils import escape as esc

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.normpath(os.path.join(AQUI, "..", ".."))
RES = os.path.join(RAIZ, "resultados")
OUT = os.path.join(RES, "artigo", "figuras") + os.sep
os.makedirs(OUT, exist_ok=True)

NOMES = {  # parameter names of the program -> names in the paper
    "Dilatação térmica do bloco": ("alpha", "Thermal expansion \u03b1"),
    "Excentricidade da carga (topo e base)": ("ecc", "Load eccentricity"),
    "Imperfeição inicial (H dividida por)": ("imp", "Initial imperfection"),
    "Fluência transitória do bloco (k_tr)": ("ktr", "Transient creep k_tr"),
    "Rigidez da alvenaria (fator na lei)": ("stiff", "Stiffness of wall law"),
    "Condutividade térmica do bloco": ("k", "Thermal conductivity k"),
    "Emissividade da superfície": ("epsf", "Surface emissivity"),
    "Calor específico do bloco": ("c", "Specific heat c"),
    "Emissividade das paredes dos furos": ("epsh", "Emissivity of holes"),
    "Resistência da alvenaria (fator na lei)": ("str", "Strength of wall law"),
    "Convecção na face exposta": ("h", "Convection h"),
}
SAIDAS = [("R", "Fire resistance (min)"), ("I", "Insulation (min)"), ("Pmax", "Cold capacity (kN)")]


def dados_morris():
    pasta = os.path.join(RES, "sensibilidade_parede_bcepq")
    linhas = {}
    for o, nome in SAIDAS:
        with open(os.path.join(pasta, f"sensibilidade_{o}.csv"), encoding="utf-8") as f:
            for r in csv.reader(f):
                if r[0] == "parâmetro":
                    continue
                p, en = NOMES.get(r[0], (r[0], r[0]))
                lo, hi = [float(x.replace(",", ".")) for x in r[3].split(" a ")] if " a " in r[3] else (0.0, 0.0)
                linhas[(o, p)] = dict(p=p, name=en, out=o, outName=nome, mu=round(float(r[2]), 1), lo=lo, hi=hi)
    ordem = [l["p"] for l in sorted((v for (o, _), v in linhas.items() if o == "R"), key=lambda v: -v["mu"])]
    rows = []
    for v in linhas.values():
        v["rank"] = ordem.index(v["p"]); rows.append(v)
    return rows


def dados_abacos():
    rows = []
    for d in sorted(glob.glob(os.path.join(RES, "artigo", "abacos", "H*_a*"))):
        pj = json.load(open(os.path.join(d, "projeto.json"), encoding="utf-8"))
        g = pj["geometria"]
        H = g["n_fiadas"] * g["altura_fiada"] / 1000.0 if g["altura_fiada"] > 1 else g["n_fiadas"] * g["altura_fiada"]
        t = g["largura"] / 1000.0 if g["largura"] > 1 else g["largura"]
        alfa = float(re.search(r"_a([0-9.]+)$", d).group(1))
        with open(os.path.join(d, "parede_tempo_carga.csv"), encoding="utf-8") as f:
            for r in list(csv.reader(f))[1:]:
                cens = r[3].strip().startswith(">")
                R = float(r[3].replace(">", "").strip())
                rows.append(dict(alpha=alfa, H=round(H, 2), Ht=round(H / t, 1), level=float(r[0]), R=R, cens=cens))
    return rows

INK, QUIET, GRID, AXIS, ACC, GREY = "#1a1a1a", "#555555", "#e4e4e4", "#8c8c8c", "#1f63a8", "#b9b9b9"
FONT = "Liberation Sans, Arial, Helvetica, sans-serif"


def svg(w, h, body):
    return (f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {w} {h}' width='{w}' height='{h}' "
            f"font-family='{FONT}' font-size='12'><rect width='{w}' height='{h}' fill='white'/>{body}</svg>")


def t(x, y, s, size=12, fill=INK, anchor="start", weight="normal", extra=""):
    return (f"<text x='{x:.1f}' y='{y:.1f}' font-size='{size}' fill='{fill}' text-anchor='{anchor}' "
            f"font-weight='{weight}' {extra}>{esc(s)}</text>")


def save(name, s, width_mm):
    p = OUT + name
    open(p + ".svg", "w", encoding="utf-8").write(s)
    if shutil.which("rsvg-convert"):
        subprocess.run(["rsvg-convert", "-w", str(int(width_mm / 25.4 * 600)), "-o", p + ".png", p + ".svg"], check=True)
        subprocess.run(["rsvg-convert", "-f", "pdf", "-o", p + ".pdf", p + ".svg"], check=True)


# ---------------- Figure 1: workflow ----------------
def box(x, y, w, h, main=False):
    if main:
        return (f"<rect x='{x}' y='{y}' width='{w}' height='{h}' rx='8' fill='{ACC}' fill-opacity='0.12' "
                f"stroke='{ACC}' stroke-width='2'/>")
    return f"<rect x='{x}' y='{y}' width='{w}' height='{h}' rx='8' fill='none' stroke='{AXIS}' stroke-width='1.25'/>"


def fig1():
    dy = -40  # no title
    b = [f"<defs><marker id='a' viewBox='0 0 10 10' refX='9' refY='5' markerWidth='6' markerHeight='6' "
         f"orient='auto-start-reverse'><path d='M0 0L10 5L0 10z' fill='{AXIS}'/></marker></defs>",
         f"<g transform='translate(0,{dy})'>",
         f"<g fill='none' stroke='{AXIS}' stroke-width='1.25'>"]
    for d in ["M110 120V152", "M110 208V240", "M196 268H228", "M356 144V176", "M356 260V228",
              "M468 366H488V202H468", "M528 252H560"]:
        b.append(f"<path d='{d}' marker-end='url(#a)'/>")
    b.append("</g>")
    for (y, name, body) in [(64, "Browser interface", "forms, log, 3D viewer"),
                            (152, "Project file", "JSON, six input groups"),
                            (240, "Validation", "checked before each run")]:
        b += [box(24, y, 172, 56), t(110, y + 24, name, 13, INK, "middle", "bold"), t(110, y + 42, body, 11.5, QUIET, "middle")]
    b.append(f"<rect x='228' y='56' width='300' height='356' rx='8' fill='none' stroke='{AXIS}' stroke-dasharray='4 4'/>")
    b.append(t(244, 80, "Analyses", 13, INK, "start", "bold"))
    for (y, name, body, main) in [(92, "2D section", "heat transfer, plane strain", False),
                                  (176, "Wall", "corotational beam, fibre section", True),
                                  (260, "3D prism", "thermal, elastic, damage, joints", False),
                                  (340, "Wall sensitivity", "Morris, runs in parallel", False)]:
        b += [box(244, y, 224, 52, main), t(356, y + 22, name, 13, INK, "middle", "bold"),
              t(356, y + 40, body, 11.5, INK if main else QUIET, "middle")]
    b += [t(364, 164, "fibre temperatures", 11.5, QUIET), t(364, 248, "prism σ–ε curves", 11.5, QUIET),
          t(506, 284, "r(k + 1) wall runs", 11.5, QUIET, "middle", extra="transform='rotate(-90 506 284)'")]
    b += [box(560, 168, 176, 168), t(576, 194, "Results folder", 13, INK, "start", "bold")]
    for k, s in enumerate(["CSV tables", "PNG figures", "VTK for ParaView", "JSON for the 3D viewer", "copy of the project"]):
        b.append(t(576, 220 + 24 * k, s, 11.5))
    b.append("</g>")
    save("Fig1_workflow", svg(760, 384, "".join(b)), 190)


# ---------------- Figure 2: Morris mu* ----------------
def fig2():
    rows = dados_morris()
    outs = ["R", "I", "Pmax"]
    n = max(r["rank"] for r in rows) + 1
    x0, pw, gap, yR0, dy = 196, 166, 18, 58, 24
    px = lambda i: x0 + i * (pw + gap)
    yr = lambda k: yR0 + k * dy
    hot = lambda r: (r["out"] == "R" and r["p"] == "alpha") or (r["out"] == "I" and r["p"] == "k")
    b = []
    for r in [r for r in rows if r["out"] == "R"]:
        if "k_tr" in r["name"]:
            b.append(f"<text x='{x0 - 12}' y='{yr(r['rank']) + 4}' font-size='12' fill='{INK}' text-anchor='end'>"
                     f"Transient creep k<tspan dy='3' font-size='9'>tr</tspan></text>")
        else:
            b.append(t(x0 - 12, yr(r["rank"]) + 4, r["name"], 12, INK, "end"))
    for i, o in enumerate(outs):
        rs = [r for r in rows if r["out"] == o]
        mx = max(r["hi"] for r in rs) * 1.18 or 1
        X = lambda v: px(i) + v / mx * pw
        big = max(r["mu"] for r in rs)
        b.append(t(px(i), 30, rs[0]["outName"], 12.5, INK, "start", "bold"))
        b.append(f"<line x1='{px(i)}' x2='{px(i)}' y1='{yR0 - 14}' y2='{yr(n - 1) + 14}' stroke='{AXIS}'/>")
        for r in rs:
            y = yr(r["rank"])
            wbar = max(0, X(r["mu"]) - px(i))
            b.append(f"<rect x='{px(i)}' y='{y - 6}' width='{wbar:.2f}' height='12' rx='2' fill='{ACC if hot(r) else GREY}'/>")
            if r["hi"] > 0:
                b.append(f"<line x1='{X(r['lo']):.2f}' x2='{X(r['hi']):.2f}' y1='{y}' y2='{y}' stroke='{INK}' stroke-width='1'/>")
            if r["mu"] >= 0.12 * big:
                lab = str(round(r["mu"])) if r["mu"] >= 10 else str(r["mu"])
                b.append(t(X(r["hi"]) + 5, y + 4, lab, 11.5))
            if r["mu"] == 0:
                b.append(t(px(i) + 5, y + 4, "0", 11.5, QUIET))
    yb = yr(n - 1) + 34
    b.append(t(x0, yb, "Bars: μ* (mean absolute elementary effect over the input range); lines: 95% bootstrap interval.", 11, QUIET))
    save("Fig2_morris", svg(760, yb + 14, "".join(b)), 190)


# ---------------- Figure 3: design chart ----------------
def fig3():
    rows = dados_abacos()
    als = sorted({r["alpha"] for r in rows}); hs = sorted({r["H"] for r in rows})
    left, pw, gap, y0p, ph = 64, 206, 22, 78, 210
    px = lambda i: left + i * (pw + gap)
    X = lambda i, v: px(i) + (v - 10) / 80 * pw
    Y = lambda v: y0p + ph - min(v, 240) / 240 * ph
    op = lambda h: 0.42 + 0.58 * hs.index(h) / (len(hs) - 1)
    ticks = [30, 60, 90, 120, 180, 240]
    b = []
    for j, h in enumerate(hs):
        r = next(q for q in rows if q["H"] == h)
        b.append(f"<line x1='{24 + j * 146}' x2='{44 + j * 146}' y1='22' y2='22' stroke='{ACC}' stroke-opacity='{op(h):.2f}' stroke-width='2.5'/>")
        b.append(t(50 + j * 146, 26, f"{h:.2f} m (H/t {r['Ht']})", 11.5))
    for tk in ticks:
        b.append(t(left - 8, Y(tk) + 4, str(tk), 11.5, QUIET, "end"))
    b.append(t(left - 8, y0p - 14, "min", 11.5, QUIET, "end"))
    for i, a in enumerate(als):
        for tk in ticks:
            b.append(f"<line x1='{px(i)}' x2='{px(i) + pw}' y1='{Y(tk):.2f}' y2='{Y(tk):.2f}' stroke='{GRID}'/>")
        b.append(f"<line x1='{px(i)}' x2='{px(i) + pw}' y1='{Y(0)}' y2='{Y(0)}' stroke='{AXIS}'/>")
        b.append(f"<text x='{px(i)}' y='{y0p - 14}' font-size='12' font-weight='bold' fill='{INK}'>α = {a:g} × 10"
                 f"<tspan dy='-5' font-size='8.5'>−6</tspan><tspan dy='5'>/°C</tspan></text>")
        for v, anc in [(10, "start"), (50, "middle"), (90, "end")]:
            b.append(t(X(i, v), Y(0) + 18, f"{v}%", 11.5, QUIET, anc))
        for h in hs:
            pts = sorted([r for r in rows if r["alpha"] == a and r["H"] == h], key=lambda r: r["level"])
            d = "".join(f"{'L' if k else 'M'}{X(i, r['level']):.2f} {Y(r['R']):.2f}" for k, r in enumerate(pts))
            b.append(f"<path d='{d}' fill='none' stroke='{ACC}' stroke-opacity='{op(h):.2f}' stroke-width='2'/>")
            for r in pts:
                xx, yy = X(i, r["level"]), Y(r["R"])
                if r["cens"]:
                    b.append(f"<path d='M{xx - 4:.2f} {yy + 2:.2f}L{xx:.2f} {yy - 5:.2f}L{xx + 4:.2f} {yy + 2:.2f}Z' fill='{ACC}' fill-opacity='{op(h):.2f}'/>")
                else:
                    b.append(f"<circle cx='{xx:.2f}' cy='{yy:.2f}' r='2.6' fill='{ACC}' fill-opacity='{op(h):.2f}'/>")
    b.append(t(left + (3 * pw + 2 * gap) / 2, Y(0) + 40, "Load level, % of the cold capacity of each wall", 11.5, QUIET, "middle"))
    b.append(t(left, Y(0) + 62, "Triangles: no failure within 240 min", 11.5, QUIET))
    save("Fig3_design_chart", svg(760, int(Y(0) + 74), "".join(b)), 190)


if __name__ == "__main__":
    fig1(); fig2(); fig3()
    print("figures written to", OUT)
