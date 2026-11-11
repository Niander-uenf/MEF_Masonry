"""Comparison with ANSYS: BCEPQ prism of the UENF research project (Sá, 2025).

The ANSYS model (Mechanical 2024 R1, ds.dat input files of 06/01/2025) is the
prism of 5 courses of one and a half units, with the real tongue-and-groove interlock, SOLID187
elements, all contact surfaces bonded (bonded, MPC), fixed base,
uniform temperature in the body (T_ref = 22 °C) and pressure on the ledge of the top
course, outside the interlock. Linear material: E = 2801 MPa (applied to the material),
ν = 0.182, α = 6e-6 /°C. The transient thermal analysis of the same project has no
specific heat, so it returns the steady state: the whole prism at the
temperature imposed on its faces (the field imported in ds.dat is 600.0 °C at all
23,070 nodes).

MEF_Masonry models the same prism as a single body (the bonded faces make the
assembly a continuous solid), of 422.03 × 140 × 301.1 mm (from the base to the loaded
ledge), with the three 84 mm holes aligned at the mean ANSYS positions
(in the ANSYS courses the holes shift by up to 7 mm) and the interlock omitted. The
"with groove" variant removes the 5 mm × 7 mm groove from the bottom face, which in
ANSYS lies outside the support, and the load goes on the ledge only (outside the projection of the
interlock). As the problem is linear, each mesh solves two unit cases
(temperature and load) and the six ANSYS cases follow by superposition.

Output: results/paper/benchmark_ansys.json and a table in the terminal."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from mef_masonry.parameters import Geometry, Ceramic
from mef_masonry.thermal import Materials
from mef_masonry import mesh as mh, mechanical as mech

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "paper")

# ---------------------------------------------------------------- ANSYS data
L, W, H = 0.42203, 0.140, 0.3011          # m: length, width, base -> top ledge
HOLES_X = (0.07368, 0.21125, 0.34850)     # m: mean of the positions in the odd and even courses
D_HOLE = 0.084
T_REF = 22.0
E, NU, ALPHA = 2801e6, 0.182, 6e-6
# groove of the bottom face (course 1: half unit 0-142.3 mm and unit 142.03-422.03 mm):
# outer rectangle of the groove in each unit (x0, x1, y0, y1) and groove width
GROOVE = [(0.0287, 0.1250, 0.0217, 0.1183), (0.14203 + 0.0215, 0.14203 + 0.2575, 0.0217, 0.1183)]
B_GROOVE, P_GROOVE = 0.00507, 0.007
# interlock of the top course (projection of the tongue, which carries no load): half unit on the left, unit on the right
INTERLOCK_TOP = [(0.0287, 0.1250, 0.0217, 0.1183), (0.14203 + 0.0215, 0.14203 + 0.2575, 0.0217, 0.1183)]
# cases of the report (T in °C, force in N) and ANSYS responses: maximum total displacement (mm)
# and maximum von Mises stress (MPa)
CASES = [
    ("600 °C", 600.0, 97806.38, 1.1151, 55.858),
    ("800 °C", 800.0, 98956.41, 1.5566, 73.889),
    ("1000 °C", 1000.0, 114405.20, 1.9745, 92.489),
    ("load of 600 °C, no temperature", T_REF, 97806.38, 0.30700, 34.694),
    ("load of 800 °C, no temperature", T_REF, 98956.41, 0.31061, 35.102),
    ("load of 1000 °C, no temperature", T_REF, 114405.20, 0.35911, 40.582),
]


def _in_rectangle(x, y, r, border=0.0):
    x0, x1, y0, y1 = r
    return (x > x0 - border) & (x < x1 + border) & (y > y0 - border) & (y < y1 + border)


def _in_groove_band(x, y, rectangles=GROOVE):
    """Band of width B_GROOVE along the inner boundary of each rectangle (groove or tongue)."""
    s = np.zeros(np.shape(x), bool)
    for r in rectangles:
        x0, x1, y0, y1 = r
        inside = _in_rectangle(x, y, r)
        core = (x > x0 + B_GROOVE) & (x < x1 - B_GROOVE) & (y > y0 + B_GROOVE) & (y < y1 - B_GROOVE)
        s |= inside & ~core
    return s


def model(h, element, groove):
    geo = Geometry(unit_length=L / 1.5, width=W, course_height=H / 5, hole_diameter=D_HOLE,
                    n_courses=5, units_per_course=1.5, joint_type="dry", joint_thickness=0.002)
    m2 = mh.section_mesh(geo, h=h, symmetry=True, bands=None, holes=[(x, W / 2) for x in HOLES_X],
                        radius=D_HOLE / 2)
    lay = [mh.Layer(0.0, P_GROOVE, "uniform", material=mh.CERAMIC),
           mh.Layer(P_GROOVE, H, "uniform", material=mh.CERAMIC)]
    m3 = mh.extrude(m2, lay, dz=h, geo=geo, kind=element)
    if groove:
        base2d = m2.quad if element == "H8" else m2.tri
        c2 = m2.xy[base2d].mean(1)
        band = _in_groove_band(c2[:, 0], c2[:, 1])
        zc = m3.xyz[m3.conn3[:, :6] if element != "H8" else m3.conn3].mean(1)[:, 2]
        strip = (zc < P_GROOVE) & band[m3.tri_of]
        m3.mat[strip & (m3.mat == mh.CERAMIC)] = mh.AIR
    return m3


def solve_(h, element, groove):
    m3 = model(h, element, groove)
    M = Materials(ceramic=Ceramic(E_unit_gross=E, net_area_ratio=1.0, nu=NU, alpha=ALPHA, T_ref=T_REF))
    T600 = np.full(m3.n_nodes, 600.0)
    region = None
    if groove:   # the ANSYS pressure covers the whole ledge, except the base of the tongue (band of ~5 mm)
        region = lambda x, y: ~_in_groove_band(x, y, INTERLOCK_TOP)
    t0 = time.time()
    R = mech.solve_3d(m3, M, T_nodes=T600, load_N=CASES[0][2], plates="steel", fixed_base=True,
                       load_region=region)
    dt = time.time() - t0
    uT, uF = R["thermal"].u, R["load"].u          # 600 °C and 97.8 kN
    sT, sF = R["thermal"].sig, R["load"].sig
    active = R["thermal"].active
    res = []
    for name, T, F, _, _ in CASES:
        a, b = (T - T_REF) / (600.0 - T_REF), F / CASES[0][2]
        u = a * uT + b * uF
        s = (a * sT + b * sF)[active]
        vm = np.sqrt(0.5 * ((s[:, 0] - s[:, 1]) ** 2 + (s[:, 1] - s[:, 2]) ** 2 + (s[:, 2] - s[:, 0]) ** 2)
                     + 3 * (s[:, 3] ** 2 + s[:, 4] ** 2 + s[:, 5] ** 2))
        umax = np.linalg.norm(u, axis=1)
        i = int(np.argmax(umax))
        res.append({"case": name, "u_max_mm": float(umax[i] * 1e3), "node": m3.xyz[i].tolist(),
                    "von_mises_max_MPa": float(vm.max() / 1e6)})
    nel = int(active.sum())
    return {"h_mm": h * 1e3, "element": element, "groove": groove, "elements": nel, "nodes": int(m3.n_nodes),
            "time_s": round(dt, 1), "cases": res, "fit": fit_scales(uT, uF)}


def fit_scales(uT, uF):
    """Factors a (thermal part) and s (load part) that bring MEF_Masonry to the six
    ANSYS values: u = a uT (T - T_ref)/578 + s uF F/97.8 kN. Shows in which part of the
    response the difference lies (a = s = 1: same response)."""
    from scipy.optimize import least_squares
    target = np.array([c[3] for c in CASES])
    kT = np.array([(c[1] - T_REF) / (600.0 - T_REF) for c in CASES])
    kF = np.array([c[2] / CASES[0][2] for c in CASES])

    def umax(p):
        a, s = p
        return np.array([np.linalg.norm(a * t * uT + s * f * uF, axis=1).max() * 1e3 for t, f in zip(kT, kF)])
    r = least_squares(lambda p: umax(p) / target - 1, [1.0, 1.0])
    fit = umax(r.x)
    return {"a_thermal": float(r.x[0]), "s_load": float(r.x[1]),
            "u_max_fitted_mm": fit.tolist(), "max_error_pct": float(np.abs(fit / target - 1).max() * 100)}


if __name__ == "__main__":
    fast = "--fast" in sys.argv
    # meshes that fit in ~2 GB with SuperLU (without pypardiso)
    meshes = [(0.010, "H8", False), (0.010, "H8", True)] if fast else \
        [(0.015, "H8", False), (0.010, "H8", False), (0.0075, "H8", False), (0.015, "W15", False),
         (0.010, "H8", True), (0.0075, "H8", True), (0.015, "W15", True), (0.010, "W15", True)]
    output = {"ansys": [{"case": c[0], "T": c[1], "F_N": c[2], "u_max_mm": c[3], "von_mises_max_MPa": c[4]}
                       for c in CASES], "mef_masonry": []}
    for h, el, groove in meshes:
        r = solve_(h, el, groove)
        output["mef_masonry"].append(r)
        print(f"\n{el} h = {h * 1e3:.1f} mm, groove = {groove}: {r['elements']} elements, {r['time_s']} s")
        for c, a in zip(r["cases"], CASES):
            d = (c["u_max_mm"] - a[3]) / a[3] * 100
            print(f"  {c['case']:36s} u_max {c['u_max_mm']:.4f} mm (ANSYS {a[3]:.4f}, {d:+.2f}%)"
                  f"  max von Mises {c['von_mises_max_MPa']:.1f} MPa (ANSYS {a[4]:.1f})", flush=True)
        aj = r["fit"]
        print(f"  fit: a = {aj['a_thermal']:.4f} (thermal), s = {aj['s_load']:.4f} (load), "
              f"maximum error {aj['max_error_pct']:.2f}%", flush=True)
    os.makedirs(HERE, exist_ok=True)
    json.dump(output, open(os.path.join(HERE, "benchmark_ansys.json"), "w"), indent=1, ensure_ascii=False)
