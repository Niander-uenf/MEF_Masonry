"""Verification of phase 4: transient creep (load-induced thermal strain).

Model: on first heating under compression, Δε_tr = Δφ / fc20 · [(1 + ν) σ⁻ - ν tr(σ⁻) I],
with φ(T) = k_tr ε_th(T) (Anderberg and Thelandersson, 1976) or tabulated.

V1  Homogeneous prism under constant stress, heated uniformly (W6, H8 and W15):
    ε_z = ε_th + σ/E + k σ/fc ε_th and ε_x = ε_th - ν σ/E - ν k σ/fc ε_th, exact.
V2  Heating with the top locked (relaxation): σ = fc/k [exp(-E k ε_th / fc) - 1], the
    solution of dσ/dε_th = -E (1 + k σ/fc); the adaptive step controls the error.
V3  Irreversibility: heating to 520 °C, cooling to 270 °C and reheating to 520 °C does not change
    ε_tr; going on to 620 °C adds only the part from 520 to 620 °C.
V4  Without compression there is no creep: under tension, ε_tr = 0.
V5  Tabulated φ(T) replaces k ε_th: ε_tr = φ(T) σ / fc.
V6  Wall (fibres): short column under constant load, heated uniformly; exact
    shortening. The coefficient of the wall with mortar is the series mean of unit and joint.
V7  BCEPQ wall under ISO 834 (k_tr = 0 by default): creep on and off give the
    same time; with a hypothetical k_tr, the time changes (recorded separately, without tolerance).
V8  Defaults and library: mortar 2.0 (Prakash et al., 2020), unit and grout 0; the
    ready-made projects pass validation.
V9  (--full) Transient test of the solid brick prism: creep of the mortar
    relieves the thermal stresses between brick and joint, and the damage during heating falls.

Usage: python examples/14_verification_creep.py [--full]
Writes results/verification_creep.csv.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mef_masonry import mesh as msh, project as pj, library as lib, curves
from mef_masonry.parameters import Geometry, CeramicDamage
from mef_masonry.thermal import Materials
from mef_masonry.nonlinear import NonlinearModel
from mef_masonry.wall import WallModel, FibreSection, LinearLaw
from mef_masonry.output import save_csv

FULL = "--full" in sys.argv
lines = []


def record(case, model, ref, tol_rel, unit_str="", tol_abs=0.0):
    error = abs(model - ref) / max(abs(ref), 1e-30)
    ok = error <= tol_rel or abs(model - ref) <= tol_abs
    tol = f"{100 * tol_rel:.2g}%" if tol_abs <= 0 else f"{tol_abs:g} {unit_str}"
    lines.append([case, round(model, 9), round(ref, 9), unit_str, f"{100 * error:.3g}%", tol,
                   "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: model {model:.6g} {unit_str}, reference {ref:.6g} ({100 * error:.3g}%)")


def condition(case, ok, a, b, unit_str, criterion):
    lines.append([case, round(a, 6), round(b, 6), unit_str, "—", criterion, "ok" if ok else "FAILED"])
    print(f"{'ok ' if ok else 'XX '} {case}: {a:.5g} and {b:.5g} {unit_str} ({criterion})")


def report(case, a, unit_str, note):
    lines.append([case, round(a, 4), "", unit_str, "—", note, "informative"])
    print(f"--  {case}: {a:.4g} {unit_str} ({note})")


# ----------------------------------------------------------------------------- test prism
E, NU, FC, ALPHA, K = 10e9, 0.2, 20e6, 8e-6, 2.0
geo = Geometry(unit_length=0.10, width=0.10, course_height=0.20, hole_diameter=0, holes_per_unit=0,
                units_per_course=1, n_courses=1, bond="stack_bond", joint_type="mortar")
A, H = geo.length * geo.width, geo.prism_height


def materials():
    M = Materials()
    c = M.ceramic
    c.E_unit_gross, c.net_area_ratio, c.nu, c.alpha, c.T_ref = E, 1.0, NU, ALPHA, 20.0
    return M


def model(m, k=K, table=None, ft=50e6):
    mod = NonlinearModel(m, materials(), CeramicDamage(ft=ft, fc=FC), None, verbose=False)
    mod.set_creep(k_unit=k, table_unit=table)
    return mod


def uniform(m, T):
    return np.full(m.n_nodes, float(T))


def heat_to(mod, m, T0, T1, n, **kw):
    mod.heat({i: uniform(m, T0 + (T1 - T0) * i / n) for i in range(1, n + 1)}, **kw)


def lateral(mod, m):
    ux, x = mod.u[0::3], m.xyz[:, 0]
    return (ux[np.isclose(x, geo.length)].mean() - ux[np.isclose(x, 0.0)].mean()) / geo.length


# ---------------------------------------------------------------------------- V1
sig = -6e6
eth = ALPHA * 500
ref_z = eth + sig / E + K * sig / FC * eth
ref_x = eth - NU * sig / E - NU * K * sig / FC * eth
for el in ("W6", "H8", "W15"):
    m = msh.prism_mesh(geo, h=0.025, dz=0.05, t_joint=0.01, symmetry=True, plates=False, capping=False, element=el)
    mod = model(m)
    mod.heat({0.0: uniform(m, 20)}, pre_load_N=-sig * A)
    heat_to(mod, m, 20, 520, 10, pre_load_N=-sig * A)
    record(f"V1 axial strain under constant load, {el}", mod.u[3 * mod.top + 2].mean() / H * 1e3,
              ref_z * 1e3, 1e-9, "‰")
    record(f"V1 lateral strain (transient Poisson = ν), {el}", lateral(mod, m) * 1e3, ref_x * 1e3, 1e-9, "‰")

m = msh.prism_mesh(geo, h=0.025, dz=0.05, t_joint=0.01, symmetry=True, plates=False, capping=False)

# ---------------------------------------------------------------------------- V2
ref2 = FC / K * (np.exp(-E * K * ALPHA * 500 / FC) - 1)
for tol_tr, tol in ((None, 1e-2), (1e-5, 1e-3)):
    mod = model(m)
    if tol_tr:
        mod.tol_tr = tol_tr
    mod.heat({0.0: uniform(m, 20)}, pre_load_N=0.0)
    heat_to(mod, m, 20, 520, 10, top_fixed=True)
    s2 = float(mod.nominal_stress()[:, 2].mean())
    record(f"V2 stress during restrained heating to 520 °C (relaxation), tol_tr = {mod.tol_tr:g}", s2 / 1e6,
              ref2 / 1e6, tol, "MPa")
report("V2 stress without creep during the same heating (-E α ΔT)", -E * ALPHA * 500 / 1e6, "MPa", "for comparison")

# ---------------------------------------------------------------------------- V3
mod = model(m)
P = 6e6 * A
mod.heat({0.0: uniform(m, 20)}, pre_load_N=P)
heat_to(mod, m, 20, 520, 5, pre_load_N=P)
e520 = float(mod.eps_tr[:, 2].mean())
heat_to(mod, m, 520, 270, 5, pre_load_N=P)
e_cool = float(mod.eps_tr[:, 2].mean())
heat_to(mod, m, 270, 520, 5, pre_load_N=P)
e_reheat = float(mod.eps_tr[:, 2].mean())
heat_to(mod, m, 520, 620, 2, pre_load_N=P)
e620 = float(mod.eps_tr[:, 2].mean())
record("V3 ε_tr a 520 °C = k σ/fc ε_th", e520 * 1e3, -K * 0.3 * ALPHA * 500 * 1e3, 1e-9, "‰")
record("V3 ε_tr after cooling to 270 °C (does not recover)", e_cool * 1e3, e520 * 1e3, 1e-12, "‰")
record("V3 ε_tr after reheating to 520 °C (does not add)", e_reheat * 1e3, e520 * 1e3, 1e-12, "‰")
record("V3 ε_tr at 620 °C (adds only from 520 to 620 °C)", e620 * 1e3, -K * 0.3 * ALPHA * 600 * 1e3, 1e-9, "‰")

# ---------------------------------------------------------------------------- V4
mod = model(m)
mod.heat({0.0: uniform(m, 20)}, pre_load_N=-1e6 * A)          # tension of 1 MPa
heat_to(mod, m, 20, 520, 5, pre_load_N=-1e6 * A)
condition("V4 under tension there is no transient creep", float(np.abs(mod.eps_tr).max()) < 1e-12,
         float(np.abs(mod.eps_tr).max()), 0.0, "", "< 1e-12 (arredondamento)")

# ---------------------------------------------------------------------------- V5
tab = [(20.0, 0.0), (300.0, 0.002), (520.0, 0.010)]
mod = model(m, k=0.0, table=tab)
mod.heat({0.0: uniform(m, 20)}, pre_load_N=P)
heat_to(mod, m, 20, 410, 6, pre_load_N=P)
phi410 = float(np.interp(410, *zip(*tab)))
record("V5 ε_tr with tabulated φ(T), at 410 °C", float(mod.eps_tr[:, 2].mean()) * 1e3, -phi410 * 0.3 * 1e3, 1e-9, "‰")

# ---------------------------------------------------------------------------- V6
t = 0.14
sec = FibreSection.rectangle(t, 20)
Ew, a_w, fc_w = 5e9, 6e-6, 10e6
c = lambda T: K * a_w * np.maximum(np.asarray(T, float) - 20, 0.0) / fc_w   # noqa: E731
mw = WallModel(1.0, sec, LinearLaw(Ew), lambda T: a_w * (np.asarray(T, float) - 20), n_el=4, T_initial=20,
                  regularise=False, creep=c)
Pw = 3e6 * t
mw.apply_load(Pw, 0.0, 0.0)
mw.fire(Pw, 0.0, 0.0, [(60.0 * i, np.full(len(sec.y), 20 + 50.0 * i)) for i in range(11)])
ref6 = a_w * 500 - 3e6 / Ew - K * a_w * 500 * 3e6 / fc_w
record("V6 wall: top displacement under constant load and heating", mw.top_displacement * 1e3,
          ref6 * 1e3, 1e-9, "mm")
p6 = pj.complete(lib.apply(pj.default_project(), "unit", "solid_brick_10x6x22"))
p6["model"].update({"analysis": "wall", "k_tr_unit": 1.5})
g6 = pj.geometry(p6)
M6 = pj.materials(p6, g6)
_, _, sec6 = pj._wall_section(p6, g6)
coef = pj._wall_creep(p6, g6, M6, sec6)
T6 = 400.0
fa = g6.joint_thickness / g6.course_height
ref_c = ((1 - fa) * 1.5 * float(M6.ceramic.thermal_strain(np.array([T6]))[0]) / (p6["model"]["fc"] * 1e6)
         + fa * 2.0 * float(M6.mortar.thermal_strain(np.array([T6]))[0]) / (p6["ceramic"]["mor_fc"] * 1e6))
record("V6 wall coefficient: unit and joint in series, at 400 °C", float(coef(np.full(len(sec6.y), T6))[0]) * 1e12,
          ref_c * 1e12, 1e-12, "×10⁻¹² /Pa")

# ---------------------------------------------------------------------------- V7
base = pj.complete(pj.apply(pj.default_project(), pj.PRESETS["BCEPQ wall 2.94 × 3.00 m — ISO 834 under load"]))
base["mesh"]["n_elements"] = 12
base["test"]["iso_duration"] = 120


def wall(**mods):
    """(fire resistance time, maximum deflection at the end in mm) of the BCEPQ wall under ISO 834."""
    import tempfile
    q = {**base, "model": {**base["model"], **mods}}
    with tempfile.TemporaryDirectory() as d:
        r = pj.run(q, d, log=lambda s: None)
    v = {name: x for name, x, _ in r["summary"]}
    return v["Fire resistance time (loss of equilibrium)"], float(v["Maximum deflection at the end"])


t_off, f_off = wall(creep=False)
t_on, f_on = wall(creep=True)
condition(f"V7 BCEPQ wall, k_tr = 0: same deflection at the end with creep on (time {t_on} min)",
         t_on == t_off and f_on == f_off, f_on, f_off, "mm", "equal")
if FULL:
    base["test"]["iso_duration"] = 240
    t0, f0 = wall(creep=False)
    for k in (1.0, 2.0):
        tk, fk = wall(creep=True, k_tr_unit=k)
        report(f"V7 BCEPQ wall, 240 min, k_tr = {k:g} (hypothetical): time {tk} min; deflection at the end", fk, "mm",
                 f"without creep: {t0} min and {f0:.1f} mm")

# ---------------------------------------------------------------------------- V8
p8 = pj.complete(pj.default_project())
condition("V8 defaults: k_tr of the mortar 2.0; unit and grout 0",
         (p8["model"]["k_tr_mortar"], p8["model"]["k_tr_unit"], p8["model"]["k_tr_grout"]) == (2.0, 0.0, 0.0),
         p8["model"]["k_tr_mortar"], p8["model"]["k_tr_unit"], "", "2,0 / 0 / 0")
bad = [k for k, v in pj.PRESETS.items() if pj.validate(pj.complete(pj.apply(pj.default_project(), v)))]
condition("V8 ready-made projects pass validation", not bad, len(pj.PRESETS) - len(bad),
         len(pj.PRESETS), "projects", "all")

# ---------------------------------------------------------------------------- V9
if FULL:
    import tempfile

    def transient_damage(creep):
        q = pj.complete(lib.apply(pj.default_project(), "unit", "solid_brick_10x6x22"))
        q["model"].update({"analysis": "nonlinear", "creep": creep})
        q["test"].update({"kind": "transient", "load": 46, "rate": 5, "T_target": 800, "dwell": 30,
                            "T_initial": 20})
        with tempfile.TemporaryDirectory() as d:
            pj.run(q, d, log=lambda s: None)
            a = np.loadtxt(os.path.join(d, "transient_test.csv"), delimiter=",", skiprows=1)
        return float(a[-1, -1])
    d_off, d_on = transient_damage(False), transient_damage(True)
    condition("V9 solid brick in the transient test: creep of the mortar reduces the damage (%)", d_on < d_off,
             d_on, d_off, "%", "with < without")

output = os.path.join(os.path.dirname(__file__), "..", "results", "verification_creep.csv")
save_csv(output, lines, ["case", "model", "reference", "unit", "error", "tolerance", "status"])
nf = sum(1 for l in lines if l[-1] == "FAILED")
nv = sum(1 for l in lines if l[-1] != "informative")
print(f"\n{nv - nf}/{nv} cases within tolerance"
      + (f" (plus {len(lines) - nv} informative result(s))" if len(lines) > nv else "")
      + f". Written to {os.path.abspath(output)}")
sys.exit(1 if nf else 0)
