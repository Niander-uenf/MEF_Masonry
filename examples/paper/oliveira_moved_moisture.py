"""Section 7.3 of the paper: a cheap test of the moisture-transport hypothesis in the walls of
Oliveira et al. (2021). The thermal model does not transport water; here the water is redistributed
by hand before the test, as if it had already migrated to the cold shell: the specific heat
peak of EN 1992-1-2 (100 to 115 °C) leaves the exposed shell and goes to the unexposed shell,
with the mass of water conserved.

  A: the water of the exposed shell moves to the unexposed shell (exposed 0%, webs 1.5%, cold 3%)
  B: all the water of the section goes to the unexposed shell (exposed and webs 0%, cold ~5%)

Peak of c(T): EN 1992-1-2 (0% -> 900, 1.5% -> 1470, 3% -> 2020 J/kg.K), extrapolated
linearly above 3%. Base: project "oliveira" (lower-limit k, 1.5%, emissivity 0.7)."""
import os, sys, json, copy
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oliveira import project, HERE
from mef_masonry import project as pj, thermal as th, mesh as msh, curves
from mef_masonry.mesh import AIR, CERAMIC, MORTAR, GROUT
from mef_masonry.thermal import Boundary


SHAPE = "en1992"   # "en1992": peak of EN 1992-1-2; "latent": all the heat of vaporisation between 100 and 110 °C


def table_c(u):
    global SHAPE
    if SHAPE == "latent":
        P = 2.257e6 * u / 100 / 10.0          # extra J/kg.K between 100 and 110 °C
        return [(20, 900), (100, 900), (100.01, 900 + P), (110, 910 + P), (110.01, 910), (200, 1000), (400, 1100), (1200, 1100)]
    peak = 900 + u / 1.5 * 570 if u <= 1.5 else 1470 + (u - 1.5) / 1.5 * 550
    return [(20, 900), (100, 900), (100.01, peak), (115, peak), (200, 1000), (400, 1100), (1200, 1100)]


def run_case(u_exp, u_web, u_cold, name):
    p = project("oliveira")
    geo = pj.geometry(p)
    M = pj.materials(p, geo)
    ts, msh_ = p["test"], p["mesh"]
    m2 = msh.section_mesh(geo, h=msh_["h"] / 1e3, symmetry=False, bands=())
    mat = th.section_material(m2, geo).copy()
    yc = m2.xy[m2.tri][:, :, 1].mean(1)
    W, tc = geo.width, p["geometry"]["face_shell"] / 1e3
    sol = mat == CERAMIC
    exp, cold = sol & (yc < tc), sol & (yc > W - tc)
    mat[exp], mat[cold] = MORTAR, GROUT          # labels only: the three zones have the same concrete
    for target in (M.mortar, M.grout):
        target.k, target.rho, target.c = M.ceramic.k, M.ceramic.rho, M.ceramic.c
        target.table_k, target.table_rho = M.ceramic.table_k, M.ceramic.table_rho
    M.mortar.table_c = table_c(u_exp)
    M.ceramic.table_c = table_c(u_web)
    M.grout.table_c = table_c(u_cold)
    d1, d2 = m2.xy[m2.tri[:, 1]] - m2.xy[m2.tri[:, 0]], m2.xy[m2.tri[:, 2]] - m2.xy[m2.tri[:, 0]]
    area = np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]) / 2
    water = sum(u * area[s].sum() for u, s in ((u_exp, exp), (u_web, sol & ~exp & ~cold), (u_cold, cold)))
    f, tf, _ = pj.thermal_programme(p)
    nexp = th.boundary_lengths_2d(m2, faces=("y1",))
    ctr = [Boundary(th.boundary_lengths_2d(m2, faces=("y0",)), f, ts["h_convection"], ts["emissivity"], "exposed"),
           Boundary(nexp, curves.constant(ts["T_initial"]), 4.0, ts["emissivity"], "unexposed")]
    face = np.nonzero(nexp > 0)[0]; w = nexp[face] / nexp[face].sum(); Ti = ts["T_initial"]
    r = th.solve_2d(m2, M, ctr, tf, dt=msh_["dt"], T_initial=Ti, mat=mat, record_every=60.0,
                       probes={"face": face}, stop_if=lambda t, T: (T[face] @ w - Ti >= 140) or (T[face].max() - Ti >= 180))
    I = r.stopped_at / 60 if r.stopped_at else None
    print(f"{name:40s} u = {u_exp:.2f} / {u_web:.2f} / {u_cold:.2f} %   relative water = {water:.4e}   insulation = {I} min", flush=True)
    return {"u_exposed": u_exp, "u_webs": u_web, "u_cold": u_cold, "water_rel": water, "I_min": I,
            "areas_mm2": {"exposed": area[exp].sum() * 1e6, "webs": area[sol & ~exp & ~cold].sum() * 1e6, "cold": area[cold].sum() * 1e6}}


if __name__ == "__main__":
    out = {}
    for shape in ("en1992", "latent"):
        SHAPE = shape
        print("peak of c(T):", shape)
        o = out[shape] = {"reference": run_case(1.5, 1.5, 1.5, "reference (1.5% throughout the section)")}
        a = o["reference"]["areas_mm2"]
        o["A"] = run_case(0.0, 1.5, 1.5 * (a["exposed"] + a["cold"]) / a["cold"], "A: water of the exposed shell in the cold one")
        o["B"] = run_case(0.0, 0.0, 1.5 * (a["exposed"] + a["webs"] + a["cold"]) / a["cold"], "B: all the water in the cold shell")
        for u in (10, 20):
            o[f"cold_{u}"] = run_case(0.0, 0.0, u, f"cold shell with {u}% (extra water)")
    json.dump(out, open(os.path.join(HERE, "oliveira_moved_moisture.json"), "w"), indent=1)
