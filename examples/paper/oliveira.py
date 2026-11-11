"""Thermal validation: walls of 3-cell concrete blocks of Oliveira et al. (2021), ISO 834."""
import sys, os, json, csv, copy
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from mef_masonry import project as pj, library as lib
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "paper")

def base():
    p = pj.default_project()
    p = lib.apply(p, "ceramic", "material_library", "concrete_nbr15220") if hasattr(lib, "apply") else p
    return p

def project(var):
    p = pj.default_project()
    m = copy.deepcopy(lib.MATERIALS["concrete_nbr15220"]["changes"])
    for g, d in m.items(): p[g].update(d)
    p["ceramic"]["material_library"] = "concrete_nbr15220"
    p["geometry"].update({"unit_library": "custom", "hole_shape": "rectangular", "unit_length": 201,
        "width": 100, "course_height": 100, "holes_per_unit": 3, "face_shell": 16, "web": 16, "end_shell": 16,
        "wall_cells": 0, "grouting": "none", "units_per_course": 1, "n_courses": 10, "joint_type": "dry",
        "joint_thickness": 2})
    p["model"].update({"analysis": "section_2d", "element": "T3"})
    p["test"].update({"kind": "iso834", "T_initial": 20, "iso_duration": 150, "h_convection": 25})
    p["mesh"].update({"h": 4, "dt": 30})
    c = p["ceramic"]
    if var in ("en_moisture", "oliveira"):
        # EN 1992-1-2: lower-limit k; c with the moisture peak u = 1.5% (1470 J/kg.K between 100 and 115 °C); rho(T)
        c["table_k"] = [[T, 1.36 - 0.136 * (T / 100) + 0.0057 * (T / 100) ** 2] for T in range(20, 1201, 100)]
        c["table_c"] = [[20, 900], [100, 900], [100.01, 1470], [115, 1470], [200, 1000], [400, 1100], [1200, 1100]]
        c["rho"] = 2300
        c["table_rho"] = [[20, 2300], [115, 2300], [200, 2300 * 0.98], [400, 2300 * 0.95], [1200, 2300 * 0.88]]
    if var == "oliveira":
        p["test"]["emissivity"] = 0.7
    return p

if __name__ == "__main__":
    out = {}
    for var in ("library", "en_moisture", "oliveira"):
        p = project(var)
        e = pj.validate(p)
        if e: print(var, e); continue
        r = pj.run(p, os.path.join(HERE, "oliveira_" + var), log=lambda s: None)
        res = {it[0]: it[1] for it in r["summary"]}
        rows = list(csv.reader(open(os.path.join(HERE, "oliveira_" + var, "temperatures_2d.csv"), encoding="utf-8")))
        out[var] = {"I": res.get("Insulation time (140/180 °C on the unexposed face)"), "ratio": res.get("Net-to-gross area ratio"),
                    "T": [[float(a), float(c)] for a, b, c in rows[1:]]}
        print(var, out[var]["I"], out[var]["ratio"], flush=True)
    json.dump(out, open(os.path.join(HERE, "oliveira.json"), "w"))
