"""Thermal validation (section 7.3 of the paper): sensitivity of the insulation prediction for the walls of
Oliveira et al. (2021) and a 100 mm solid section (comparison with the EI 90 of EN 1992-1-2)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oliveira import project, HERE
from mef_masonry import project as pj

def run_case(name, p):
    r = pj.run(p, os.path.join(HERE, "oliveira_" + name), log=lambda s: None)
    res = {it[0]: it[1] for it in r["summary"]}
    print(f"{name:28s} insulation = {res.get('Insulation time (140/180 °C on the unexposed face)')} min", flush=True)

if __name__ == "__main__":
    for var in ("library", "en_moisture"):
        p = project(var); p["geometry"].update(holes_per_unit=0, hole_shape="circular", hole_diameter=0)
        run_case("solid_" + var, p)
    p = project("oliveira"); p["ceramic"]["table_c"] = [[20, 900], [100, 900], [100.01, 2020], [115, 2020], [200, 1000], [400, 1100], [1200, 1100]]
    run_case("moisture_3pct", p)
    p = project("oliveira"); p["test"]["hole_emissivity"] = 0.7; run_case("eps_holes_07", p)
    p = project("oliveira"); p["mesh"]["h"] = 2; run_case("mesh_2mm", p)
    p = project("oliveira"); p["geometry"].update(holes_per_unit=0, hole_shape="circular", hole_diameter=0)
    p["ceramic"].update(rho=2300 * 0.5365, table_rho=[], k=1.36 * 0.5365, table_k=[]); run_case("equivalent_homogeneous", p)
