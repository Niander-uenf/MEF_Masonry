"""Compact export for the visual post-processor (WebGL) of the app.

Writes a JSON with the mesh (coordinates in mm, connectivity, material) and
a list of states; each state carries nodal or element fields and,
optionally, the nodal displacement (mm) for the deformed shape.
Missing values (for example, damage outside the ceramic) are written as null.
"""
import json
from typing import Dict, List, Optional

import numpy as np

from .mesh import Mesh2D, Mesh3D, NAMES

# id: (name, unit, location, suggested colour map)
FIELDS = {
    "T": ("Temperature", "°C", "node", "inferno"),
    "umag": ("Total displacement", "mm", "node", "viridis"),
    "uz": ("Vertical displacement", "mm", "node", "diverging"),
    "s1": ("Maximum principal stress", "MPa", "element", "diverging"),
    "s3": ("Minimum principal stress", "MPa", "element", "diverging"),
    "szz": ("Vertical stress", "MPa", "element", "diverging"),
    "d": ("Effective damage", "", "element", "damage"),
    "dt": ("Tensile damage", "", "element", "damage"),
    "dc": ("Compressive damage", "", "element", "damage"),
    "gap": ("Initial gap of the joints", "mm", "element", "viridis"),
    "material": ("Material", "", "element", "categories"),
}


def _list(v, fmt="{:.4g}"):
    v = np.asarray(v, float).ravel()
    return "[" + ",".join("null" if not np.isfinite(x) else fmt.format(x) for x in v) + "]"


def _int_list(v):
    return "[" + ",".join(str(int(x)) for x in np.asarray(v).ravel()) + "]"


def export(path: str, m, states: List[Dict], description: str = "", material=None):
    """states: [{"label": str, "fields": {id: array}, "u": (n, 3) in m or None}]
    The W15 wedge is written as W6 (corner nodes only): the edges are straight and the
    viewer paints the nodal fields at the corners."""
    faces = None
    sub = None
    if isinstance(m, Mesh3D):
        con = m.conn3
        if m.kind == "W15":
            sub = np.unique(m.conn3[:, :6])
            new = -np.ones(m.n_nodes, int); new[sub] = np.arange(len(sub))
            con = new[m.conn3[:, :6]]
            xyz = m.xyz[sub] * 1e3
        else:
            xyz = m.xyz * 1e3
        kind, npe = "3d", con.shape[1]
        if m.kind == "H8":
            from .elements3d import kind as element_type
            faces = element_type("H8").faces_vis
        symmetry, wid = m.m2d.symmetry, m.m2d.model_width * 1e3
        mat = m.mat
    else:
        xy = m.xy[:, ::-1] if getattr(m, "vertical", False) else m.xy   # vertical section: thickness × height
        xyz, con, kind, npe = np.c_[xy, np.zeros(m.n_nodes)] * 1e3, m.tri, "2d", 3
        symmetry, wid = m.symmetry, m.model_width * 1e3
        mat = np.where(m.mark == 0, 0, 1) if material is None else np.asarray(material)
    used = sorted({k for e in states for k in e["fields"]})
    fields = [{"id": k, "name": FIELDS.get(k, (k,))[0], "unit": FIELDS.get(k, ("", ""))[1],
               "local": FIELDS.get(k, ("", "", "node"))[2], "cmap": FIELDS.get(k, ("", "", "", "viridis"))[3]}
              for k in used]
    nn = m.n_nodes

    def nodal(v):
        v = np.asarray(v)
        return v[sub] if (sub is not None and len(v) == nn) else v
    parts = ['{"kind":"%s","description":%s,"symmetry":%s,"model_width":%.6g,"nodes_per_element":%d,'
              % (kind, json.dumps(description, ensure_ascii=False), "true" if symmetry else "false", wid, npe)
              + ('"faces":%s,' % json.dumps(faces) if faces else ""),
              '"material_names":' + json.dumps({int(k): v for k, v in NAMES.items()}, ensure_ascii=False) + ",",
              '"nodes":' + _list(xyz, "{:.6g}") + ",",
              '"elements":' + _int_list(con) + ",",
              '"material":' + _int_list(mat) + ",",
              '"fields":' + json.dumps(fields, ensure_ascii=False) + ',"states":[']
    for i, e in enumerate(states):
        s = '{"label":%s,"values":{' % json.dumps(e["label"], ensure_ascii=False)
        s += ",".join('"%s":%s' % (k, _list(nodal(v))) for k, v in e["fields"].items())
        s += '},"u":' + ("null" if e.get("u") is None else _list(nodal(np.asarray(e["u"]).reshape(-1, 3)) * 1e3, "{:.5g}"))
        s += "}" + ("," if i < len(states) - 1 else "")
        parts.append(s)
    parts.append("]}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(parts))


def u_fields(u):
    u = np.asarray(u, float).reshape(-1, 3)
    return {"umag": np.linalg.norm(u, axis=1) * 1e3, "uz": u[:, 2] * 1e3}


def nonlinear_states(mod, load_max: int = 24) -> List[Dict]:
    """States of the non-linear model: heating fields and loading steps
    (sampled, always including the peak and the last one)."""
    inst = mod.instantaneous or []
    ht = [s for s in inst if s["phase"] == "heating"]
    cg = [s for s in inst if s["phase"] == "load"]
    sel = list(range(len(cg)))
    if len(cg) > load_max:
        ip = int(np.argmax([s["load_N"] for s in cg]))
        sel = sorted(set(np.linspace(0, len(cg) - 1, load_max).round().astype(int)) | {ip, len(cg) - 1})
    ne = len(mod.m.conn3)
    T_last = mod.T

    def per_element(v):
        out = np.full(ne, np.nan)
        out[mod.id_glob[mod.ic]] = v
        return out
    sty = []
    for s in ht:
        T = s["T"] if s["T"] is not None else T_last
        sty.append({"label": f"heating, t = {s['label'] / 60:.0f} min",
                    "fields": {"T": T, **u_fields(s["u"]), "dt": per_element(s["dt"]), "dc": per_element(s["dc"]),
                               "d": per_element(s["d"])}, "u": s["u"]})
    ip = int(np.argmax([s["load_N"] for s in cg])) if cg else -1
    for k in sel:
        s = cg[k]
        r = f"load {s['load_N'] / 1e3:.1f} kN, δ = {s['label'] * 1e3:.3f} mm" + (" (peak)" if k == ip else "")
        sty.append({"label": r, "fields": {"T": T_last, **u_fields(s["u"]), "dt": per_element(s["dt"]),
                                           "dc": per_element(s["dc"]), "d": per_element(s["d"])}, "u": s["u"]})
    return sty
