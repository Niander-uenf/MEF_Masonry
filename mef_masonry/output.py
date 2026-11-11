"""Outputs: VTK files (open in ParaView), CSV and figures."""
import csv
from typing import Dict, Optional

import numpy as np

from .mesh import Mesh2D, Mesh3D


def save_vtk(path: str, m, node_data: Optional[Dict[str, np.ndarray]] = None,
               elem_data: Optional[Dict[str, np.ndarray]] = None):
    """Writes mesh and results in legacy VTK (ASCII), readable by ParaView.
    node_data/elem_data: {name: array (n,) scalar or (n, 3) vector}."""
    if isinstance(m, Mesh3D):
        pts = m.xyz
        tp = m.elem
        cell = m.conn3[:, tp.vtk_order]   # VTK order (W6: VTK_WEDGE, H8: VTK_HEXAHEDRON, W15: VTK_QUADRATIC_WEDGE)
        kind = tp.vtk
    else:
        pts = np.c_[m.xy[:, ::-1] if getattr(m, "vertical", False) else m.xy, np.zeros(m.n_nodes)]
        cell = m.tri
        kind = 5
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\nmef_masonry\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(pts)} double\n")
        np.savetxt(f, pts * 1000.0, fmt="%.6g")  # mm
        nn = cell.shape[1]
        f.write(f"CELLS {len(cell)} {len(cell) * (nn + 1)}\n")
        np.savetxt(f, np.c_[np.full(len(cell), nn), cell], fmt="%d")
        f.write(f"CELL_TYPES {len(cell)}\n")
        np.savetxt(f, np.full(len(cell), kind), fmt="%d")
        for section, data, n in (("POINT_DATA", node_data, len(pts)), ("CELL_DATA", elem_data, len(cell))):
            if not data:
                continue
            f.write(f"{section} {n}\n")
            for name, v in data.items():
                v = np.nan_to_num(np.asarray(v, float), nan=0.0)
                if v.ndim == 1:
                    f.write(f"SCALARS {name} double 1\nLOOKUP_TABLE default\n")
                    np.savetxt(f, v, fmt="%.6g")
                else:
                    f.write(f"VECTORS {name} double\n")
                    np.savetxt(f, v, fmt="%.6g")


def save_csv(path: str, lines, header=None):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if header:
            w.writerow(header)
        w.writerows(lines)


def section_figure(m2: Mesh2D, elem_values, title, path, unit_str="", mask=None, cmap="inferno"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    v = np.asarray(elem_values, float)
    tri = m2.tri if mask is None else m2.tri[mask]
    v = v if mask is None else v[mask]
    if getattr(m2, "vertical", False):      # vertical section: thickness horizontal, height vertical
        fig, ax = plt.subplots(figsize=(5.2, 5.6))
        tp = ax.tripcolor(m2.xy[:, 1] * 1000, m2.xy[:, 0] * 1000, tri, facecolors=v, cmap=cmap)
        ax.set_aspect("equal"); ax.set_xlabel("y, thickness (mm; y = 0 is the exposed face)")
        ax.set_ylabel("height in the course (mm)")
    else:
        fig, ax = plt.subplots(figsize=(10, 3.2 if not m2.symmetry else 2.2))
        tp = ax.tripcolor(m2.xy[:, 0] * 1000, m2.xy[:, 1] * 1000, tri, facecolors=v, cmap=cmap)
        ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    ax.set_title(title, fontsize=10)
    fig.colorbar(tp, ax=ax, label=unit_str, shrink=0.8)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
