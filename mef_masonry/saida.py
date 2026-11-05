"""Saídas: arquivos VTK (abrir no ParaView), CSV e figuras."""
import csv
from typing import Dict, Optional

import numpy as np

from .malha import Malha2D, Malha3D


def salvar_vtk(caminho: str, m, dados_nos: Optional[Dict[str, np.ndarray]] = None,
               dados_elem: Optional[Dict[str, np.ndarray]] = None):
    """Grava malha e resultados em VTK legado (ASCII), legível pelo ParaView.
    dados_nos/dados_elem: {nome: array (n,) escalar ou (n, 3) vetor}."""
    if isinstance(m, Malha3D):
        pts = m.xyz
        tp = m.elem
        cel = m.cun[:, tp.vtk_ordem]   # ordem do VTK (W6: VTK_WEDGE, H8: VTK_HEXAHEDRON, W15: VTK_QUADRATIC_WEDGE)
        tipo = tp.vtk
    else:
        pts = np.c_[m.xy[:, ::-1] if getattr(m, "vertical", False) else m.xy, np.zeros(m.n_nos)]
        cel = m.tri
        tipo = 5
    with open(caminho, "w") as f:
        f.write("# vtk DataFile Version 3.0\nmef_masonry\nASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(pts)} double\n")
        np.savetxt(f, pts * 1000.0, fmt="%.6g")  # mm
        nn = cel.shape[1]
        f.write(f"CELLS {len(cel)} {len(cel) * (nn + 1)}\n")
        np.savetxt(f, np.c_[np.full(len(cel), nn), cel], fmt="%d")
        f.write(f"CELL_TYPES {len(cel)}\n")
        np.savetxt(f, np.full(len(cel), tipo), fmt="%d")
        for secao, dados, n in (("POINT_DATA", dados_nos, len(pts)), ("CELL_DATA", dados_elem, len(cel))):
            if not dados:
                continue
            f.write(f"{secao} {n}\n")
            for nome, v in dados.items():
                v = np.nan_to_num(np.asarray(v, float), nan=0.0)
                if v.ndim == 1:
                    f.write(f"SCALARS {nome} double 1\nLOOKUP_TABLE default\n")
                    np.savetxt(f, v, fmt="%.6g")
                else:
                    f.write(f"VECTORS {nome} double\n")
                    np.savetxt(f, v, fmt="%.6g")


def salvar_csv(caminho: str, linhas, cabecalho=None):
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if cabecalho:
            w.writerow(cabecalho)
        w.writerows(linhas)


def figura_secao(m2: Malha2D, valores_elem, titulo, caminho, unidade="", mascara=None, cmap="inferno"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    v = np.asarray(valores_elem, float)
    tri = m2.tri if mascara is None else m2.tri[mascara]
    v = v if mascara is None else v[mascara]
    if getattr(m2, "vertical", False):      # seção vertical: espessura na horizontal, altura na vertical
        fig, ax = plt.subplots(figsize=(5.2, 5.6))
        tp = ax.tripcolor(m2.xy[:, 1] * 1000, m2.xy[:, 0] * 1000, tri, facecolors=v, cmap=cmap)
        ax.set_aspect("equal"); ax.set_xlabel("y, espessura (mm; y = 0 é a face exposta)")
        ax.set_ylabel("altura na fiada (mm)")
    else:
        fig, ax = plt.subplots(figsize=(10, 3.2 if not m2.simetria else 2.2))
        tp = ax.tripcolor(m2.xy[:, 0] * 1000, m2.xy[:, 1] * 1000, tri, facecolors=v, cmap=cmap)
        ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    ax.set_title(titulo, fontsize=10)
    fig.colorbar(tp, ax=ax, label=unidade, shrink=0.8)
    fig.tight_layout(); fig.savefig(caminho, dpi=130); plt.close(fig)
