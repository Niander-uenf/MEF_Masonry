"""Exportação compacta para o pós-processador visual (WebGL) do app.

Grava um JSON com a malha (coordenadas em mm, conectividade, material) e
uma lista de estados; cada estado traz campos nodais ou por elemento e,
opcionalmente, o deslocamento nodal (mm) para a forma deformada.
Valores ausentes (por exemplo, dano fora da cerâmica) vão como null.
"""
import json
from typing import Dict, List, Optional

import numpy as np

from .malha import Malha2D, Malha3D, NOMES

# id: (nome, unidade, local, mapa de cores sugerido)
CAMPOS = {
    "T": ("Temperatura", "°C", "no", "inferno"),
    "umag": ("Deslocamento total", "mm", "no", "viridis"),
    "uz": ("Deslocamento vertical", "mm", "no", "divergente"),
    "s1": ("Tensão principal máxima", "MPa", "elemento", "divergente"),
    "s3": ("Tensão principal mínima", "MPa", "elemento", "divergente"),
    "szz": ("Tensão vertical", "MPa", "elemento", "divergente"),
    "d": ("Dano efetivo", "", "elemento", "dano"),
    "dt": ("Dano de tração", "", "elemento", "dano"),
    "dc": ("Dano de compressão", "", "elemento", "dano"),
    "folga": ("Folga inicial das juntas", "mm", "elemento", "viridis"),
    "material": ("Material", "", "elemento", "categorias"),
}


def _lista(v, fmt="{:.4g}"):
    v = np.asarray(v, float).ravel()
    return "[" + ",".join("null" if not np.isfinite(x) else fmt.format(x) for x in v) + "]"


def _lista_int(v):
    return "[" + ",".join(str(int(x)) for x in np.asarray(v).ravel()) + "]"


def exportar(caminho: str, m, estados: List[Dict], descricao: str = "", material=None):
    """estados: [{"rotulo": str, "campos": {id: array}, "u": (n, 3) em m ou None}]
    A cunha W15 vai como W6 (só os nós de canto): as arestas são retas e o
    visualizador pinta os campos nodais nos cantos."""
    faces = None
    sub = None
    if isinstance(m, Malha3D):
        con = m.cun
        if m.tipo == "W15":
            sub = np.unique(m.cun[:, :6])
            novo = -np.ones(m.n_nos, int); novo[sub] = np.arange(len(sub))
            con = novo[m.cun[:, :6]]
            xyz = m.xyz[sub] * 1e3
        else:
            xyz = m.xyz * 1e3
        tipo, npe = "3d", con.shape[1]
        if m.tipo == "H8":
            from .elementos3d import tipo as tipo_elemento
            faces = tipo_elemento("H8").faces_vis
        simetria, larg = m.m2d.simetria, m.m2d.largura_modelo * 1e3
        mat = m.mat
    else:
        xy = m.xy[:, ::-1] if getattr(m, "vertical", False) else m.xy   # seção vertical: espessura × altura
        xyz, con, tipo, npe = np.c_[xy, np.zeros(m.n_nos)] * 1e3, m.tri, "2d", 3
        simetria, larg = m.simetria, m.largura_modelo * 1e3
        mat = np.where(m.marca == 0, 0, 1) if material is None else np.asarray(material)
    usados = sorted({k for e in estados for k in e["campos"]})
    campos = [{"id": k, "nome": CAMPOS.get(k, (k,))[0], "unidade": CAMPOS.get(k, ("", ""))[1],
               "local": CAMPOS.get(k, ("", "", "no"))[2], "mapa": CAMPOS.get(k, ("", "", "", "viridis"))[3]}
              for k in usados]
    nn = m.n_nos

    def nodal(v):
        v = np.asarray(v)
        return v[sub] if (sub is not None and len(v) == nn) else v
    partes = ['{"tipo":"%s","descricao":%s,"simetria":%s,"largura_modelo":%.6g,"nos_por_elemento":%d,'
              % (tipo, json.dumps(descricao, ensure_ascii=False), "true" if simetria else "false", larg, npe)
              + ('"faces":%s,' % json.dumps(faces) if faces else ""),
              '"nomes_materiais":' + json.dumps({int(k): v for k, v in NOMES.items()}, ensure_ascii=False) + ",",
              '"nos":' + _lista(xyz, "{:.6g}") + ",",
              '"elementos":' + _lista_int(con) + ",",
              '"material":' + _lista_int(mat) + ",",
              '"campos":' + json.dumps(campos, ensure_ascii=False) + ',"estados":[']
    for i, e in enumerate(estados):
        s = '{"rotulo":%s,"valores":{' % json.dumps(e["rotulo"], ensure_ascii=False)
        s += ",".join('"%s":%s' % (k, _lista(nodal(v))) for k, v in e["campos"].items())
        s += '},"u":' + ("null" if e.get("u") is None else _lista(nodal(np.asarray(e["u"]).reshape(-1, 3)) * 1e3, "{:.5g}"))
        s += "}" + ("," if i < len(estados) - 1 else "")
        partes.append(s)
    partes.append("]}")
    with open(caminho, "w", encoding="utf-8") as f:
        f.write("".join(partes))


def campos_u(u):
    u = np.asarray(u, float).reshape(-1, 3)
    return {"umag": np.linalg.norm(u, axis=1) * 1e3, "uz": u[:, 2] * 1e3}


def estados_nao_linear(mod, max_carga: int = 24) -> List[Dict]:
    """Estados do modelo não linear: campos do aquecimento e passos da carga
    (amostrados, sempre com o pico e o último)."""
    inst = mod.instantaneos or []
    aq = [s for s in inst if s["fase"] == "aquecimento"]
    cg = [s for s in inst if s["fase"] == "carga"]
    sel = list(range(len(cg)))
    if len(cg) > max_carga:
        ip = int(np.argmax([s["carga_N"] for s in cg]))
        sel = sorted(set(np.linspace(0, len(cg) - 1, max_carga).round().astype(int)) | {ip, len(cg) - 1})
    ne = len(mod.m.cun)
    T_ult = mod.T

    def por_elemento(v):
        out = np.full(ne, np.nan)
        out[mod.id_glob[mod.ic]] = v
        return out
    est = []
    for s in aq:
        T = s["T"] if s["T"] is not None else T_ult
        est.append({"rotulo": f"aquecimento, t = {s['rotulo'] / 60:.0f} min",
                    "campos": {"T": T, **campos_u(s["u"]), "dt": por_elemento(s["dt"]), "dc": por_elemento(s["dc"]),
                               "d": por_elemento(s["d"])}, "u": s["u"]})
    ip = int(np.argmax([s["carga_N"] for s in cg])) if cg else -1
    for k in sel:
        s = cg[k]
        r = f"carga {s['carga_N'] / 1e3:.1f} kN, δ = {s['rotulo'] * 1e3:.3f} mm" + (" (pico)" if k == ip else "")
        est.append({"rotulo": r, "campos": {"T": T_ult, **campos_u(s["u"]), "dt": por_elemento(s["dt"]),
                                           "dc": por_elemento(s["dc"]), "d": por_elemento(s["d"])}, "u": s["u"]})
    return est
