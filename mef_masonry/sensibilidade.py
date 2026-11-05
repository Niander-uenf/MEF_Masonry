"""Análise de sensibilidade global pelo método de Morris (efeitos elementares).

Morris (1991), com a medida μ* de Campolongo, Cariboni e Saltelli (2007):
cada trajetória parte de um ponto da grade de p níveis no hipercubo [0, 1]^k e
muda um parâmetro de cada vez, em ordem aleatória, por Δ = p / [2 (p - 1)].
O efeito elementar do parâmetro i é EE_i = [y(x + Δ e_i) - y(x)] / (±Δ), na
escala da faixa inteira do parâmetro: μ* ≈ quanto a saída muda, em média, quando
o parâmetro percorre a sua faixa; σ mede a não linearidade e as interações.

O estudo aplicado é o da parede no ISO 834 (projeto.avaliar_parede): para cada
ponto, o projeto é montado com os valores dos parâmetros e roda sem figuras.
"""
import copy
import os
import time
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Parâmetros incertos da parede (faixas padrão e fontes na interface)
# ---------------------------------------------------------------------------
# tipo "rel": faixa em ± % do valor do projeto; "abs": mínimo e máximo absolutos
PARAMETROS = [
    {"chave": "k", "nome": "Condutividade térmica do bloco", "curto": "k", "tipo": "rel", "unidade": "%",
     "padrao": (20.0,)},
    {"chave": "c", "nome": "Calor específico do bloco", "curto": "c", "tipo": "rel", "unidade": "%",
     "padrao": (15.0,)},
    {"chave": "eps_furos", "nome": "Emissividade das paredes dos furos", "curto": "ε furos", "tipo": "abs",
     "unidade": "", "padrao": (0.70, 0.95)},
    {"chave": "h", "nome": "Convecção na face exposta", "curto": "h", "tipo": "abs", "unidade": "W/m²·K",
     "padrao": (15.0, 35.0)},
    {"chave": "eps", "nome": "Emissividade da superfície", "curto": "ε face", "tipo": "abs", "unidade": "",
     "padrao": (0.60, 0.95)},
    {"chave": "alfa", "nome": "Dilatação térmica do bloco", "curto": "α", "tipo": "abs", "unidade": "×10⁻⁶/°C",
     "padrao": (4.0, 8.0)},
    {"chave": "ktr", "nome": "Fluência transitória do bloco (k_tr)", "curto": "k_tr", "tipo": "abs", "unidade": "",
     "padrao": (0.0, 2.35)},
    {"chave": "res", "nome": "Resistência da alvenaria (fator na lei)", "curto": "f_res", "tipo": "rel",
     "unidade": "%", "padrao": (15.0,)},
    {"chave": "rig", "nome": "Rigidez da alvenaria (fator na lei)", "curto": "f_rig", "tipo": "rel", "unidade": "%",
     "padrao": (25.0,)},
    {"chave": "imp", "nome": "Imperfeição inicial (H dividida por)", "curto": "imperf.", "tipo": "abs", "unidade": "",
     "padrao": (1000.0, 250.0)},
    {"chave": "exc", "nome": "Excentricidade da carga (topo e base)", "curto": "excentr.", "tipo": "abs",
     "unidade": "mm", "padrao": (-10.0, 10.0)},
]
POR_CHAVE = {d["chave"]: d for d in PARAMETROS}
SAIDAS = [("R", "tempo de resistência", "min"), ("I", "tempo de isolamento", "min"),
          ("Pmax", "capacidade a frio", "kN")]


def ativos(sens: Dict) -> List[Dict]:
    """Parâmetros marcados na página Sensibilidade, com a faixa (lo, hi) informada."""
    out = []
    for d in PARAMETROS:
        ch = d["chave"]
        if not sens.get(f"var_{ch}", True):
            continue
        if d["tipo"] == "rel":
            f = sens.get(f"faixa_{ch}", d["padrao"][0])
            out.append({**d, "lo": -f, "hi": f})
        else:
            out.append({**d, "lo": sens.get(f"min_{ch}", d["padrao"][0]), "hi": sens.get(f"max_{ch}", d["padrao"][1])})
    return out


def valor(par: Dict, x: float) -> float:
    """Valor físico do parâmetro na posição x ∈ [0, 1] da faixa. Na imperfeição, a
    amplitude H/n varia linearmente (n vai de lo a hi)."""
    lo, hi = par["lo"], par["hi"]
    if par["chave"] == "imp":
        a = (1 - x) / lo + x / hi
        return 1.0 / a
    return lo + (hi - lo) * x


def aplicar(base: Dict, pars: Sequence[Dict], x: Sequence[float]) -> Dict:
    """Projeto com os parâmetros na posição x (cópia)."""
    q = copy.deepcopy(base)
    c, en, mo = q["ceramica"], q["ensaio"], q["modelo"]
    for par, xi in zip(pars, x):
        v = valor(par, float(xi))
        ch = par["chave"]
        if ch in ("k", "c"):
            f = 1 + v / 100
            c[ch] = c[ch] * f
            tab = f"tabela_{ch}"
            if c.get(tab):
                c[tab] = [[T, val * f] for T, val in c[tab]]
        elif ch == "eps_furos":
            en["emissividade_furos"] = v
        elif ch == "h":
            en["h_conveccao"] = v
        elif ch == "eps":
            en["emissividade"] = v
        elif ch == "alfa":
            if c.get("tabela_eps_termica") and c.get("alfa", 0) > 0:
                f = v / c["alfa"]
                c["tabela_eps_termica"] = [[T, e * f] for T, e in c["tabela_eps_termica"]]
            c["alfa"] = v
        elif ch == "ktr":
            mo["fluencia"] = True
            mo["k_tr_bloco"] = v
            mo["tabela_fluencia_bloco"] = []
        elif ch == "res":
            mo["fator_resistencia_parede"] = mo.get("fator_resistencia_parede", 1.0) * (1 + v / 100)
        elif ch == "rig":
            mo["fator_rigidez_parede"] = mo.get("fator_rigidez_parede", 1.0) * (1 + v / 100)
        elif ch == "imp":
            en["imperfeicao"] = v
        elif ch == "exc":
            en["excentricidade_topo"] = en["excentricidade_base"] = v
    return q


# ---------------------------------------------------------------------------
# Plano de Morris e efeitos elementares
# ---------------------------------------------------------------------------
def trajetorias(k: int, r: int, p: int = 4, semente: int = 1) -> np.ndarray:
    """r trajetórias de k + 1 pontos no hipercubo [0, 1]^k (grade de p níveis)."""
    rng = np.random.default_rng(semente)
    niveis = np.arange(p) / (p - 1)
    delta = p / (2 * (p - 1))
    X = np.empty((r, k + 1, k))
    for j in range(r):
        x = rng.choice(niveis, size=k)
        X[j, 0] = x
        for s, i in enumerate(rng.permutation(k)):
            x = x.copy()
            x[i] = x[i] + delta if x[i] + delta <= 1 + 1e-12 else x[i] - delta
            X[j, s + 1] = x
    return X


def efeitos(X: np.ndarray, Y: np.ndarray, n_boot: int = 1000, semente: int = 7) -> Dict[str, np.ndarray]:
    """μ*, μ, σ e o intervalo de 95% de μ* (bootstrap nas trajetórias) para a
    saída Y (r, k + 1). Pares com NaN ficam de fora."""
    r, n, k = X.shape
    EE = np.full((r, k), np.nan)
    for j in range(r):
        for s in range(k):
            dx = X[j, s + 1] - X[j, s]
            i = int(np.argmax(np.abs(dx)))
            EE[j, i] = (Y[j, s + 1] - Y[j, s]) / dx[i]
    a = np.abs(EE)
    mu_s = np.nanmean(a, axis=0)
    mu = np.nanmean(EE, axis=0)
    sig = np.array([np.nanstd(EE[:, i], ddof=1) if np.sum(np.isfinite(EE[:, i])) > 1 else np.nan
                    for i in range(k)])
    rng = np.random.default_rng(semente)
    bs = np.array([np.nanmean(a[rng.integers(0, r, r)], axis=0) for _ in range(n_boot)])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    return {"EE": EE, "mu_estrela": mu_s, "mu": mu, "sigma": sig, "ic_lo": lo, "ic_hi": hi,
            "validos": np.sum(np.isfinite(EE), axis=0)}


# ---------------------------------------------------------------------------
# Execução das rodadas (em paralelo quando possível)
# ---------------------------------------------------------------------------
def _avaliar(args):
    import contextlib
    import io
    from .projeto import avaliar_parede
    q, lei = args
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return avaliar_parede(q, lei)
    except Exception as e:  # noqa: BLE001 - uma rodada que falha vira NaN
        return {"erro": f"{type(e).__name__}: {e}"}


def processos_padrao():
    """Núcleos menos um (a interface continua respondendo); com 2 núcleos, os 2."""
    n = os.cpu_count() or 1
    return n if n <= 2 else n - 1


_UMA_LINHA = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def rodar(projetos: Sequence[Dict], lei, processos: int = 0, log: Callable[[str], None] = print) -> List[Dict]:
    """Avalia os projetos (paralelo com processos > 1; 0 = automático)."""
    n = len(projetos)
    procs = processos or processos_padrao()
    res: List[Optional[Dict]] = [None] * n
    t0 = time.time()

    def progresso(i_feitas):
        dt = time.time() - t0
        falta = dt / i_feitas * (n - i_feitas)
        log(f"  {i_feitas}/{n} rodadas ({dt / 60:.1f} min; faltam cerca de {falta / 60:.0f} min)")

    pool = None
    if procs > 1 and n > 1:
        # processos novos ("spawn", como no Windows): o servidor tem outras linhas de execução. Cada
        # processo usa uma linha só do BLAS, para não disputar os núcleos; a variável vale até os
        # processos nascerem, nos submit
        antes = {k: os.environ.get(k) for k in _UMA_LINHA}
        os.environ.update({k: "1" for k in _UMA_LINHA})
        try:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor, as_completed
            pool = ProcessPoolExecutor(max_workers=min(procs, n), mp_context=mp.get_context("spawn"))
            futs = {pool.submit(_avaliar, (q, lei)): i for i, q in enumerate(projetos)}
        except Exception as e:  # noqa: BLE001
            log(f"  Sem processos em paralelo ({e}); rodando em sequência.")
            if pool is not None:
                pool.shutdown(wait=False, cancel_futures=True)
            pool = None
        finally:
            for k, v in antes.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
    if pool is not None:
        try:
            feitas = 0
            for f in as_completed(futs):
                try:
                    res[futs[f]] = f.result()
                except Exception as e:  # noqa: BLE001 - processo que caiu
                    res[futs[f]] = {"erro": f"{type(e).__name__}: {e}"}
                feitas += 1
                if feitas % max(1, n // 20) == 0 or feitas == n:
                    progresso(feitas)
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
    else:
        for i, q in enumerate(projetos):
            res[i] = _avaliar((q, lei))
            if (i + 1) % max(1, n // 20) == 0 or i + 1 == n:
                progresso(i + 1)
    return res
