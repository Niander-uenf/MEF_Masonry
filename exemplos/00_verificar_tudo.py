"""Roda todas as verificações do programa e resume o resultado.

Cada script compara o modelo com soluções exatas, com o 3D ou com ensaios,
e grava o seu CSV em resultados/. Este roteiro chama todos, na ordem, e grava
resultados/resumo_verificacoes.csv com uma linha por script.

Uso:
    python exemplos/00_verificar_tudo.py              (rápido, cerca de 7 min)
    python exemplos/00_verificar_tudo.py --completo   (com os casos longos, cerca de 45 min)

Sai com código 1 se algum caso falhar.
"""
import os
import re
import subprocess
import sys
import time

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, ".."))
from mef_masonry.saida import salvar_csv  # noqa: E402

COMPLETO = "--completo" in sys.argv

ROTEIROS = [
    ("01", "01_verificacao.py", [], "Térmica e termoelástica linear (W6)"),
    ("01 H8", "01_verificacao.py", ["--elemento", "H8"], "Os mesmos casos com o hexaedro H8"),
    ("01 W15", "01_verificacao.py", ["--elemento", "W15"], "Os mesmos casos com a cunha quadrática W15"),
    ("06", "06_verificacao_nao_linear.py", [], "Não linear: dano e juntas com atrito (W6)"),
    ("06 H8", "06_verificacao_nao_linear.py", ["--elemento", "H8"], "Não linear com H8"),
    ("06 W15", "06_verificacao_nao_linear.py", ["--elemento", "W15"], "Não linear com W15"),
    ("10", "10_verificacao_parede.py", [] if COMPLETO else ["--rapido"], "Modo Parede"),
    ("11", "11_verificacao_materiais_argamassa.py", [], "c(T), junta de argamassa e biblioteca"),
    ("12", "12_verificacao_celulas_graute.py", ["--completo"] if COMPLETO else [], "Células retangulares e graute"),
    ("13", "13_verificacao_furos_horizontais.py", ["--completo"] if COMPLETO else [], "Furos horizontais"),
    ("14", "14_verificacao_fluencia.py", ["--completo"] if COMPLETO else [], "Fluência transitória"),
    ("15", "15_verificacao_sensibilidade.py", ["--completo"] if COMPLETO else [], "Sensibilidade (Morris)"),
]

PADROES = [re.compile(r"(\d+) de (\d+) verificações OK"), re.compile(r"(\d+)/(\d+) casos dentro da tolerância")]


def resumo_saida(texto):
    for p in PADROES:
        achados = p.findall(texto)
        if achados:
            ok, total = achados[-1]
            return int(ok), int(total)
    return None, None


linhas = []
falhou = False
t_ini = time.time()
print(f"Verificações do mef_masonry ({'completas' if COMPLETO else 'rápidas'})\n")
for rot, arq, args, desc in ROTEIROS:
    t0 = time.time()
    proc = subprocess.run([sys.executable, os.path.join(AQUI, arq)] + args, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", cwd=os.path.join(AQUI, ".."))
    dt = time.time() - t0
    ok, total = resumo_saida(proc.stdout)
    situacao = "ok" if (proc.returncode == 0 and ok is not None and ok == total) else "FALHOU"
    falhou |= situacao != "ok"
    linhas.append([rot, desc, f"{ok}/{total}" if ok is not None else "—", round(dt, 1), situacao])
    print(f"{'ok ' if situacao == 'ok' else 'XX '} {rot:7s} {desc:45s} {ok}/{total} casos  ({dt:5.0f} s)")
    if situacao != "ok":
        print("    últimas linhas da saída:")
        for l in (proc.stdout + proc.stderr).strip().splitlines()[-8:]:
            print("    " + l)

saida = os.path.join(AQUI, "..", "resultados", "resumo_verificacoes.csv")
salvar_csv(saida, linhas, ["roteiro", "o que verifica", "casos ok", "tempo (s)", "situação"])
n_casos = sum(int(l[2].split("/")[1]) for l in linhas if "/" in l[2])
print(f"\n{sum(l[-1] == 'ok' for l in linhas)}/{len(linhas)} roteiros sem falhas ({n_casos} casos) em "
      f"{(time.time() - t_ini) / 60:.1f} min. Gravado em {os.path.abspath(saida)}")
sys.exit(1 if falhou else 0)
