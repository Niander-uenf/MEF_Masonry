"""Servidor local do programa (só biblioteca padrão do Python).

Uso: python app.py [--porta 8765] [--pasta resultados_app] [--sem-navegador]
Abre http://127.0.0.1:<porta>/ no navegador. Uma análise roda por vez, numa
linha de execução separada; a página acompanha o registro e mostra os
resultados e o pós-processador.
"""
import argparse
import base64
import contextlib
import datetime
import json
import mimetypes
import os
import re
import socket
import threading
import time
import traceback
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs, unquote

from .. import projeto as pj

WEB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
PASTA = os.path.abspath("resultados_app")


class Tarefa:
    def __init__(self, p, pasta):
        self.id = uuid.uuid4().hex[:10]
        self.projeto, self.pasta = p, pasta
        self.estado = "na fila"
        self.log = []
        self.inicio = time.time()
        self.fim = None
        self.erro = None
        self.cancelar = threading.Event()
        self._parcial = ""

    # objeto tipo arquivo para redirecionar print() do solver
    def write(self, s):
        if self.cancelar.is_set():
            raise pj.Cancelado()
        self._parcial += s
        while "\n" in self._parcial:
            linha, self._parcial = self._parcial.split("\n", 1)
            if linha.strip():
                self.log.append(linha.rstrip())
        return len(s)

    def flush(self):
        pass

    def linha(self, s):
        self.write(str(s) + "\n")

    def info(self, desde=0):
        return {"id": self.id, "estado": self.estado, "log": self.log[desde:], "n": len(self.log),
                "tempo_s": round((self.fim or time.time()) - self.inicio, 1), "erro": self.erro,
                "pasta": os.path.basename(self.pasta)}


TAREFAS = {}
ATUAL = {"tarefa": None}
TRAVA = threading.Lock()


def _rodar(t: Tarefa):
    t.estado = "rodando"
    try:
        with contextlib.redirect_stdout(t):
            pj.executar(t.projeto, t.pasta, log=t.linha)
        t.estado = "concluída"
    except pj.Cancelado:
        t.estado = "cancelada"
        t.log.append("Análise cancelada.")
    except Exception as e:  # noqa
        t.estado = "erro"
        t.erro = str(e)
        t.log.extend(traceback.format_exc().splitlines()[-6:])
    finally:
        t.fim = time.time()


def _slug(s):
    s = re.sub(r"[^\w\-]+", "_", s.strip(), flags=re.UNICODE).strip("_")
    return s[:40] or "projeto"


def _historico():
    out = []
    if not os.path.isdir(PASTA):
        return out
    for nome in sorted(os.listdir(PASTA), reverse=True):
        d = os.path.join(PASTA, nome)
        f = os.path.join(d, "resultados.json")
        pf = os.path.join(d, "projeto.json")
        if not os.path.isdir(d) or not os.path.exists(pf):
            continue
        try:
            p = json.load(open(pf, encoding="utf-8"))
        except Exception:  # noqa
            continue
        out.append({"pasta": nome, "nome": p.get("nome", nome), "analise": p.get("modelo", {}).get("analise"),
                    "ensaio": p.get("ensaio", {}).get("tipo"), "concluida": os.path.exists(f),
                    "visual": os.path.exists(os.path.join(d, "visual.json"))})
    return out


class Manipulador(BaseHTTPRequestHandler):
    server_version = "mef_masonry/" + pj.VERSAO

    def log_message(self, fmt, *args):  # silencioso
        pass

    # -- utilitários
    def _json(self, obj, codigo=200):
        dados = json.dumps(obj, ensure_ascii=False, default=float).encode("utf-8")
        self.send_response(codigo)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(dados)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(dados)

    def _arquivo(self, caminho, baixar=False):
        if not os.path.isfile(caminho):
            return self._json({"erro": "arquivo não encontrado"}, 404)
        tipo = mimetypes.guess_type(caminho)[0] or "application/octet-stream"
        if caminho.endswith(".js"):
            tipo = "text/javascript"
        tam = os.path.getsize(caminho)
        self.send_response(200)
        self.send_header("Content-Type", tipo + ("; charset=utf-8" if tipo.startswith("text") or tipo.endswith("json") else ""))
        self.send_header("Content-Length", str(tam))
        self.send_header("Cache-Control", "no-store")
        if baixar:
            self.send_header("Content-Disposition", f'attachment; filename="{os.path.basename(caminho)}"')
        self.end_headers()
        with open(caminho, "rb") as f:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                self.wfile.write(b)

    def _corpo(self):
        n = int(self.headers.get("Content-Length", 0) or 0)
        return json.loads(self.rfile.read(n).decode("utf-8")) if n else {}

    def _pasta_resultado(self, nome):
        nome = os.path.basename(unquote(nome))
        d = os.path.join(PASTA, nome)
        return d if os.path.isdir(d) else None

    # -- rotas
    def do_GET(self):
        u = urlparse(self.path)
        partes = [unquote(x) for x in u.path.strip("/").split("/") if x]
        q = parse_qs(u.query)
        try:
            if not partes:
                return self._arquivo(os.path.join(WEB, "index.html"))
            if partes[0] == "web" and len(partes) == 2:
                return self._arquivo(os.path.join(WEB, os.path.basename(partes[1])))
            if partes[:2] == ["api", "esquema"]:
                padrao = pj.projeto_padrao()
                return self._json({"esquema": pj.ESQUEMA, "padrao": padrao, "versao": pj.VERSAO,
                                   "predefinidos": {k: pj.aplicar(padrao, v) for k, v in pj.PREDEFINIDOS.items()},
                                   "biblioteca": pj._bib.para_interface(), "pasta": PASTA})
            if partes[:2] == ["api", "tarefa"] and len(partes) == 3:
                t = TAREFAS.get(partes[2])
                if not t:
                    return self._json({"erro": "tarefa desconhecida"}, 404)
                return self._json(t.info(int(q.get("desde", ["0"])[0])))
            if partes[:2] == ["api", "atual"]:
                t = ATUAL["tarefa"]
                return self._json(t.info() if t else {})
            if partes[:2] == ["api", "historico"]:
                return self._json({"execucoes": _historico(), "pasta": PASTA})
            if partes[:2] == ["api", "resultado"] and len(partes) == 3:
                d = self._pasta_resultado(partes[2])
                if not d:
                    return self._json({"erro": "resultado não encontrado"}, 404)
                f = os.path.join(d, "resultados.json")
                if not os.path.exists(f):
                    return self._json({"erro": "análise sem resultados (interrompida?)"}, 404)
                return self._json(json.load(open(f, encoding="utf-8")))
            if partes[0] == "arquivos" and len(partes) == 3:
                d = self._pasta_resultado(partes[1])
                if not d:
                    return self._json({"erro": "pasta não encontrada"}, 404)
                return self._arquivo(os.path.join(d, os.path.basename(partes[2])), baixar="baixar" in q)
            return self._json({"erro": "rota desconhecida"}, 404)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception as e:  # noqa
            return self._json({"erro": str(e)}, 500)

    def do_POST(self):
        u = urlparse(self.path)
        partes = [x for x in u.path.strip("/").split("/") if x]
        try:
            if partes[:2] == ["api", "validar"]:
                p = pj.completar(self._corpo())
                return self._json({"erros": pj.validar(p)})
            if partes[:2] == ["api", "malha"]:
                p = pj.completar(self._corpo())
                erros = [e for e in pj.validar(p) if "Geometria" in e or "furo" in e.lower() or "bloco" in e.lower()]
                if erros:
                    return self._json({"erros": erros})
                png, info = pj.figura_malha(p)
                return self._json({"png": base64.b64encode(png).decode(), "resumo": info, "erros": []})
            if partes[:2] == ["api", "executar"]:
                p = pj.completar(self._corpo())
                erros = pj.validar(p)
                if erros:
                    return self._json({"erros": erros}, 400)
                with TRAVA:
                    t0 = ATUAL["tarefa"]
                    if t0 and t0.estado in ("na fila", "rodando"):
                        return self._json({"erros": ["Já há uma análise rodando; aguarde ou cancele."]}, 409)
                    nome = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + _slug(p.get("nome", "projeto"))
                    t = Tarefa(p, os.path.join(PASTA, nome))
                    TAREFAS[t.id] = t
                    ATUAL["tarefa"] = t
                threading.Thread(target=_rodar, args=(t,), daemon=True).start()
                return self._json({"id": t.id, "pasta": nome})
            if partes[:2] == ["api", "cancelar"] and len(partes) == 3:
                t = TAREFAS.get(partes[2])
                if t:
                    t.cancelar.set()
                return self._json({"ok": True})
            return self._json({"erro": "rota desconhecida"}, 404)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception as e:  # noqa
            return self._json({"erros": [str(e)]}, 500)


def _porta_livre(preferida):
    for p in [preferida] + list(range(preferida + 1, preferida + 50)):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue
    return 0


def principal(argv=None):
    global PASTA
    ap = argparse.ArgumentParser(description="Programa MEF para prismas de alvenaria seca em temperatura")
    ap.add_argument("--porta", type=int, default=8765)
    ap.add_argument("--pasta", default="resultados_app", help="pasta onde ficam as análises")
    ap.add_argument("--sem-navegador", action="store_true")
    a = ap.parse_args(argv)
    PASTA = os.path.abspath(a.pasta)
    os.makedirs(PASTA, exist_ok=True)
    porta = _porta_livre(a.porta)
    srv = ThreadingHTTPServer(("127.0.0.1", porta), Manipulador)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"mef_masonry {pj.VERSAO}: abra {url} no navegador (Ctrl+C encerra). Resultados em {PASTA}", flush=True)
    if not a.sem_navegador:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("Encerrado.")
    finally:
        srv.server_close()


if __name__ == "__main__":
    principal()
