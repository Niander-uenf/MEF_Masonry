"""Local server of the program (Python standard library only).

Usage: python app.py [--port 8765] [--folder results_app] [--no-browser]
Opens http://127.0.0.1:<port>/ in the browser. One analysis runs at a time, in a
separate thread; the page follows the log and shows the
results and the post-processor.
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

from .. import project as pj, __version__

WEB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
FOLDER = os.path.abspath("results_app")


class Task:
    def __init__(self, p, folder):
        self.id = uuid.uuid4().hex[:10]
        self.project, self.folder = p, folder
        self.state = "queued"
        self.log = []
        self.begin = time.time()
        self.end = None
        self.error = None
        self.cancel = threading.Event()
        self._partial = ""

    # file-like object to redirect print() of the solver
    def write(self, s):
        if self.cancel.is_set():
            raise pj.Cancelled()
        self._partial += s
        while "\n" in self._partial:
            line, self._partial = self._partial.split("\n", 1)
            if line.strip():
                self.log.append(line.rstrip())
        return len(s)

    def flush(self):
        pass

    def line(self, s):
        self.write(str(s) + "\n")

    def info(self, since=0):
        return {"id": self.id, "state": self.state, "log": self.log[since:], "n": len(self.log),
                "time_s": round((self.end or time.time()) - self.begin, 1), "error": self.error,
                "folder": os.path.basename(self.folder)}


TASKS = {}
CURRENT = {"task": None}
LOCK = threading.Lock()


def _run(t: Task):
    t.state = "running"
    try:
        with contextlib.redirect_stdout(t):
            pj.run(t.project, t.folder, log=t.line)
        t.state = "done"
    except pj.Cancelled:
        t.state = "cancelled"
        t.log.append("Analysis cancelled.")
    except Exception as e:  # noqa
        t.state = "error"
        t.error = str(e)
        t.log.extend(traceback.format_exc().splitlines()[-6:])
    finally:
        t.end = time.time()


def _slug(s):
    s = re.sub(r"[^\w\-]+", "_", s.strip(), flags=re.UNICODE).strip("_")
    return s[:40] or "project"


def _history():
    out = []
    if not os.path.isdir(FOLDER):
        return out
    for name in sorted(os.listdir(FOLDER), reverse=True):
        d = os.path.join(FOLDER, name)
        f = os.path.join(d, "results.json")
        pf = os.path.join(d, "project.json")
        old = not os.path.exists(pf)          # run of version 1.x: only its inputs can be reused
        if old:
            pf = os.path.join(d, "projeto.json")
        if not os.path.isdir(d) or not os.path.exists(pf):
            continue
        try:
            p = pj.complete(json.load(open(pf, encoding="utf-8")))
        except Exception:  # noqa
            continue
        out.append({"folder": name, "name": p.get("name", name), "analysis": p.get("model", {}).get("analysis"),
                    "test": p.get("test", {}).get("kind"), "done": os.path.exists(f) and not old,
                    "visual": os.path.exists(os.path.join(d, "visual.json")) and not old})
    return out


class Handler(BaseHTTPRequestHandler):
    server_version = "mef_masonry/" + __version__

    def log_message(self, fmt, *args):  # silent
        pass

    # -- utilities
    def _json(self, obj, code=200):
        data = json.dumps(obj, ensure_ascii=False, default=float).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def _file(self, path, download=False):
        if not os.path.isfile(path):
            return self._json({"error": "file not found"}, 404)
        kind = mimetypes.guess_type(path)[0] or "application/octet-stream"
        if path.endswith(".js"):
            kind = "text/javascript"
        nbytes = os.path.getsize(path)
        self.send_response(200)
        self.send_header("Content-Type", kind + ("; charset=utf-8" if kind.startswith("text") or kind.endswith("json") else ""))
        self.send_header("Content-Length", str(nbytes))
        self.send_header("Cache-Control", "no-store")
        if download:
            self.send_header("Content-Disposition", f'attachment; filename="{os.path.basename(path)}"')
        self.end_headers()
        with open(path, "rb") as f:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                self.wfile.write(b)

    def _body(self):
        n = int(self.headers.get("Content-Length", 0) or 0)
        return json.loads(self.rfile.read(n).decode("utf-8")) if n else {}

    def _result_folder(self, name):
        name = os.path.basename(unquote(name))
        d = os.path.join(FOLDER, name)
        return d if os.path.isdir(d) else None

    # -- routes
    def do_GET(self):
        u = urlparse(self.path)
        parts = [unquote(x) for x in u.path.strip("/").split("/") if x]
        q = parse_qs(u.query)
        try:
            if not parts:
                return self._file(os.path.join(WEB, "index.html"))
            if parts[0] == "web" and len(parts) == 2:
                return self._file(os.path.join(WEB, os.path.basename(parts[1])))
            if parts[:2] == ["api", "schema"]:
                default = pj.default_project()
                return self._json({"schema": pj.SCHEMA, "default": default, "version": __version__,
                                   "presets": {k: pj.apply(default, v) for k, v in pj.PRESETS.items()},
                                   "library": pj._lib.for_interface(), "folder": FOLDER})
            if parts[:2] == ["api", "task"] and len(parts) == 3:
                t = TASKS.get(parts[2])
                if not t:
                    return self._json({"error": "unknown task"}, 404)
                return self._json(t.info(int(q.get("since", ["0"])[0])))
            if parts[:2] == ["api", "current"]:
                t = CURRENT["task"]
                return self._json(t.info() if t else {})
            if parts[:2] == ["api", "history"]:
                return self._json({"runs": _history(), "folder": FOLDER})
            if parts[:2] == ["api", "result"] and len(parts) == 3:
                d = self._result_folder(parts[2])
                if not d:
                    return self._json({"error": "result not found"}, 404)
                f = os.path.join(d, "results.json")
                if not os.path.exists(f):
                    return self._json({"error": "analysis without results (interrupted?)"}, 404)
                return self._json(json.load(open(f, encoding="utf-8")))
            if parts[0] == "files" and len(parts) == 3:
                d = self._result_folder(parts[1])
                if not d:
                    return self._json({"error": "folder not found"}, 404)
                return self._file(os.path.join(d, os.path.basename(parts[2])), download="download" in q)
            return self._json({"error": "unknown route"}, 404)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception as e:  # noqa
            return self._json({"error": str(e)}, 500)

    def do_POST(self):
        u = urlparse(self.path)
        parts = [x for x in u.path.strip("/").split("/") if x]
        try:
            if parts[:2] == ["api", "complete"]:      # project file (also of version 1.x) -> current format
                return self._json({"project": pj.complete(self._body())})
            if parts[:2] == ["api", "validate"]:
                p = pj.complete(self._body())
                return self._json({"errors": pj.validate(p)})
            if parts[:2] == ["api", "mesh"]:
                p = pj.complete(self._body())
                errors = [e for e in pj.validate(p) if "Geometry" in e or "hole" in e.lower() or "unit" in e.lower()]
                if errors:
                    return self._json({"errors": errors})
                png, info = pj.mesh_figure(p)
                return self._json({"png": base64.b64encode(png).decode(), "summary": info, "errors": []})
            if parts[:2] == ["api", "run"]:
                p = pj.complete(self._body())
                errors = pj.validate(p)
                if errors:
                    return self._json({"errors": errors}, 400)
                with LOCK:
                    t0 = CURRENT["task"]
                    if t0 and t0.state in ("queued", "running"):
                        return self._json({"errors": ["An analysis is already running; wait or cancel it."]}, 409)
                    name = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + _slug(p.get("name", "project"))
                    t = Task(p, os.path.join(FOLDER, name))
                    TASKS[t.id] = t
                    CURRENT["task"] = t
                threading.Thread(target=_run, args=(t,), daemon=True).start()
                return self._json({"id": t.id, "folder": name})
            if parts[:2] == ["api", "cancel"] and len(parts) == 3:
                t = TASKS.get(parts[2])
                if t:
                    t.cancel.set()
                return self._json({"ok": True})
            return self._json({"error": "unknown route"}, 404)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception as e:  # noqa
            return self._json({"errors": [str(e)]}, 500)


def _free_port(preferred):
    for p in [preferred] + list(range(preferred + 1, preferred + 50)):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue
    return 0


def main(argv=None):
    global FOLDER
    ap = argparse.ArgumentParser(description="FEM program for masonry in fire")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--folder", default="results_app", help="folder where the analyses are stored")
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args(argv)
    FOLDER = os.path.abspath(a.folder)
    os.makedirs(FOLDER, exist_ok=True)
    port = _free_port(a.port)
    srv = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"mef_masonry {__version__}: open {url} in the browser (Ctrl+C quits). Results in {FOLDER}", flush=True)
    if not a.no_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("Stopped.")
    finally:
        srv.server_close()


if __name__ == "__main__":
    main()
