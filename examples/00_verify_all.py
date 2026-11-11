"""Runs all the verifications of the program and summarises the outcome.

Each script compares the model with exact solutions, with the 3D model or with tests,
and writes its CSV to results/. This script calls all of them, in order, and writes
results/verification_summary.csv with one line per script.

Usage:
    python examples/00_verify_all.py           (fast, about 7 min)
    python examples/00_verify_all.py --full    (with the long cases, about 45 min)

Exits with code 1 if any case fails.
"""
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from mef_masonry.output import save_csv  # noqa: E402

FULL = "--full" in sys.argv

SCRIPTS = [
    ("01", "01_verification.py", [], "Thermal and linear thermoelastic (W6)"),
    ("01 H8", "01_verification.py", ["--element", "H8"], "The same cases with the H8 hexahedron"),
    ("01 W15", "01_verification.py", ["--element", "W15"], "The same cases with the quadratic W15 wedge"),
    ("06", "06_verification_nonlinear.py", [], "Non-linear: damage and joints with friction (W6)"),
    ("06 H8", "06_verification_nonlinear.py", ["--element", "H8"], "Non-linear with H8"),
    ("06 W15", "06_verification_nonlinear.py", ["--element", "W15"], "Non-linear with W15"),
    ("10", "10_verification_wall.py", [] if FULL else ["--fast"], "Wall mode"),
    ("11", "11_verification_mortar_materials.py", [], "c(T), mortar joint and library"),
    ("12", "12_verification_cells_grout.py", ["--full"] if FULL else [], "Rectangular cells and grout"),
    ("13", "13_verification_horizontal_holes.py", ["--full"] if FULL else [], "Horizontal holes"),
    ("14", "14_verification_creep.py", ["--full"] if FULL else [], "Transient creep"),
    ("15", "15_verification_sensitivity.py", ["--full"] if FULL else [], "Sensitivity (Morris)"),
]

PATTERNS = [re.compile(r"(\d+) of (\d+) checks OK"), re.compile(r"(\d+)/(\d+) cases within tolerance")]


def output_summary(txt_):
    for p in PATTERNS:
        found = p.findall(txt_)
        if found:
            ok, total = found[-1]
            return int(ok), int(total)
    return None, None


lines = []
failed = False
t_ini = time.time()
print(f"Verifications of mef_masonry ({'full' if FULL else 'fast'})\n")
for lbl, fname, args, desc in SCRIPTS:
    t0 = time.time()
    proc = subprocess.run([sys.executable, os.path.join(HERE, fname)] + args, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", cwd=os.path.join(HERE, ".."))
    dt = time.time() - t0
    ok, total = output_summary(proc.stdout)
    status = "ok" if (proc.returncode == 0 and ok is not None and ok == total) else "FAILED"
    failed |= status != "ok"
    lines.append([lbl, desc, f"{ok}/{total}" if ok is not None else "—", round(dt, 1), status])
    print(f"{'ok ' if status == 'ok' else 'XX '} {lbl:7s} {desc:45s} {ok}/{total} cases  ({dt:5.0f} s)")
    if status != "ok":
        print("    last lines of the output:")
        for l in (proc.stdout + proc.stderr).strip().splitlines()[-8:]:
            print("    " + l)

output = os.path.join(HERE, "..", "results", "verification_summary.csv")
save_csv(output, lines, ["script", "what it verifies", "cases ok", "time (s)", "status"])
n_cases = sum(int(l[2].split("/")[1]) for l in lines if "/" in l[2])
print(f"\n{sum(l[-1] == 'ok' for l in lines)}/{len(lines)} scripts without failures ({n_cases} cases) in "
      f"{(time.time() - t_ini) / 60:.1f} min. Written to {os.path.abspath(output)}")
sys.exit(1 if failed else 0)
