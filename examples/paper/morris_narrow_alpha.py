"""Section 8.1 of the paper: the Morris screening of the BCEPQ wall repeated with the α range
narrowed from 4 to 8 to 5 to 7 ×10⁻⁶/°C (the other ten inputs, the trajectories and the
seed equal to those of the screening in the paper). Shows how much the first place of α depends on the
range assigned to it. ~40 min on two cores."""
import os, sys, json
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, ROOT)
from mef_masonry import project as pj

p = json.load(open(os.path.join(ROOT, "results", "sensitivity_wall_bcepq", "project.json")))
p["sensitivity"].update(min_alpha=5.0, max_alpha=7.0)
p["name"] = "BCEPQ wall — Morris with α from 5 to 7 ×10⁻⁶/°C"
pj.run(p, os.path.join(ROOT, "results", "paper", "morris_alpha_5_7"))
