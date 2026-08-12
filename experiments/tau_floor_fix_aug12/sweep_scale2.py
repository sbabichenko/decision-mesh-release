#!/usr/bin/env python3
"""Scale-prior and honest-v_i sweep on the SMM design.

Arms cross {prior: eb (per-depth MMLE+borrow), ebscale (smooth scale)} with
{v_i: stratum (default), smooth POINT_VAR}. off/depth as anchors.
"""
import csv, json, os, shutil, subprocess, time
from collections import defaultdict

ROOT = "/home/user/decision-mesh-release"
WORK = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(WORK, "dm_snapshot7")
PV = os.path.join(WORK, "pointvar_smooth_smm.txt")
SMM = os.path.join(ROOT, "data", "smm_design_202606.csv")

BASE = {"DMESH_SPLIT":"1","DMESH_SPLIT_MODE":"local","DMESH_PL_SOLVER":"1",
        "DMESH_SIGMA0_MIN":"1.0","DMESH_SIGMA0_MAX":"6","DMESH_PI0_MIN":"0.05",
        "DMESH_TAU_FLOOR_LOGIT_SD":"0.005",
        "DMESH_ONELAW":"0.1428,0.2415,0.0286,0.0038;8.5,8.1,5.4,0.0"}
ARMS = {
    "off":              {},
    "eb":               {"DMESH_FINAL_RESHRINK":"1","DMESH_RESHRINK_LADDER":"eb"},
    "eb_pv":            {"DMESH_FINAL_RESHRINK":"1","DMESH_RESHRINK_LADDER":"eb","DMESH_POINT_VAR":PV},
    "ebscale":          {"DMESH_FINAL_RESHRINK":"1","DMESH_RESHRINK_LADDER":"ebscale"},
    "ebscale_free":     {"DMESH_FINAL_RESHRINK":"1","DMESH_RESHRINK_LADDER":"ebscale","DMESH_EB_FREE":"1"},
    "ebscale_free_pv":  {"DMESH_FINAL_RESHRINK":"1","DMESH_RESHRINK_LADDER":"ebscale","DMESH_EB_FREE":"1","DMESH_POINT_VAR":PV},
}
SEEDS = [7] + list(range(101, 111))

def run(seed, arm):
    outdir = os.path.join(WORK, "runs5", f"s{seed}_{arm}")
    shutil.rmtree(outdir, ignore_errors=True); os.makedirs(outdir)
    env = dict(os.environ); env.update(BASE); env.update(ARMS[arm])
    env["DMESH_DATA"] = SMM; env["DMESH_DUMP"] = os.path.join(outdir, "run")
    t0 = time.time()
    p = subprocess.run([BIN,"0","0.10","1","24",str(seed)], env=env,
                       capture_output=True, text=True, cwd=outdir)
    rec = {"seed":seed,"arm":arm,"exit":p.returncode,"wall_s":round(time.time()-t0,1)}
    for l in p.stdout.splitlines():
        if l.startswith("HELDOUT:"):
            rec["heldout_dev"] = float(l.split("|")[0].split("deviance/pool")[1])
    shutil.rmtree(outdir, ignore_errors=True)
    return rec

results = []
for seed in SEEDS:
    for arm in ARMS:
        r = run(seed, arm); results.append(r)
        print(f"s{seed:<4} {arm:12s} exit={r['exit']} dev={r.get('heldout_dev',float('nan')):.4f} wall={r['wall_s']}s", flush=True)
json.dump(results, open(os.path.join(WORK,"results_scale2.json"),"w"), indent=1)
print("wrote results_scale2.json")
