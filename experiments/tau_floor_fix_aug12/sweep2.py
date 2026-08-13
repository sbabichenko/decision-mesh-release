#!/usr/bin/env python3
"""Ladder-variant sweep: depth (pooled, fixed floor) vs cells vs twogroups vs off."""
import csv, json, os, shutil, subprocess, sys, time
from collections import defaultdict

ROOT = "/home/user/decision-mesh-release"
WORK = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(WORK, "dm_snapshot2")

SMM_ENV = {
    "DMESH_SPLIT": "1", "DMESH_SPLIT_MODE": "local", "DMESH_PL_SOLVER": "1",
    "DMESH_SIGMA0_MIN": "1.0", "DMESH_SIGMA0_MAX": "6", "DMESH_PI0_MIN": "0.05",
    "DMESH_TAU_FLOOR_LOGIT_SD": "0.005",
    "DMESH_ONELAW": "0.1428,0.2415,0.0286,0.0038;8.5,8.1,5.4,0.0",
}
APRIL_ENV = {
    "DMESH_SPLIT": "1", "DMESH_SPLIT_MODE": "local", "DMESH_PL_SOLVER": "1",
    "DMESH_HIER_FIT": "2", "DMESH_HEIGHT_LIMIT_LOGIT": "12",
    "DMESH_IRLS_STEP_CAP": "4", "DMESH_HIER_MAP": "1", "DMESH_QUAD_PRIOR": "1",
    "DMESH_SIGMA0_MIN": "0.3", "DMESH_SIGMA0_MAX": "6", "DMESH_PI0_MIN": "0.05",
    "DMESH_TAU_FLOOR_LOGIT_SD": "0.005", "DMESH_STAGE_COUNT": "3",
    "DMESH_MAX_GATE_ROUNDS": "80", "DMESH_PARENT_VAR_SCALE": "0.25",
    "DMESH_GATE_EXTRA_VAR": "0.35", "DMESH_B3_USE_SCORE": "1",
    "DMESH_B3_STAGE0_EXACT": "0",
}

ARMS = {
    "off":       {},
    "depth":     {"DMESH_FINAL_RESHRINK": "1"},
    "cells":     {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "cells"},
    "cells_rawkeff": {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "cells",
                      "DMESH_TAU_CELL_RAWKEFF": "1"},
    "twogroups": {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "twogroups"},
    "twogroups_free": {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "twogroups",
                       "DMESH_RESHRINK_2G_FREE": "1"},
}

def cases():
    smm = os.path.join(ROOT, "data", "smm_design_202606.csv")
    april = os.path.join(ROOT, "data", "ginnie_design.csv")
    flatnull = os.path.join(ROOT, "tests", "data", "flat_p78_clean_seed0.csv")
    for seed in [7] + list(range(101, 111)):
        yield f"smm_s{seed}", smm, SMM_ENV, seed
    for seed in (7, 101, 103):
        yield f"april_s{seed}", april, APRIL_ENV, seed
    for seed in (7, 103):
        yield f"aprilnull_s{seed}", flatnull, APRIL_ENV, seed

def run_one(case_id, data, env_extra, seed, arm):
    outdir = os.path.join(WORK, "runs2", f"{case_id}_{arm}")
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    env = dict(os.environ)
    env.update(env_extra)
    env.update(ARMS[arm])
    env["DMESH_DATA"] = data
    env["DMESH_DUMP"] = os.path.join(outdir, "run")
    t0 = time.time()
    p = subprocess.run([BIN, "0", "0.10", "1", "24", str(seed)], env=env,
                       capture_output=True, text=True, cwd=outdir)
    rec = {"case": case_id, "arm": arm, "seed": seed, "exit": p.returncode,
           "wall_s": round(time.time() - t0, 1)}
    for line in p.stdout.splitlines():
        if line.startswith("HELDOUT:"):
            parts = line.split("|")
            rec["heldout_dev"] = float(parts[0].split("deviance/pool")[1])
            rec["x2_info"] = float(parts[1].split("X2/info")[1])
    rec["diag"] = [l.strip() for l in p.stderr.splitlines()
                   if "[reshrink]" in l or "[ladder" in l]
    ladder = defaultdict(lambda: [0, 0.0])
    path = os.path.join(outdir, "run_hier_vertices.csv")
    if os.path.exists(path):
        with open(path) as f:
            for r in csv.DictReader(f):
                try:
                    if (r["active"] == "1" and r["gate_admitted"] == "1" and
                            r["free_coefficient"] == "1" and r["parent0"]):
                        d = int(r["depth"])
                        ladder[d][0] += 1
                        ladder[d][1] = max(ladder[d][1], float(r["lambda_v"]))
                except (KeyError, ValueError):
                    pass
    rec["ladder"] = {str(d): {"n": e[0], "lam_max": e[1]}
                     for d, e in sorted(ladder.items())}
    rec["admitted"] = sum(e[0] for e in ladder.values())
    shutil.rmtree(outdir, ignore_errors=True)
    return rec

def main():
    results = []
    for case_id, data, env, seed in cases():
        for arm in ARMS:
            rec = run_one(case_id, data, env, seed, arm)
            results.append(rec)
            print(f"{case_id:16s} {arm:10s} exit={rec['exit']} "
                  f"dev={rec.get('heldout_dev', float('nan')):.4f} wall={rec['wall_s']}s",
                  flush=True)
    with open(os.path.join(WORK, "results_ladder.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("wrote results_ladder.json")

if __name__ == "__main__":
    main()
