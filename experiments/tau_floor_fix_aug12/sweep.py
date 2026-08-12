#!/usr/bin/env python3
"""A/B sweep for the final-reshrink tau-floor fix.

Arms per (dataset, seed):
  legacy : DMESH_FINAL_RESHRINK=1 + DMESH_RESHRINK_NO_FLOOR=1 (historical behavior)
  fixed  : DMESH_FINAL_RESHRINK=1 (floor retained; the fix)
  off    : no reshrink (anchor)

Captures: held-out deviance, X2/info, admitted-coefficient count, per-depth
lambda_v max, [reshrink] tau lines, run time. Dumps are deleted after summary.
"""
import csv, json, os, shutil, subprocess, sys, time
from collections import defaultdict

ROOT = "/home/user/decision-mesh-release"
BIN = os.path.join(ROOT, "build-release", "decision_mesh")
WORK = os.path.dirname(os.path.abspath(__file__))

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

CASES = []  # (case_id, data_path, env, q, seed)

def add_cases():
    smm = os.path.join(ROOT, "data", "smm_design_202606.csv")
    april = os.path.join(ROOT, "data", "ginnie_design.csv")
    flatnull = os.path.join(ROOT, "tests", "data", "flat_p78_clean_seed0.csv")
    for seed in [7] + list(range(101, 111)):
        CASES.append((f"smm_s{seed}", smm, SMM_ENV, "0.10", seed))
    for seed in [7] + list(range(101, 106)):
        CASES.append((f"april_s{seed}", april, APRIL_ENV, "0.10", seed))
    for seed in [7] + list(range(101, 105)):
        CASES.append((f"aprilnull_s{seed}", flatnull, APRIL_ENV, "0.10", seed))
    for nseed in range(3):
        p = os.path.join(WORK, f"smmnull_{nseed}.csv")
        if os.path.exists(p):
            for seed in (7, 101):
                CASES.append((f"smmnull{nseed}_s{seed}", p, SMM_ENV, "0.10", seed))

def parse_stdout(text):
    out = {}
    for line in text.splitlines():
        if line.startswith("HELDOUT:"):
            parts = line.split("|")
            out["heldout_dev"] = float(parts[0].split("deviance/pool")[1])
            out["x2_info"] = float(parts[1].split("X2/info")[1])
        if "vertices" in line and "admitted" in line:
            out.setdefault("vertices_line", line.strip())
    return out

def lambda_ladder(dump_prefix):
    path = dump_prefix + "_hier_vertices.csv"
    ladder, admitted = defaultdict(lambda: [0, 0.0, 0.0]), 0
    if not os.path.exists(path):
        return {}, -1
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                if r["active"] != "1" or r["gate_admitted"] != "1":
                    continue
                if r["free_coefficient"] != "1" or not r["parent0"]:
                    continue
                d, lam = int(r["depth"]), float(r["lambda_v"])
            except (KeyError, ValueError):
                continue
            admitted += 1
            e = ladder[d]
            e[0] += 1
            e[1] = max(e[1], lam)
            e[2] += lam
    return {d: {"n": e[0], "lam_max": e[1], "lam_mean": e[2] / e[0]}
            for d, e in ladder.items()}, admitted

def run_one(case_id, data, env_extra, q, seed, arm):
    tag = f"{case_id}_{arm}"
    outdir = os.path.join(WORK, "runs", tag)
    shutil.rmtree(outdir, ignore_errors=True)
    os.makedirs(outdir)
    env = dict(os.environ)
    env.update(env_extra)
    env["DMESH_DATA"] = data
    env["DMESH_DUMP"] = os.path.join(outdir, "run")
    if arm in ("legacy", "fixed"):
        env["DMESH_FINAL_RESHRINK"] = "1"
    if arm == "legacy":
        env["DMESH_RESHRINK_NO_FLOOR"] = "1"
    t0 = time.time()
    p = subprocess.run([BIN, "0", q, "1", "24", str(seed)], env=env,
                       capture_output=True, text=True, cwd=outdir)
    wall = time.time() - t0
    rec = {"case": case_id, "arm": arm, "seed": seed, "exit": p.returncode,
           "wall_s": round(wall, 1)}
    rec.update(parse_stdout(p.stdout))
    rec["reshrink_lines"] = [l for l in p.stderr.splitlines() if "[reshrink]" in l and "->" in l]
    ladder, admitted = lambda_ladder(os.path.join(outdir, "run"))
    rec["admitted"] = admitted
    rec["ladder"] = {str(k): v for k, v in sorted(ladder.items())}
    shutil.rmtree(outdir, ignore_errors=True)
    return rec

def main():
    add_cases()
    only = sys.argv[1] if len(sys.argv) > 1 else ""
    results = []
    for case_id, data, env, q, seed in CASES:
        if only and not case_id.startswith(only):
            continue
        for arm in ("legacy", "fixed", "off"):
            rec = run_one(case_id, data, env, q, seed, arm)
            results.append(rec)
            dev = rec.get("heldout_dev", float("nan"))
            print(f"{case_id:16s} {arm:6s} exit={rec['exit']} dev={dev:.4f} "
                  f"adm={rec['admitted']} wall={rec['wall_s']}s", flush=True)
    outname = os.path.join(WORK, f"results_{only or 'all'}.json")
    with open(outname, "w") as f:
        json.dump(results, f, indent=1)
    print("wrote", outname)

if __name__ == "__main__":
    main()
