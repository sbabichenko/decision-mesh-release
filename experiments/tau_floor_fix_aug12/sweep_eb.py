#!/usr/bin/env python3
"""EB modes sweep; merges eb/ebcell/eb_free into results_ladder.json."""
import json, os, sys
sys.argv = [sys.argv[0]]
import importlib.util
spec = importlib.util.spec_from_file_location("sweep2", os.path.join(os.path.dirname(os.path.abspath(__file__)), "sweep2.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
m.BIN = os.path.join(m.WORK, "dm_snapshot4")
m.ARMS = {
    "eb":      {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "eb"},
    "ebcell":  {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "ebcell"},
    "eb_free": {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "eb", "DMESH_EB_FREE": "1"},
}
results = []
for case_id, data, env, seed in m.cases():
    for arm in m.ARMS:
        rec = m.run_one(case_id, data, env, seed, arm)
        results.append(rec)
        print(f"{case_id:16s} {arm:8s} exit={rec['exit']} dev={rec.get('heldout_dev', float('nan')):.4f} wall={rec['wall_s']}s", flush=True)
existing = json.load(open(os.path.join(m.WORK, "results_ladder.json")))
existing = [r for r in existing if r.get("arm") not in m.ARMS]  # idempotent
json.dump(existing + results, open(os.path.join(m.WORK, "results_ladder.json"), "w"), indent=1)
print("merged EB arms into results_ladder.json")
