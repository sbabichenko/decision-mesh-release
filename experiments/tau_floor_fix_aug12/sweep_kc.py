#!/usr/bin/env python3
"""keffcont-only sweep; appends rows to results_ladder.json for merge."""
import json, os, sys
sys.argv = [sys.argv[0]]
import importlib.util
spec = importlib.util.spec_from_file_location("sweep2", os.path.join(os.path.dirname(os.path.abspath(__file__)), "sweep2.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
m.BIN = os.path.join(m.WORK, "dm_snapshot3")
m.ARMS = {"keffcont": {"DMESH_FINAL_RESHRINK": "1", "DMESH_RESHRINK_LADDER": "keffcont"}}
results = []
for case_id, data, env, seed in m.cases():
    rec = m.run_one(case_id, data, env, seed, "keffcont")
    results.append(rec)
    print(f"{case_id:16s} keffcont exit={rec['exit']} dev={rec.get('heldout_dev', float('nan')):.4f} wall={rec['wall_s']}s", flush=True)
existing = json.load(open(os.path.join(m.WORK, "results_ladder.json")))
json.dump(existing + results, open(os.path.join(m.WORK, "results_ladder.json"), "w"), indent=1)
print("merged keffcont into results_ladder.json")
