#!/usr/bin/env python3
"""Summarize ladder-variant sweep."""
import json, math, sys
from collections import defaultdict

recs = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "results_ladder.json"))
ARMS = ["off", "depth", "cells_rawkeff", "cells", "keffcont", "twogroups", "twogroups_free"]
by_case = defaultdict(dict)
for r in recs:
    by_case[r["case"]][r["arm"]] = r

groups = defaultdict(list)
for case in sorted(by_case):
    groups[case.split("_s")[0]].append(case)

def mean_sd(xs):
    xs = [x for x in xs if x == x]
    if not xs:
        return float("nan"), float("nan")
    m = sum(xs) / len(xs)
    v = sum((x - m) ** 2 for x in xs) / max(len(xs) - 1, 1)
    return m, math.sqrt(v)

for g, cases in groups.items():
    print(f"\n=== {g} (held-out deviance/pool) ===")
    print(f"{'case':>15} | " + " ".join(f"{a:>10}" for a in ARMS) +
          " | best arm")
    cols = defaultdict(list)
    for case in cases:
        vals = {a: by_case[case].get(a, {}).get("heldout_dev", float("nan")) for a in ARMS}
        best = min((v, a) for a, v in vals.items() if v == v)[1]
        print(f"{case:>15} | " + " ".join(f"{vals[a]:>10.4f}" for a in ARMS) +
              f" | {best}")
        for a in ARMS:
            cols[a].append(vals[a])
    line = f"{'MEAN (sd)':>15} | "
    for a in ARMS:
        m, s = mean_sd(cols[a])
        line += f"{m:.4f}({s:.3f}) "
    print(line)

print("\n=== two-groups EM diagnostics (SMM) ===")
for case in sorted(by_case):
    if not case.startswith("smm"):
        continue
    tg = by_case[case].get("twogroups", {})
    for l in tg.get("diag", []):
        if "ladder2g" in l:
            print(f"{case:>12}: {l}")
