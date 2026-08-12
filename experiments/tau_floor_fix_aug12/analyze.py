#!/usr/bin/env python3
"""Summarize A/B sweep results: fixed vs legacy vs off."""
import json, math, sys
from collections import defaultdict

recs = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "results_all.json"))
by_case = defaultdict(dict)
for r in recs:
    by_case[r["case"]][r["arm"]] = r

def group_of(case):
    return case.split("_s")[0].rstrip("0123456789") if case.startswith("smmnull") else case.split("_s")[0]

groups = defaultdict(list)
for case in sorted(by_case):
    groups[group_of(case)].append(case)

def mean_sd(xs):
    xs = [x for x in xs if x == x]
    if not xs:
        return float("nan"), float("nan")
    m = sum(xs) / len(xs)
    v = sum((x - m) ** 2 for x in xs) / max(len(xs) - 1, 1)
    return m, math.sqrt(v)

print("=" * 100)
print("Per-run table: held-out deviance (fixed vs legacy vs off), admitted count, max lambda_v")
print("=" * 100)
for g, cases in groups.items():
    print(f"\n--- {g} ---")
    print(f"{'case':>15} | {'dev legacy':>10} {'dev fixed':>10} {'dev off':>10} | {'d(fix-leg)':>10} | "
          f"{'adm L/F':>9} | {'maxlam L':>9} {'maxlam F':>9} | identical?")
    deltas, devL, devF, devO = [], [], [], []
    for case in cases:
        arms = by_case[case]
        dl = arms.get("legacy", {}).get("heldout_dev", float("nan"))
        df = arms.get("fixed", {}).get("heldout_dev", float("nan"))
        do = arms.get("off", {}).get("heldout_dev", float("nan"))
        al = arms.get("legacy", {}).get("admitted", -1)
        af = arms.get("fixed", {}).get("admitted", -1)
        ml = max((v["lam_max"] for v in arms.get("legacy", {}).get("ladder", {}).values()), default=0)
        mf = max((v["lam_max"] for v in arms.get("fixed", {}).get("ladder", {}).values()), default=0)
        ident = "YES" if arms.get("legacy", {}).get("ladder") == arms.get("fixed", {}).get("ladder") \
                and abs(dl - df) < 1e-12 else "no"
        print(f"{case:>15} | {dl:>10.4f} {df:>10.4f} {do:>10.4f} | {df-dl:>+10.4f} | "
              f"{al:>4}/{af:<4} | {ml:>9.3g} {mf:>9.3g} | {ident}")
        deltas.append(df - dl); devL.append(dl); devF.append(df); devO.append(do)
    mL, sL = mean_sd(devL); mF, sF = mean_sd(devF); mO, sO = mean_sd(devO)
    mD, sD = mean_sd(deltas)
    print(f"{'MEAN':>15} | {mL:>10.4f} {mF:>10.4f} {mO:>10.4f} | {mD:>+10.4f} | "
          f"(sd legacy {sL:.4f}, fixed {sF:.4f}, delta {sD:.4f})")

print()
print("=" * 100)
print("Runs where the floor bound (legacy != fixed): reshrink taus and lambda ladders")
print("=" * 100)
for case in sorted(by_case):
    arms = by_case[case]
    l, f = arms.get("legacy", {}), arms.get("fixed", {})
    if l.get("ladder") == f.get("ladder"):
        continue
    print(f"\n### {case}")
    for line in l.get("reshrink_lines", []):
        print("  legacy:", line.strip())
    for arm_name, rec in (("legacy", l), ("fixed", f)):
        lad = rec.get("ladder", {})
        tops = {d: (v["n"], f'{v["lam_max"]:.3g}') for d, v in lad.items()}
        print(f"  {arm_name:6s} dev={rec.get('heldout_dev'):.4f} adm={rec.get('admitted')} "
              f"ladder(n,max_lam)={tops}")
