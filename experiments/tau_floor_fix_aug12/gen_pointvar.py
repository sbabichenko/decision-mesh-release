#!/usr/bin/env python3
"""Honest heterogeneous per-pool variance v_i, smooth in pool size.

The shipped one-law prices pool-effect variance in 4 coarse size strata
(0.1428/0.2415/0.0286/0.0038 for n0 <128/<512/<4096/>=4096), a step
function with an 8.5x cliff at n0=512. That discretization is an artifact:
a pool at n0=511 and one at n0=513 get wildly different v_i. Replace the
step with a smooth log-log interpolation of the stratum variances against
the stratum geometric-mean sizes, so v_i varies continuously with n0 and
small pools keep their large (prior) variance. Pure covariate (size) model:
no conditioning on the pool's own residual, so no surface-error absorption.

Writes one v_i per data row, in file order, for DMESH_POINT_VAR.
"""
import csv, math, sys

CSV = "/home/user/decision-mesh-release/data/smm_design_202606.csv"
OUT = sys.argv[1] if len(sys.argv) > 1 else "pointvar_smooth_smm.txt"

# (geometric-mean size of stratum, stratum s2)
ANCHORS = [
    (math.sqrt(25 * 128),  0.1428),   # <128
    (math.sqrt(128 * 512), 0.2415),   # 128-512
    (math.sqrt(512 * 4096),0.0286),   # 512-4096
    (8192.0,               0.0038),   # >=4096
]
lx = [math.log(a[0]) for a in ANCHORS]
ly = [math.log(a[1]) for a in ANCHORS]

def s2_of(n0):
    x = math.log(max(n0, 1.0))
    if x <= lx[0]:  return math.exp(ly[0])
    if x >= lx[-1]: return math.exp(ly[-1])
    for i in range(len(lx) - 1):
        if lx[i] <= x <= lx[i + 1]:
            t = (x - lx[i]) / (lx[i + 1] - lx[i])
            return math.exp(ly[i] + t * (ly[i + 1] - ly[i]))
    return math.exp(ly[-1])

rows = list(csv.DictReader(open(CSV)))
with open(OUT, "w") as f:
    for r in rows:
        f.write(f"{s2_of(float(r['n0'])):.6f}\n")
print(f"wrote {OUT}: {len(rows)} rows")
# show the smooth curve at a few sizes vs the step
for n in (30, 100, 127, 129, 300, 511, 513, 1000, 4096, 20000):
    print(f"  n0={n:6d}  smooth s2={s2_of(n):.4f}")
