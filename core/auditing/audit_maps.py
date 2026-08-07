#!/usr/bin/env python3
"""Whole-surface audit maps.

Usage: python3 audit_maps.py <pools.csv> <dump_prefix> <law_json> <out.png>

(a) dominant ladder depth: which hierarchy level carries the largest
    absolute contribution to the estimate at each location;
(b) kernel n_eff of the estimate across the surface (log10);
(c) max single-pool leave-one-out influence per cell, centi-logits
    (the mega-risk map).
Issuer-share map: pending an issuer column in the pool csv (re-parse).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from audit_lib import MeshAudit

pools_csv, prefix, law_json, out = sys.argv[1:5]
A = MeshAudit(pools_csv, prefix, law_json)

GX, GY = 72, 46
wq = np.linspace(A.wala.min() + 1, A.wala.max() - 1, GX)
cq = np.linspace(A.wac.min() + 0.05, A.wac.max() - 0.05, GY)
dom = np.full((GX, GY), np.nan)
neff = np.full((GX, GY), np.nan)
for ix, wv in enumerate(wq):
    for iy, cv in enumerate(cq):
        try:
            rows, _ = A.ladder(wv, cv)
        except ValueError:
            continue
        by = {}
        for hr, dep, psi, s, contrib in rows:
            by[dep] = by.get(dep, 0.0) + contrib
        if by:
            dom[ix, iy] = max(by, key=lambda d: abs(by[d]))
        kern = A.kernel(wv, cv)
        neff[ix, iy] = A.kernel_stats(kern)["n_eff"]

lev = A.leverage()
infl = A.self_influence_map(lev)
BX, BY = 60, 42
num, xe, ye = np.histogram2d(A.wala, A.wac, bins=[BX, BY], weights=infl)
cnt, _, _ = np.histogram2d(A.wala, A.wac, bins=[BX, BY])
mx = np.zeros_like(num)
# max per cell
ib = np.clip(np.digitize(A.wala, xe) - 1, 0, BX - 1)
jb = np.clip(np.digitize(A.wac, ye) - 1, 0, BY - 1)
for i in np.where(infl > 0)[0]:
    mx[ib[i], jb[i]] = max(mx[ib[i], jb[i]], infl[i])
mx[cnt == 0] = np.nan

fig, ax = plt.subplots(1, 3, figsize=(18, 4.6), dpi=170)
ext = [A.wala.min(), A.wala.max(), A.wac.min(), A.wac.max()]
im0 = ax[0].imshow(dom.T, origin="lower", extent=ext, aspect="auto",
                   cmap="plasma")
fig.colorbar(im0, ax=ax[0], shrink=0.85).set_label("dominant depth")
ax[0].set_title("(a) which scale carries the estimate")
neff = np.where(neff > 1, neff, np.nan)
im1 = ax[1].imshow(np.log10(neff).T, origin="lower", extent=ext,
                   aspect="auto", cmap="viridis")
fig.colorbar(im1, ax=ax[1], shrink=0.85).set_label("log10 n_eff")
ax[1].set_title("(b) how many pools power it")
im2 = ax[2].imshow(np.minimum(mx, 30).T, origin="lower", extent=ext,
                   aspect="auto", cmap="inferno")
fig.colorbar(im2, ax=ax[2], shrink=0.85).set_label("centi-logits (capped 30)")
ax[2].set_title("(c) max single-pool influence (mega risk)")
for a in ax:
    a.set_xlabel("WALA"); a.set_ylabel("WAC")
fig.suptitle("Surface-wide attribution audit (conditional on selection; see README.md)")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(out, bbox_inches="tight")
print(f"sole-support pools (leverage > 0.99): {A.n_sole_support}")
print(f"dominant-depth range {np.nanmin(dom):.0f}-{np.nanmax(dom):.0f}; "
      f"n_eff range {np.nanmin(neff):.0f}-{np.nanmax(neff):.0f}; "
      f"max single-pool influence {np.nanmax(mx):.1f} centi")
print(f"wrote {out}")
