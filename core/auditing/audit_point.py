#!/usr/bin/env python3
"""Audit one query point: where did this number come from?

Usage: python3 audit_point.py <pools.csv> <dump_prefix> <law_json> \
                              <wala> <wac> <out.png>

Left: the hierarchy waterfall (base rate -> depth-grouped shrunken
surpluses -> final estimate; exact, verified against the surface).
Right: the effective kernel (per-pool influence weights, binned), with
n_eff, stratum shares, top contributors, and variance attribution.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from audit_lib import MeshAudit, STRATA_EDGES

pools_csv, prefix, law_json = sys.argv[1:4]
wq, cq = float(sys.argv[4]), float(sys.argv[5])
out = sys.argv[6]

A = MeshAudit(pools_csv, prefix, law_json)
rows, eta0 = A.ladder(wq, cq)
tot, direct = A.ladder_check(wq, cq)
assert abs(tot - direct) < 1e-9, "ladder does not reconstruct the surface"
kern = A.kernel(wq, cq)
ks = A.kernel_stats(kern)
va = A.variance_attribution(wq, cq)
pred, direct = A.kernel_check(wq, cq)
recon_gap = abs(pred - direct) * 100

# depth-grouped waterfall
by_depth = {}
for hr, dep, psi, s, contrib in rows:
    by_depth[dep] = by_depth.get(dep, 0.0) + contrib
deps = sorted(by_depth)
fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2), dpi=170,
                       gridspec_kw={"width_ratios": [1.0, 1.15]})
cum = eta0
xs, labels = [], []
ax[0].bar(0, eta0, color="#666666")
labels.append("base"); xs.append(0)
for j, dep in enumerate(deps, start=1):
    c = by_depth[dep]
    ax[0].bar(j, c, bottom=cum, color="#27496d" if c >= 0 else "#8c3b5e")
    cum += c
    xs.append(j); labels.append(f"d{dep}")
ax[0].bar(len(deps) + 1, cum, color="#2e8b6a")
xs.append(len(deps) + 1); labels.append("estimate")
ax[0].set_xticks(xs); ax[0].set_xticklabels(labels, fontsize=8)
ax[0].axhline(0, color="k", lw=0.7)
ax[0].set_ylabel("log-odds")
ax[0].set_title(f"waterfall at WALA {wq:.0f}, WAC {cq:.2f}: "
                f"estimate {cum:+.3f} logits (verified)")

# kernel heatmap
BX, BY = 60, 42
H, xe, ye = np.histogram2d(A.wala, A.wac, bins=[BX, BY], weights=np.abs(kern))
im = ax[1].imshow(np.log10(np.maximum(H.T, 1e-9)), origin="lower",
                  extent=[xe[0], xe[-1], ye[0], ye[-1]], aspect="auto",
                  cmap="magma")
ax[1].plot([wq], [cq], "o", ms=8, mfc="none", mec="#2e8b6a", mew=2)
fig.colorbar(im, ax=ax[1], shrink=0.85).set_label("log10 |influence|")
sh = ks["stratum_abs_shares"]
ax[1].set_title(f"effective kernel: n_eff {ks['n_eff']:.0f} pools | "
                f"|w| shares <128/512/4096/mega: "
                f"{sh[0]:.2f}/{sh[1]:.2f}/{sh[2]:.2f}/{sh[3]:.2f}")
ax[1].set_xlabel("WALA"); ax[1].set_ylabel("WAC")
txt = (f"predictive sd {va['sd_logits']:.3f} logits; variance shares by "
       f"stratum {'/'.join(f'{s:.2f}' for s in va['stratum_shares'])}; "
       f"kernel reconstruction gap {recon_gap:.1f} centi (see README caveat 3)")
fig.text(0.5, 0.005, txt, ha="center", fontsize=9, color="#444444")
fig.suptitle("Estimate attribution (exact, conditional on selection; see auditing/README.md)")
fig.tight_layout(rect=[0, 0.03, 1, 0.97])
fig.savefig(out, bbox_inches="tight")

top = np.argsort(-np.abs(kern))[:10]
print(f"estimate {cum:+.4f} logits = base {eta0:+.3f} "
      + " ".join(f"d{d}:{by_depth[d]:+.3f}" for d in deps))
print(f"n_eff {ks['n_eff']:.0f}; top pools by |influence|:")
for i in top:
    print(f"  wala {A.wala[i]:6.1f} wac {A.wac[i]:.3f} n {A.n[i]:8.0f} "
          f"kernel {kern[i]:+.5f}")
print(f"kernel reconstruction gap {recon_gap:.1f} centi-logits")
print(f"wrote {out}")
