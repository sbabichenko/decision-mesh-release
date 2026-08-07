#!/usr/bin/env python3
"""Price every admitted feature: held-out value of each certified surplus.

Usage: python3 audit_ablation.py <pools.csv> <dump_prefix> <law_json> <out.csv>

For each gate-admitted vertex, the surface is ablated by zeroing its
shrunken surplus (exact, linear, no refit) and the change in held-out
law-weighted quadratic score is recorded. Positive value = the feature
helps held-out fit; near-zero = admitted but predictively idle;
negative = the feature hurts held-out fit (flag for review).
"""
import sys
import numpy as np
from audit_lib import MeshAudit

pools_csv, prefix, law_json, out = sys.argv[1:5]
A = MeshAudit(pools_csv, prefix, law_json)
tab = A.ablation_table()
tab.to_csv(out, index=False)
if len(tab):
    v = tab.heldout_value.values
    print(f"admitted features priced: {len(tab)}")
    print(f"value quantiles (held-out score units): "
          f"p10 {np.percentile(v,10):+.1f}  median {np.median(v):+.1f}  "
          f"p90 {np.percentile(v,90):+.1f}")
    print(f"negative-value features: {(v < 0).sum()} "
          f"({(v < 0).mean()*100:.0f}%) - review list at top of {out}")
    print(tab.head(5).to_string(index=False))
else:
    print("no admitted interior features in this dump")
