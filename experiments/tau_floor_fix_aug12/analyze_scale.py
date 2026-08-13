#!/usr/bin/env python3
import json, math
from collections import defaultdict
recs = json.load(open("results_scale2.json"))
ARMS = ["off","eb","eb_pv","ebscale","ebscale_free","ebscale_free_pv"]
by = defaultdict(dict)
for r in recs: by[r["seed"]][r["arm"]] = r.get("heldout_dev", float("nan"))
def ms(xs):
    xs=[x for x in xs if x==x]; m=sum(xs)/len(xs)
    return m, math.sqrt(sum((x-m)**2 for x in xs)/max(len(xs)-1,1))
print(f"{'seed':>5} | " + " ".join(f"{a:>11}" for a in ARMS) + " | best")
cols=defaultdict(list)
for s in sorted(by):
    v={a:by[s].get(a,float('nan')) for a in ARMS}
    best=min((x,a) for a,x in v.items() if x==x)[1]
    print(f"{s:>5} | " + " ".join(f"{v[a]:>11.4f}" for a in ARMS) + f" | {best}")
    for a in ARMS: cols[a].append(v[a])
print("-"*90)
line=f"{'MEAN':>5} | "
for a in ARMS:
    m,sd=ms(cols[a]); line+=f"{m:.4f}     "
print(line)
line=f"{'sd':>5} | "
for a in ARMS:
    m,sd=ms(cols[a]); line+=f"{sd:.4f}     "
print(line)
# win counts vs off and vs eb
for ref in ("off","eb"):
    print(f"\nwins vs {ref}:")
    for a in ARMS:
        if a==ref: continue
        w=sum(by[s][a]<by[s][ref]-1e-9 for s in by if a in by[s] and ref in by[s])
        t=sum(abs(by[s][a]-by[s][ref])<=1e-9 for s in by if a in by[s] and ref in by[s])
        n=len(by)
        d=sum(by[s][ref]-by[s][a] for s in by)/n
        print(f"  {a:12s} {w}/{n} win, {t} tie, mean delta {d:+.4f}")
