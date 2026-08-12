#!/usr/bin/env python3
"""Robustness + validation for the shrinkage-cell prototype, across seeds.

1. Validation: does the offline per-depth ML tau reproduce the C++ stored
   shrinkage per depth?  (Sanity that the offline model is faithful.)
2. Held-out-vertices EB marginal log-lik for several exchangeability units,
   averaged over seeds:
     A  per-depth                     (current)
     B  depth x Keff fixed bins
     C  B + partial pooling
     Q  depth x Keff tertiles
     R  continuous log-tau ~ a + b*depth + c*log Keff (global regression)
   plus a trimmed variant (drop top 5% |draw|) to check heavy-tail sensitivity.
"""
import sys, math, random, collections, statistics
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from eb_cells_prototype import (load, tau2_ml, group, marg_ll, keff_bin,
                                scheme_depth, scheme_cells)

DUMPS = sys.argv[1:] or [
    "results/permutation_gate/real_ginnie_unfiltered/baseline_q10__seed7/run_hier_vertices.csv",
    "/tmp/reshrink_repro/seed11_hier_vertices.csv",
    "/tmp/reshrink_repro/seed13_hier_vertices.csv",
    "/tmp/reshrink_repro/seed21_hier_vertices.csv",
]

def scheme_tertiles(fit):
    ks = sorted(r['keff'] for r in fit)
    if len(ks) < 6:
        return scheme_depth(fit)
    q1, q2 = ks[len(ks)//3], ks[2*len(ks)//3]
    def b(k): return 0 if k < q1 else (1 if k < q2 else 2)
    by_cell = group(fit, lambda r: (r['depth'], b(r['keff'])))
    by_depth = group(fit, lambda r: r['depth'])
    dt = {d: (tau2_ml(rs) or 0.0) for d, rs in by_depth.items()}
    glob = tau2_ml(fit) or 0.0
    ct = {}
    for c, rs in by_cell.items():
        n = len(rs); w = n/(n+6.0)
        tc = tau2_ml(rs); td = dt.get(c[0], glob)
        ct[c] = td if tc is None else math.exp(w*math.log(max(tc,1e-9))+(1-w)*math.log(max(td,1e-9)))
    return lambda r: ct.get((r['depth'], b(r['keff'])), dt.get(r['depth'], glob))

def scheme_regress(fit):
    """log tau_v^2 = a + b*depth + c*log(Keff), one global surface, fit by
    a crude grid/normal-equations on per-vertex 'target' log(max(draw^2-s,eps))."""
    xs, ys = [], []
    for r in fit:
        t = math.log(max(r['draw']**2 - r['s'], 1e-6))
        xs.append((1.0, float(r['depth']), math.log(max(r['keff'], 1e-4))))
        ys.append(t)
    # normal equations (3x3)
    import itertools
    A = [[0.0]*3 for _ in range(3)]; g = [0.0]*3
    for x, y in zip(xs, ys):
        for i in range(3):
            g[i] += x[i]*y
            for j in range(3):
                A[i][j] += x[i]*x[j]
    beta = solve3(A, g)
    if beta is None:
        return scheme_depth(fit)
    def tau_of(r):
        lx = (1.0, float(r['depth']), math.log(max(r['keff'], 1e-4)))
        return math.exp(sum(b*v for b, v in zip(beta, lx)))
    return tau_of

def solve3(A, b):
    import copy
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(3):
        p = max(range(i, 3), key=lambda r: abs(M[r][i]))
        if abs(M[p][i]) < 1e-12: return None
        M[i], M[p] = M[p], M[i]
        piv = M[i][i]
        M[i] = [v/piv for v in M[i]]
        for r in range(3):
            if r != i:
                f = M[r][i]
                M[r] = [a-f*c for a, c in zip(M[r], M[i])]
    return [M[i][3] for i in range(3)]

def cv(rows, build, reps=200, seed=1, trim=0.0):
    rng = random.Random(seed); tot = 0.0
    for _ in range(reps):
        idx = list(range(len(rows))); rng.shuffle(idx)
        h = len(idx)//2
        fit = [rows[i] for i in idx[:h]]; sco = [rows[i] for i in idx[h:]]
        if trim > 0:
            sco = sorted(sco, key=lambda r: abs(r['draw']))[:int(len(sco)*(1-trim))]
        tot += marg_ll(sco, build(fit))
    return tot/reps

def validate(rows):
    print("  validation: offline per-depth ML tau vs C++ stored shrinkage")
    tau = scheme_depth(rows)
    print(f"    {'depth':>5} {'n':>3} {'S_offline_med':>14} {'S_stored_med':>13}")
    for d, rs in sorted(group(rows, lambda r: r['depth']).items()):
        so = sorted(1.0/(1.0+(1.0/max(tau(r),1e-12))*r['s']) for r in rs)
        ss = sorted(r['S_stored'] for r in rs)
        print(f"    {d:5d} {len(rs):3d} {so[len(so)//2]:14.4g} {ss[len(ss)//2]:13.4g}")

SCHEMES = [
    ("A per-depth",         scheme_depth),
    ("B fixed bins",        lambda f: scheme_cells(f, pool=False)),
    ("C bins+pool",         lambda f: scheme_cells(f, pool=True)),
    ("Q tertiles+pool",     scheme_tertiles),
    ("R log-tau regress",   scheme_regress),
]

def main():
    print("VALIDATION (first dump):")
    validate(load(DUMPS[0]))
    print("\nHELD-OUT EB marginal log-lik, delta vs A (per-depth); mean over 200 splits\n")
    print(f"  {'dump/seed':>10} {'n':>4} | " + " ".join(f"{s[0].split()[0]:>7}" for s in SCHEMES)
          + " | " + " ".join(f"{s[0].split()[0]+'*':>7}" for s in SCHEMES[1:]) + "   (* = 5%-trimmed)")
    agg = collections.defaultdict(float); aggt = collections.defaultdict(float); nd = 0
    for path in DUMPS:
        rows = load(path)
        if len(rows) < 20: continue
        nd += 1
        base = cv(rows, scheme_depth)
        cells = [cv(rows, b) - base for _, b in SCHEMES]
        based = base  # trimmed baseline
        baset = cv(rows, scheme_depth, trim=0.05)
        cellst = [cv(rows, b, trim=0.05) - baset for _, b in SCHEMES[1:]]
        tag = path.split('/')[-2] if 'seed' not in path.split('/')[-1] else path.split('/')[-1].split('_')[0]
        print(f"  {tag[:10]:>10} {len(rows):4d} | " +
              " ".join(f"{c:+7.2f}" for c in cells) + " | " +
              " ".join(f"{c:+7.2f}" for c in cellst))
        for (name, _), c in zip(SCHEMES, cells): agg[name] += c
        for (name, _), c in zip(SCHEMES[1:], cellst): aggt[name] += c
    print("\n  MEAN delta vs A across dumps:")
    for name, _ in SCHEMES:
        extra = f"   trimmed {aggt[name]/nd:+.2f}" if name in aggt else ""
        print(f"    {name:20s} {agg[name]/nd:+8.2f}{extra}")
    print("\n  (delta > 0 means the richer unit beats plain per-depth out of sample)")

if __name__ == "__main__":
    main()
