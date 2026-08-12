#!/usr/bin/env python3
"""Offline prototype: hierarchical shrinkage beyond pure depth.

Reads a shipped run_hier_vertices.csv dump and compares the shrinkage prior
under three exchangeability units:

  (A) per-depth tau               -- the current production rule
  (B) per-(depth x Keff-bin) tau  -- evidence-geometry cells
  (C) (B) with partial pooling    -- cell tau shrunk toward its depth mean

Everything operates on the *raw* admission surpluses recovered from the dump
(undoing the pooled shrinkage the C++ stored), so the prior is fit on the
draws its own model assumes -- delta_raw ~ N(0, tau^2 + s_v) -- rather than on
already-shrunk coefficients.

Honest scoring: admitted vertices are randomly split into fit/score halves
many times; each scheme's tau is estimated on the fit half and its EB marginal
log-likelihood  sum log N(delta_raw; 0, tau_cell^2 + s_v)  is evaluated on the
held-out score half.  A better exchangeability unit predicts the held-out
coefficients' dispersion better.  This controls for the extra parameters the
cell schemes spend.

Units: sigma_sq, lambda_v, tau^2 are all per-centilogit^2 (the dump's native
"height x100" convention); lambda*s and shrinkage factors are dimensionless.
"""
import csv, sys, math, random, collections

def load(path):
    out = []
    for r in csv.DictReader(open(path)):
        if r['active'] != '1' or r['free_coefficient'] != '1':
            continue
        if r['gate_admitted'] != '1' or r['parent0'] == '-1':
            continue
        try:
            s = float(r['sigma_sq']); lam = float(r['lambda_v'])
            dp = float(r['delta_pooled']); keff = float(r['current_keff'])
            d = int(r['depth'])
        except (ValueError, KeyError):
            continue
        if not (s > 0 and s < 1e299):        # need a usable data variance
            continue
        # Recover the raw surplus: delta_pooled = S * delta_raw with
        # S = 1/(1 + lambda*s) the MAP shrinkage factor toward the parent.
        S = 1.0 / (1.0 + lam * s)
        draw = dp / S if S > 0 else dp
        out.append(dict(id=int(r['id']), depth=d, s=s, keff=keff,
                        draw=draw, lam_stored=lam, S_stored=S))
    return out

# ---- tau estimators -------------------------------------------------------

def tau2_ml(rows, floor=1e-6):
    """Marginal MLE of tau^2 for draws ~ N(0, tau^2 + s): 1-D bisection on the
    score equation sum (draw^2 - (tau2+s)) / (tau2+s)^2 = 0."""
    if len(rows) < 2:
        return None
    def score(t2):
        return sum((r['draw']**2 - (t2 + r['s'])) / (t2 + r['s'])**2 for r in rows)
    lo, hi = 0.0, 1.0
    if score(lo) <= 0:                       # no slab evidence beyond noise
        return 0.0
    while score(hi) > 0 and hi < 1e12:
        hi *= 2
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if score(mid) > 0: lo = mid
        else: hi = mid
    return max(0.5 * (lo + hi), floor)

def keff_bin(keff):
    if keff < 0.1:  return 0        # near-void corridor
    if keff < 1.0:  return 1        # thin
    if keff < 10.0: return 2        # workhorse
    return 3                        # broad

# ---- schemes: return a dict vertex-id -> tau^2 ----------------------------

def scheme_depth(fit):
    tau = {}
    for d, rs in group(fit, lambda r: r['depth']).items():
        tau[d] = tau2_ml(rs)
    glob = tau2_ml(fit) or 0.0
    return lambda r: _borrow_depth(tau, r['depth'], glob)

def scheme_cells(fit, pool=False):
    by_cell = group(fit, lambda r: (r['depth'], keff_bin(r['keff'])))
    by_depth = group(fit, lambda r: r['depth'])
    depth_tau = {d: (tau2_ml(rs) or 0.0) for d, rs in by_depth.items()}
    glob = tau2_ml(fit) or 0.0
    cell_tau = {}
    for c, rs in by_cell.items():
        d = c[0]
        t_cell = tau2_ml(rs)
        t_depth = depth_tau.get(d, glob)
        if t_cell is None:
            cell_tau[c] = t_depth
        elif pool:
            # Partial pooling in log-variance toward the depth mean; strength
            # w = n/(n+kappa) so small cells lean on the depth, large cells
            # speak for themselves.  kappa is one prototype knob.
            n = len(rs); kappa = 6.0
            w = n / (n + kappa)
            lt = w * math.log(max(t_cell, 1e-9)) + (1 - w) * math.log(max(t_depth, 1e-9))
            cell_tau[c] = math.exp(lt)
        else:
            cell_tau[c] = t_cell
    return lambda r: cell_tau.get((r['depth'], keff_bin(r['keff'])),
                                  depth_tau.get(r['depth'], glob))

def _borrow_depth(tau, d, glob):
    if tau.get(d) is not None:
        return tau[d]
    known = [k for k, v in tau.items() if v is not None]
    if not known:
        return glob
    nd = min(known, key=lambda k: abs(k - d))
    return tau[nd]

def group(rows, key):
    g = collections.defaultdict(list)
    for r in rows:
        g[key(r)].append(r)
    return g

# ---- scoring --------------------------------------------------------------

def marg_ll(rows, tau_of):
    """EB marginal log-likelihood of held-out draws under a fitted tau map."""
    ll = 0.0
    for r in rows:
        v = (tau_of(r) or 0.0) + r['s']
        ll += -0.5 * (math.log(2 * math.pi * v) + r['draw']**2 / v)
    return ll

def cv(rows, build, reps=200, seed=1):
    rng = random.Random(seed)
    tot = 0.0
    for _ in range(reps):
        idx = list(range(len(rows)))
        rng.shuffle(idx)
        half = len(idx) // 2
        fit = [rows[i] for i in idx[:half]]
        sco = [rows[i] for i in idx[half:]]
        tot += marg_ll(sco, build(fit))
    return tot / reps

# ---- main -----------------------------------------------------------------

def main(path):
    rows = load(path)
    print(f"loaded {len(rows)} admitted free coefficients from\n  {path}\n")

    print("Exchangeability problem -- shrinkage vs evidence WITHIN each depth")
    print("(current rule gives every vertex at a depth the same tau):")
    print(f"  {'depth':>5} {'n':>3}  {'Keff min/med/max':>22}  {'stored shrink min/med/max':>26}")
    for d, rs in sorted(group(rows, lambda r: r['depth']).items()):
        k = sorted(r['keff'] for r in rs)
        sh = sorted(r['S_stored'] for r in rs)
        md = lambda a: a[len(a)//2]
        print(f"  {d:5d} {len(rs):3d}  {k[0]:6.3g}/{md(k):6.3g}/{k[-1]:6.3g}  "
              f"{sh[0]:7.3g}/{md(sh):7.3g}/{sh[-1]:7.3g}")

    print("\nHeld-out-vertices EB marginal log-likelihood (higher = better prior;")
    print("mean over 200 random 50/50 splits of the admitted coefficients):")
    schemes = [
        ("A  per-depth (current unit)",       scheme_depth),
        ("B  depth x Keff cells",             lambda f: scheme_cells(f, pool=False)),
        ("C  cells + partial pooling",        lambda f: scheme_cells(f, pool=True)),
    ]
    base = None
    for name, build in schemes:
        val = cv(rows, build)
        if base is None: base = val
        print(f"  {name:32s} {val:12.3f}   (delta vs A: {val-base:+8.3f})")

    # Decoupling: fit each scheme on ALL vertices and show shrinkage vs Keff
    # for the DEEP well-supported vertices (depth>=8), the ones the depth rule
    # over-shrinks.
    print("\nDeep well-supported vertices (depth>=8, Keff>=3): shrinkage toward parent")
    print("  under A (per-depth) vs C (cells+pool). S->1 means fully pinned.")
    tauA = scheme_depth(rows)
    tauC = scheme_cells(rows, pool=True)
    deep = sorted([r for r in rows if r['depth'] >= 8 and r['keff'] >= 3.0],
                  key=lambda r: -r['keff'])
    print(f"  {'id':>6} {'depth':>5} {'Keff':>7} {'s':>8}  {'S_A':>7} {'S_C':>7}  {'freed?':>7}")
    for r in deep[:12]:
        SA = 1.0 / (1.0 + (1.0/max(tauA(r),1e-12)) * r['s'])
        SC = 1.0 / (1.0 + (1.0/max(tauC(r),1e-12)) * r['s'])
        freed = "yes" if SC < SA - 0.02 else ""
        print(f"  {r['id']:6d} {r['depth']:5d} {r['keff']:7.3g} {r['s']:8.3g}  "
              f"{SA:7.3g} {SC:7.3g}  {freed:>7}")

if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else \
        "results/permutation_gate/real_ginnie_unfiltered/baseline_q10__seed7/run_hier_vertices.csv"
    main(p)
