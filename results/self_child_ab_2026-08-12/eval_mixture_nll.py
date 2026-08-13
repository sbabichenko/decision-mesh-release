#!/usr/bin/env python3
"""Held-out mixture NLL against one fixed reference law (the marginal
currency of sec:aug6-release / sec:aug12).

For each held-out pool: -log integral Binom(k | n, expit(eta + u)) dF(u),
with eta = logit(p_hat) from the run's obs dump and F the FIXED reference
law: per size stratum (n<128 / <512 / <4096 / rest) a two-component
Gaussian scale mixture matched to the DMESH_ONELAW variance s2 and excess
kurtosis kurt via the same construction as the gate null:
  w = 1.8/(kurt+3), dd = sqrt(kurt/(3 w (1-w))), va = 1-w*dd, vb = 1+(1-w)*dd
(kurt=0 collapses to a single Gaussian). Every arm is evaluated against the
same law so differences reflect the surface alone. Integration: 61-node
Gauss-Hermite per component; binomial log-pmf with constants dropped
(identical across arms at fixed data, so ranking-invariant) -- kept for
comparability with reported absolute values (log C(n,k) included).

Usage: eval_mixture_nll.py <obs_csv> [more obs_csvs...]
Prints one line per file: path, m, mean_nll.
"""
import csv, math, sys
import numpy as np

S2   = [0.1428, 0.2415, 0.0286, 0.0038]
KURT = [8.5, 8.1, 5.4, 0.0]

def components(s2, kurt):
    if kurt <= 1e-9:
        return [(1.0, s2)]
    w = 1.8 / (kurt + 3.0)
    dd = math.sqrt(kurt / (3.0 * w * (1.0 - w)))
    va, vb = 1.0 - w * dd, 1.0 + (1.0 - w) * dd
    return [(1.0 - w, va * s2), (w, vb * s2)]

NODES, WEIGHTS = np.polynomial.hermite_e.hermegauss(61)
LOGW = np.log(WEIGHTS / math.sqrt(2.0 * math.pi))

def pool_nll(eta, n, k, s2, kurt):
    best = []
    lgc = math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    for wc, var in components(s2, kurt):
        u = NODES * math.sqrt(var)
        e = eta + u
        # stable log expit / log(1-expit)
        lp = -np.logaddexp(0.0, -e)
        lq = -np.logaddexp(0.0, e)
        ll = lgc + k * lp + (n - k) * lq + LOGW + math.log(wc)
        best.append(ll)
    allv = np.concatenate(best)
    m = allv.max()
    return -(m + math.log(np.exp(allv - m).sum()))

def eval_file(path):
    tot = 0.0
    m = 0
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["is_train"] == "1":
                continue
            n = float(r["n"]); k = float(r["k"])
            p = min(max(float(r["p_hat"]), 1e-12), 1 - 1e-12)
            eta = math.log(p / (1 - p))
            b = 0 if n < 128 else 1 if n < 512 else 2 if n < 4096 else 3
            tot += pool_nll(eta, n, k, S2[b], KURT[b])
            m += 1
    return m, tot / m

if __name__ == "__main__":
    for path in sys.argv[1:]:
        m, nll = eval_file(path)
        print(f"{path}  m={m}  mean_nll={nll:.4f}")
