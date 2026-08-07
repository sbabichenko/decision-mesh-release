#!/usr/bin/env python3
"""Standing evaluation harness for the clampless line.

Usage: python3 eval_harness.py <obs_dump.csv> <real_pools_all.csv> <law_shared.json>
Reports: held-out mixture NLL (all / ramp / mid / seasoned), pooled-rate
base profile WAC 5.75-6.25 (never empirical-logit dots at the edge), and
fitted base value. Obs dumps are pid-aligned to the data csv (dump
wala/wac columns are internal-scale; do not use them).
"""
import sys, json
import numpy as np, pandas as pd
from scipy.special import expit, gammaln

def evaluate(obs_csv, data_csv, law_json, label=""):
    o = pd.read_csv(obs_csv)
    if 'pid' in o.columns: o = o.sort_values('pid').reset_index(drop=True)
    d = pd.read_csv(data_csv)
    n = d.n0.values.astype(float); k = d.k_term.values.astype(float)
    wala, wac = d.wala.values, d.wac.values
    held = ~o['is_train'].astype(bool).values if 'is_train' in o.columns else (np.arange(len(d)) % 2 == 1)
    law = json.load(open(law_json))['law']
    xs, ws = np.polynomial.hermite_e.hermegauss(20); lw = np.log(ws/ws.sum())
    const = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
    def loglik(kk,nn,eta):
        p = np.clip(expit(eta),1e-14,1-1e-14); return kk*np.log(p)+(nn-kk)*np.log1p(-p)
    def comp(kk,nn,gg,s2):
        uh = np.zeros(len(kk))
        for _ in range(60):
            p = expit(gg+uh); uh = np.clip(uh-np.clip((nn*p-kk+uh/s2)/(nn*p*(1-p)+1/s2),-2,2),-12,12)
        p = expit(gg+uh); s = 1/np.sqrt(nn*p*(1-p)+1/s2)
        T = np.empty((20,len(kk)))
        for j in range(20):
            uu = uh+s*xs[j]; T[j] = lw[j]+0.5*xs[j]**2+loglik(kk,nn,gg+uu)-0.5*uu**2/s2
        m = T.max(0)
        return m+np.log(np.exp(np.clip(T-m,-700,0)).sum(0))+np.log(s)-0.5*np.log(s2)
    pcol = 'p_hat' if 'p_hat' in o.columns else 'phat'
    ph = np.clip(o[pcol].values,1e-12,1-1e-12)
    g = np.log(ph/(1-ph))
    out = np.full(len(d), np.nan)
    for (lo,hi),(pi,s0,sL,tot) in zip([(25,128),(128,512),(512,4096),(4096,10**9)], law):
        te = held&(n>=lo)&(n<hi)
        a = np.log(pi)+comp(k[te],n[te],g[te],sL); b = np.log(1-pi)+comp(k[te],n[te],g[te],s0)
        mx = np.maximum(a,b); out[te] = -(mx+np.log(np.exp(a-mx)+np.exp(b-mx))+const[te])
    ramp = held&(wala<24)&(wac>=5.5)&(wac<7.0)
    mid  = held&(wala>=60)&(wala<200)
    seas = held&(wala>=200)&(wala<300)&(wac>=5.5)&(wac<7.0)
    res = dict(nll_all=float(np.nanmean(out[held])), nll_ramp=float(np.nanmean(out[ramp])),
               nll_mid=float(np.nanmean(out[mid])), nll_seasoned=float(np.nanmean(out[seas])))
    print(f"{label or obs_csv}: NLL all {res['nll_all']:.4f} | ramp {res['nll_ramp']:.3f} | mid {res['nll_mid']:.3f} | seasoned {res['nll_seasoned']:.3f}")
    # Diagnostic standard v2: NO empirical logits, and NO fixed-effects
    # weighting. The band reference is the marginal MLE of a constant eta
    # under the mixture law (megas discounted by their own pool effect,
    # exactly as the model prescribes); the fitted band average uses the
    # law-tempered information weight w = 1/(1/(n p_b (1-p_b)) + s2_strat)
    # at the bin's reference rate. Per-pool references are k/n with the
    # likelihood value.
    def law_of(nn):
        si = np.digitize(nn, [128,512,4096])
        return si
    def band_mle(m2):
        kk, nn = k[m2], n[m2]
        si = law_of(nn)
        grid = np.linspace(-10, 3, 40)
        def negll(e):
            tot = 0.0
            for s_ in range(4):
                ms = si==s_
                if ms.sum()==0: continue
                pi_, s0_, sL_, _ = law[s_]
                a = np.log(pi_)+comp(kk[ms],nn[ms],np.full(ms.sum(),e),sL_)
                b = np.log(1-pi_)+comp(kk[ms],nn[ms],np.full(ms.sum(),e),s0_)
                mx = np.maximum(a,b); tot += (mx+np.log(np.exp(a-mx)+np.exp(b-mx))).sum()
            return -tot
        v = [negll(e) for e in grid]
        e0 = grid[int(np.argmin(v))]
        fine = np.linspace(e0-0.4, e0+0.4, 17)
        v2 = [negll(e) for e in fine]
        return fine[int(np.argmin(v2))]
    mm = held&(wac>=5.75)&(wac<6.25)
    print("  base profile (law-marginal MLE dot | fit, law-weighted):")
    for a in range(0,9,2):
        m2 = mm&(wala>=a)&(wala<a+2)
        if m2.sum() < 5: continue
        dot = band_mle(m2)
        p_b = 1/(1+np.exp(-dot))
        s2b = np.array([law[s_][1]**2 for s_ in law_of(n[m2])])  # s0^2 per pool's stratum
        wlaw = 1/(1/np.maximum(n[m2]*p_b*(1-p_b),1e-9) + s2b)
        print(f"    WALA {a}-{a+2}: {dot:+.2f} | {np.average(g[m2],weights=wlaw):+.2f}")
    return res

if __name__ == "__main__":
    evaluate(sys.argv[1], sys.argv[2], sys.argv[3])
