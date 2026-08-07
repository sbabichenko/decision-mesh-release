#!/usr/bin/env python3
"""t_split acceptance battery: 20 locally balanced seeded holdouts.

Usage: python3 t_split.py <binary> <real_pools_all.csv> <law_shared.json>
Pass/fail: 0/20 with NLL > 3.2; report mean/sd/worst.
Extra env for the runs can be set in the caller's environment.
"""
import numpy as np, pandas as pd, json, subprocess, os, sys, tempfile
from scipy.special import expit, gammaln
B, DATA, LAW = sys.argv[1], sys.argv[2], sys.argv[3]
law = json.load(open(LAW))['law']
xs, ws = np.polynomial.hermite_e.hermegauss(20); lw = np.log(ws/ws.sum())
def loglik(kk,nn,eta):
    p = np.clip(expit(eta),1e-14,1-1e-14); return kk*np.log(p)+(nn-kk)*np.log1p(-p)
def comp(kk,nn,gg,s2):
    uh = np.zeros(len(kk))
    for _ in range(45):
        p = expit(gg+uh); uh = np.clip(uh-np.clip((nn*p-kk+uh/s2)/(nn*p*(1-p)+1/s2),-2,2),-12,12)
    p = expit(gg+uh); s = 1/np.sqrt(nn*p*(1-p)+1/s2)
    T = np.empty((20,len(kk)))
    for j in range(20):
        uu = uh+s*xs[j]; T[j] = lw[j]+0.5*xs[j]**2+loglik(kk,nn,gg+uu)-0.5*uu**2/s2
    m = T.max(0)
    return m+np.log(np.exp(np.clip(T-m,-700,0)).sum(0))+np.log(s)-0.5*np.log(s2)
def nll_held(g, n, k, held):
    const = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
    out = np.full(len(n), np.nan)
    for (lo,hi),(pi,s0,sL,_t) in zip([(25,128),(128,512),(512,4096),(4096,10**9)], law):
        te = held&(n>=lo)&(n<hi)
        a = np.log(pi)+comp(k[te],n[te],g[te],sL); b = np.log(1-pi)+comp(k[te],n[te],g[te],s0)
        mx = np.maximum(a,b); out[te] = -(mx+np.log(np.exp(a-mx)+np.exp(b-mx))+const[te])
    return float(np.nanmean(out[held]))
d0 = pd.read_csv(DATA)
seeds = [0,101,202,303,404,505,606,707,808,909,11,22,33,44,55,66,77,88,99,110]
vals = []
tmp = tempfile.mkdtemp()
for seed in seeds:
    d = d0
    env = dict(os.environ, DMESH_DATA=os.path.abspath(DATA), DMESH_SPLIT="1",
               DMESH_SPLIT_MODE="local", DMESH_SPLIT_BLOCK_SIZE="24",
               DMESH_DUMP=f"{tmp}/r{seed}")
    subprocess.run([B,"0","0.10","1","24",str(seed)], env=env, capture_output=True)
    n = d.n0.values.astype(float); k = d.k_term.values.astype(float)
    o = pd.read_csv(f'{tmp}/r{seed}_obs.csv').sort_values('pid').reset_index(drop=True)
    held = ~o.is_train.astype(bool).values
    ph = np.clip(o.p_hat.values,1e-12,1-1e-12)
    v = nll_held(np.log(ph/(1-ph)), n, k, held)
    vals.append(v)
    print(f"seed {seed:>4}: {v:.4f}{'  <-- ERUPTION' if v>3.2 else ''}", flush=True)
V = np.array(vals); er = V>3.2
print(f"\nt_split: eruptions {int(er.sum())}/20 | mean {V.mean():.4f} | sd {V.std():.4f} | worst {V.max():.4f}")
print("PASS" if er.sum()==0 else "FAIL")
