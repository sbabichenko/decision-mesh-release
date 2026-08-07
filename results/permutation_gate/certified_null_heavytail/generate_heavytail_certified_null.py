#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path('/mnt/data/lindsey_delivery/lindsey_quantitative_delivery_v2/experiments/corrected_homogeneous_certified_null')
OUT=Path('/mnt/data/heavytail_certified_null_inputs_2026-08-01')
OUT.mkdir(exist_ok=True)

def expit(x):
    x=np.asarray(x,float)
    out=np.empty_like(x)
    pos=x>=0
    out[pos]=1/(1+np.exp(-x[pos]))
    ex=np.exp(x[~pos]); out[~pos]=ex/(1+ex)
    return out

def plateau_cos(z,ol,il,ih,oh):
    z=np.asarray(z); a=np.zeros_like(z,float); a[(z>=il)&(z<=ih)]=1
    left=(z>ol)&(z<il); a[left]=.5*(1-np.cos(np.pi*(z[left]-ol)/(il-ol)))
    right=(z>ih)&(z<oh); a[right]=.5*(1+np.cos(np.pi*(z[right]-ih)/(oh-ih)))
    return a

def mixture_effect(rng,size,total_var=.35,wide_weight=.10,variance_ratio=10.0):
    core_var=total_var/((1-wide_weight)+wide_weight*variance_ratio)
    wide_var=variance_ratio*core_var
    wide=rng.random(size)<wide_weight
    sd=np.where(wide,np.sqrt(wide_var),np.sqrt(core_var))
    return rng.normal(0,sd),wide,core_var,wide_var

df=pd.read_csv(ROOT/'instruments/ginnie_certified_null_v3.csv')
beta=np.load(ROOT/'instruments/null_truth_beta.npy')
wala=df.wala.to_numpy(float); wac=df.wac.to_numpy(float); n0=df.n0.to_numpy(float)
x=2*wala/359-1; y=2*(wac-1.831)/(9-1.831)-1
X=np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y,x*x*y,x*y*y,x*x*y*y]); poly=X@beta
W=(wala>=80)&(wala<=200)&(wac>=3)&(wac<=6); xc=(wala-140)/60; yc=(wac-4.5)/1.5
Q=np.column_stack([np.ones(len(x)),xc,yc,xc*xc,xc*yc,yc*yc]); qb=np.linalg.solve(Q[W].T@(n0[W,None]*Q[W]),Q[W].T@(n0[W]*poly[W])); quad=Q@qb
blend=plateau_cos(wala,80,105,175,200)*plateau_cos(wac,3,3.5,5.5,6); eta=(1-blend)*poly+blend*quad
np.savez(OUT/'hcn_heavytail_truth.npz',eta=eta,polynomial=poly,quadratic=quad,blend=blend,quad_beta=qb,wala=wala,wac=wac,n0=n0)
ledger=[]
for seed in range(5):
    rng=np.random.default_rng(np.random.SeedSequence([20260801, 8675309, seed]))
    u,wide,core_var,wide_var=mixture_effect(rng,len(df))
    k=rng.binomial(n0.astype(np.int64),expit(eta+u))
    out=df[['wala','wac','n0']].copy(); out['k_term']=k
    path=OUT/f'hcn_heavytail_seed{seed}.csv'; out.to_csv(path,index=False)
    np.savez(OUT/f'hcn_heavytail_latent_seed{seed}.npz',u=u,wide=wide)
    ledger.append(dict(seed=seed,rows=len(out),total_exposure=int(n0.sum()),pooled_rate=float(k.sum()/n0.sum()),wide_count=int(wide.sum()),wide_share=float(wide.mean()),sample_mean_u=float(u.mean()),sample_var_u=float(u.var()),core_var=core_var,wide_var=wide_var,total_target_var=.35,wide_weight=.10,variance_ratio=10.0))
pd.DataFrame(ledger).to_csv(OUT/'generation_ledger.csv',index=False)
print(pd.DataFrame(ledger).to_string(index=False))
