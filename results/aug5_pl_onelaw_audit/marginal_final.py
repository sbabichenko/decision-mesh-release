#!/usr/bin/env python3
"""Law-marginal final height solve at frozen topology (joint (eta,u) MAP under
the converged mixture law + EB prior). Topology, gate, EB lambdas all frozen
from a recommended-config dump; only free-vertex heights move; conformity
vertices track parent interpolation. Held-out pools carry zero weight and no u.
Usage: marginal_final.py <prefix> -> writes <prefix>_marginal_obs.csv, prints diag."""
import sys, numpy as np, pandas as pd
from scipy.special import expit

prefix = sys.argv[1]
d = pd.read_csv('smm_design_202606.csv')
n = d.n0.values.astype(float); k = d.k_term.values.astype(float)
o = pd.read_csv(f'{prefix}_obs.csv').sort_values('pid').reset_index(drop=True)
hv = pd.read_csv(f'{prefix}_hier_vertices.csv')
mesh = pd.read_csv(f'{prefix}_mesh.csv')
train = o.is_train.astype(bool).values
SUF=''
if len(sys.argv)>3 and sys.argv[3]=='all':
    train = np.ones(len(o),bool)   # post-hoc oracle: refit heights on ALL pools
    SUF='_all'
BASE = -4.7571; S = 0.01   # eta = BASE + S*h (centilogit heights), verified R2=1

# ---- geometry: pool -> (triangle vertex ids, barycentric) ----
act = hv[hv.active==1]
key = {(round(v.x,9),round(v.y,9)): int(v.id) for _,v in act.iterrows()}
H = {int(r.id): float(r.height) for _,r in hv.iterrows()}          # all vertices (incl. retired parents)
tris=[]; bary_ids=np.full((len(o),3),-1); bary_w=np.zeros((len(o),3))
ux=o.wala.values/60; uy=(o.wac.values-2.6)/5
A=[]
for _,t in mesh.iterrows():
    ids=[key.get((round(cx,9),round(cy,9))) for cx,cy in [(t.x0,t.y0),(t.x1,t.y1),(t.x2,t.y2)]]
    if None in ids: continue
    A.append((ids,(t.x0,t.y0,t.x1,t.y1,t.x2,t.y2)))
assigned=np.zeros(len(o),bool)
for ids,(x0,y0,x1,y1,x2,y2) in A:
    den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2)
    l0=((y1-y2)*(ux-x2)+(x2-x1)*(uy-y2))/den
    l1=((y2-y0)*(ux-x2)+(x0-x2)*(uy-y2))/den
    l2=1-l0-l1
    m=(l0>=-1e-9)&(l1>=-1e-9)&(l2>=-1e-9)&(~assigned)
    bary_ids[m]=ids; bary_w[m]=np.column_stack([l0,l1,l2])[m]; assigned|=m
assert assigned.all()

free = act[(act.free_coefficient==1)].copy()
conf = act[(act.free_coefficient==0)].sort_values('depth')
free_ids = free.id.astype(int).tolist()
lam = {int(r.id): float(r.lambda_v) for _,r in free.iterrows()}
parents = {int(r.id):(int(r.parent0),int(r.parent1)) for _,r in act.iterrows()}
# per-free-vertex pool lists with hat values
pools_of = {vid: ([],[]) for vid in free_ids}
for i in range(len(o)):
    for c in range(3):
        vid=int(bary_ids[i,c])
        if vid in pools_of: pools_of[vid][0].append(i); pools_of[vid][1].append(bary_w[i,c])
pools_of = {vid:(np.array(ii,dtype=int),np.array(bb,dtype=float)) for vid,(ii,bb) in pools_of.items()}

# ---- converged law ----
law=[(0.0374,1.362,0.1428),(0.0763,1.624,0.2415),(0.0180,0.313,0.0286),(0.0097,0.004,0.0038)]
strata=[(25,128),(128,512),(512,4096),(4096,10**9)]
s0v=np.zeros(len(o)); sLv=np.zeros(len(o)); piv=np.zeros(len(o))
for (lo,hi),(pi,sL,tot) in zip(strata,law):
    m=(n>=lo)&(n<hi); s0=(tot-pi*sL)/(1-pi)
    s0v[m]=max(s0,1e-6); sLv[m]=sL; piv[m]=pi

def eta_of(H):
    h=np.array([[H[int(bary_ids[i,c])] for c in range(3)] for i in range(len(o))])
    return BASE+S*(bary_w*h).sum(1)
def refresh_conf(H):
    for _,r in conf.iterrows():
        p0,p1=parents[int(r.id)]
        if p0>=0 and p1>=0: H[int(r.id)]=0.5*(H[p0]+H[p1])

# posterior mode + responsibility + Laplace component marginals (vectorized)
def u_step(eta):
    def mode(s2):
        uh=np.zeros(len(o))
        for _ in range(40):
            p=expit(eta+uh); uh=np.clip(uh-np.clip((n*p-k+uh/s2)/(n*p*(1-p)+1/s2),-2,2),-12,12)
        p=expit(eta+uh); c=n*p*(1-p)+1/s2
        ll=k*np.log(np.clip(p,1e-14,1))+(n-k)*np.log(np.clip(1-p,1e-14,1))-0.5*uh*uh/s2-0.5*np.log(s2*c)
        return uh,ll
    u0,l0=mode(s0v); uL,lL=mode(sLv)
    a=np.log(piv)+lL; b=np.log(1-piv)+l0
    mx=np.maximum(a,b); g=np.exp(a-mx)/(np.exp(a-mx)+np.exp(b-mx))
    return g*uL+(1-g)*u0, g

MAXIT=int(sys.argv[2]) if len(sys.argv)>2 else 10
best_ll=None; best_H=None; bad=0
Hc=dict(H); refresh_conf(Hc)
u=np.zeros(len(o)); gam=np.zeros(len(o))
for it in range(MAXIT):
    eta=eta_of(Hc)
    u_new,gam=u_step(eta); u_new[~train]=0.0
    if True:
        def _mll(eta):
            def mode(s2):
                uh=np.zeros(len(o))
                for _ in range(40):
                    p=expit(eta+uh); uh=np.clip(uh-np.clip((n*p-k+uh/s2)/(n*p*(1-p)+1/s2),-2,2),-12,12)
                p=expit(eta+uh); c=n*p*(1-p)+1/s2
                return k*np.log(np.clip(p,1e-14,1))+(n-k)*np.log(np.clip(1-p,1e-14,1))-0.5*uh*uh/s2-0.5*np.log(s2*c)
            l0=mode(s0v); lL=mode(sLv)
            a=np.log(piv)+lL; b=np.log(1-piv)+l0
            mx=np.maximum(a,b)
            return (mx+np.log(np.exp(a-mx)+np.exp(b-mx)))[train].sum()
        cur=_mll(eta)
        if it==0 or cur>best_ll: best_ll=cur; best_H=dict(Hc); bad=0
        else: bad+=1
        if 'VERBOSE' in dir(): print(f'  it{it}: train marginal ll {cur:.1f} | v1894 p {100*expit(BASE+S*Hc.get(1894,0)):.1f}%', flush=True)
        if bad>=2:
            Hc=best_H; break
    u=u_new                                      # envelope score needs the current mode
    # h-step: one damped tempered-Newton sweep (marginal ascent)
    for sweep in range(1):
        for vid in free_ids:
            ii,bb=pools_of[vid]
            m=train[ii]
            if not m.any(): continue
            i2,b2=ii[m],bb[m]
            eta_loc=BASE+S*np.array([sum(Hc[int(bary_ids[j,c])]*bary_w[j,c] for c in range(3)) for j in i2])
            p=expit(eta_loc+u[i2])
            w=n[i2]*p*(1-p)
            s2m=gam[i2]*sLv[i2]+(1-gam[i2])*s0v[i2]
            c_i=1.0/(1.0+w*s2m)                     # Laplace marginal curvature temper
            # envelope theorem: marginal score is the UNtempered score at the
            # posterior mode; only the curvature carries the 1/(1+w*s2) temper.
            score=(b2*(k[i2]-n[i2]*p)).sum()*S                  # d/dh, h in centilogits
            info=(c_i*w*b2*b2).sum()*S*S
            p0,p1=parents[vid]
            pint=0.5*(Hc[p0]+Hc[p1]) if p0>=0 and p1>=0 else Hc[vid]
            g=score-lam[vid]*(Hc[vid]-pint)
            hcurv=info+lam[vid]
            Hc[vid]=Hc[vid]+np.clip(0.5*g/max(hcurv,1e-12),-30,30)
        refresh_conf(Hc)
if best_H is not None: Hc=best_H
eta=eta_of(Hc)
pd.DataFrame({'id':list(Hc.keys()),'height':[Hc[i] for i in Hc]}).to_csv(f'{prefix}_marginal{SUF}_heights.csv',index=False)
out=o.copy(); out['p_hat']=expit(eta)
out.to_csv(f'{prefix}_marginal{SUF}_obs.csv',index=False)
for vid in [1894,1893,1831]:
    if vid in Hc and vid in H:
        print(f'v{vid}: {BASE+S*H[vid]:+.2f} -> {BASE+S*Hc[vid]:+.2f} log-odds (p {100*expit(BASE+S*H[vid]):.1f}% -> {100*expit(BASE+S*Hc[vid]):.1f}%)')
i0=o.index[o.pid==33681][0] if (o.pid==33681).any() else None
if i0 is not None: print('canary 33681 p:', round(float(expit(eta[i0])),4))
