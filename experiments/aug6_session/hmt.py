import pandas as pd, numpy as np
from scipy.special import expit, gammaln
from scipy.optimize import minimize

def build(hier):
    hv=hier
    hid=hv.set_index('id')
    free=hv[(hv.active==1)&(hv.free_coefficient==1)].copy()
    # delta-hat vs edge-parent midpoint; roots (depth 0 / no parents) excluded from shrinkage
    rows=[]
    for _,r in free.iterrows():
        if r.parent0<0 or r.depth==0: continue
        pm=0.5*(hid.loc[int(r.parent0),'height']+hid.loc[int(r.parent1),'height'])
        rows.append(dict(id=int(r.name) if False else int(r.id),depth=int(r.depth),x=r.x,y=r.y,
                         dh=(r.height-pm),xtx=max(r.current_xTx,1e-8),keff=r.current_keff,
                         p0=int(r.parent0),p1=int(r.parent1)))
    V=pd.DataFrame(rows).set_index('id')
    # tree parent: deeper edge parent, walked up through non-coefficient vertices
    isfree=set(V.index)
    def treepar(vid):
        r=hid.loc[vid]
        cands=[]
        for pid in [int(r.parent0),int(r.parent1)]:
            if pid>=0: cands.append(pid)
        if not cands: return -1
        pid=max(cands,key=lambda q:hid.loc[q,'depth'])
        seen=0
        while pid not in isfree and seen<20:
            rr=hid.loc[pid]
            if rr.parent0<0: return -1
            pid=max([int(rr.parent0),int(rr.parent1)],key=lambda q:hid.loc[q,'depth'])
            seen+=1
        return pid if pid in isfree else -1
    V['pa']=[treepar(v) for v in V.index]
    return V

def fit_hmt(V,iters=25,lam=None):
    Told=V.xtx.values/(V.xtx.values+lam) if lam is not None else np.ones(len(V))
    dh=V.dh.values/np.maximum(Told,0.25)   # invert prior shrinkage, capped
    u=dh*np.sqrt(V.xtx.values)
    phi=np.percentile(u**2,25)             # init low; refined in EM from the null component
    se2=phi/V.xtx.values
    Xc=np.column_stack([np.ones(len(V)),(V.depth.values-V.depth.mean()),
                        (np.log(V.xtx.values)-np.log(V.xtx.values).mean())])
    beta=np.array([np.log(max(np.var(dh)-se2.mean(),phi*0.1)),0.0,0.0])
    pi0,r0,r1=0.2,0.05,0.7
    ch={v:[] for v in V.index}
    for v,pa in V.pa.items():
        if pa>=0 and pa in ch: ch[pa].append(v)
    roots=[v for v,pa in V.pa.items() if pa<0 or pa not in ch]
    order=[]
    stack=list(roots)
    while stack:
        v=stack.pop(); order.append(v); stack+=ch[v]
    idx={v:i for i,v in enumerate(V.index)}
    def estep(beta,pi0,r0,r1):
        tau2=np.exp(Xc@beta)
        L0=np.exp(-0.5*dh**2/se2)/np.sqrt(se2)
        L1=np.exp(-0.5*dh**2/(se2+tau2))/np.sqrt(se2+tau2)
        Tr=np.array([[1-r0,r0],[1-r1,r1]])
        up={}; msg={}
        for v in reversed(order):
            i=idx[v]
            b=np.array([L0[i],L1[i]])
            for c in ch[v]: b=b*msg[c]
            b=b/max(b.sum(),1e-300)
            up[v]=b
            msg[v]=Tr@b
        q=np.zeros(len(V)); down={}; pair=np.zeros((2,2))
        for v in order:
            i=idx[v]
            if V.pa.loc[v]<0 or V.pa.loc[v] not in idx:
                d=np.array([1-pi0,pi0])
            else:
                u_=V.pa.loc[v]
                du=down[u_]*up[u_]/np.maximum(msg[v],1e-300)
                du=du/max(du.sum(),1e-300)
                d=du@Tr
                J=np.outer(du,np.ones(2))*Tr*np.outer(np.ones(2),up[v])
                pair+=J/max(J.sum(),1e-300)
            down[v]=d
            bel=d*up[v]; q[i]=bel[1]/max(bel.sum(),1e-300)
        return q,pair,tau2
    for it in range(iters):
        q,pair,tau2=estep(beta,pi0,r0,r1)
        w0=(1-q)
        if w0.sum()>3:
            phi=np.clip((w0*dh**2*V.xtx.values).sum()/w0.sum(),phi*0.2,phi*5) if it<5 else (w0*dh**2*V.xtx.values).sum()/w0.sum()
            se2=phi/V.xtx.values
        rootq=[q[idx[v]] for v in roots]
        pi0=np.clip(np.mean(rootq),1e-3,0.9)
        r0=np.clip(pair[0,1]/max(pair[0].sum(),1e-9),1e-3,0.5)
        r1=np.clip(pair[1,1]/max(pair[1].sum(),1e-9),0.5,0.999)
        def negm(b):
            t2=np.exp(np.clip(Xc@b,-12,8))
            return -(q*(-0.5*np.log(se2+t2)-0.5*dh**2/(se2+t2))).sum()+0.5*(b@b)/25
        beta=minimize(negm,beta,method='L-BFGS-B').x
    q,pair,tau2=estep(beta,pi0,r0,r1)
    shr=q*(tau2/(tau2+se2))
    return dict(dh_raw=dh,q=q,shrink=shr,beta=beta,pi0=pi0,r0=r0,r1=r1,phi=phi,tau2=tau2,V=V)

def reshrunk_heights(hier,V,shr,dh_used=None):
    hid=hier.set_index('id').copy()
    newh=hid.height.copy()
    base=dh_used if dh_used is not None else V.dh.values
    dh=dict(zip(V.index,base*shr))
    for vid in hid[hid.active==1].sort_values('depth').index:
        r=hid.loc[vid]
        if r.parent0<0: continue
        pm=0.5*(newh.loc[int(r.parent0)]+newh.loc[int(r.parent1)])
        if vid in dh: newh.loc[vid]=pm+dh[vid]
        elif r.free_coefficient==1: newh.loc[vid]=pm+(r.height-0.5*(hid.loc[int(r.parent0),'height']+hid.loc[int(r.parent1),'height']))
        else: newh.loc[vid]=pm
    return newh
