import pandas as pd, numpy as np, os, glob, json
BASE='/mnt/data/true_null_cert_compare'
REAL_RANGE=(0.0,359.0,1.831,9.0)
rows=[]
for ddir in sorted(glob.glob(BASE+'/*_*')):
    if not os.path.isdir(ddir) or not os.path.exists(ddir+'/run_hier_vertices.csv'): continue
    name=os.path.basename(ddir)
    h=pd.read_csv(ddir+'/run_hier_vertices.csv'); m=pd.read_csv(ddir+'/run_mesh.csv'); o=pd.read_csv(ddir+'/run_obs.csv')
    active=h.active.eq(1); free=h.free_coefficient.eq(1); gate=h.gate_admitted.eq(1)
    w=np.maximum(o['info'].to_numpy(float),1e-12); g=o.g_hat.to_numpy(float)
    mu=np.average(g,weights=w); grms=np.sqrt(np.average((g-mu)**2,weights=w))
    # triangle area depth and certified window
    ar=.5*np.abs((m.x1-m.x0)*(m.y2-m.y0)-(m.x2-m.x0)*(m.y1-m.y0))
    dep=np.rint(np.log2(.5/np.maximum(ar,1e-12))).astype(int)
    cx=(m.x0+m.x1+m.x2)/3*359
    cy=(m.y0+m.y1+m.y2)/3*(9-1.831)+1.831
    win=(cx>=105)&(cx<=175)&(cy>=3.5)&(cy<=5.5)
    deep=win&(dep>=8)
    # obs window quadratic residual
    ow=(o.wala.between(105,175)&o.wac.between(3.5,5.5)).to_numpy()
    qrms=np.nan
    if ow.sum()>20:
        x=(o.loc[ow,'wala'].to_numpy()-140)/35; y=(o.loc[ow,'wac'].to_numpy()-4.5)/1.0
        X=np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y]); ww=w[ow]; gg=g[ow]
        A=X.T@(ww[:,None]*X); b=X.T@(ww*gg); beta=np.linalg.lstsq(A,b,rcond=None)[0]
        qrms=np.sqrt(np.average((gg-X@beta)**2,weights=ww))
    # candidate counts in window based vertex coords
    vx=h.x*359; vy=h.y*(9-1.831)+1.831
    vwin=vx.between(105,175)&vy.between(3.5,5.5)
    rows.append(dict(run=name,vertices=len(h),faces=len(m),active=int(active.sum()),gate_free=int((active&free&gate).sum()),gate_in_window=int((active&free&gate&vwin).sum()),gate_deep_in_window=int((active&free&gate&vwin&h.depth.ge(8)).sum()),max_hdepth=int(h.loc[active,'depth'].max()),global_spatial_rms_centi=grms,window_nonquad_rms_centi=qrms,window_faces=int(win.sum()),window_deep_faces=int(deep.sum()),max_face_depth=int(dep.max()),runtime=float(open(ddir+'/time.txt').read().strip()) if os.path.exists(ddir+'/time.txt') else np.nan))
res=pd.DataFrame(rows)
print(res.to_string(index=False,float_format=lambda x:f'{x:.3f}'))
res.to_csv(BASE+'/comparison_metrics.csv',index=False)
