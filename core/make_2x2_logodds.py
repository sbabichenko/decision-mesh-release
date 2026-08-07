#!/usr/bin/env python3
"""Standard 2x2 experiment figure for Decision Mesh runs (matches the June 2026
template): observed probability / fitted probability / held-out error / mesh.
Coordinates are recovered from the data csv (row-aligned with the obs dump),
avoiding the sim-convention scaling in the dump. Usage:
  make_2x2.py <data.csv> <dump_prefix> <title> <out.png>"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import tri as mtri

data_csv, prefix, title, out = sys.argv[1:5]
d = pd.read_csv(data_csv)
obs = pd.read_csv(f"{prefix}_obs.csv")
assert len(d) == len(obs), "row alignment"
wala, wac = d.wala.values, d.wac.values
n, k, phat = d.n0.values.astype(float), d.iloc[:,3].values.astype(float), obs.p_hat.values
def lo_(p): return np.log(np.clip(p,1e-6,1-1e-6)/np.clip(1-p,1e-6,1-1e-6))
MINMASS = 300.0  # loans per cell below which a cell is masked (kills small-sample speckle)
lo_w, hi_w = 0.0, wala.max()*1.02; lo_c, hi_c = wac.min()-0.2, wac.max()+0.2
BX, BY = 72, 56
def binned(vnum, vden):
    hn,_,_ = np.histogram2d(wala, wac, bins=[BX,BY], range=[[lo_w,hi_w],[lo_c,hi_c]], weights=vnum)
    hd,_,_ = np.histogram2d(wala, wac, bins=[BX,BY], range=[[lo_w,hi_w],[lo_c,hi_c]], weights=vden)
    with np.errstate(invalid="ignore"): r = hn/hd
    return np.ma.masked_invalid(r)
ext = [lo_w, hi_w, lo_c, hi_c]

fig, ax = plt.subplots(2, 2, figsize=(19.5, 12), dpi=300)
mass,_,_ = np.histogram2d(wala, wac, bins=[BX,BY], range=[[lo_w,hi_w],[lo_c,hi_c]], weights=n)
obs_surf = np.ma.filled(lo_(binned(k, n)), np.nan); obs_surf[mass < MINMASS] = np.nan
fit_surf = np.ma.filled(lo_(binned(n*phat, n)), np.nan); fit_surf[mass < MINMASS] = np.nan
v0o, v1o = np.percentile(obs_surf[np.isfinite(obs_surf)], [2, 98]); v0o, v1o = np.floor(v0o*2)/2, np.ceil(v1o*2)/2
v0f, v1f = np.percentile(fit_surf[np.isfinite(fit_surf)], [1, 99]); v0f, v1f = np.floor(v0f*10)/10, np.ceil(v1f*10)/10
im0 = ax[0,0].imshow(obs_surf.T, origin="lower", extent=ext, aspect="auto", cmap="viridis", vmin=v0o, vmax=v1o)
ax[0,0].set_title("Observed cumulative termination (log-odds)")

# The actual finite element function: gouraud-shaded piecewise-linear
# surface on the fitted triangulation itself, no rasterization.
import matplotlib.tri as mtri
mesh = pd.read_csv(f"{prefix}_mesh.csv")
eta0 = np.log(k.sum()/n.sum()/(1-k.sum()/n.sum()))
corners = np.vstack([mesh[[f"x{i}", f"y{i}", f"h{i}"]].values for i in (0,1,2)])
xy, inv = np.unique(np.round(corners[:, :2], 9), axis=0, return_inverse=True)
hv = np.zeros(len(xy)); hv[inv] = corners[:, 2]          # per-vertex heights
tris = inv.reshape(3, len(mesh)).T                        # triangle -> vertex ids
vx = lo_w + xy[:, 0] * (hi_w - lo_w)                      # unit square -> data coords
vy = lo_c + xy[:, 1] * (hi_c - lo_c)
veta = eta0 + hv / 100.0
tri = mtri.Triangulation(vx, vy, triangles=tris)
im1 = ax[0,1].tripcolor(tri, veta, shading="gouraud", cmap="viridis", vmin=v0f, vmax=v1f)
ax[0,1].triplot(tri, color="white", lw=0.15, alpha=0.5)
ax[0,1].set_xlim(lo_w, hi_w); ax[0,1].set_ylim(lo_c, hi_c)
ax[0,1].set_title("Adaptive-mesh fitted: the P1 element function on its triangulation")
cb = fig.colorbar(im0, ax=ax[0,0], shrink=0.9, pad=0.02); cb.set_label("observed log-odds")
cb1 = fig.colorbar(im1, ax=ax[0,1], shrink=0.9, pad=0.02); cb1.set_label("fitted log-odds")

held = ~obs['is_train'].astype(bool).values if 'is_train' in obs.columns else (np.arange(len(d)) % 2 == 1)
# per-pool residual in log-odds and law-aware precision, cell-aggregated
ph_e = np.clip((k+0.5)/(n+1), 1e-6, 1-1e-6)
r_i = np.log(ph_e/(1-ph_e)) - lo_(phat)
l2n = np.log2(np.maximum(n, 32))
s2_i = np.interp(l2n, [5.5, 7.6, 10.0, 13.5], [0.60, 0.46, 0.47, 0.33])   # stratum totals
w_i = np.where(held, 1.0/(1.0/((n+1)*ph_e*(1-ph_e)) + s2_i), 0.0)
num,_,_ = np.histogram2d(wala, wac, bins=[BX,BY], range=[[lo_w,hi_w],[lo_c,hi_c]], weights=w_i*r_i)
den,_,_ = np.histogram2d(wala, wac, bins=[BX,BY], range=[[lo_w,hi_w],[lo_c,hi_c]], weights=w_i)
zmap = np.where(den > 0, num/np.maximum(den,1e-12)*np.sqrt(den), np.nan)
zmap[mass < MINMASS] = np.nan
im2 = ax[1,0].imshow(zmap.T, origin="lower", extent=ext, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax[1,0].set_title("Held-out miss, standardized by the law (cell z)")
cb2 = fig.colorbar(im2, ax=ax[1,0], shrink=0.9, pad=0.02); cb2.set_label("z (fast +, slow -)")

mesh = pd.read_csv(f"{prefix}_mesh.csv")
pts = {}; tris = []
def vid(x, y):
    key = (round(x,12), round(y,12))
    if key not in pts: pts[key] = len(pts)
    return pts[key]
for _, r in mesh.iterrows():
    tris.append([vid(r.x0,r.y0), vid(r.x1,r.y1), vid(r.x2,r.y2)])
xy = np.array(list(pts.keys()))
mx = wala.min() + xy[:,0]*(wala.max()-wala.min())
my = wac.min() + xy[:,1]*(wac.max()-wac.min())
# Bottom right: refinement depth as face fill, the hierarchical
# (bisection-tree) topology as parent edges, and the shrinkage state as
# vertex dots sized by the post-shrinkage hierarchical surplus.
hv_df = pd.read_csv(f"{prefix}_hier_vertices.csv")
hv_df = hv_df[hv_df.active == 1].reset_index(drop=True)
hx = lo_w + hv_df.x.values * (hi_w - lo_w); hy = lo_c + hv_df.y.values * (hi_c - lo_c)
pos = {int(i): (px, py) for i, px, py in zip(hv_df.id, hx, hy)}
hgt = {int(i): h for i, h in zip(hv_df.id, hv_df.new_height)}
face_depth = np.round(np.log2(0.5 / np.maximum(
    0.5*np.abs((mesh.x1-mesh.x0)*(mesh.y2-mesh.y0)-(mesh.x2-mesh.x0)*(mesh.y1-mesh.y0)), 1e-12))).astype(int)
im3 = ax[1,1].tripcolor(tri, facecolors=face_depth.astype(float), cmap="plasma",
                        vmin=face_depth.min(), vmax=face_depth.max(), alpha=0.28)
cb3 = fig.colorbar(im3, ax=ax[1,1], shrink=0.9, pad=0.02); cb3.set_label("refinement depth (bisections)")
segs = []
for vid, p0, p1 in zip(hv_df.id, hv_df.parent0, hv_df.parent1):
    if p0 >= 0 and int(p0) in pos: segs.append([pos[int(vid)], pos[int(p0)]])
    if p1 >= 0 and int(p1) in pos: segs.append([pos[int(vid)], pos[int(p1)]])
from matplotlib.collections import LineCollection
ax[1,1].add_collection(LineCollection(segs, colors="#14213d", linewidths=0.35, alpha=0.55))
# dot size decreasing in depth: coarse hierarchy levels read as landmarks,
# the dense fine-scale regions stay uncluttered
dsize = 26.0 * 0.70 ** hv_df.depth.values + 0.4
ax[1,1].scatter(hx, hy, s=dsize, c="#14213d", edgecolors="white", linewidths=0.2, alpha=0.9)
ax[1,1].set_xlim(lo_w, hi_w); ax[1,1].set_ylim(lo_c, hi_c)
ax[1,1].set_xlabel("WALA (months)"); ax[1,1].set_ylabel("WAC (%)")
ax[1,1].set_title("Shrinkage topology: bisection-tree edges, vertex size = coarseness, muted depth fill")
for a in ax.flat: a.set_xlabel("WALA (months)"); a.set_ylabel("WAC (%)")
fig.suptitle(title, fontsize=17)
fig.savefig(out, bbox_inches="tight")
print(f"wrote {out}: {len(mesh)} faces, {len(xy)} vertices")
