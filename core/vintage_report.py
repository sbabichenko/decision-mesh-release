#!/usr/bin/env python3
"""Layer-1 monthly vintage report for the Decision Mesh overlay.

Usage:  python3 vintage_report.py <monthlySFPS_YYYYMM.txt|pools.csv> <outdir>

Pipeline: parse -> composite q30 bootstrap fit -> mixture EM per stratum
-> one-law q30 production fit -> metrics (held-out deviance, mixture NLL,
PIT tail occupancy) -> lump census -> one-page PDF. Saves state JSON so the
next vintage's page shows deltas. Requires the decision_mesh binary and
make_2x2.py; paths configurable below.
"""
import sys, os, json, subprocess, tempfile, datetime
import numpy as np, pandas as pd
from scipy.special import expit, gammaln
from scipy.stats import binom

BIN   = os.environ.get("DMESH_BIN", "/home/claude/decision_mesh_latest/build/decision_mesh")
MK2X2 = os.environ.get("MK2X2", "/home/claude/realdata/make_2x2.py")
CORE  = {"DMESH_SIGMA0_MIN":"0.3","DMESH_SIGMA0_MAX":"6","DMESH_PI0_MIN":"0.05",
         "DMESH_TAU_FLOOR_LOGIT_SD":"0.005","DMESH_SPLIT":"1"}
STRATA = [(25,128),(128,512),(512,4096),(4096,10**9)]
XS, WS = np.polynomial.hermite_e.hermegauss(20); LW = np.log(WS/WS.sum())

def parse_raw(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\r\n").split("|")
            if len(p) < 20 or p[0] != "PS" or p[4] != "SF": continue
            try:
                rows.append((float(p[8]), int(p[13]), float(p[15]) if p[15] else np.nan,
                             float(p[17]) if p[17] else np.nan, float(p[19]) if p[19] else np.nan))
            except ValueError: pass
    df = pd.DataFrame(rows, columns=["orig","n_now","aols","wac","wala"])
    df["n0"] = np.round(df.orig/df.aols)
    df = df[np.isfinite(df.n0)&np.isfinite(df.wac)&np.isfinite(df.wala)]
    df = df[(df.n0>=25)&(df.n_now<=df.n0)&(df.wac>=1.5)&(df.wac<=9)&(df.wala>=0)&(df.wala<=360)]
    df["k_term"] = (df.n0-df.n_now).astype(int)
    return df[["wala","wac","n0","k_term"]].reset_index(drop=True)

def run_fit(csv, outdir, tag, extra_env, q):
    env = dict(os.environ); env.update(CORE); env.update(extra_env)
    env["DMESH_DATA"] = csv; env["DMESH_DUMP"] = os.path.join(outdir, tag)
    env["DMESH_FINAL_RESHRINK"] = "1"
    out = subprocess.run([BIN,"0",q,"1","24","7"], env=env, capture_output=True, text=True).stdout
    dev = [l for l in out.split("\n") if "deviance/pool" in l]
    vtx = [l for l in out.split("\n") if "vertices" in l]
    return (float(dev[-1].split("deviance/pool")[-1].split()[0]) if dev else np.nan,
            int(vtx[-1].split("vertices")[-1].split()[0].rstrip(",")) if vtx else -1)

def loglik(k,n,eta):
    p = np.clip(expit(eta),1e-14,1-1e-14); return k*np.log(p)+(n-k)*np.log1p(-p)

def component(k,n,g,s2, want_u=False):
    uh = np.zeros(len(k))
    for _ in range(60):
        p = expit(g+uh); step = (n*p-k+uh/s2)/(n*p*(1-p)+1/s2)
        uh = np.clip(uh-np.clip(step,-2,2),-12,12)
    p = expit(g+uh); c = n*p*(1-p)+1/s2; s = 1/np.sqrt(c)
    T = np.empty((20,len(k))); U = np.empty((20,len(k)))
    for j in range(20):
        uu = uh+s*XS[j]; T[j] = LW[j]+0.5*XS[j]**2+loglik(k,n,g+uu)-0.5*uu**2/s2; U[j]=uu
    m = T.max(0); E = np.exp(np.clip(T-m,-700,0)); Z = E.sum(0)
    lm = m+np.log(Z)+np.log(s)-0.5*np.log(s2)
    return (lm, (E*U).sum(0)/Z, (E*U*U).sum(0)/Z) if want_u else (lm, (E*U*U).sum(0)/Z)

def em(k,n,g, iters=20):
    pi,s0,sL = 0.05,0.10,2.0
    for _ in range(iters):
        lm0,E0 = component(k,n,g,s0); lmL,EL = component(k,n,g,sL)
        a = np.log(pi)+lmL; b = np.log(1-pi)+lm0; mx = np.maximum(a,b)
        gam = np.exp(a-mx)/(np.exp(a-mx)+np.exp(b-mx))
        pi = float(np.clip(gam.mean(),1e-5,0.5)); s0 = float(np.sum((1-gam)*E0)/np.sum(1-gam))
        sL = float(max(np.sum(gam*EL)/max(np.sum(gam),1e-9), s0))
    return pi,s0,sL

def mix_cdf(kq, nn, gg, pi, s0, sL):
    F = np.zeros(len(kq))
    for wgt, s2 in [(1-pi, s0), (pi, sL)]:
        s = np.sqrt(max(s2,1e-6))
        for j in range(20):
            F += wgt*(WS[j]/WS.sum())*binom.cdf(kq, nn, expit(gg + s*XS[j]))
    return F

def main(inp, outdir, vintage=None):
    os.makedirs(outdir, exist_ok=True)
    vintage = vintage or "".join(c for c in os.path.basename(inp) if c.isdigit())[:6] or "unknown"
    csv = os.path.join(outdir, f"pools_{vintage}.csv")
    df = pd.read_csv(inp) if inp.endswith(".csv") else parse_raw(inp)
    df.to_csv(csv, index=False)
    n = df.n0.values.astype(float); k = df.k_term.values.astype(float)
    even = np.arange(len(df))%2==0

    # 1) bootstrap surface, 2) law, 3) production fit
    run_fit(csv, outdir, "boot", {"DMESH_GATE_EXTRA_VAR":"0.35"}, "0.30")
    gb = pd.read_csv(os.path.join(outdir,"boot_obs.csv")).p_hat.values
    gb = np.log(np.clip(gb,1e-9,1-1e-9)/np.clip(1-gb,1e-9,1-1e-9))
    law, kurts = [], []
    for lo,hi in STRATA:
        m = (n>=lo)&(n<hi)&even
        pi,s0,sL = em(k[m], n[m], gb[m])
        tot = (1-pi)*s0+pi*sL
        law.append((pi,s0,sL,tot))
        kurts.append(3*((1-pi)*s0**2+pi*sL**2)/tot**2-3)
    LAW = ",".join(f"{t[3]:.3f}" for t in law) + ";" + ",".join(f"{ku:.2f}" for ku in kurts)
    dev, vtx = run_fit(csv, outdir, "prod", {"DMESH_ONELAW":LAW}, "0.30")

    g = pd.read_csv(os.path.join(outdir,"prod_obs.csv")).p_hat.values
    g = np.log(np.clip(g,1e-9,1-1e-9)/np.clip(1-g,1e-9,1-1e-9))
    # metrics + census on the production surface
    const = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
    nll, cnt, pit, census = 0.0, 0, [], []
    rng = np.random.default_rng(1)
    for (lo,hi),(pi,s0,sL,tot) in zip(STRATA, law):
        m = (n>=lo)&(n<hi); tr, te = m&even, m&~even
        pi2,s02,sL2 = em(k[tr], n[tr], g[tr], iters=15)
        lm0,_ = component(k[te],n[te],g[te],s02); lmL,_ = component(k[te],n[te],g[te],sL2)
        a = np.log(pi2)+lmL; b = np.log(1-pi2)+lm0; mx = np.maximum(a,b)
        nll -= (mx+np.log(np.exp(a-mx)+np.exp(b-mx)) + const[te]).sum(); cnt += te.sum()
        Fk  = mix_cdf(k[te], n[te], g[te], pi2, s02, sL2)
        Fk1 = mix_cdf(k[te]-1, n[te], g[te], pi2, s02, sL2)
        pit.append(Fk1 + rng.uniform(size=te.sum())*(Fk-Fk1))
        _, Eu, _ = component(k[m], n[m], g[m], sL2, want_u=True)
        lm0m,_ = component(k[m],n[m],g[m],s02); lmLm,_ = component(k[m],n[m],g[m],sL2)
        am = np.log(pi2)+lmLm; bm = np.log(1-pi2)+lm0m; mxm = np.maximum(am,bm)
        gam = np.exp(am-mxm)/(np.exp(am-mxm)+np.exp(bm-mxm))
        lum = gam>0.5
        census.append((int(lum.sum()), float(Eu[lum].mean()) if lum.any() else 0.0,
                       float(np.mean(Eu[lum]>0)) if lum.any() else 0.0))
    pit = np.concatenate(pit)
    metrics = {"vintage": vintage, "pools": int(len(df)), "vertices": vtx,
               "heldout_dev": round(dev,3), "mix_nll": round(nll/cnt,4),
               "pit_tail": round(float(np.mean((pit<0.05)|(pit>0.95))),4),
               "law": [[round(x,4) for x in t] for t in law],
               "kurts": [round(x,1) for x in kurts],
               "census": census}
    # 2x2
    subprocess.run(["python3", MK2X2, csv, os.path.join(outdir,"prod"),
                    f"Ginnie {vintage} - one-law gate q=0.30 - held-out adaptive Decision Mesh",
                    os.path.join(outdir,"panel_2x2.png")], check=True)
    # prior state for deltas
    state_path = os.path.join(outdir, "vintage_state.json")
    prior = json.load(open(state_path)) if os.path.exists(state_path) else None
    render(outdir, metrics, prior)
    json.dump(metrics, open(state_path,"w"), indent=1)
    print("report:", os.path.join(outdir, f"vintage_report_{vintage}.pdf"))

def render(outdir, M, prior):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, matplotlib.image as mpimg
    fig = plt.figure(figsize=(15.5, 8.7))
    fig.suptitle(f"Decision Mesh vintage report - Ginnie Mae {M['vintage']}", fontsize=15, y=0.985)
    axA = fig.add_axes([0.015, 0.05, 0.60, 0.88]); axA.axis("off")
    axA.imshow(mpimg.imread(os.path.join(outdir,"panel_2x2.png")))
    def delta(key, cur, fmt="+.3f"):
        if not prior or key not in prior: return ""
        return f"  ({cur-prior[key]:{fmt}} vs prior)"
    ax = fig.add_axes([0.635, 0.60, 0.35, 0.31]); ax.axis("off")
    lines = [f"pools {M['pools']:,}    mesh vertices {M['vertices']:,}",
             f"held-out deviance/pool   {M['heldout_dev']:.3f}{delta('heldout_dev',M['heldout_dev'])}",
             f"held-out mixture NLL     {M['mix_nll']:.4f}{delta('mix_nll',M['mix_nll'],'+.4f')}",
             f"PIT tail occupancy       {M['pit_tail']*100:.1f}%  (nominal 10.0%)",
             "", "certification ledger (standing benchmarks, current code):",
             "  composite gate: 2 (v3), 7 (v4) certified-false deep faces"]
    ax.text(0, 1, "Metrics\n" + "\n".join(lines), va="top", family="monospace", fontsize=9.5)
    ax = fig.add_axes([0.635, 0.30, 0.35, 0.27]); ax.axis("off")
    rows = ["stratum      pi     s0^2   sL^2   total  kurt"]
    for (lo,hi), t, ku in zip(STRATA, M["law"], M["kurts"]):
        lab = f"{lo}-{'inf' if hi>10**8 else hi}"
        rows.append(f"{lab:<11}{t[0]:6.3f} {t[1]:6.3f} {t[2]:6.2f} {t[3]:6.3f} {ku:5.1f}")
    ax.text(0, 1, "Pool-effect law (mixture EM, train half)\n" + "\n".join(rows),
            va="top", family="monospace", fontsize=9.5)
    ax = fig.add_axes([0.635, 0.05, 0.35, 0.22]); ax.axis("off")
    rows = ["stratum      lumps   mean u   share fast"]
    for (lo,hi), (c, mu, fr) in zip(STRATA, M["census"]):
        lab = f"{lo}-{'inf' if hi>10**8 else hi}"
        rows.append(f"{lab:<11}{c:6d}  {mu:+7.2f}   {fr*100:5.0f}%")
    note = "" if not prior else f"\nlump total {sum(c for c,_,_ in M['census'])} vs prior {sum(c for c,_,_ in prior['census'])}"
    ax.text(0, 1, "Lump census (responsibility > 0.5)\n" + "\n".join(rows) + note,
            va="top", family="monospace", fontsize=9.5)
    fig.text(0.015, 0.012, "one-law gate q=0.30; law re-estimated per vintage; "
             "guard suite and ledger per EXPERIMENT_NOTES.md", fontsize=7.5, color="#555555")
    fig.savefig(os.path.join(outdir, f"vintage_report_{M['vintage']}.pdf"), bbox_inches="tight")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
