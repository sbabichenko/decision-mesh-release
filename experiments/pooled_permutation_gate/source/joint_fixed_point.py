#!/usr/bin/env python3
"""Joint fixed point: surface and pool-effect law estimated to mutual
consistency, one law object feeding gate, working weights, posteriors,
and metric alike.

Usage: python3 joint_fixed_point.py <pools.csv> <outdir> [iters]

Per iteration: fit (gate + weights consume the current law) -> mixture EM
on train-half residuals -> new law + per-pool gamma-weighted variances +
posteriors -> held-out mixture NLL evaluated BY REFERENCE with the same
law. Convergence tracked on law parameters, mesh size, held-out metrics.
"""
import sys, os, json, subprocess
import numpy as np, pandas as pd
from scipy.special import expit, gammaln

BIN = os.environ.get("DMESH_BIN", "/home/claude/decision_mesh_latest/build/decision_mesh")
CORE = {"DMESH_SIGMA0_MIN":"0.3","DMESH_SIGMA0_MAX":"6","DMESH_PI0_MIN":"0.05",
        "DMESH_TAU_FLOOR_LOGIT_SD":"0.005","DMESH_SPLIT":"1","DMESH_FINAL_RESHRINK":"1"}
STRATA = [(25,128),(128,512),(512,4096),(4096,10**9)]
XS, WS = np.polynomial.hermite_e.hermegauss(20); LW = np.log(WS/WS.sum())

def loglik(k,n,eta):
    p = np.clip(expit(eta),1e-14,1-1e-14); return k*np.log(p)+(n-k)*np.log1p(-p)

def component(k,n,g,s2):
    uh = np.zeros(len(k))
    for _ in range(60):
        p = expit(g+uh); step = (n*p-k+uh/s2)/(n*p*(1-p)+1/s2)
        uh = np.clip(uh-np.clip(step,-2,2),-12,12)
    p = expit(g+uh); c = n*p*(1-p)+1/s2; s = 1/np.sqrt(c)
    T = np.empty((20,len(k))); U = np.empty((20,len(k)))
    for j in range(20):
        uu = uh+s*XS[j]; T[j] = LW[j]+0.5*XS[j]**2+loglik(k,n,g+uu)-0.5*uu**2/s2; U[j]=uu
    m = T.max(0); E = np.exp(np.clip(T-m,-700,0)); Z = E.sum(0)
    return m+np.log(Z)+np.log(s)-0.5*np.log(s2), (E*U).sum(0)/Z, (E*U*U).sum(0)/Z

def em(k,n,g, iters=20):
    pi,s0,sL = 0.05,0.10,2.0
    for _ in range(iters):
        lm0,_,E0 = component(k,n,g,s0); lmL,_,EL = component(k,n,g,sL)
        a = np.log(pi)+lmL; b = np.log(1-pi)+lm0; mx = np.maximum(a,b)
        gam = np.exp(a-mx)/(np.exp(a-mx)+np.exp(b-mx))
        pi = float(np.clip(gam.mean(),1e-5,0.5)); s0 = float(np.sum((1-gam)*E0)/np.sum(1-gam))
        sL = float(max(np.sum(gam*EL)/max(np.sum(gam),1e-9), s0))
    return pi,s0,sL

def fit(csv, outdir, tag, law_str, point_var, law_weights):
    env = dict(os.environ); env.update(CORE)
    env["DMESH_DATA"] = csv; env["DMESH_DUMP"] = os.path.join(outdir, tag)
    if law_str: env["DMESH_ONELAW"] = law_str
    else: env["DMESH_GATE_EXTRA_VAR"] = "0.35"
    if point_var: env["DMESH_POINT_VAR"] = point_var
    if law_weights: env["DMESH_LAW_WEIGHTS"] = "1"
    out = subprocess.run([BIN,"0","0.30","1","24","7"], env=env, capture_output=True, text=True).stdout
    dev = float([l for l in out.split("\n") if "deviance/pool" in l][-1].split("deviance/pool")[-1].split()[0].rstrip(","))
    vtx = int([l for l in out.split("\n") if "vertices" in l][-1].split("vertices")[-1].split()[0].rstrip(","))
    p = pd.read_csv(os.path.join(outdir, tag+"_obs.csv")).p_hat.values
    return dev, vtx, np.log(np.clip(p,1e-9,1-1e-9)/np.clip(1-p,1e-9,1-1e-9))

def main(csv, outdir, iters=4):
    os.makedirs(outdir, exist_ok=True)
    d = pd.read_csv(csv); n = d.n0.values.astype(float); k = d.k_term.values.astype(float)
    even = np.arange(len(d))%2==0
    const = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
    law_str, pv_path, prev_law = None, None, None
    DAMP = 0.5  # geometric damping on law updates; per-pool variances feed the
                # gate only (working weights use stratum totals) because the
                # gamma->weight->residual loop is self-exciting (measured:
                # undamped oscillation, it-2 blowup)
    print(f"{'it':>3} {'vertices':>9} {'heldout dev':>12} {'mix NLL (ref)':>14}  law totals per stratum")
    for it in range(iters+1):
        tag = f"it{it}"
        dev, vtx, g = fit(csv, outdir, tag, law_str, pv_path, law_weights=(it>0))
        law, pointvar, nll, cnt = [], np.zeros(len(d)), 0.0, 0
        post = {"gamma": np.zeros(len(d)), "u_mix": np.zeros(len(d))}
        for lo,hi in STRATA:
            m = (n>=lo)&(n<hi); tr, te = m&even, m&~even
            pi,s0,sL = em(k[tr], n[tr], g[tr])
            if prev_law is not None:
                po,so,slo,_ = prev_law[len(law)]
                pi = DAMP*pi + (1-DAMP)*po; s0 = DAMP*s0 + (1-DAMP)*so; sL = DAMP*sL + (1-DAMP)*slo
            law.append((pi,s0,sL,(1-pi)*s0+pi*sL))
            # posteriors and gamma-weighted variances for ALL pools in stratum
            lm0,u0,_ = component(k[m],n[m],g[m],s0); lmL,uL,_ = component(k[m],n[m],g[m],sL)
            a = np.log(pi)+lmL; b = np.log(1-pi)+lm0; mx = np.maximum(a,b)
            gam = np.exp(a-mx)/(np.exp(a-mx)+np.exp(b-mx))
            tot = (1-pi)*s0 + pi*sL
            pointvar[m] = 0.5*(gam*sL + (1-gam)*s0) + 0.5*tot  # damped per-pool
            post["gamma"][m] = gam; post["u_mix"][m] = gam*uL + (1-gam)*u0
            # metric BY REFERENCE: same law, held-out half
            teh = te[m]
            nll -= ((mx + np.log(np.exp(a-mx)+np.exp(b-mx)))[teh] + const[te & m]).sum()
            cnt += int(teh.sum())
        kurts = [3*((1-p_)*s0_**2+p_*sL_**2)/t_**2-3 for p_,s0_,sL_,t_ in law]
        prev_law = law
        law_str = ",".join(f"{t[3]:.4f}" for t in law) + ";" + ",".join(f"{ku:.2f}" for ku in kurts)
        np.savetxt(os.path.join(outdir, f"pointvar_it{it}.txt"), pointvar, fmt="%.6f")
        pv_path = None  # per-pool weight feedback is self-exciting (measured twice);
                        # all consumers use the stratum law until the C++ separates
                        # gate and weight arrays (POINT_VAR hook reserved for gate-only)
        totals = "/".join(f"{t[3]:.3f}" for t in law)
        print(f"{it:>3} {vtx:>9} {dev:>12.3f} {nll/cnt:>14.4f}  {totals}")
        json.dump({"law": law, "kurts": kurts, "iteration": it,
                   "vertices": vtx, "heldout_dev": dev, "mix_nll": round(nll/cnt,4)},
                  open(os.path.join(outdir, "law_shared.json"), "w"), indent=1)
        pd.DataFrame({"gamma": post["gamma"], "u_mix": post["u_mix"],
                      "point_var": pointvar}).to_csv(os.path.join(outdir,"posteriors.csv"), index=False)
    print("shared law object: law_shared.json | posteriors: posteriors.csv")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv)>3 else 4)
