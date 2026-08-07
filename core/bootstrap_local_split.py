#!/usr/bin/env python3
"""Run repeated locally balanced Ginnie holdouts.

Usage:
  bootstrap_local_split.py <binary> <data.csv> <law.json> <output_dir> [workers=4]

The default seed battery is 101..120. These are repeated seeded holdouts, not
an i.i.d. row bootstrap: each run keeps whole pools together and regenerates the
local exposure-balanced split inside the C++ loader.
"""
from __future__ import annotations
import concurrent.futures as cf
import json, os, re, subprocess, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit, gammaln

BIN, DATA, LAW, OUT = map(str, sys.argv[1:5])
WORKERS = int(sys.argv[5]) if len(sys.argv) > 5 else 4
SEEDS = list(range(101, 121))
OUTP = Path(OUT)
OUTP.mkdir(parents=True, exist_ok=True)
law = json.load(open(LAW))["law"]
xs, ws = np.polynomial.hermite_e.hermegauss(20)
lw = np.log(ws / ws.sum())
d = pd.read_csv(DATA)
n = d.n0.values.astype(float)
k = d.k_term.values.astype(float)
wala, wac = d.wala.values.astype(float), d.wac.values.astype(float)
const = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def loglik(kk, nn, eta):
    p = np.clip(expit(eta), 1e-14, 1 - 1e-14)
    return kk * np.log(p) + (nn - kk) * np.log1p(-p)


def comp(kk, nn, gg, s2):
    uh = np.zeros(len(kk))
    for _ in range(60):
        p = expit(gg + uh)
        uh = np.clip(
            uh - np.clip((nn * p - kk + uh / s2) / (nn * p * (1 - p) + 1 / s2), -2, 2),
            -12,
            12,
        )
    p = expit(gg + uh)
    sc = 1 / np.sqrt(nn * p * (1 - p) + 1 / s2)
    T = np.empty((20, len(kk)))
    for j in range(20):
        uu = uh + sc * xs[j]
        T[j] = lw[j] + 0.5 * xs[j] ** 2 + loglik(kk, nn, gg + uu) - 0.5 * uu**2 / s2
    m = T.max(0)
    return m + np.log(np.exp(np.clip(T - m, -700, 0)).sum(0)) + np.log(sc) - 0.5 * np.log(s2)


def evaluate(prefix: Path):
    o = pd.read_csv(str(prefix) + "_obs.csv").sort_values("pid").reset_index(drop=True)
    held = ~o.is_train.astype(bool).values
    ph = np.clip(o.p_hat.values, 1e-12, 1 - 1e-12)
    g = np.log(ph / (1 - ph))
    out = np.full(len(d), np.nan)
    for (lo, hi), (pi, s0, sL, _tot) in zip(
        [(25, 128), (128, 512), (512, 4096), (4096, 10**9)], law
    ):
        te = held & (n >= lo) & (n < hi)
        a = np.log(pi) + comp(k[te], n[te], g[te], sL)
        b = np.log(1 - pi) + comp(k[te], n[te], g[te], s0)
        mx = np.maximum(a, b)
        out[te] = -(mx + np.log(np.exp(a - mx) + np.exp(b - mx)) + const[te])
    masks = {
        "nll_all": held,
        "nll_ramp": held & (wala < 24) & (wac >= 5.5) & (wac < 7.0),
        "nll_mid": held & (wala >= 60) & (wala < 200),
        "nll_seasoned": held & (wala >= 200) & (wala < 300) & (wac >= 5.5) & (wac < 7.0),
    }
    res = {q: float(np.nanmean(out[m])) for q, m in masks.items()}
    low = held & (wala < 12) & (wac >= 5.5) & (wac < 6.25)
    info = n * ph * (1 - ph)
    z = (k - n * ph) / np.sqrt(np.maximum(info, 1e-12))
    res["low_wala_mean_z"] = (
        float(np.average(z[low], weights=np.sqrt(np.maximum(info[low], 1e-12))))
        if low.any()
        else np.nan
    )
    res["low_wala_max_abs_z"] = float(np.max(np.abs(z[low]))) if low.any() else np.nan
    res["heldout_rows"] = int(held.sum())
    return res


def env_for(prefix: Path):
    e = dict(os.environ)
    e.update(
        {
            "DMESH_DATA": os.path.abspath(DATA),
            "DMESH_SPLIT": "1",
            "DMESH_SPLIT_MODE": "local",
            "DMESH_SPLIT_BLOCK_SIZE": "24",
            "DMESH_PL_SOLVER": "1",
            "DMESH_HIER_FIT": "2",
            "DMESH_HEIGHT_LIMIT_LOGIT": "12",
            "DMESH_IRLS_STEP_CAP": "4",
            "DMESH_HIER_MAP": "1",
            "DMESH_SIGMA0_MIN": "0.3",
            "DMESH_SIGMA0_MAX": "6",
            "DMESH_PI0_MIN": "0.05",
            "DMESH_TAU_FLOOR_LOGIT_SD": "0.005",
            "DMESH_DUMP": str(prefix),
            "DMESH_STAGE_COUNT": "3",
            "DMESH_MAX_GATE_ROUNDS": "80",
            "DMESH_PARENT_VAR_SCALE": "0.25",
            "DMESH_GATE_EXTRA_VAR": "0.35",
            "DMESH_B3_USE_SCORE": "1",
            "DMESH_B3_STAGE0_EXACT": "0",
            "DMESH_CANDIDATE_BH": "1",
        }
    )
    for key in [
        "DMESH_STOP_AFTER_ADAPTATION",
        "DMESH_HIST_NULL",
        "DMESH_HIST_TOTAL_VAR",
        "DMESH_HIST_MIX_SHAPE",
        "DMESH_B3_FIXED_EMPIRICAL_NULL",
        "DMESH_ONELAW",
        "DMESH_POINT_VAR",
        "DMESH_LAW_WEIGHTS",
    ]:
        e.pop(key, None)
    return e


def run(seed):
    rd = OUTP / f"seed_{seed}"
    rd.mkdir(exist_ok=True)
    prefix = rd / "run"
    t = time.time()
    if Path(str(prefix) + "_obs.csv").exists() and (rd / "stdout.log").exists() and (rd / "stderr.log").exists():
        class Completed:
            returncode = 0
            stdout = (rd / "stdout.log").read_text()
            stderr = (rd / "stderr.log").read_text()
        cp = Completed()
    else:
        cp = subprocess.run(
            [BIN, "0", "0.10", "1", "24", str(seed)],
            env=env_for(prefix),
            capture_output=True,
            text=True,
        )
        (rd / "stdout.log").write_text(cp.stdout)
        (rd / "stderr.log").write_text(cp.stderr)
        if cp.returncode:
            return {"seed": seed, "error": f"exit {cp.returncode}", "wall_seconds": time.time() - t}
    r = evaluate(prefix)
    r.update(seed=seed, evaluation_wall_seconds=time.time() - t)
    m = pd.read_csv(str(prefix) + "_mesh.csv")
    hv = pd.read_csv(str(prefix) + "_hier_vertices.csv")
    r["faces"] = len(m)
    r["active_vertices"] = int((hv.active == 1).sum())
    r["gate_admitted_vertices"] = int(((hv.active == 1) & (hv.gate_admitted == 1)).sum())
    tm = re.search(r"TOTAL ([0-9.]+) s", cp.stdout)
    if tm:
        r["internal_runtime_seconds"] = float(tm.group(1))
    mm = re.search(r"HELDOUT: mean deviance/pool ([0-9.]+) \| X2/info ([0-9.]+)", cp.stdout)
    if mm:
        r["heldout_deviance_per_pool"] = float(mm.group(1))
        r["x2_per_information"] = float(mm.group(2))
    sm = re.search(
        r"\[split\] mode local: train (\d+) rows / ([0-9.]+) exposure / rate ([0-9.]+); held (\d+) rows / ([0-9.]+) exposure / rate ([0-9.]+)",
        cp.stderr,
    )
    if sm:
        r.update(
            train_rows=int(sm.group(1)),
            train_exposure=float(sm.group(2)),
            train_rate=float(sm.group(3)),
            held_exposure=float(sm.group(5)),
            held_rate=float(sm.group(6)),
        )
    return r


rows = []
with cf.ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = {ex.submit(run, s): s for s in SEEDS}
    for fut in cf.as_completed(futs):
        r = fut.result()
        rows.append(r)
        if "error" in r:
            print(f"seed {r['seed']}: {r['error']}", flush=True)
        else:
            print(
                f"seed {r['seed']}: NLL {r['nll_all']:.4f}, ramp {r['nll_ramp']:.3f}, faces {r['faces']}",
                flush=True,
            )
R = pd.DataFrame(rows).sort_values("seed")
R.to_csv(OUTP / "bootstrap_20seed_metrics.csv", index=False)
num = R.select_dtypes(include=[np.number]).drop(columns=["seed"], errors="ignore")
S = pd.DataFrame(
    {"mean": num.mean(), "sd": num.std(ddof=1), "min": num.min(), "median": num.median(), "max": num.max()}
)
S.to_csv(OUTP / "bootstrap_20seed_summary.csv")
q = R.nll_all.quantile([0.025, 0.5, 0.975])
report = f"""# Twenty-seed locally balanced Ginnie holdouts

Seeds: 101--120. Whole pools are kept intact. Each split uses recursive spatial blocks of at most 24 pools and exposure-matched pairing within each block.

- Marginal NLL: mean **{R.nll_all.mean():.4f}**, SD **{R.nll_all.std(ddof=1):.4f}**, range **[{R.nll_all.min():.4f}, {R.nll_all.max():.4f}]**.
- Descriptive 2.5/50/97.5 percentiles: **{q.loc[0.025]:.4f} / {q.loc[0.5]:.4f} / {q.loc[0.975]:.4f}**.
- Ramp NLL: mean **{R.nll_ramp.mean():.4f}**, SD **{R.nll_ramp.std(ddof=1):.4f}**.
- Faces: mean **{R.faces.mean():.1f}**, range **[{R.faces.min()}, {R.faces.max()}]**.
- Active vertices: mean **{R.active_vertices.mean():.1f}**, range **[{R.active_vertices.min()}, {R.active_vertices.max()}]**.
- Low-WALA/WAC 5.5--6.25 weighted mean Pearson z: mean **{R.low_wala_mean_z.mean():+.3f}**, range **[{R.low_wala_mean_z.min():+.3f}, {R.low_wala_mean_z.max():+.3f}]**.
- Topology-collapse runs with fewer than 500 faces: **{int((R.faces < 500).sum())}/20**. They remain in the primary summary.

The noncollapsed main cluster (faces >= 500) has mean NLL **{R.loc[R.faces >= 500, 'nll_all'].mean():.4f}** and SD **{R.loc[R.faces >= 500, 'nll_all'].std(ddof=1):.4f}**. These are repeated seeded holdouts, not an iid nonparametric bootstrap; the spread measures split sensitivity.
"""
(OUTP / "BOOTSTRAP_REPORT.md").write_text(report)
print(report)
