#!/usr/bin/env python3
"""Analyze Phase-B2 exact signed-root gain dumps.

The script treats the exact score as a ranking coordinate, not as an assumed
N(0,1) pivot.  It checks objective identities and summarizes where it differs
from the legacy working-scale score.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def safe_spearman(a: pd.Series, b: pd.Series) -> float:
    if len(a) < 3 or a.nunique() < 2 or b.nunique() < 2:
        return float("nan")
    return float(spearmanr(a, b).statistic)


def top_overlap(g: pd.DataFrame, frac: float = 0.10) -> float:
    n = max(1, int(np.ceil(frac * len(g))))
    a = set(g.nlargest(n, "legacy_abs_rank_input").index)
    b = set(g.nlargest(n, "exact_abs_rank_input").index)
    return len(a & b) / n


def analyze(path: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    d = pd.read_csv(path)
    root_err = np.abs(d.signed_root_gain**2 - 2.0 * d.gain_nats)
    decomp_err = np.abs(d.data_gain_nats + d.penalty_gain_nats - d.gain_nats)
    nonzero = d.exact_displacement_units.abs() > 1e-8
    sign_mismatch = (
        np.sign(d.loc[nonzero, "signed_root_gain"])
        != np.sign(d.loc[nonzero, "exact_displacement_units"])
    )
    sign_disagree = np.sign(d.legacy_z) != np.sign(d.signed_root_gain)

    groups = []
    for (stage, rnd), g in d.groupby(["stage", "round"], sort=True):
        groups.append({
            "stage": int(stage),
            "round": int(rnd),
            "n": int(len(g)),
            "legacy_mean": float(g.legacy_z.mean()),
            "legacy_sd": float(g.legacy_z.std(ddof=1)),
            "exact_mean": float(g.signed_root_gain.mean()),
            "exact_sd": float(g.signed_root_gain.std(ddof=1)),
            "abs_spearman": safe_spearman(g.legacy_abs_rank_input, g.exact_abs_rank_input),
            "top10_overlap": top_overlap(g),
            "sign_disagree_fraction": float((np.sign(g.legacy_z) != np.sign(g.signed_root_gain)).mean()),
            "median_keff": float(g.keff.median()),
            "median_min_p": float(g.min_p_null.median()),
        })
    by_round = pd.DataFrame(groups)

    bins = pd.cut(
        d.min_p_null,
        [-np.inf, 1e-6, 1e-4, 1e-2, 0.1, np.inf],
        labels=["<=1e-6", "1e-6..1e-4", "1e-4..1e-2", "1e-2..0.1", ">0.1"],
    )
    regime = (
        d.assign(p_regime=bins)
        .groupby("p_regime", observed=True)
        .agg(
            n=("id", "size"),
            legacy_abs_median=("legacy_abs_rank_input", "median"),
            exact_abs_median=("exact_abs_rank_input", "median"),
            exact_abs_p95=("exact_abs_rank_input", lambda x: x.quantile(0.95)),
            sign_disagree_fraction=("legacy_z", lambda x: float((np.sign(x) != np.sign(d.loc[x.index, "signed_root_gain"])).mean())),
            median_keff=("keff", "median"),
        )
        .reset_index()
    )

    penalty_fraction = np.divide(
        d.penalty_gain_nats,
        d.gain_nats,
        out=np.zeros(len(d), dtype=float),
        where=d.gain_nats.to_numpy() > 1e-12,
    )
    nearest = d.iloc[np.argmin(np.hypot(d.x.to_numpy(), d.y.to_numpy() - 0.3750871809178406))]
    out = {
        "file": str(path),
        "n_candidates": int(len(d)),
        "identity": {
            "max_signed_root_error": float(root_err.max()),
            "max_gain_decomposition_error": float(decomp_err.max()),
            "negative_gain_count": int((d.gain_nats < -1e-10).sum()),
            "signed_direction_mismatch_count": int(sign_mismatch.sum()),
        },
        "comparison": {
            "signed_correlation": float(d[["legacy_z", "signed_root_gain"]].corr().iloc[0, 1]),
            "absolute_spearman": safe_spearman(d.legacy_abs_rank_input, d.exact_abs_rank_input),
            "sign_disagreement_count": int(sign_disagree.sum()),
            "sign_disagreement_fraction": float(sign_disagree.mean()),
            "median_round_top10_overlap": float(by_round.top10_overlap.median()),
            "median_round_abs_spearman": float(by_round.abs_spearman.median()),
        },
        "exact_score_quantiles": {
            str(q): float(d.signed_root_gain.quantile(q))
            for q in [0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999]
        },
        "penalty": {
            "median_penalty_share": float(np.median(penalty_fraction)),
            "p05_penalty_share": float(np.quantile(penalty_fraction, 0.05)),
            "p95_penalty_share": float(np.quantile(penalty_fraction, 0.95)),
            "data_gain_negative_count": int((d.data_gain_nats < -1e-10).sum()),
            "penalty_only_direction_count": int(((d.data_gain_nats < 0) & (d.gain_nats > 1e-12)).sum()),
        },
        "regime_straddle_nearest": {
            k: (int(nearest[k]) if k in {"stage", "round", "id", "depth"} else float(nearest[k]))
            for k in ["stage", "round", "id", "x", "y", "depth", "legacy_z",
                      "signed_root_gain", "null_units", "prior_units", "optimum_units",
                      "gain_nats", "data_gain_nats", "penalty_gain_nats", "min_p_null",
                      "max_p_null", "keff"]
        },
    }
    return out, by_round, regime


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("scores", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary, by_round, regime = analyze(args.scores)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    by_round.to_csv(args.out / "by_round.csv", index=False)
    regime.to_csv(args.out / "by_probability_regime.csv", index=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
