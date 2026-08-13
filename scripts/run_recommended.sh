#!/usr/bin/env bash
# The writeup's "Recommended configuration" (DecisionMesh writeup, Aug 5
# 2026 revision, front matter): exact solver, one-law gate with the
# converged per-stratum law, sigma0 >= 1 central matching. Binomial-
# information weights throughout; no law-tempered fitting weights
# (DMESH_LAW_WEIGHTS is convicted of double-charging and retired).
#
# Target dataset: the June 2026 SMM design (data/smm_design_202606.csv,
# 33,740 pools, n0 >= 25, pooled rate 0.852%). Measured properties at
# this configuration (ten seeds, locally balanced splits): held-out
# deviance 1.151, held-out mixture NLL 0.9242, realized null size
# 0.2-0.3 false admissions per run at q = .10.
#
# Structural environment follows scripts/run_ginnie_benchmark.sh with
# the writeup's recommended overrides: DMESH_ONELAW replaces the scalar
# DMESH_GATE_EXTRA_VAR pricing, and DMESH_SIGMA0_MIN is 1.0 (the
# documented central-matching constraint), not the benchmark script's
# drifted 0.3.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-"$ROOT/build-release/decision_mesh"}
SEED=${2:-7}
OUT=${3:-"$ROOT/output/recommended_seed${SEED}"}
DATA=${DMESH_BENCH_DATA:-"$ROOT/data/smm_design_202606.csv"}
mkdir -p "$OUT"

export DMESH_DATA="$DATA" DMESH_SPLIT=1 DMESH_SPLIT_MODE=${DMESH_SPLIT_MODE:-local}
export DMESH_PL_SOLVER=1 DMESH_HIER_FIT=2 DMESH_HIER_MAP=1 DMESH_QUAD_PRIOR=1
export DMESH_HEIGHT_LIMIT_LOGIT=12 DMESH_IRLS_STEP_CAP=4
export DMESH_TAU_FLOOR_LOGIT_SD=0.005
export DMESH_STAGE_COUNT=3 DMESH_MAX_GATE_ROUNDS=80
export DMESH_PARENT_VAR_SCALE=0.25
export DMESH_B3_USE_SCORE=1 DMESH_B3_STAGE0_EXACT=0
# --- the recommended trio ---
export DMESH_ONELAW="0.1428,0.2415,0.0286,0.0038;8.5,8.1,5.4,0.0"
export DMESH_SIGMA0_MIN=1.0
export DMESH_SIGMA0_MAX=6 DMESH_PI0_MIN=0.05
unset DMESH_GATE_EXTRA_VAR DMESH_LAW_WEIGHTS DMESH_LAW_WEIGHTS_FINAL 2>/dev/null || true
unset DMESH_CANDIDATE_BH DMESH_POINT_VAR 2>/dev/null || true
export DMESH_DUMP="$OUT/run"

/usr/bin/time -f '%e' -o "$OUT/wall_seconds.txt" \
  "$BIN" 0 0.10 1 24 "$SEED" >"$OUT/stdout.log" 2>"$OUT/stderr.log"
echo "done: $OUT"
