#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/data/perm_gate_impl/DecisionMesh_complete_release_2026-08-01
BIN="$ROOT/build-perm/decision_mesh"
DATA="$ROOT/data/ginnie_design.csv"
OUTROOT=/mnt/data/real_ginnie_perm_q15_multiseed_2026-08-01
mkdir -p "$OUTROOT"
run_one() {
  local mode=$1 seed=$2
  local out="$OUTROOT/${mode}__seed${seed}"
  mkdir -p "$out"
  if [[ -s "$out/run_mesh.csv" && -s "$out/stdout.log" ]]; then return; fi
  export DMESH_DATA="$DATA"
  export DMESH_SPLIT=1 DMESH_SPLIT_MODE=local
  export DMESH_PL_SOLVER=1 DMESH_HIER_FIT=2
  export DMESH_HEIGHT_LIMIT_LOGIT=12 DMESH_IRLS_STEP_CAP=4 DMESH_HIER_MAP=1
  export DMESH_SIGMA0_MIN=0.3 DMESH_SIGMA0_MAX=6 DMESH_PI0_MIN=0.05
  export DMESH_TAU_FLOOR_LOGIT_SD=0.005 DMESH_DUMP="$out/run"
  export DMESH_STAGE_COUNT=3 DMESH_MAX_GATE_ROUNDS=80
  export DMESH_PARENT_VAR_SCALE=0.25
  export DMESH_GATE_EXTRA_VAR=0.35 DMESH_B3_USE_SCORE=1 DMESH_B3_STAGE0_EXACT=0
  unset DMESH_STOP_AFTER_ADAPTATION DMESH_HIST_NULL DMESH_HIST_TOTAL_VAR
  unset DMESH_HIST_MIX_SHAPE DMESH_B3_FIXED_EMPIRICAL_NULL DMESH_CANDIDATE_BH
  unset DMESH_PERM_GATE DMESH_PERM_BH DMESH_PERM_DRAWS DMESH_PERM_MIN_KEFF
  unset DMESH_MIN_CURRENT_KEFF DMESH_PERM_ALLOW_PANEL
  local q=0.10
  if [[ "$mode" == "perm_bh_k4_q15" ]]; then
    q=0.15
    export DMESH_PERM_GATE=1 DMESH_PERM_BH=1 DMESH_PERM_DRAWS=32 DMESH_PERM_MIN_KEFF=4
  fi
  /usr/bin/time -f '%e' -o "$out/wall_seconds.txt" \
    "$BIN" 0 "$q" 1 24 "$seed" >"$out/stdout.log" 2>"$out/stderr.log"
}
export -f run_one
export ROOT BIN DATA OUTROOT
jobs=$(mktemp)
for seed in $(seq 1 10); do
  echo "baseline_q10 $seed" >> "$jobs"
  echo "perm_bh_k4_q15 $seed" >> "$jobs"
done
cat "$jobs" | xargs -n2 -P6 bash -lc 'run_one "$0" "$1"'
rm -f "$jobs"
