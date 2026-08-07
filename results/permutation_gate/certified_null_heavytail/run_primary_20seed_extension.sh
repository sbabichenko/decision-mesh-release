#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/data/perm_gate_impl/DecisionMesh_complete_release_2026-08-01
BIN=$ROOT/build-perm/decision_mesh
DATA_ROOT=/mnt/data/heavytail_certified_null_inputs_2026-08-01
OUT_ROOT=/mnt/data/perm_gate_heavytail_certified_null_2026-08-01
run_one() {
  local seed=$1 mode=$2 q=$3 k=$4
  local data="$DATA_ROOT/hcn_heavytail_seed${seed}.csv" out="$OUT_ROOT/${mode}__seed${seed}"
  mkdir -p "$out"; [[ -f "$out/run_mesh.csv" ]] && return
  export DMESH_DATA="$data" DMESH_SPLIT=1 DMESH_SPLIT_MODE=local DMESH_PL_SOLVER=1 DMESH_HIER_FIT=2 DMESH_HEIGHT_LIMIT_LOGIT=12 DMESH_IRLS_STEP_CAP=4 DMESH_HIER_MAP=1 DMESH_SIGMA0_MIN=0.3 DMESH_SIGMA0_MAX=6 DMESH_PI0_MIN=0.05 DMESH_TAU_FLOOR_LOGIT_SD=0.005 DMESH_DUMP="$out/run" DMESH_STAGE_COUNT=3 DMESH_MAX_GATE_ROUNDS=80 DMESH_PARENT_VAR_SCALE=0.25 DMESH_GATE_EXTRA_VAR=0.35 DMESH_B3_USE_SCORE=1 DMESH_B3_STAGE0_EXACT=0
  unset DMESH_STOP_AFTER_ADAPTATION DMESH_HIST_NULL DMESH_HIST_TOTAL_VAR DMESH_HIST_MIX_SHAPE DMESH_B3_FIXED_EMPIRICAL_NULL DMESH_CANDIDATE_BH DMESH_PERM_GATE DMESH_PERM_BH DMESH_PERM_DRAWS DMESH_PERM_MIN_KEFF DMESH_MIN_CURRENT_KEFF
  if [[ "$mode" == perm* ]]; then export DMESH_PERM_GATE=1 DMESH_PERM_BH=1 DMESH_PERM_DRAWS=32 DMESH_PERM_MIN_KEFF="$k"; fi
  /usr/bin/time -f '%e' -o "$out/time.txt" "$BIN" 0 "$q" 1 24 7 >"$out/stdout.log" 2>"$out/stderr.log"
}
export -f run_one; export ROOT BIN DATA_ROOT OUT_ROOT
jobs=$(mktemp)
for seed in $(seq 5 19); do
 echo "$seed baseline_q10 0.10 0.5" >>$jobs
 echo "$seed perm_bh_k05_q10 0.10 0.5" >>$jobs
 echo "$seed perm_bh_k4_q10 0.10 4" >>$jobs
 echo "$seed perm_bh_k4_q20 0.20 4" >>$jobs
done
cat $jobs | xargs -n4 -P8 bash -lc 'run_one "$@"' _
rm -f $jobs
