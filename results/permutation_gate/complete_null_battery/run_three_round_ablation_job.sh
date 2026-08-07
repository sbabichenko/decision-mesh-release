#!/usr/bin/env bash
set -euo pipefail
BIN=$1; ANALYZE=$2; RUNS=$3; condition=$4; seed=$5; mode=$6; variance=$7; data=$8
OUT="$RUNS/${condition}__${mode}__seed${seed}"; mkdir -p "$OUT"
export DMESH_DATA="$data" DMESH_SPLIT=1 DMESH_PL_SOLVER=1 DMESH_HIER_FIT=2
export DMESH_HEIGHT_LIMIT_LOGIT=12 DMESH_IRLS_STEP_CAP=4 DMESH_HIER_MAP=1
export DMESH_SIGMA0_MIN=0.3 DMESH_SIGMA0_MAX=6 DMESH_PI0_MIN=0.05
export DMESH_TAU_FLOOR_LOGIT_SD=0.005 DMESH_DUMP="$OUT/run"
export DMESH_STAGE_COUNT=3 DMESH_MAX_GATE_ROUNDS=3 DMESH_STOP_AFTER_ADAPTATION=1
export DMESH_PARENT_VAR_SCALE=0.25 DMESH_GATE_EXTRA_VAR="$variance" DMESH_B3_USE_SCORE=1 DMESH_B3_STAGE0_EXACT=0
unset DMESH_HIST_NULL DMESH_HIST_TOTAL_VAR DMESH_HIST_MIX_SHAPE DMESH_B3_FIXED_EMPIRICAL_NULL DMESH_CANDIDATE_BH
unset DMESH_PERM_GATE DMESH_PERM_BH DMESH_PERM_DRAWS DMESH_PERM_MIN_KEFF DMESH_MIN_CURRENT_KEFF
if [[ "$mode" == baseline_k4 ]]; then
 export DMESH_MIN_CURRENT_KEFF=4
elif [[ "$mode" == perm_bh_k05 ]]; then
 export DMESH_PERM_GATE=1 DMESH_PERM_BH=1 DMESH_PERM_DRAWS=32 DMESH_PERM_MIN_KEFF=0.5
fi
/usr/bin/time -f '%e' -o "$OUT/time.txt" "$BIN" 0 0.10 1 24 7 >"$OUT/stdout.log" 2>"$OUT/stderr.log"
python3 "$ANALYZE" "$OUT" --condition "$condition" --seed "$seed" --mode "$mode" --stages 3 --rounds 3 --parent-scale 0.25 > "$OUT/metrics_stdout.json"
rm -f "$OUT/run_obs.csv" "$OUT/run_pools.csv" "$OUT/run_surface.csv"
