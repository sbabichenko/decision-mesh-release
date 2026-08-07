#!/usr/bin/env bash
set -euo pipefail
BIN=${DMESH_BIN:-./build/decision_mesh}
DATA=${1:?usage: $0 null.csv output-dir}
OUT=${2:?usage: $0 null.csv output-dir}
mkdir -p "$OUT"
export DMESH_DATA="$DATA"
export DMESH_SPLIT=1
export DMESH_PL_SOLVER=1
export DMESH_HIER_FIT=2
export DMESH_HEIGHT_LIMIT_LOGIT=12
export DMESH_IRLS_STEP_CAP=4
export DMESH_HIER_MAP=1
export DMESH_SIGMA0_MIN=0.3
export DMESH_SIGMA0_MAX=6
export DMESH_PI0_MIN=0.05
export DMESH_TAU_FLOOR_LOGIT_SD=0.005
export DMESH_HIST_NULL=1
export DMESH_HIST_TOP=8
export DMESH_HIST_TOTAL_VAR='0.7052106170,0.4648783763,0.4790827059,0.3806532231'
export DMESH_HIST_MIX_SHAPE='0.2710449345,5.85195;0.1583151119,5.923;0.0688630065,7.7869;0.1039544004,91.617'
export DMESH_DUMP="$OUT/run"
export DMESH_B2_SCORE_DUMP="$OUT/scores.csv"
"$BIN" 0 0.10 1 24 7 >"$OUT/stdout.log" 2>"$OUT/stderr.log"
