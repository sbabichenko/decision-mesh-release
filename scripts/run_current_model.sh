#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-"$ROOT/build-release/decision_mesh"}
SEED=${2:-7}
OUT=${3:-"$ROOT/output/current_seed${SEED}"}
mkdir -p "$OUT"

export DMESH_PL_SOLVER=1
export DMESH_SPLIT=1
export DMESH_SPLIT_MODE=${DMESH_SPLIT_MODE:-local}
export DMESH_HIER_FIT=2
export DMESH_HIER_MAP=1
export DMESH_QUAD_PRIOR=1
export DMESH_DUMP="$OUT/run"

"$BIN" 0 0.10 1 24 "$SEED" >"$OUT/stdout.log" 2>"$OUT/stderr.log"
