#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-"$ROOT/build-release/decision_mesh"}
OUT=${2:-"$ROOT/bootstrap-output/local-20seed"}
WORKERS=${3:-4}
python3 "$ROOT/core/bootstrap_local_split.py" \
  "$BIN" "$ROOT/data/ginnie_design.csv" \
  "$ROOT/core/auditing/law_shared.json" "$OUT" "$WORKERS"
