#!/usr/bin/env bash
# Regenerate the per-seed vertex dumps used by robust.py, from the shipped
# Ginnie design under the current recommended configuration.  Outputs go to
# $OUTDIR (default /tmp/reshrink_repro), matching robust.py's default paths.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
BIN=${1:-"$ROOT/build-release/decision_mesh"}
OUTDIR=${2:-/tmp/reshrink_repro}
mkdir -p "$OUTDIR"
for s in 11 13 21; do
  DMESH_PL_SOLVER=1 DMESH_SPLIT=1 DMESH_SPLIT_MODE=local DMESH_HIER_FIT=2 \
  DMESH_HIER_MAP=1 DMESH_QUAD_PRIOR=1 DMESH_DUMP="$OUTDIR/seed${s}" \
    "$BIN" 0 0.10 1 24 "$s" >/dev/null 2>&1
  echo "seed $s -> $OUTDIR/seed${s}_hier_vertices.csv"
done
