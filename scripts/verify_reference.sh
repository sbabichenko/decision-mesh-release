#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN=${1:-"$ROOT/build-release/decision_mesh"}
TMP=${2:-"$ROOT/verification-output"}
rm -rf "$TMP"; mkdir -p "$TMP"
for MODE in bh lindsey; do
  "$ROOT/scripts/run_ginnie_benchmark.sh" "$BIN" "$MODE" "$TMP/$MODE"
  (cd "$TMP/$MODE" && sha256sum -c "$ROOT/reference/$MODE/SHA256SUMS.txt" --quiet --status) || {
      echo "$MODE reference hash mismatch" >&2
      exit 1
    }
  echo "$MODE reference outputs verified"
done
