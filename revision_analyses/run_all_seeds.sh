#!/usr/bin/env bash
# Background pipeline for the Monte-Carlo convergence study (idempotent: finished stages are skipped).
# Usage: PY=<python> bash revision_analyses/run_all_seeds.sh <out_dir> <published_configs.npz>
# (DIM_REPO is optional; by default the repository containing this folder is used)
set -euo pipefail
OUT="$1"; PUB="$2"
HERE="$(cd "$(dirname "$0")" && pwd)"
export MPLBACKEND=Agg
"$PY" -W error "$HERE/r3_mc_runs.py" evaluate --out "$OUT" --configs "$PUB"
for s in ${SEEDS:-1 2 3 4 5}; do
  "$PY" -W error "$HERE/r3_mc_runs.py" sample --out "$OUT" --seed "$s" --n 10000000
  "$PY" -W error "$HERE/r3_mc_runs.py" evaluate --out "$OUT" --configs "$OUT/mc_seed${s}_10000000_configs.npz"
done
echo ALL_DONE
