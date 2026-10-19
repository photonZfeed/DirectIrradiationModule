#!/usr/bin/env bash
# Background pipeline for the Monte-Carlo convergence study (idempotent: finished stages are skipped).
# Seed 42 reproduces the published 5 million-sample set in its first 5 million distinct configurations.
# Usage: PY=<python> [SEEDS="42 1 2 3"] bash revision_analyses/run_all_seeds.sh <out_dir>
# (DIM_REPO is optional; by default the repository containing this folder is used)
# A sampler needs about 2 GB of memory at the end of a run; do not run more than two at a time on 8 GB.
set -euo pipefail
OUT="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
export MPLBACKEND=Agg
for s in ${SEEDS:-42 1 2 3}; do
  "$PY" -W error "$HERE/r3_mc_runs.py" sample --out "$OUT" --seed "$s" --n 10000000
  "$PY" -W error "$HERE/r3_mc_runs.py" evaluate --out "$OUT" --configs "$OUT/mc_seed${s}_10000000_configs.npz"
done
echo ALL_DONE
