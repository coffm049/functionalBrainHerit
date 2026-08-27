#!/usr/bin/env bash
# Runner for all brain visualizations (FC + SA).
#   cd <project-root> && bash summary/02_run_brain_viz.sh
# Order: 01_compare_mash_twin.R must have been run first to create mash_twin_wide.csv
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV="$ROOT/.venv"

source "$VENV/bin/activate"
export MPLBACKEND=Agg

# FC — Gordon (352) and ProbaConns (80)
python "$ROOT/summary/02a_brain_gordon.py"
python "$ROOT/summary/02b_brain_probaconns.py"

# SA — 17 networks, w/ and w/o total surface
python "$ROOT/summary/02c_brain_sa.py"
python "$ROOT/summary/02d_brain_sa_total.py"

echo "Done. Outputs in $ROOT/results/summary/brain/"
