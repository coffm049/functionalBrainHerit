#!/usr/bin/env bash
# Simple runner for the two per-atlas brain visualizations.
#   cd <project-root> && bash summary/08_run_brain_viz.sh
# Order: 01_compare_mash_twin.R must have been run first to create mash_twin_wide.csv
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV="$ROOT/.venv"

source "$VENV/bin/activate"
export MPLBACKEND=Agg

# Gordon (352 parcels, network dictionary via shortDicionary.csv) — 02
python "$ROOT/summary/02_brain_gordon.py"

# ProbaConns (80 parcels, networks from dlabel label table) — 03
python "$ROOT/summary/03_brain_probaconns.py"

echo "Done. Outputs in $ROOT/results/summary/brain/"
