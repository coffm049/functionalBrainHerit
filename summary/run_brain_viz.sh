#!/usr/bin/env bash
# Simple runner for the two per-atlas brain visualizations.
#   cd <project-root> && bash summary/run_brain_viz.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV="$ROOT/.venv"

source "$VENV/bin/activate"
export MPLBACKEND=Agg

# Gordon (352 parcels, network dictionary via shortDicionary.csv)
python "$ROOT/summary/brain_gordon.py"

# ProbaConns (80 parcels, networks from dlabel label table)
python "$ROOT/summary/brain_probaconns.py"

echo "Done. Outputs in $ROOT/results/summary/brain/"
