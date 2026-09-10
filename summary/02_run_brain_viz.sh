#!/usr/bin/env bash
# Runner for all brain visualizations (FC + SA) + FC distance (requires workbench for geodesic).
#   cd <project-root> && bash summary/02_run_brain_viz.sh
# Order: 01_compare_mash_twin.R must have been run first to create mash_twin_wide.csv
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV="$ROOT/.venv"

# Make wb_command available for geodesic distance (04f) if workbench module exists
if command -v module >/dev/null 2>&1; then
  module load workbench 2>/dev/null || module load connectome-workbench 2>/dev/null || true
fi
# Also try common absolute paths
export PATH="/usr/local/workbench/bin:/opt/workbench/bin:$PATH"

source "$VENV/bin/activate"
export MPLBACKEND=Agg

# FC — Gordon (352) and ProbaConns (80)
python "$ROOT/summary/02a_brain_gordon.py"
python "$ROOT/summary/02b_brain_probaconns.py"

# SA — 14 networks (excl. 4/6/17), w/ and w/o total surface
python "$ROOT/summary/02c_brain_sa.py"
python "$ROOT/summary/02d_brain_sa_total.py"

# FC heritability vs distance (Euclidean + geodesic via wb_command, per edge, Twin + AdjHE-RE, simple scatter)
python "$ROOT/summary/04f_FC_distance.py" || echo "04f failed (wb_command may be missing, but Euclidean scatter still produced if possible)"

echo "Done. Outputs in $ROOT/results/summary/brain/ and $ROOT/results/summary/plots/fc_distance*"
