#!/usr/bin/env bash
# Run the brain-space heritability visualizations (twin vs AdjHE) for both
# atlases, inside the project-root Python venv (.venv).
#
#   cd <project-root> && bash summary/run_brain_viz.sh
#   (optional) conda activate gdc   # if R + python live in the 'gdc' env
#
# Before first use, create the venv (once) from the conda-base python:
#   ~/miniconda3/bin/python -m venv .venv
#   .venv/bin/pip install "numpy>=1.24" "pandas>=2.0" "nibabel>=5.0" \
#                          "nilearn>=0.10" "plotly>=5.0"
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WIDE="$ROOT/results/summary/mash_twin_wide.csv"

# Regenerate the summary table if it is missing (needs R + tidyverse on PATH).
if [ ! -f "$WIDE" ]; then
  echo "WIDE summary not found at $WIDE -> running compare_mash_twin.R"
  Rscript "$ROOT/summary/compare_mash_twin.R"
fi

VENV="$ROOT/.venv"

# The project .venv is self-contained (its own python + nilearn), so we do NOT
# `conda activate base` here: under `set -u` that crashed on the inherited gdc
# env's deactivate.d scripts (CONDA_BACKUP_* unbound). Rscript is on PATH.
CONDA_ROOT="$HOME/miniconda3"

if [ ! -f "$VENV/bin/activate" ]; then
  echo "ERROR: venv not found at $VENV" >&2
  echo "Create it with: $CONDA_ROOT/bin/python -m venv $VENV && \\" >&2
  echo "                 $VENV/bin/pip install \"numpy>=1.24\" \"pandas>=2.0\" \\" >&2
  echo "                 \"nibabel>=5.0\" \"nilearn>=0.10\" \"plotly>=5.0\"" >&2
  exit 1
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"
export MPLBACKEND=Agg   # headless HPC: no X display needed for PNG output

# ---- EDIT THESE: paths on the HPC ----------------------------------------
# WIDE is defined above (and auto-generated from compare_mash_twin.R if missing)
DLABEL_GORDON="/users/4/coffm049/papers/brainTemplates/Gordon.networks.32k_fs_LR.dlabel.nii"
DLABEL_PROBA="/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii"
SURF_L="/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii"
SURF_R="/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii"
WB_CMD="wb_command"          # external HCP Workbench binary (for MNI centroids)
CENTROIDS_GORDON=""          # optional Nx3 CSV; leave "" to derive via wb_command
CENTROIDS_PROBA=""
NODE_PCT=90
EDGE_PCT=90
# --------------------------------------------------------------------------

run_atlas () {
  local atlas="$1" n="$2" dlabel="$3" centroids="$4"
  local out="$ROOT/results/summary/brain/$atlas"
  mkdir -p "$out"
  local cflag=()
  if [ -n "$centroids" ]; then cflag=(--centroids "$centroids"); fi
  python "$ROOT/summary/brain_connectome.py" \
    --wide "$WIDE" --atlas "$atlas" --n "$n" \
    --dlabel "$dlabel" --surf-l "$SURF_L" --surf-r "$SURF_R" \
    --wb-command "$WB_CMD" \
    --methods twin AdjHE --node-pct "$NODE_PCT" --edge-pct "$EDGE_PCT" \
    --outdir "$out" "${cflag[@]}"
}

run_atlas gordon 352 "$DLABEL_GORDON" "$CENTROIDS_GORDON"
run_atlas probaConns 80 "$DLABEL_PROBA" "$CENTROIDS_PROBA"

echo "Done. Outputs in $ROOT/results/summary/brain/"
