#!/usr/bin/env bash
# Run the brain-space heritability visualizations (twin vs AdjHE) for both
# atlases, inside the project-root Python venv (.venv).
#
#   cd <project-root> && bash summary/run_brain_viz.sh
#
# Before first use, create the venv (once) from the conda-base python:
#   ~/miniconda3/bin/python -m venv .venv
#   .venv/bin/pip install "numpy>=1.24" "pandas>=2.0" "nibabel>=5.0" "nilearn>=0.10"
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WIDE="$ROOT/results/summary/mash_twin_wide.csv"

VENV="$ROOT/.venv"
CONDA_ROOT="$HOME/miniconda3"

if [ ! -f "$VENV/bin/activate" ]; then
  echo "ERROR: venv not found at $VENV" >&2
  echo "Create it with: $CONDA_ROOT/bin/python -m venv $VENV && \\" >&2
  echo "                 $VENV/bin/pip install \"numpy>=1.24\" \"pandas>=2.0\" \\" >&2
  echo "                 \"nibabel>=5.0\" \"nilearn>=0.10\"" >&2
  exit 1
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"
export MPLBACKEND=Agg   # headless HPC: no X display needed for PNG output

# ---- EDIT THESE: paths on the HPC ----------------------------------------
DLABEL_GORDON="/users/4/coffm049/papers/brainTemplates/Gordon.networks.32k_fs_LR.dlabel.nii"
DLABEL_PROBA="/projects/standard/faird/shared/data/Probabilitic_network_ROIs_small_package/ABCD/combined_clusters/combined_clusters_thresh0.8.dlabel.nii"
SURF_L="/users/4/coffm049/papers/brainTemplates/Conte69.L.inflated.32k_fs_LR.surf.gii"
SURF_R="/users/4/coffm049/papers/brainTemplates/Conte69.R.inflated.32k_fs_LR.surf.gii"
NODE_PCT=90
# --------------------------------------------------------------------------

run_atlas () {
  local atlas="$1" n="$2" dlabel="$3"
  local out="$ROOT/results/summary/brain/$atlas"
  mkdir -p "$out"
  python "$ROOT/summary/brain_connectome.py" \
    --wide "$WIDE" --atlas "$atlas" --n "$n" \
    --dlabel "$dlabel" --surf-l "$SURF_L" --surf-r "$SURF_R" \
    --node-pct "$NODE_PCT" --write-dscalar \
    --outdir "$out"
}

run_atlas gordon 352 "$DLABEL_GORDON"
run_atlas probaConns 80 "$DLABEL_PROBA"

echo "Done. Outputs in $ROOT/results/summary/brain/"
