#!/bin/bash
# Create filtered_ids.csv from GRM IDs (2 columns: FID IID, no header)
GRM_ID="/projects/standard/rando149/coffm049/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/no_rels/no_rels.grm.id"
OUT="/projects/standard/rando149/coffm049/filtered_ids.csv"

if [ ! -f "$GRM_ID" ]; then
    echo "ERROR: GRM ID file not found: $GRM_ID"
    exit 1
fi

# Copy first two columns (FID, IID) with no header
awk '{print $1, $2}' "$GRM_ID" > "$OUT"
echo "Created $OUT with $(wc -l < "$OUT") lines"
head -5 "$OUT"