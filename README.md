# functionalBrainHerit

ABCD heritability pipeline (SNP heritability via MASH + twin ACE estimates via `mets::twinlm`)
for three phenotype sets:

- **SA** — 17 PFN-network cortical surface-area phenotypes
- **FCs / gordon** — 61,776 Gordon-parcellation functional-connectivity edges (`pconns.parquet`)
- **FCs / probaConns** — 3,160 probabilistic-connectome FC edges (`probaConns.parquet`)

All IDs use the unified `sub-<id>` scheme. MASH outputs land under
`results/`; twin estimates under `results/*/twinEsts/`; all SLURM stdout/stderr
under `logs/` (job names: `SA_MASH`, `SA_twin`, `gordon_twin`,
`probaConns_twin`, `FCgordon_<METHOD>_<KIND>`, `FCproba_<METHOD>_<KIND>`).

## Prerequisites (run once on the HPC)

```bash
cd /users/4/coffm049/papers/functionalBrainHerit
git pull
mkdir -p logs

# SA MASH requires FID+IID in the phenotype CSV (MASH needs both columns).
# Add FID = IID to TotalCorticalRepresentation_ByPFN_ABCD.csv if not already done:
python - <<'PY'
import pandas as pd
f = "/projects/standard/rando149/coffm049/ABCD/Results/02_Phenotypes/TotalCorticalRepresentation_ByPFN_ABCD.csv"
d = pd.read_csv(f)
if "FID" not in d.columns:
    d.insert(0, "FID", d["IID"])
    d.to_csv(f, index=False)
PY
```

## Jobs

### 1. SA — SNP heritability (MASH)
```bash
sbatch SA/Est.SLURM
```
Runs `GCTA_wo_total.json` and `AdjHE_RE_wo_total.json` (17 `network_surfarea*`
phenotypes). Output: `results/SA/{GCTA_wo_total,AdjHE_RE_wo_total}/`.
*Depends on the FID-column prerequisite above.*

### 2. SA — twin ACE estimates
```bash
sbatch SA/twinEsts/estimate.SLURM
```
Runs `estimate.R` (with `totalNetworkSurface` covariate → `herit_w_total.Rds`)
and `estimateWOtot.R` (without → `herit_wo_total.Rds`).
Output: `results/SA/twinEsts/`.

### 3. FCs / gordon — SNP heritability (MASH)
```bash
bash FCs/gordon/submit.sh <Method> <FE|RE>      # e.g. AdjHE FE / AdjHE RE
```
Generates per-chunk configs and submits a SLURM array (`run.slurm`).
Output: `results/FCs/gordon/`.

### 4. FCs / gordon — twin ACE estimates
```bash
sbatch FCs/gordon/twinEsts.SLURM
```
Array `0-116` (117 chunks × 528 edges). Output: `results/FCs/gordon/herit_<i>.Rds`.

### 5. FCs / probaConns — SNP heritability (MASH)
```bash
bash FCs/probaConns/submit.sh <Method> <FE|RE>
```
Generates per-chunk configs and submits a SLURM array (`run.slurm`).
Output: `results/FCs/probaConns/`.

### 6. FCs / probaConns — twin ACE estimates
```bash
sbatch FCs/probaConns/twinEsts.SLURM
```
Array `0-157` (158 chunks × 20 edges). Output: `results/FCs/probaConns/herit_<i>.Rds`.

## Notes
- Twin scripts wrap each `twinlm` in `tryCatch`; a failed edge is stored as a
  `twinlm_error` entry (job still completes) rather than killing the chunk.
- After twin runs, verify: RDS counts (gordon 117, probaConns 158, SA 2), scan for
  `twinlm_error` entries, and confirm `zyg` MZ/DZ pairs exist per edge.
