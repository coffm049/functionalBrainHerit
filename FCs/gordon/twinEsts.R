# srun -N 1 --mem=32gb -t 1:00:00 -p interactive --pty bash
# module load R/4.4.0-openblas-rocky8
# conda activate twinEst
library(arrow)
library(tidyverse)
library(mets)

args = commandArgs(trailingOnly = TRUE)
iteration = as.numeric(args[1])

CHUNK = 528
# gordon: 61776 edges -> 117 chunks of 528 (o0..o61775, 0-indexed)
phenoNames <- paste0("o", (iteration * CHUNK) : ((iteration + 1) * CHUNK - 1))

# 1) Load covariates + filtered IDs + FID mapping ONCE
cat("Loading covariates and ID maps...\n")
filtered_ids <- read_table("/projects/standard/rando149/coffm049/filtered_ids.tsv",
                           col_names = c("FID", "IID"), show_col_types = FALSE) %>%
  select(IID)

IDs <- read_table("/projects/standard/rando149/coffm049/ABCD/Results/IDs/IDs.txt",
                  col_names = c("FID", "IID"))

covars <- read_csv("/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv", show_col_types = FALSE) %>%
  select(FID, IID, age, female, site_id_l, household.income, high.educ, genetic_zygosity_status_1) %>%
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case = TRUE) ~ "MZ",
    grepl("di",   genetic_zygosity_status_1, ignore.case = TRUE) ~ "DZ",
    .default = NA
  )) %>%
  select(-genetic_zygosity_status_1) %>%
  drop_na() %>%
  inner_join(IDs, by = c("FID", "IID")) %>%
  inner_join(filtered_ids, by = "IID") %>%
  # Filter to FIDs with at least 2 members (twin pairs)
  add_count(FID) %>%
  filter(n >= 2) %>%
  select(-n)

cat("Covariate base: ", nrow(covars), " individuals, ", length(unique(covars$FID)), " families\n")

# 2) Loop through phenotypes ONE AT A TIME
results <- list()
for (phenoName in phenoNames) {
  cat("  Processing ", phenoName, "...\n")
  
  # Read ONLY this phenotype column + IID
  pheno <- read_parquet(
    "/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/pconns.parquet",
    col_select = c("IID", phenoName)
  ) %>% distinct(IID, .keep_all = TRUE) %>%
    inner_join(filtered_ids, by = "IID")
  
  # Join with covariates
  d <- left_join(covars, pheno, by = c("FID", "IID")) %>%
    drop_na(all_of(phenoName))
  
  if (nrow(d) < 50) {
    cat("    Skipping ", phenoName, ": only ", nrow(d), " obs\n")
    results[[phenoName]] <- tibble(
      phenotype = phenoName,
      herit = list(structure(list(error = "insufficient data"), class = "twinlm_error"))
    )
    next
  }
  
  # Estimate
  est <- tryCatch(
    summary(twinlm(
      as.formula(paste(phenoName, "~ site_id_l + age + female + household.income + high.educ")),
      data = as.data.frame(d), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace"
    )),
    error = function(e) structure(list(error = conditionMessage(e)), class = "twinlm_error")
  )
  
  results[[phenoName]] <- tibble(
    phenotype = phenoName,
    herit = list(est)
  )
}

# 3) Combine and save
out_df <- bind_rows(results)
out <- paste0("/users/4/coffm049/papers/functionalBrainHerit/results/FCs/gordon/herit_", iteration, ".Rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(out_df, out)
cat("Saved ", out, "\n")