# srun -N 1 --mem=32gb -t 1:00:00 -p interactive --pty bash
# module load R/4.4.0-openblas-rocky8
# conda activate twinEst
library(arrow)
library(tidyverse)
library(mets)

args = commandArgs(trailingOnly = TRUE)

to_num <- function(x) suppressWarnings(as.numeric(x))
iteration <- NULL
if (length(args) >= 1 && !is.na(to_num(args[1]))) {
  iteration <- to_num(args[1])
} else {
  task <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  if (task != "" && !is.na(to_num(task))) iteration <- to_num(task)
}
if (is.null(iteration) || iteration < 0 || iteration > 157) {
  stop("No valid SLURM_ARRAY_TASK_ID / iteration argument provided (expected 0-157)")
}
cat("iteration =", iteration, "\n")

CHUNK = 20
# probaConns: 3160 edges -> 158 chunks of 20 (o0..o3159, 0-indexed)
phenoNames <- paste0("o", (iteration * CHUNK) : ((iteration + 1) * CHUNK - 1))

# 1) Load covariates + filtered IDs + FID mapping ONCE
cat("Loading covariates and ID maps...\n")
filtered_ids <- read_table("/projects/standard/rando149/coffm049/filtered_ids_broad.tsv",
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
    "/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet",
    col_select = c("IID", phenoName)
  ) %>% distinct(IID, .keep_all = TRUE) %>%
    inner_join(filtered_ids, by = "IID")
  
  # Join with covariates
  d <- left_join(covars, pheno, by = "IID") %>%
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
out <- paste0("/users/4/coffm049/papers/functionalBrainHerit/results/FCs/probaConns/herit_", iteration, ".Rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(out_df, out)
cat("Saved ", out, "\n")