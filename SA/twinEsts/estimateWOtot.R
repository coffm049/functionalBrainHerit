# conda activate twinEst
library(tidyverse)
library(mets)

# 17 PFN network surface-area phenotypes (current pipeline)
phenoNames <- paste0("network_surfarea", 1:17)

# 1) Load covariates + filtered IDs + FID mapping ONCE
cat("Loading covariates and ID maps...\n")
filtered_ids <- read_csv("/projects/standard/rando149/coffm049/filtered_ids.csv", col_names = c("IID"), show_col_types = FALSE)

IDs <- read_table("/projects/standard/rando149/coffm049/ABCD/Results/IDs/IDs.txt",
                  col_names = c("FID", "IID"))

pheno <- read_csv("/projects/standard/rando149/coffm049/ABCD/Results/02_Phenotypes/TotalCorticalRepresentation_ByPFN_ABCD.csv", show_col_types = FALSE) %>%
  select(IID, all_of(phenoNames)) %>%
  inner_join(filtered_ids, by = "IID") %>%
  left_join(IDs, by = "IID") %>% distinct()

covars <- read_csv("/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv", show_col_types = FALSE) %>%
  select(FID, IID, age, female, site_id_l, household.income, high.educ, genetic_zygosity_status_1) %>%
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case = TRUE) ~ "MZ",
    grepl("di",   genetic_zygosity_status_1, ignore.case = TRUE) ~ "DZ",
    .default = NA
  )) %>%
  select(-genetic_zygosity_status_1) %>%
  drop_na() %>%
  inner_join(pheno, by = c("FID", "IID")) %>%
  # Filter to FIDs with at least 2 members (twin pairs)
  add_count(FID) %>%
  filter(n >= 2) %>%
  select(-n)

cat("Covariate base: ", nrow(covars), " individuals, ", length(unique(covars$FID)), " families\n")

# 2) Loop through phenotypes ONE AT A TIME
results <- list()
for (phenoName in phenoNames) {
  cat("  Processing ", phenoName, "...\n")
  
  d <- covars %>% drop_na(all_of(phenoName))
  
  if (nrow(d) < 50) {
    cat("    Skipping ", phenoName, ": only ", nrow(d), " obs\n")
    results[[phenoName]] <- tibble(
      phenotype = phenoName,
      herit = list(structure(list(error = "insufficient data"), class = "twinlm_error"))
    )
    next
  }
  
  # Estimate - NO totalNetworkSurface covariate (wo_total)
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
out <- "/standard/projects/coffm049/papers/functionalBrainHerit/results/SA/twinEsts/herit_wo_total.Rds"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(out_df, out)
cat("Saved ", out, "\n")