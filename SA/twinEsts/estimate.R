library(tidyverse)
library(mets)

# 17 PFN network surface-area phenotypes (current pipeline)
pheno <- read_csv("/projects/standard/rando149/coffm049/ABCD/Results/02_Phenotypes/TotalCorticalRepresentation_ByPFN_ABCD.csv") %>%
  select(IID, network_surfarea1:network_surfarea17)

# Attach family ID (IDs.txt has no header) and total surface area
IDs <- read_table("/projects/standard/rando149/coffm049/ABCD/Results/IDs/IDs.txt",
                  col_names = c("FID", "IID"))
pheno <- left_join(pheno, IDs, by = "IID") %>% distinct() %>%
  mutate(totalNetworkSurface = rowSums(select(., starts_with("network_surfarea"))))

df <- read_csv("/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv") %>%
  select(FID, IID, age, female, site_id_l, household.income, high.educ, genetic_zygosity_status_1) %>%
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case = TRUE) ~ "MZ",
    grepl("di",   genetic_zygosity_status_1, ignore.case = TRUE) ~ "DZ",
    .default = NA
  )) %>%
  select(-genetic_zygosity_status_1) %>%
  drop_na() %>%
  left_join(pheno, by = c("FID", "IID")) %>%
  drop_na() %>%
  pivot_longer(cols = network_surfarea1:network_surfarea17, names_to = "phenotype") %>%
  nest(data = -phenotype) %>%
  mutate(herit = map(data, function(d) {
    tryCatch(
      summary(twinlm(value ~ site_id_l + age + female + household.income + high.educ + totalNetworkSurface,
                     data = as.data.frame(d), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace")),
      error = function(e) structure(list(error = conditionMessage(e)), class = "twinlm_error"))
  }, .progress = TRUE))

out <- "/users/4/coffm049/papers/functionalBrainHerit/results/SA/twinEsts/herit_w_total.Rds"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(df, out)
