library(arrow)
library(tidyverse)
library(mets)

phenoNames <- c("o1")


pheno <- read_parquet(
  # "/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet",
  "/scratch.global/coffm049/fctwinCheck/probaConns.parquet",,
  col_select = c("IID", phenoNames)
) %>% distinct(IID, .keep_all = TRUE)

# Family IDs for twin pairing (IDs.txt has no header)
IDs <- read_table(#"/projects/standard/rando149/coffm049/ABCD/Results/IDs/IDs.txt",
      "/scratch.global/coffm049/fctwinCheck/IDs.txt",
                  col_names = c("FID", "IID"))
pheno <- left_join(pheno, IDs, by = "IID") %>% distinct()

df <- read_csv(#"/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv"
        "/scratch.global/coffm049/fctwinCheck/Covars2.csv"
               ) %>%
  select(FID, IID, age, female, site_id_l, household.income, high.educ, genetic_zygosity_status_1) %>%
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case = TRUE) ~ "MZ",
    grepl("di",   genetic_zygosity_status_1, ignore.case = TRUE) ~ "DZ",
    .default = NA
  )) %>%
  select(-genetic_zygosity_status_1) %>%
  drop_na() %>%
  left_join(pheno, by = c("FID", "IID")) %>%
  drop_na()

mod <- summary(twinlm(value ~ site_id_l + age + female + household.income + high.educ, data = as.data.frame(df), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace"))

library(arrow)
library(arrow)
