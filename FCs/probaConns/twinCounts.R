# srun -N 1 --mem=32gb -t 1:00:00 -p interactive --pty bash
# module load R/4.4.0-openblas-rocky8
library(arrow)
library(tidyverse)
library(mets)

args = commandArgs(trailingOnly = TRUE)
iteration = as.numeric(args[2])

CHUNK = 20
# probaConns: 3160 edges -> 158 chunks of 20 (o0..o3159, 0-indexed)
phenoNames <- paste0("o", (iteration * CHUNK) : ((iteration + 1) * CHUNK - 1))

pheno <- read_parquet(
  "/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet",
  col_select = c("IID", phenoNames)
) %>% distinct(IID, .keep_all = TRUE)

# Family IDs for twin pairing (IDs.txt has no header)
IDs <- read_table("/projects/standard/rando149/coffm049/ABCD/Results/IDs/IDs.txt",
                  col_names = c("FID", "IID"))
pheno <- left_join(pheno, IDs, by = "IID") %>% distinct()

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
  pivot_longer(cols = starts_with("o"), names_to = "phenotype") %>%
  nest(data = -phenotype) %>%
  mutate(herit = map(data, function(d) {
    tryCatch(
      summary(twinlm(value ~ site_id_l + age + female + household.income + high.educ,
                     data = as.data.frame(d), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace")),
      error = function(e) structure(list(error = conditionMessage(e)), class = "twinlm_error"))
  }))

out <- paste0("/users/4/coffm049/papers/functionalBrainHerit/results/FCs/probaConns/herit_", iteration, ".Rds")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(df, out)
