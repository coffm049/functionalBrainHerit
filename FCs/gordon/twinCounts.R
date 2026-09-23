# srun -N 1 --mem=32gb -t 1:00:00 -p interactive --pty bash
# module load R/4.4.0-openblas-rocky8
library(arrow)
library(tidyverse)
library(mets)


pheno <- read_parquet(
  "/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/pconns.parquet",
  col_select = c("IID", "o1")
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
  #left_join(pheno, by = c("FID", "IID")) %>%
  drop_na()

df |>
  reframe(.by = FID, n= n()) |>
  arrange(n) |>
  filter(n>1) |>
  distinct(FID) |>
  dim()

herit <- readRDS("herit_w_total.Rds")
table(herit$zyg)
