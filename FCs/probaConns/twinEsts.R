# srun -N 1  --mem=32gb -t 1:00:00 -p interactive --pty bash 
# module load R/4.4.0-openblas-rocky8
#options(readr.show_progress = TRUE)
library(arrow)
library(tidyverse)
library(mets)

args = commandArgs(trailingOnly=TRUE)
iteration = as.numeric(args[2])
# 3160 phenos
# 158 evenly spaced parallel threads of 20 blocks
phenoNames <- paste0("o", (1 + iteration * 20) : ((iteration+ 1) * 20))

pheno <- read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet") %>%
  mutate(
      IID = gsub('^(.{4})(.*)$', '\\1_\\2', IID)
  ) %>%
  distinct(IID, .keep_all = T) %>%
  select(IID, all_of(phenoNames))
IDs = read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/proba.files", col_names = c("IID")) %>%
  mutate(
    IID = basename(IID),
    # grab the 9-20 characers
    IID = substr(IID, 5, 19),
    # insert a _ 5 position
    IID = paste0(substr(IID, 1, 4), "_", substr(IID, 5, 19))
) %>% 
  left_join(read_table("~/ABCD/Results/IDs/IDs.txt", col_names=c("FID", "IID")), by = "IID") %>%
  distinct()

pheno = drop_na(left_join(pheno, IDs, by = c("IID")))

df <- read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv") %>%
  # Select the relevant columns
  select(FID, IID, age, female, site_id_l, household.income, high.educ, anthro_height_calc, genetic_zygosity_status_1) %>% 
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case=T) ~ "MZ",
    grepl("di", genetic_zygosity_status_1, ignore.case=T) ~ "DZ",
    .default = NA
  )) %>%
    #anthro_height_calc = as.numeric(anthro_height_calc)) %>% 
    as.data.frame() %>%
    select(-genetic_zygosity_status_1) %>%
    drop_na() %>% 
    left_join(pheno, by = c("FID", "IID")) %>%
    drop_na() %>%
    pivot_longer(cols = starts_with("o"), names_to = "phenotype") %>%
    nest(data = -phenotype) %>%
    mutate(herit = map(data, 
    ~summary(twinlm(value ~ site_id_l + age + female + household.income + high.educ, data = as.data.frame(.), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" ))))

saveRDS(df, paste0("herit_", iteration, ".Rds"))

