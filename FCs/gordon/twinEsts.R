# srun -N 1  --mem=32gb -t 1:00:00 -p interactive --pty bash 
# module load R/4.4.0-openblas-rocky8
#options(readr.show_progress = TRUE)
library(arrow)
library(tidyverse)
library(mets)
library(txtplot)
args = commandArgs(trailingOnly=TRUE)
iteration = as.numeric(args[2])
phenoNames <- paste0("o", (1 + iteration * 528) : ((iteration+ 1) * 528))

pheno <- read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/pconns.parquet") %>%
  mutate(
      IID = gsub('^(.{4})(.*)$', '\\1_\\2', IID)
  ) %>%
  distinct(IID, .keep_all = T) %>%
  select(IID, all_of(phenoNames))

# 61776 phenos
# 117 evenly spaced parallel threads of 528 blocks
IDs = read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/FC.files", col_names = c("IID")) %>%
  left_join(read_table("~/ABCD/Results/IDs/IDs.txt", col_names=c("FID", "IID")), by = "IID") %>%
  distinct()
pheno = drop_na(left_join(pheno, IDs, by = c("IID")))

# dim(pheno) 1876 x 61777

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

