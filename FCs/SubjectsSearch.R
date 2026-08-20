# module load R/4.4.0-openblas-rocky8

library(tidyverse)
library(arrow)
#pheno = read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/pconns.csv", n_max= Inf, col_select = 1, col_names=F)
#pheno = t(read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/pconns.parquet", n_max= 2, col_names = T))

FCfiles = read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/FC.files", col_names = c("IID")) %>%
    mutate(
        IID = stringr::str_remove(IID, "^NDAR_")
    )
SAs = read_table("/panfs/jay/groups/31/rando149/coffm049/Temppheno.csv") %>% 
    select(IID, FID) %>%
    mutate(
        IID = stringr::str_remove(IID, "^NDAR_")
    )
phenoSubjects = data.frame(IID = rownames(as.data.frame(t(read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/data/probaConns/part.0.parquet"))))) %>%
    mutate(
        IID = stringr::str_remove(IID, "^NDAR")) %>%
    filter(!grepl("_null_", IID, fixed = TRUE))


twinIDs <- read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv") %>%
  # Select the relevant columns
  select(FID, IID, age, female, site_id_l, household.income, high.educ, genetic_zygosity_status_1, rel_ingroup_order) %>% 
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case=T) ~ "MZ",
    grepl("di", genetic_zygosity_status_1, ignore.case=T) ~ "DZ",
    .default = NA
  )) %>%
    select(-genetic_zygosity_status_1) %>%
    mutate(
        IID = stringr::str_remove(IID, "^NDAR_"))


# Find which twin subjects were in SAs
# See which of those are in gordon
# See if those files exist
