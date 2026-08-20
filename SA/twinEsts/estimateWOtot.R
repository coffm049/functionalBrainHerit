library(tidyverse)
library(mets)
library(txtplot)

# Load data
pheno = read_table("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/allTopoPhenos.csv") %>%
    select(-c(network_surfarea4, network_surfarea6))
#txtdensity(pheno$network_surfarea1)


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
    pivot_longer(cols = c(anthro_height_calc, network_surfarea1:network_surfarea16), names_to = "phenotype") %>%
    nest(data = -phenotype) %>%
    mutate(herit = map(data, 
    ~summary(twinlm(value ~ site_id_l + age + female + household.income + high.educ, data = as.data.frame(.), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" )), .progress = TRUE))

saveRDS(df, "/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/SA/twinEsts/herit_wo_total.Rds")
