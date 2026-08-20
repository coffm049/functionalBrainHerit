library(tidyverse)
library(mets)
library(txtplot)

# Load data
pheno = read_table("~/Temppheno.csv") %>%
    select(-c(network_surfarea4, network_surfarea6))
#txtdensity(pheno$network_surfarea1)


df <- read_csv("/projects/standard/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv") %>%
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
    mutate(totalNetworkSurface = network_surfarea1 + network_surfarea2 + network_surfarea3 +network_surfarea5 +network_surfarea7 +network_surfarea8 +network_surfarea9 +network_surfarea10 +network_surfarea11 +network_surfarea12 +network_surfarea13 +network_surfarea14 +network_surfarea15+network_surfarea16) %>%
    pivot_longer(cols = c(anthro_height_calc, network_surfarea1:network_surfarea16), names_to = "phenotype") %>%
    nest(-phenotype) %>%
    mutate(herit = map(data, 
    ~summary(twinlm(value ~ site_id_l + age + female + household.income + high.educ + totalNetworkSurface, data = as.data.frame(.), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" )), .progress = TRUE))

saveRDS(df, "/users/4/coffm049/papers/functionalBrainHerit/results/SA/twinEsts/herit_w_total.Rds")
