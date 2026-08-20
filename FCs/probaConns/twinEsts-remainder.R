# srun -N 1  --mem=32gb -t 1:00:00 -p interactive --pty bash 
# module load R/4.4.0-openblas-rocky8
#options(readr.show_progress = TRUE)
library(arrow)
library(tidyverse)
library(mets)

args = commandArgs(trailingOnly=TRUE)
iteration = as.numeric(args[2])

pheno <- read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet") %>%
  mutate(
      IID = gsub('^(.{4})(.*)$', '\\1_\\2', IID)
  )

# 3160 phenos
# 158 evenly spaced parallel threads of 20 blocks
IDs = read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/FC.files", col_names = c("IID"))

pheno = left_join(pheno, IDs, by = c("IID"))

df <- read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv") %>%
  # Select the relevant columns
  select(FID, IID, age, female, site_id_l, household.income, high.educ, genetic_zygosity_status_1, rel_ingroup_order) %>% 
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case=T) ~ "MZ",
    grepl("di", genetic_zygosity_status_1, ignore.case=T) ~ "DZ",
    .default = NA
  )) %>%
    select(-genetic_zygosity_status_1) %>%
    left_join(pheno, by = c("IID"))
    # args 0 - 116
phenoNames <- paste0("o", (1 + iteration * 20) : (((iteration+ 1) * 20) - 1))
phenoNames <- c("o0", phenoNames)
results <- data.frame("phenotype" = c(), "herit" = c(), "variance" = c(), "lowerCI" = c(), "upperCI" = c())
for (ph in phenoNames) {
  print(ph)
  result = tryCatch({
      tempResult <- twinlm(as.formula(paste0(as.character(ph), "~ site_id_l + age + female + household.income + high.educ")), data = as.data.frame(df), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" )
      data.frame("phenotype" = ph, "herit" = summary(tempResult)$heritability[1, 1], "variance" = summary(tempResult)$heritability[1, 2] ** 2, "lowerCI" = summary(tempResult)$heritability[1, 3],
        "upperCI" =  summary(tempResult)$heritability[1, 4])
    }, warning = function(w) {
      tempResult <- twinlm(as.formula(paste0(as.character(ph), "~ site_id_l + age + female + household.income + high.educ")), data = as.data.frame(df), DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" )
       data.frame("phenotype" = ph, "herit" = summary(tempResult)$heritability[1, 1], "variance" = summary(tempResult)$heritability[1, 2] ** 2, "lowerCI" = summary(tempResult)$heritability[1, 3],
        "upperCI" =  summary(tempResult)$heritability[1, 4])
    }, error = function(e) {
      data.frame("phenotype" = ph, "herit" = NA, "variance" = NA, "lowerCI" = NA, "upperCI" = NA)
    })
  results <- rbind(results, result)
}


results %>%
    write_csv(paste0("/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/03_heritability/Topography/FCs/probaConntwinEsts", "_", 
    as.character(args[2]), ".csv"))


