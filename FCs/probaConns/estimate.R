library(arrow)
library(tidyverse)
library(mets)

# Load data
pheno = t(read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/data/probaConns/part.0.parquet")) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column(var = "IID") %>%
    mutate(
        IID = gsub('^(.{4})(.*)$', '\\1_\\2', IID)
    )
df <- read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv") %>%
  # Select the relevant columns
  select(FID, IID, age, female, site_id_l, household.income, high.educ, genetic_zygosity_status_1, rel_ingroup_order) %>% 
  mutate(zyg = case_when(
    grepl("mono", genetic_zygosity_status_1, ignore.case=T) ~ "MZ",
    grepl("di", genetic_zygosity_status_1, ignore.case=T) ~ "DZ",
    .default = NA
  )) %>%
    select(-genetic_zygosity_status_1) %>%
    drop_na() %>%
    left_join(pheno, by = "IID") %>%
    drop_na() %>%
    pivot_longer(cols = c(starts_with("V")), names_to = "phenotype") %>%
    mutate(phenotype = as.factor(phenotype)) %>%
    nest(data = -phenotype) 

df %>%
    mutate(herit = map(data, 
    ~as.data.frame(data) %>% twinlm(value ~ site_id_l + age + female + household.income + high.educ, data = ., DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" )))
    unnest(herit) 

test <- as.data.frame(df$data[[1]])
twinlm(age ~ female, data = test, DZ = "DZ", zyg = "zyg", id = "FID", type = "ace", twinnum = "rel_ingroup_order")


twinlm(value ~ site_id_l + age + female + household.income + high.educ, data = test2, DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" )
twinlm(value ~ site_id_l + age + female + household.income + high.educ, data = test3, DZ = "DZ", zyg = "zyg", id = "FID", type = "ace" )


df %>%
    select(phenotype, herit) %>%
    write_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/03_heritability/Topography/SA/twinEsts.csv")



# Assuming test3 is correctly defined as a data frame from df$data[[1]]
test3 <- as.data.frame(df$data[[1]])
test3 <- as.data.frame(df$data[[1]])

# Verify column names in test3
colnames(test3)

# Check if 'value' column exists in test3
"value" %in% colnames(test3)

# Perform the heritability analysis using twinlm
twinlm(value ~ age, data = test3, zyg = "zyg", DZ = "DZ", type = "ace", id = "FID", twinnum = "rel_ingroup_order")
twinlm(value ~ age, data = test3, zyg = "zyg", DZ = "DZ", type = "ace", id = "FID")

# Access and print the heritability estimate
print(summary(result)$heritability[1, 1])

run_twinlm <- function(data) {
  tryCatch(
    {
      result <- twinlm(value ~ age, data = drop_na(as.data.frame(data)), zyg = "zyg", DZ = "DZ", type = "ace", id = "FID", twinnum = "rel_ingroup_order")
      summary(result)$heritability[1, 1]
    },
    error = function(err) {
      # Return NA if no observations or other error occurs
      NA
    }
  )
}


results <- df |>
  split(df$phenotype) |>
  map(\(data) run_twinlm(data))