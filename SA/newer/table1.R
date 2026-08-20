library(tidyverse)
library(table1)
library(kableExtra)
library(arrow)

json_data <- jsonlite::fromJSON("GCTAbigN.json")
phenojson <- jsonlite::fromJSON("/home/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Topography/FCs/probaConns/temp/GCTA.FE.0.json")
pheno <- read_table(json_data$pheno)
covar <- read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/Covars2.csv") %>%
  # Select the relevant columns
  select(FID, IID, Age = age, Female= female, Site= site_id_l, `Household income` = household.income, `Higher ed` = high.educ, genetic_zygosity_status_1, `Race/Ethnicity` = race_ethnicity) %>%
  mutate(IID = str_remove(IID, "_"))

gordon <- data.frame("IID" = names(read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/data/probaConns/part.0.parquet"))) %>%
  mutate(gordon = T)
proba <- read_parquet("/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/probaConns.parquet", col_select = "IID") %>%
  mutate(proba = T)

full_join(gordon, proba) %>%
  left_join(covar, by = c("IID")) %>%
  as.data.frame() %>%
  select(-starts_with("network"), -starts_with("genus")) %>%
  slice_head(n=2, by = FID)  %>%
    mutate(
        `Higher ed` = factor(`Higher ed`, levels = c("< HS Diploma", "HS Diploma/GED", "Some College", "Bachelor", "Post Graduate Degree")),
        `Household income` = factor(`Household income`, levels = c("[<50K]", "[>=50K & <100K]", "[>=100K]")),
        Female = as.logical(ifelse(Female == "yes", 1, 0)),
        gordon = ifelse(is.na(gordon), F, T),
        proba = ifelse(is.na(proba), F, T),
        Imaging = case_when(
            gordon & proba ~ "Gordon + Proba",
            proba  ~ "Proba",
    ),
      Relative = case_when(
        grepl("mono", genetic_zygosity_status_1) ~ "MZ",
        grepl("di", genetic_zygosity_status_1) ~ "DZ",
        grepl("sib", genetic_zygosity_status_1) ~ "Sib",
      .default = "Unrelated"
    )
  )  %>%
  mutate(.by = FID,
    Relative = case_when(
      sum(Relative == "DZ") ==2 ~ "DZ",
      sum(Relative == "MZ") ==2 ~ "MZ",
      .default = "Unrelated"
    )
  ) %>%
  table1(~ Female + Age + `Higher ed` + `Household income` +  `Race/Ethnicity` + Imaging | Relative, data = .,
  render.missing=NULL) %>%
  saveRDS("table1.Rds")

pheno %>%
  left_join(ids, by = c("FID", "IID")) %>%
  left_join(covar, by = c("FID", "IID")) %>%
  as.data.frame() %>%
  mutate(
      `Race/Ethnicity` = factor(ifelse(is.na(`Race/Ethnicity`), "Missing", `Race/Ethnicity`)),
      `Higher ed` = factor(`Higher ed`, levels = c("< HS Diploma", "HS Diploma/GED", "Some College", "Bachelor", "Post Graduate Degree")),
      `Household income` = factor(`Household income`, levels = c("[<50K]", "[>=50K & <100K]", "[>=100K]")),
      Female = as.logical(ifelse(Female == "yes", 1, 0)),
      Relative = case_when(
        grepl("mono", genetic_zygosity_status_1) ~ "MZ",
        grepl("di", genetic_zygosity_status_1) ~ "DZ",
        grepl("sib", genetic_zygosity_status_1) ~ "Sib",
      .default = "Unrelated"
    )
  )  %>%
  select(-starts_with("network"), -starts_with("genus")) %>%
  slice_head(n=2, by = FID)  %>%
  table1(~ Female + Age + `Higher ed` + `Household income` +  `Race/Ethnicity` | Relative, data = .) %>%
  saveRDS("table1.Rds")



