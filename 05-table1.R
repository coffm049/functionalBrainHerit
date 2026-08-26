library(tidyverse)
library(table1)



avoid <- read_csv("data/derivatives_subjects_to_avoid_minus_dwi.csv", col_names = c("IID", "year")) %>%
    filter(grepl("baseline", year)) %>%
    mutate(IID = paste0(str_sub(IID, 1, 4), "_", str_sub(IID, 5, -1)))

# goopd
# geno <- read_table("data/no_rels.grm.id", col_names = c("FID", "IID"))
geno <- read_table("data/full.grm.id", col_names = c("FID", "IID"))

# good
probas <- read_csv("data/proba.files", col_names = "file") %>%
    mutate(
        IID = str_extract(file, "sub-.*_ses"),
        IID = str_remove(IID, "sub-"),
        IID = str_remove(IID, "_ses"),
        # add _ between 4 and 5th charcater and the rest of the string using paste
        IID = paste0(str_sub(IID, 1, 4), "_", str_sub(IID, 5, -1)),
    )

# good
fcFiles <- read_csv("data/FC.files", col_names = "IID")

df <- read_csv("data/deitend.csv") %>%
    distinct() %>%
    filter(IID %in% fcFiles$IID | IID %in% probas$IID) %>%
    filter(IID %in% geno$IID) %>%
    filter(!IID %in% avoid$IID) %>%
    # inner_join(geno, by = c("FID", "IID")) %>%
    # inner_join(probas, by = c("IID")) %>%
    # inner_join(fcFiles, by = "IID") %>%
    # only include individuals who had scans
    mutate(
        race_ethnicity = factor(ifelse(is.na(race_ethnicity), "Missing", race_ethnicity)),
        high.educ = factor(high.educ, levels = c("< HS Diploma", "HS Diploma/GED", "Some College", "Bachelor", "Post Graduate Degree")),
        household.income = factor(household.income, levels = c("[<50K]", "[>=50K & <100K]", "[>=100K]", "Bachelor", "Post Graduate Degree")),
        genetic_zygosity_status_1 = ifelse(is.na(genetic_zygosity_status_1), "None", genetic_zygosity_status_1),
        genetic_zygosity_status_1 = factor(genetic_zygosity_status_1, levels = c("None", "siblings", "dizygotic", "monozygotic")),
        female = as.logical(ifelse(female == "yes", 1, 0))
    ) %>%
    rename(
        Age = age,
        Female = female,
        `Household income` = household.income,
        `Parent education` = high.educ,
        `Race/ethnicity` = race_ethnicity,
        `MRI manufacturer` = mri_info_manufacturer,
        Relatives = genetic_zygosity_status_1
    ) %>%
    distinct()

# make sure only 1 from each family is included
df %>%
    slice_head(n = 1, by = FID) %>%
    table1(~ Female + Age + `Parent education` + `Household income` + `MRI manufacturer` + Relatives + `Race/ethnicity`, data = .)




# how mnay twin pairs
df %>%
    filter(n() >= 2, .by = FID) %>%
    slice_head(n = 2, by = FID) %>%
    filter(Relatives != "None") %>%
    filter(Relatives != "siblings") %>%
    table1(~ Female + Age + `Parent education` + `Household income` + `MRI manufacturer` + `Race/ethnicity` | Relatives, data = .)
df %>%
    slice_head(n = 2, by = FID) %>%
    filter(n() == 2, .by = FID) %>%
    table1(~ Female + Age + `Parent education` + `Household income` + `MRI manufacturer` + `Race/ethnicity` + Relatives | site_id_l, data = .)
