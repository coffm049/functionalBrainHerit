library(tidyverse)

df <- read_csv("data/allEstimatesLabeled.csv")
# kendall tau

cMethod <- df %>%
  filter(h2 != 0) %>%
  pivot_wider(id_cols = c(pheno, phenoClass), names_from = Type, values_from = h2) %>%
  select(pheno, GCTA, twin, phenoClass) %>%
  # head() %>%
  reframe(c = cor(GCTA, twin, method = "kendall", use = "pairwise.complete.obs"), .by = phenoClass)

# System level
cgMethod <- df %>%
  filter(h2 != 0) %>%
  reframe(h2 = median(h2, na.rm = T), .by = c(phenoClass, connection, Type)) %>%
  pivot_wider(id_cols = c(connection, phenoClass), names_from = Type, values_from = h2) %>%
  select(connection, GCTA, twin, phenoClass) %>%
  # head() %>%
  reframe(c = cor(GCTA, twin, method = "kendall", use = "pairwise.complete.obs"), .by = phenoClass)

# system level methyod correlation
cgPheno <- df %>%
  filter(h2 != 0) %>%
  reframe(h2 = median(h2, na.rm = T), .by = c(phenoClass, connection, Type)) %>%
  pivot_wider(id_cols = c(connection, Type), names_from = phenoClass, values_from = h2) %>%
  select(connection, Gordon, PFN, Type) %>%
  # head() %>%
  reframe(c = cor(Gordon, PFN, method = "kendall", use = "pairwise.complete.obs"), .by = Type)

df %>%
  filter(h2 > 0) %>%
  reframe(m = median(h2, na.rm = T), .by = c(phenoClass, Type))
