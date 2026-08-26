library(tidyverse)
df <- read_csv("data/allEstimatesLabeled.csv") %>%
  filter(Type %in% c("GCTA", "twin")) %>%
  mutate(
    testStat = ifelse(Type == "GCTA",
      h2**2 / var.h2,
      (h2 / ((upperCI - h2) / qnorm(0.975)))**2
    ),
    # for snp estimates
    pval = 1 - (0.5 + pchisq(testStat, df = 1) / 2),
    h2 = ifelse(is.na(h2), 0, h2),
    pval = ifelse(is.na(pval), 1, pval),
  ) %>%
  mutate(
    .by = c(phenoClass, Type),
    adj.pval = p.adjust(pval, method = "bonferroni"),
    # calculate pvalue following mixture of chis
    signif = case_when(
      adj.pval < 0.05 ~ TRUE,
      adj.pval >= 0.05 ~ FALSE,
      is.na(adj.pval) & (testStat == 0) ~ FALSE,
      is.na(adj.pval) & (testStat > 0) ~ 0.5 + pchisq(testStat, df = 1) / 2 > 0.95,
      is.na(adj.pval) & (lowerCI > 0) ~ TRUE,
      .default = FALSE
    ),
    signif = as.logical(signif)
  )

df %>%
  reframe(
    across(c("h2", "testStat", "pval", "adj.pval"), function(x) mean(is.na(x))),
    .by = Type
  )


df %>%
  reframe(fivenum(h2), .by = c(phenoClass, Type), Position = c("Minimum", "First Quartile", "Median", "Third Quartile", "Maximum")) %>%
  pivot_wider(id_cols = c(phenoClass, Type), names_from = Position, values_from = `fivenum(h2)`)
# snpdf %>%
#   reframe(pval = fivenum(pval, na.rm = T), .by = c(phenoClass, Type), Position = c("Minimum", "First Quartile", "Median", "Third Quartile", "Maximum")) %>%
#   pivot_wider(id_cols = c(phenoClass, Type), names_from = Position, values_from = pval)
df %>%
  write_csv("data/allEstimates2.csv")
