library(tidyverse)
PFNgrid <- expand.grid(
  "phenoClass" = "PFN",
  "Type" = c("GCTA", "twin"),
  "pheno" = 0:3159
)
gordongrid <- expand.grid(
  "phenoClass" = "Gordon",
  "Type" = c("GCTA", "twin"),
  "pheno" = 0:61775
)
fullGrid <- rbind(PFNgrid, gordongrid)

ests <- list.files("../results/FCs100824/Topography/FCs", pattern = ".*.csv$", full.names = T)
ests <- ests[!grepl(".log", ests)]
snpests <- ests[!grepl("twin.*.csv", ests)]
snpests <- snpests[grepl("gcta", snpests, ignore.case=T)]
# snpests <- snpests[!grepl("*HEreg*", snpests)]
# remove from snpests if G not in column names given each element it is a csv
# snpEstsG <- sapply(snpests, function(x) "G" %in% names(read_csv(x)))
# snpests <- snpests[snpEstsG]
# deleted one empty file GCTA 1
# hiddenproba <- list.files("../results/FCs100824/FCs/probaConns", pattern = ".*.csv", full.names = T)
# hiddenEstsG <- sapply(hiddenproba, function(x) "G" %in% names(read_csv(x)))
# hiddenproba <- hiddenproba[hiddenEstsG]

PFNsnp <- ests[grepl("proba", snpests, ignore.case= T)]
test <- read_csv(PFNsnp[1])

snpdf <- c(snpests, twinEsts) %>%
  # map read and join the csv files
  map(
    .progress = TRUE,
    function(x) {
      read_csv(
        x,
        # enforce h2 as a numeric column
        col_types = cols(
          pheno = col_character(),
          h2 = col_double(),
          herit = col_double(),
          `var(h2)` = col_double(),
          G = col_double(),
          E = col_double(),
          pval = col_double(),
          PCs = col_double(),
          time = col_double(),
          mem = col_double(),
          lowerCI = col_double(),
          upperCI = col_double()
        ),
      ) %>%
        mutate(
          Type = x,
          # extract the number from Type
          number = as.integer(str_extract(basename(Type), "[[:digit:]]+")),
          Type = basename(Type),
          Type = str_extract(Type, "HEreg|GCTA|twin|AdjHE\\.FE|AdjHE\\.RE"),
          phenoClass = ifelse(grepl("proba", x, ignore.case = T), "PFN", "Gordon")
        )
    }
  ) %>%
  reduce(bind_rows) %>%
  arrange(phenoClass, Type) %>%
  # rename(h2 = herit, var.h2 = variance) %>%
  mutate(
    # Type = factor(ifelse(Type == "AdjHE.R", "AdjHE.RE", Type)),
    h2 = coalesce(h2, herit),
    h2 = ifelse(h2 > 0.97, 0, h2),
    h2 = ifelse(h2 < 0, 0, h2),
    phenotype = as.integer(str_extract(phenotype, "[[:digit:]]+")),
    pheno = as.integer(str_extract(pheno, "[[:digit:]]+")),
    pheno = coalesce(pheno, phenotype),
    var.h2 = coalesce(`var(h2)`, variance),
    testStat = ifelse(Type == "GCTA",
      h2**2 / var.h2,
      (h2 /((upperCI - h2) / qnorm(0.975)))**2),
    # for snp estimates
    pval = 1 - (0.5 + pchisq(testStat, df = 1) / 2),
  ) %>%
  filter((PCs != 0) | is.na(PCs)) %>%
  select(-c(phenotype, herit, `var(h2)`, variance, time, mem, Covariates, PCs)) %>%
  distinct() %>%
  reframe(across(where(is.numeric), function(x) mean(x, na.rm = T)), across(where(is.logical), function(x) max(x, na.rm = T)), .by = c(Type, phenoClass, pheno)) %>%
  right_join(fullGrid, by = c("phenoClass", "Type", "pheno")) %>%
  #mutate(
  #  h2 = ifelse(is.na(h2), 0, h2),
  #  pval = ifelse(is.na(pval), 1, pval),
  #) %>%
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

snpdf %>% 
  reframe(
    across(c("h2", "testStat", "pval", "adj.pval"),
      function(x) mean(is.na(x))), .by = Type)


snpdf %>%
  reframe(fivenum(h2), .by = c(phenoClass, Type), Position = c("Minimum", "First Quartile", "Median", "Third Quartile", "Maximum")) %>%
  pivot_wider(id_cols = c(phenoClass, Type), names_from = Position, values_from = `fivenum(h2)`)
# snpdf %>%
#   reframe(pval = fivenum(pval, na.rm = T), .by = c(phenoClass, Type), Position = c("Minimum", "First Quartile", "Median", "Third Quartile", "Maximum")) %>%
#   pivot_wider(id_cols = c(phenoClass, Type), names_from = Position, values_from = pval)
snpdf %>%
  write_csv("data/allEstimates.csv")
