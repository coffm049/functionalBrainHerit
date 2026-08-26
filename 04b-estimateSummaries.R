library(tidyverse)
library(gt)
library(gridExtra)
library(gtsummary)
library(cluster)

df <- read_csv("data/allEstimates.csv") %>%
  #  mutate(sampleSize = (MZ + DZ) * 2) %>%
  # fisher z transform h2
  mutate(z = atanh(h2))

df %>%
    reframe(q = quantile(h2, c(0.05, 0.95), na.rm = T), .by = phenoClass)

df %>%
  ggplot(aes(x = h2, fill = Type)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~phenoClass)



# Mean and IQR
df %>%
  filter(h2<0.9) %>%
  reframe(h2 = fivenum(h2), .by = c(phenoClass, Type), Position = c("Minimum", "First Quartile", "Median", "Third Quartile", "Maximum"), missing = mean(is.na(h2))) %>%
  pivot_wider(id_cols = c(phenoClass, Type), names_from = Position, values_from = h2) %>%
  mutate(IQR = `Third Quartile` - `First Quartile`) %>%
  select(phenoClass, Type, Median, IQR, Maximum) %>%
  filter(Type %in% c("GCTA", "twin"))

# number of signif
df %>%
  mutate(signif = as.logical(signif)) %>%
  reframe(pSignif = mean(signif, na.rm = T), .by = c(phenoClass, Type)) %>%
  filter(Type %in% c("GCTA", "twin"))




# summarize per ROI
temp <- df %>%
  # %% confusion matrix for significance
  mutate(signif = as.logical(signif)) %>%
  pivot_wider(id_cols = c(phenoClass, pheno), names_from = Type, values_from = signif) %>%
  temp() %>%
  count(phenoClass, GCTA, twin) %>%
  drop_na() %>%
  reframe(
    .by = phenoClass,
    TP = sum(ifelse((GCTA == T) & (twin == T), n, 0)),
    FP = sum(ifelse((GCTA == T) & (twin == F), n, 0)),
    FN = sum(ifelse((GCTA == F) & (twin == T), n, 0)),
  ) %>%
  mutate(
    prec = TP / (TP + FP),
    recall = TP / (TP + FN),
    f1 = 2 * ((prec * recall) / (prec + recall))
  )

prec <-
  recall <-
  #


  temp %>%
  filter(phenoClass == "PFN") %>%
  drop_na(GCTA, twin) %>%
  select(SNP = GCTA, Twin = twin) %>%
  tbl_cross(row = SNP, percent = "cell", col = Twin)

temp %>%
  filter(phenoClass == "Gordon") %>%
  drop_na(GCTA, twin) %>%
  select(SNP = GCTA, Twin = twin) %>%
  tbl_cross(row = SNP, percent = "cell", col = Twin)


# Order of networks from Sanju
networks <- c(
  "1" = "DMN", "2" = "VIS", "3" = "FP", "5" = "DAN", "7" = "VAN",
  "8" = "SAL", "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD", "13" = "Tpole",
  "14" = "MTL", "15" = "PMN", "16" = "PON"
)

sas <- list.files("../results/FCs100824/Topography/SA/") %>%
  # map read and join the csv files
  map(function(x) {
    read_csv(paste0("../results/FCs100824/Topography/SA/", x)) %>%
      mutate(Type = x)
  }) %>%
  reduce(bind_rows) %>%
  select(h2, pheno, Type) %>%
  mutate(
    Type = factor(
      case_when(
        Type == "twinEsts.csv" ~ "Twin",
        .default = Type
      )
    ),
    pheno = parse_number(pheno),
    pheno = factor(networks[as.character(pheno)])
  ) %>%
  filter(!grepl("anthro", pheno, ignore.case = T)) %>%
  drop_na(pheno) %>%
  filter(!grepl("Naive", Type, ignore.case = T)) %>%
  mutate(phenoClass = "SA")

df %>%
  mutate(pheno = as.character(pheno)) %>%
  bind_rows(sas) %>%
  select(phenoClass, h2, Type) %>%
  mutate(Type = tolower(Type)) %>%
  filter(Type %in% c("twin", "gcta")) %>%
  mutate(Type = ifelse(Type == "twin", "Twin", "SNP")) %>%
  ggplot(aes(x = phenoClass, y = h2, group = interaction(phenoClass, Type), fill = Type)) +
  geom_boxplot() +
  theme_minimal()
ggsave("h2AcrossScales.png", dpi = 300)

?remove.packages()

df %>%
  mutate(pheno = as.character(pheno)) %>%
  bind_rows(sas) %>%
  select(phenoClass, h2, Type) %>%
  mutate(Type = tolower(Type)) %>%
  filter(Type %in% c("twin", "gcta")) %>%
  mutate(Type = ifelse(Type == "twin", "Twin", "SNP"), index = 1:n()) %>%
  pivot_wider(id_cols = c(Type, index), values_from = h2, names_from = phenoClass) %>%
  select(-index) %>%
  tbl_summary(by = Type, missing = "no", statistic = list(all_continuous() ~ "{median} ({p25}-{p75})")) %>%
  as_gt() %>%
  gt::gtsave(file = "h2AcrossScalesSummary.html")
# temp %>%
#   drop_na(HEreg, twin) %>%
#   tbl_cross(row = HEreg, percent = "cell", col = twin)
# temp %>%
#   drop_na(AdjHE.RE, twin) %>%
#   tbl_cross(row = AdjHE.RE, percent = "cell", col = twin)
# temp %>%
#   drop_na(HEreg, GCTA) %>%
#   tbl_cross(row = GCTA, percent = "cell", col = HEreg)
