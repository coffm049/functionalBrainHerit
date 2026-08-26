library(tidyverse)
library(cowplot)
library(ggsci)
scaleFUN <- function(x) sprintf("%.2f", x)
theme_set(theme_bw())
options(
  ggplot2.discrete.colour = ggsci::scale_colour_d3,
  ggplot2.discrete.fill = ggsci::scale_fill_d3,
  axis.text = element_text(face = "bold")
)


# Order of networks from Sanju
networks <- c(
  "1" = "DMN", "2" = "VIS", "3" = "FP", "5" = "DAN", "7" = "VAN",
  "8" = "SAL", "9" = "CO", "10" = "SMD", "11" = "SML", "12" = "AUD", "13" = "Tpole",
  "14" = "MTL", "15" = "PMN", "16" = "PON"
)

df <- list.files("../results/FCs100824/Topography/SA/") %>%
  # map read and join the csv files
  map(function(x) {
    read_csv(paste0("../results/FCs100824/Topography/SA/", x)) %>%
      mutate(Type = x)
  }) %>%
  reduce(bind_rows) %>%
  #  select(h2, pheno, Type) %>%
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
  filter(Type %in% c("GCTA", "Twin")) %>%
  mutate(h2 = ifelse(is.na(h2), 0, h2)) %>%
  mutate(h2 = ifelse(h2 < 0, 0, h2)) %>%
#  filter(!grepl("anthro", pheno, ignore.case = T)) %>%
  drop_na(pheno) %>%
  filter(!grepl("Naive", Type, ignore.case = T)) # %>%
#  mutate(
#    testStat = h2**2 / `var(h2)`,
#    # calculate pvalue following mixture of chis
#    signif = case_when(
#      pval < 0.05 ~ TRUE,
#      pval >= 0.05 ~ FALSE,
#      is.na(pval) & (testStat == 0) ~ FALSE,
#      is.na(pval) & (testStat > 0) ~ 0.5 + pchisq(testStat, df = 1) / 2 > 0.95,
#      is.na(pval) & (lowerCI > 0) ~ TRUE,
#      .default = FALSE
#    ),
#    signif = as.logical(signif)
#  )

newOrder <- df %>%
  reframe(m = max(h2, na.rm = T), .by = pheno) %>%
  arrange(m) %>%
  pull(pheno)




df %>%
  arrange(h2) %>%
  mutate(pheno = factor(pheno, levels = newOrder)) %>%
  filter(Type %in% c("Twin", "GCTA")) %>%
  mutate(Type = ifelse(Type == "GCTA", "SNP", "Twin")) %>%
  ggplot(aes(y = pheno, x = h2, fill = Type, group = Type)) +
  geom_bar(stat = "identity", position = "dodge2") +
  #xlim(0, 0.6) +
  theme_minimal() +
  scale_color_npg() +
  scale_fill_npg() +
  labs(x = expression(paste("Heritability (", h^2, ")")), y = "") +
  theme(
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.line = element_line(linewidth = 0.5),
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 24, face = "bold")
  )




# USE npg theme
ggsave2("figures/topo/SAs.png", dpi = 300)



# g2 <- df %>%
#   ggplot(aes(y = pheno, x = h2, fill = Type, group = Type)) +
#   geom_bar(stat = "identity", position = "dodge2") +
#   # remove the legend
#   theme_bw() +
#   scale_color_npg() +
#   scale_fill_npg() +
#   theme(
#     legend.position = "none",
#     text = element_text(size = 8)
#   ) +
#   labs(x = "", y = "") +
#   xlim(0, 0.1)
# # scale_x_continuous(labels = scaleFUN, limits = c(0, 0.1))
# combined_plot <- g1 + annotation_custom(ggplotGrob(g2),
#   xmin = 0.5, xmax = 1
# )
#
# ggsave2("SAs.png", plot = combined_plot, dpi = 300)
sas <- list.files("../results/SA110624", full.names = T)[1:2] %>%
  # map read and join the csv files
  map(function(x) {
    read_csv(x) %>%
      mutate(
        Type = basename(x),
        PC = rep(c(10, 30), 30),
        phenoClass = str_extract(Type, "genus"),
        Type = str_extract(Type, "AdjHE_RE|GCTA")
      ) %>%
      filter(PC == 30) %>%
      filter(!grepl("network_surfarea", pheno))
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
    pheno = rep(networks, 2)
  ) %>%
  mutate(
    phenoClass = "Genus",
  ) %>%
  filter(Type == "GCTA") %>%
  mutate(Type = "SNP")

twinSA <- readRDS("../results/SA110624/SA_genus_twin.Rds")$herit %>%
  map(function(x) {
    x[["heritability"]]
  }) %>%
  reduce(rbind) %>%
  as.data.frame() %>%
  mutate(
    phenoClass = c("height", rep(c("genus", "sa"), each = 14)),
    pheno = c("height", rep(networks, 2)),
    Type = "Twin",
    signif = `2.5%` > 0,
  ) %>%
  rename(h2 = Estimate) %>%
  filter(phenoClass == "genus")
rownames(twinSA) <- NULL

df <- bind_rows(sas, twinSA)
df %>%
  mutate(
    pheno = factor(pheno, levels = newOrder),
    Type = str_replace(Type, "GCTA", "SNP"),
  ) %>%
  filter(pheno != "height") %>%
  mutate(h2 = ifelse(is.na(h2), 0, h2)) %>%
  ggplot(aes(x = h2, y = pheno, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge2") +
  xlim(0, 1) +
  theme_minimal() +
  scale_color_npg() +
  scale_fill_npg() +
  labs(x = expression(paste("Heritability (", h^2, ")")), y = "") +
  theme(
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.line = element_line(linewidth = 0.5),
  )
ggsave("figures/topo/genusH2.png", dpi = 300)



phenos <- read_table("allTopoPhenos.csv") %>%
  left_join(read_csv("deitend.csv")) %>%
  select(race_ethnicity, contains("genus"), contains("network_surfarea")) %>%
  pivot_longer(
    cols = -race_ethnicity, names_pattern = "([[:alpha:]]+_*[[:alpha:]]+)(\\d{1,2})",
    names_to = c("phenoClass", "phenoNumber")
  ) %>%
  mutate(
    phenoNumber = as.numeric(phenoNumber),
    phenoNumber = case_when(
      (phenoClass == "network_surfarea") & (phenoNumber >= 7) ~ phenoNumber - 2,
      (phenoClass == "network_surfarea") & (phenoNumber >= 5) ~ phenoNumber - 1,
      .default = phenoNumber
    ),
    phenoNumber = factor(phenoNumber)
  ) %>%
  drop_na(race_ethnicity)

phenos %>%
  ggplot(aes(x = phenoNumber, y = value, fill = race_ethnicity)) +
  geom_boxplot(position = "dodge2") +
  facet_wrap(~phenoClass, scales = "free", nrow = 2) +
  theme_minimal() +
  ggsci::scale_fill_npg()
ggsave("topoVeth.png", dpi = 300)
