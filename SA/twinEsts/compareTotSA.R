library(tidyverse)

with <- readRDS("herit_w_total.Rds") %>%
  mutate(h2 = map(herit, function(x) x$heritability[1, "Estimate"])) %>%
  unnest(h2) %>%
  filter()
wo <- readRDS("herit.Rds") %>%
  mutate(h2 = map(herit, function(x) x$heritability[1, "Estimate"])) %>%
  unnest(h2)

cor(with$h2, wo$h2)



# snp
with <- read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/03_heritability/Topography/SA/GCTA")
wo <- read_csv("/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/03_heritability/Topography/SA/GCTA_wo_total")


cor(with$h2, wo$h2, use = "pairwise.complete.obs")
