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
with <- read_csv("/users/4/coffm049/papers/functionalBrainHerit/results/SA/GCTA")
wo <- read_csv("/users/4/coffm049/papers/functionalBrainHerit/results/SA/GCTA_wo_total")


cor(with$h2, wo$h2, use = "pairwise.complete.obs")
