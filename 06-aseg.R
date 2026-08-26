library(tidyverse)
twinVol <- readRDS("data/asegTwinHerit.Rds") %>%
  select(-data)  %>%
  mutate(
    h2 = map(herit, function(x) {
      x["heritability"][[1]][1]
    }),
    lower = map(herit, function(x) {
      x["heritability"][[1]][3]
    }),
  ) %>% 
  unnest(h2, lower) %>%
  mutate(
    phenoClass = "aseg",
    Type = "Twin",
  ) %>%
  mutate(
    signif = lower > 0
  ) %>%
  select(phenotype, h2, lower, signif)

write_csv(twinVol, "formattedData/twinVol.csv")


