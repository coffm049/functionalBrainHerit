library(tidyverse)


fs = dir("twinEstResults", full.names= T, pattern = "*.Rds") %>%
  map_dfr(readRDS) %>%
  select(herit)  %>%
  mutate(
    A = map(herit, function(x) x$coef[1,1]),
    C = map(herit, function(x) x$coef[2,1]),
    E = map(herit, function(x) x$coef[3,1]),
  )  %>%
  unnest(A, C, E) %>%
  mutate(h2 = A / (A + C + E))

fs %>%
  select(g, common, e) %>%
  cor()

cor(fs$g, fs$common)


mean(fs$common)
sd(fs$common)


p <- fs %>%
  arrange(desc(A), desc(C), desc(E)) %>%
  mutate(trait = 1:nrow(.)) %>%
  pivot_longer(c(A,C, E)) %>%
  rename(Source = name, Proportion = value) %>%
  ggplot(aes(x = trait, y = Proportion, fill = Source)) +
  geom_area() +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

ggsave("~/GordoncompositionPlot.png")


fs %>%
  reframe(range = quantile(h2, c(0.05, 0.95)))

write_csv(fs, "~/gordonComps.csv")

