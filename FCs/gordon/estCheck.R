library(tidyverse)
library(txtplot)

df <- readRDS("herit_1.Rds")


# [ ] extract p, G, E, SE
# [ ] plot distributions



df$herit[[1]]$all %>% as.data.frame() %>% rownames_to_column("partition") %>% pivot_wider(names_from = partition, values_from = Estimate)
df$herit[[1]]$heritability
df$herit[[1]]
# loop over the above code and form a dataframe

df2 <- map(df$herit, \(ddf) {
    ddf$all %>% as.data.frame() %>% rownames_to_column("partition") %>% pivot_wider(names_from = partition, values_from = Estimate)
}) %>%
  reduce(rbind) %>% 
  as.data.frame()


txtdensity(na.omit(df2$A))
txtdensity(na.omit(df2$C))
txtdensity(na.omit(df2$E))
txtdensity(na.omit(df2$`Broad-sense heritability`))
txtdensity(na.omit(df2$Std.Err))
txtdensity(na.omit(df2$`2.5%`))

df2 %>% arrange(desc(`Broad-sense heritability`)) %>% head()

df2 %>%
    reframe(mean(`2.5%` > 0))


df <- read_csv("/users/4/coffm049/papers/functionalBrainHerit/results/FCs/gordon/pconns.GCTA.FE.3.csv")

txtdensity(na.omit(as.numeric(df$h2)))
txtdensity(na.omit(as.numeric(df$G)))
txtdensity(na.omit(as.numeric(df$E)))
txtdensity(na.omit(as.numeric(df$pval)))
