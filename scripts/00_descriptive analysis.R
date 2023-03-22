# 00_descriptive statistics

# traits only RDS
datat <- readRDS("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/synurbic and traits only.rds")

# Few family plots
ggplot(datat, aes(fam, fill = Synurbic)) + 
  #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
  geom_bar(position = position_dodge()) +
  coord_flip() +
  theme_minimal()

ggplot(datat, aes(fam, fill = Synurbic)) + 
  #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
  geom_bar(position = "Fill") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_minimal()

ggplot(datat, aes(fam, fill = Synurbic)) + 
  #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
  geom_bar() +
  coord_flip() +
  theme_minimal()


## DB adding family-level sample size and some ideas

## fix family
library(tools)
datat$family = toTitleCase(tolower(datat$fam))

## tabulate
n_fam = data.frame(table(datat$family))
names(n_fam) = c("family","sample")

## derive coverage
tab = table(datat$family, is.na(datat$Synurbic))
tab = as.data.frame.matrix(tab)
tab$frac = tab$`FALSE` / (tab$`FALSE` + tab$`TRUE`)
tab$family = rownames(tab)
tab$missing = tab$`TRUE`

## merge
n_fam = merge(n_fam, tab[c("family", "frac", "missing")], by = "family")
rm(tab)

## reorder
n_fam = n_fam[order(n_fam$missing, n_fam$sample, decreasing=T),]
n_fam$family = factor(n_fam$family, levels = as.character(n_fam$family))

## set levels in datat
datat$family = factor(datat$family, levels = levels(n_fam$family))

## plot
ggplot(datat, aes(family)) + 
  geom_bar(aes(fill = Synurbic), position = "Fill") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "coverage") +
  theme(legend.position = "top") +
  geom_text(data = n_fam,
            aes(x = family,
                y = 1.05,
                label = sample), size = 3)





