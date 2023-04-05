# 00_descriptive statistics

# clear
rm(list=ls()) 
graphics.off()

# pachages
library(tidyverse)
library(tools)

# traits only RDS
datat <- readRDS("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/synurbic and traits only.rds")

## simple family tweak for now
#datat$fam = ifelse(datat$gen == "Miniopterus", "MINIOPTERIDAE", datat$fam)

# filter out the additional extinct bats in IUCN
names <- datat[!is.na(datat$category) & datat$category == "EX", ]$species

# remove the extinct bats
datat <- filter(datat, !species %in% names)
rm(names)

# uppercase
datat$family = toTitleCase(tolower(datat$fam))

# calculate families
fam_tab <- datat %>% 
  group_by(family) %>%
  count() %>%
  rename(sample = n)

# pull missing data
fam_na <- datat %>% 
  group_by(family) %>% 
  count(is.na(Synurbic)) %>% 
  pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
  replace(is.na(.), 0) %>%
  rename(nonmissing = "FALSE", 
         missing = "TRUE")

# merge together by family 
fam_cov <- merge(fam_tab, fam_na, by = "family", all.x = TRUE)

fam_cov <- fam_cov %>% 
  mutate(frac = nonmissing/sample) %>%
  select(-nonmissing)

# add total label
df <- data.frame(x = 20, y = 1.05, text = "Total spp.") 

png("/Users/brianabetke/Desktop/coverage.png",width=10,height=8,units="in",res=600)
ggplot(datat, aes(reorder(family, is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
  theme(legend.position = "top") +
  scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  geom_text(data = fam_cov,
            aes(x = family,
                y = 1.05,
                label = sample), size = 3) +
  geom_text(data = df, aes(x = x, y = y, label = text), size = 3, vjust = "inward")
dev.off()

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
n_fam = n_fam[order(n_fam$frac, decreasing=F),]
n_fam$family = factor(n_fam$family, levels = as.character(n_fam$family))

## set levels in datat
datat$family = factor(datat$family, levels = levels(n_fam$family))

## plot
ggplot(datat, aes(family)) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "coverage") +
  theme(legend.position = "top") +
  geom_text(data = n_fam,
            aes(x = family,
                y = 1.05,
                label = sample), size = 3)

## Category - change colors and match the family break down, stat_count works but legend position didn't change
# IUCN data 
datat %>% 
  group_by(category) %>%
  count() %>%
  rename(sample = n)

# datat %>% 
#   group_by(!is.na(category)) %>% 
#   count(is.na(Synurbic)) %>% 
#   pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
#   replace(is.na(.), 0) %>%
#   rename(nonmissing = "FALSE", 
#          missing = "TRUE")

png("/Users/brianabetke/Desktop/IUCNcoverage.png",width=10,height=8,units="in",res=600)
ggplot(datat, aes(reorder(category, is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  coord_flip() +
  stat_count(geom = "text", 
             aes(label = after_stat(count)),
             position=position_fill(vjust=1.05), colour="black") +
  labs(x = NULL, y = "coverage") +
  theme(legend.position = "top") +
  theme_minimal() 
dev.off()

library(DataExplorer)
library(gtsummary)

datat %>% select(category, Synurbic) %>% tbl_summary(by = category)
datat %>% select(Synurbic, Indomalayan:Nearctic) %>% tbl_summary(by = Synurbic)

# population trends
ggplot(datat, aes(population_trend)) +
  geom_bar(aes(fill = Synurbic), position = "Fill") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "coverage") +
  theme(legend.position = "top")

