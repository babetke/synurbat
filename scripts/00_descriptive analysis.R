# 00_descriptive statistics

# clear
rm(list=ls()) 
graphics.off()

# pachages
library(tidyverse)
library(tools)
library(patchwork)

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
df <- data.frame(x = 20, y = 1.05, text = "no. spp.") 

png("/Users/brianabetke/Desktop/coverage.png",width=10,height=8,units="in",res=600)
fam_gg <- ggplot(datat, aes(reorder(family, is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  geom_text(data = fam_cov,
            aes(x = family,
                y = 1.05,
                label = sample), size = 1.5) +
  geom_text(data = df, aes(x = x, y = y, label = text), size = 1.5, vjust = "inward")
dev.off()

fam_gg <- fam_gg + theme(axis.text.y=element_text(size=6),
                       axis.text.x=element_text(size=6),
                       axis.title=element_text(size=7))


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
iucn_gg <- ggplot(datat, aes(reorder(category, is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  coord_flip() +
  stat_count(geom = "text", 
             aes(label = after_stat(count)),
             position=position_fill(vjust=1.05), colour="black") +
  labs(x = NULL, y = "coverage") +
  theme(legend.position = "top") +
  theme_minimal() 
dev.off()

cat_tab <- datat %>% group_by(factor(category)) %>% count() %>% rename("category"="factor(category)")
df3 <- data.frame(x = 7.75, y = 1.05, text = "no. spp.") 
factor(datat$category, exclude = NULL) -> datat$category

iucn_gg <- ggplot(datat, aes(reorder(category, is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  geom_text(data = cat_tab,
            aes(x = category,
                y = 1.05,
                label = n), size = 1.25) +
  geom_text(data = df3, aes(x = x, y = y, label = text), size = 1.25, vjust = "inward")

iucn_gg <- iucn_gg + theme(axis.text.y=element_text(size=6),
                           axis.text.x=element_text(size=6),
                           axis.title=element_text(size=7))

# Biogeographical realms
br_tab <- datat %>%
  select(species,biogeographical_realm) %>% # same results without selecting by the way....
  separate_rows(biogeographical_realm, sep = ", ") %>%
  group_by(biogeographical_realm) %>%
  count() %>%
  rename(sample = n)

br_na <- datat %>%
  select(species,biogeographical_realm,Synurbic) %>%
  separate_rows(biogeographical_realm, sep = ", ") %>%
  group_by(biogeographical_realm) %>%
  count(is.na(Synurbic)) %>%
  pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
  rename(nonmissing = "FALSE", 
         missing = "TRUE")

# merge together by family 
br_cov <- merge(br_tab, br_na, by = "biogeographical_realm", all.x = TRUE)

br_cov <- br_cov %>% 
  mutate(frac = nonmissing/sample) %>%
  select(-nonmissing)

# add total label
df2 <- data.frame(x = 8.75, y = 1.05, text = "no. spp.") 

br_gg <- separate_rows(datat, biogeographical_realm, sep = ", ") %>%
ggplot(aes(reorder(biogeographical_realm, is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  geom_text(data = br_cov,
            aes(x = biogeographical_realm,
                y = 1.05,
                label = sample), size = 1.25) +
  geom_text(data = df2, aes(x = x, y = y, label = text), size = 1.25, vjust = "inward")

br_gg <- br_gg + theme(axis.text.y=element_text(size=6),
                           axis.text.x=element_text(size=6),
                           axis.title=element_text(size=7))

# multiplot, need to have A-C lables for pub
fig_1 <- (fam_gg + plot_layout(widths = c(1.75, 1.25), guide = "collect") & 
            theme(legend.position = "bottom", legend.text = element_text(size = 6), 
                  legend.title = element_text(size=6), legend.key.size = unit(0.50, 'cm'), plot.tag = element_text(size = 8))) + 
  (br_gg / iucn_gg & theme(legend.position = "none", plot.tag = element_text(size = 8)))


png("/Users/brianabetke/Desktop/fig1.png",width=7,height=5,units="in",res=300)
fig_1 + plot_annotation(tag_levels = 'A')
dev.off()


library(ggpubr)


# Table/data summaries
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


# attempting to balloon plot...
f <- data.frame(table(datat$category,datat$Synurbic))
names(f) <- c("category", "Synurbic", "Freq")

ggballoonplot(n, fill = "Synurbic", size = "Freq", show.label = TRUE)

n = f[order(f$Freq, decreasing=F),]

d <- datat %>% 
  group_by(category) %>% 
  count(Synurbic) %>%
  arrange(n)
  
ggplot(f, aes(y = reorder(category,Freq) , x = Synurbic , size = Freq)) +
  geom_point()



