# 00_descriptive statistics
# babetke@utexas.edu

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

# # filter out the additional extinct bats in IUCN -> moved to data cleaning file.
# names <- datat[!is.na(datat$category) & datat$category == "EX", ]$species
# 
# # remove the extinct bats
# datat <- filter(datat, !species %in% names)
# rm(names)

###### Alternative 1: plots with values next to labels and no no. spp. col
# families
# uppercase
datat$family <- toTitleCase(tolower(datat$fam))

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

# calculate coverage
fam_cov <- fam_cov %>% 
  mutate(frac = nonmissing/sample) %>%
  select(-nonmissing)

# add values and family names together for plot, sort in descending order of coverage
label_fam <- fam_cov %>% 
  mutate(with(fam_cov, paste(family, "(", sample, ")"))) %>% 
  rename() %>%
  arrange(frac) 

# plot family coverage
fam_gg <- ggplot(datat, aes(reorder(family, -is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Percent of Species", fill = "Anthropogenic Roosting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("no","yes","unkown"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_x_discrete(labels = label_fam$`with(fam_cov, paste(family, "(", sample, ")"))`)

# format text
fam_gg <- fam_gg + theme(axis.text.y=element_text(size=6),
                         axis.text.x=element_text(size=6),
                         axis.title=element_text(size=7))

# percentages of synurbic
datat %>% 
  group_by(family) %>%
  count(Synurbic) %>%
  rename(sample = n) %>%
  pivot_wider(names_from = Synurbic, values_from = sample) %>%
  ungroup() %>%
  rename(no = `0`, yes = `1`, missing = `NA`) %>%
  mutate(across(everything(), replace_na, 0)) -> fam_syn
  
# merge with family tab
syn_cov <- merge(fam_tab, fam_syn, by = "family")

syn_cov <- syn_cov %>% 
  mutate(synprop = yes/sample) %>%
  arrange(desc(synprop))

## geographic realms
# separate br categories
gr_dat <- datat %>%
  select(species,biogeographical_realm,Synurbic) %>%
  separate_rows(biogeographical_realm, sep = ", ")

# change NA to no data
gr_dat$fbr <- ifelse(is.na(gr_dat$biogeographical_realm),"No Data", gr_dat$biogeographical_realm) 

# count 
br_tab <- gr_dat %>%
  group_by(fbr) %>%
  count() %>%
  rename(sample = n)

# missing
br_na <- gr_dat %>%
  group_by(fbr) %>%
  count(is.na(Synurbic)) %>%
  pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
  rename(nonmissing = "FALSE", 
         missing = "TRUE")

# merge together by brtest 
br_cov <- merge(br_tab, br_na, by = "fbr", all.x = TRUE)

# calculate coverage
br_cov <- br_cov %>% 
  mutate(frac = nonmissing/sample) %>%
  select(-nonmissing)

# create axi tick labels
label_gr <- br_cov %>% 
  mutate(with(br_cov, paste(fbr, "(", sample, ")"))) %>% 
  rename() %>%
  arrange(frac) 

# plot coverage
br_gg <- ggplot(gr_dat, aes(reorder(fbr, -is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Percent of Species", fill = "Anthropogenic Roosting") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("no","yes","unknown"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_x_discrete(labels = label_gr$`with(br_cov, paste(fbr, "(", sample, ")"))`)

br_gg <- br_gg + theme(axis.text.y=element_text(size=6),
                       axis.text.x=element_text(size=6),
                       axis.title=element_text(size=7))

# percentages of synurbic
gr_dat %>% 
  group_by(fbr) %>% 
  count(Synurbic) %>% 
  pivot_wider(names_from = Synurbic, values_from = n) %>%
  ungroup() %>%
  rename(no = `0`, yes = `1`, missing = `NA`) -> syn_gr

# merge with family tab
syn_gr <- merge(br_tab, syn_gr, by = "fbr")

syn_gr <- syn_gr %>% 
  mutate(synprop = yes/sample) %>%
  arrange(desc(synprop))

## IUCN
# have values recoded for IUCN
datat$fcat <- recode(datat$category, "CR" = "Critically Endangered",
                     "DD" =  "Data Deficient", "EN" = "Endangered",
                     "LC" = "Least Concern", "NT" =  "Near Threatened",
                     "VU" = "Vulnerable")

# to get Not Evaluated category
datat$fcat <- as.character(datat$fcat)
datat$fcat <- ifelse(is.na(datat$fcat),"No Data", datat$fcat) 

# table
cat_tab <- datat %>% 
  group_by(factor(fcat)) %>% 
  count() %>% 
  rename(sample = n, fcat = `factor(fcat)`) 

# pull missing
cat_na <- datat %>% 
  group_by(fcat) %>% 
  count(is.na(Synurbic)) %>% 
  pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
  rename(nonmissing = "FALSE", missing = "TRUE")

# merge together by new category
cat_cov <- merge(cat_tab, cat_na, by = "fcat", all.x = TRUE)

# calculate coverage
cat_cov <- cat_cov %>% 
  mutate(frac = nonmissing/sample) %>%
  select(-nonmissing)

# create axis tick labels
label_cat <- cat_cov %>% 
  mutate(with(cat_cov, paste(fcat, "(", sample, ")"))) %>% 
  rename() %>%
  arrange(frac)

# plot coverage
iucn_gg <- ggplot(datat, aes(reorder(fcat, -is.na(Synurbic)))) + 
  geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = "Percent of Species", fill = "Anthropogenic Roosting") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(labels = c("no","yes","unknown"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_x_discrete(labels = label_cat$`with(cat_cov, paste(fcat, "(", sample, ")"))`)

iucn_gg <- iucn_gg + theme(axis.text.y=element_text(size=6),
                           axis.text.x=element_text(size=6),
                           axis.title=element_text(size=7))

# multiplot
fig1 <- (fam_gg + plot_layout(widths = c(1.75, 1.25), guide = "collect") & 
            theme(legend.position = "bottom", legend.text = element_text(size = 6), 
                  legend.title = element_text(size=6), legend.key.size = unit(0.50, 'cm'), plot.tag = element_text(size = 8))) + 
  (br_gg / iucn_gg & theme(legend.position = "none", plot.tag = element_text(size = 8)))

png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/figure 1.png",width=7,height=4,units="in",res=300)
fig1 + plot_annotation(tag_levels = 'A')
dev.off()

# percentages of synurbic
datat %>% 
  group_by(fcat) %>%
  count(Synurbic) %>%
  rename(sample = n) %>%
  pivot_wider(names_from = Synurbic, values_from = sample) %>%
  ungroup() %>%
  rename(no = `0`, yes = `1`, missing = `NA`) -> fcat_syn

# merge with family tab
syn_iucn <- merge(cat_tab, fcat_syn, by = "fcat")

syn_iucn <- syn_iucn %>% 
  mutate(synprop = yes/sample) %>%
  arrange(desc(synprop))

# # population trends
# ggplot(datat, aes(population_trend)) +
#   geom_bar(aes(fill = Synurbic), position = "Fill") +
#   scale_y_continuous(labels = scales::percent) +
#   coord_flip() +
#   theme_minimal() +
#   labs(x = NULL, y = "coverage") +
#   theme(legend.position = "top")

# ## DB adding family-level sample size and some ideas
# 
# ## fix family
# library(tools)
# datat$family = toTitleCase(tolower(datat$fam))
# 
# ## tabulate
# n_fam = data.frame(table(datat$family))
# names(n_fam) = c("family","sample")
# 
# ## derive coverage
# tab = table(datat$family, is.na(datat$Synurbic))
# tab = as.data.frame.matrix(tab)
# tab$frac = tab$`FALSE` / (tab$`FALSE` + tab$`TRUE`)
# tab$family = rownames(tab)
# tab$missing = tab$`TRUE`
# 
# ## merge
# n_fam = merge(n_fam, tab[c("family", "frac", "missing")], by = "family")
# rm(tab)
# 
# ## reorder
# n_fam = n_fam[order(n_fam$frac, decreasing=F),]
# n_fam$family = factor(n_fam$family, levels = as.character(n_fam$family))
# 
# ## set levels in datat
# datat$family = factor(datat$family, levels = levels(n_fam$family))
# 
# ## plot
# ggplot(datat, aes(family)) +
#   geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
#   scale_y_continuous(labels = scales::percent) +
#   coord_flip() +
#   theme_minimal() +
#   labs(x = NULL, y = "coverage") +
#   theme(legend.position = "top") +
#   geom_text(data = n_fam,
#             aes(x = family,
#                 y = 1.05,
#                 label = sample), size = 3)

# ## Initial plots
# ## family plot
# # uppercase
# datat$family <- toTitleCase(tolower(datat$fam))
# 
# # calculate families
# fam_tab <- datat %>% 
#   group_by(family) %>%
#   count() %>%
#   rename(sample = n)
# 
# # pull missing data
# fam_na <- datat %>% 
#   group_by(family) %>% 
#   count(is.na(Synurbic)) %>% 
#   pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
#   replace(is.na(.), 0) %>%
#   rename(nonmissing = "FALSE", 
#          missing = "TRUE")
# 
# # merge together by family 
# fam_cov <- merge(fam_tab, fam_na, by = "family", all.x = TRUE)
# 
# fam_cov <- fam_cov %>% 
#   mutate(frac = nonmissing/sample) %>%
#   select(-nonmissing)
# 
# # add total label
# df <- data.frame(x = 20, y = 1.05, text = "no. spp.") 
# 
# fam_gg <- ggplot(datat, aes(reorder(family, is.na(Synurbic)))) + 
#   geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
#   scale_y_continuous(labels = scales::percent) +
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
#   geom_text(data = fam_cov,
#             aes(x = family,
#                 y = 1.05,
#                 label = sample), size = 1.5) +
#   geom_text(data = df, aes(x = x, y = y, label = text), size = 1.5, vjust = "inward")
# 
# fam_gg <- fam_gg + theme(axis.text.y=element_text(size=6),
#                          axis.text.x=element_text(size=6),
#                          axis.title=element_text(size=7))
# 
# ## Category - change colors and match the family break down, stat_count works but legend position didn't change
# # IUCN data 
# 
# # have values recoded for IUCN
# datat$fcat <- recode(datat$category, "CR" = "Critically Endangered",
#                      "DD" =  "Data Deficient", "EN" = "Endangered",
#                      "LC" = "Least Concern", "NT" =  "Near Threatened",
#                      "VU" = "Vulnerable")
# 
# # to get Not Evaluated category
# datat$fcat <- as.character(datat$fcat)
# datat$fcat <- ifelse(is.na(datat$fcat),"No Data", datat$fcat) 
# 
# # table
# cat_tab <- datat %>% 
#   group_by(factor(fcat)) %>% 
#   count() %>% 
#   rename(sample = n)
# 
# # pull missing
# cat_na <- datat %>% 
#   group_by(fcat) %>% 
#   count(is.na(Synurbic)) %>% 
#   pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
#   rename(nonmissing = "FALSE", missing = "TRUE")
# 
# # merge together by family 
# cat_cov <- merge(cat_tab, cat_na, by = "fcat", all.x = TRUE)
# 
# # calculate coverage
# cat_cov <- cat_cov %>% 
#   mutate(frac = nonmissing/sample) %>%
#   select(-nonmissing)
# 
# # label for species numbers
# df3 <- data.frame(x = 7.75, y = 1.05, text = "no. spp.") 
# factor(datat$fcat, exclude = NULL) -> datat$fcat
# 
# iucn_gg <- ggplot(datat, aes(reorder(fcat, is.na(Synurbic)))) + 
#   geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
#   scale_y_continuous(labels = scales::percent) +
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
#   theme(legend.position = "none",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
#   geom_text(data = cat_tab,
#             aes(x = fcat,
#                 y = 1.05,
#                 label = sample), size = 1.25) +
#   geom_text(data = df3, aes(x = x, y = y, label = text), size = 1.25, vjust = "inward")
# 
# iucn_gg <- iucn_gg + theme(axis.text.y=element_text(size=6),
#                            axis.text.x=element_text(size=6),
#                            axis.title=element_text(size=7))
# 
# ## Biogeographical realms
# # separate br categories
# gr_dat <- datat %>%
#   select(species,biogeographical_realm,Synurbic) %>%
#   separate_rows(biogeographical_realm, sep = ", ")
# 
# # change to no data
# gr_dat$brtest <- ifelse(is.na(gr_dat$biogeographical_realm),"No Data", gr_dat$biogeographical_realm) 
# 
# # count 
# br_tab <- gr_dat %>%
#   group_by(brtest) %>%
#   count() %>%
#   rename(sample = n)
# 
# # missing
# br_na <- gr_dat %>%
#   group_by(brtest) %>%
#   count(is.na(Synurbic)) %>%
#   pivot_wider(names_from = "is.na(Synurbic)", values_from = n) %>%
#   rename(nonmissing = "FALSE", 
#          missing = "TRUE")
# 
# # merge together by family 
# br_cov <- merge(br_tab, br_na, by = "brtest", all.x = TRUE)
# 
# br_cov <- br_cov %>% 
#   mutate(frac = nonmissing/sample) %>%
#   select(-nonmissing)
# 
# # add total label
# df2 <- data.frame(x = 8.75, y = 1.05, text = "no. spp.") 
# 
# # br_gg <- separate_rows(datat, biogeographical_realm, sep = ", ") %>%
# # ggplot(aes(reorder(biogeographical_realm, is.na(Synurbic)))) + 
# #   geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
# #   scale_y_continuous(labels = scales::percent) +
# #   coord_flip() +
# #   theme_bw() +
# #   labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
# #   theme(legend.position = "none",
# #         panel.grid.major = element_blank(),
# #         panel.grid.minor = element_blank()) +
# #   scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
# #   geom_text(data = br_cov,
# #             aes(x = brtest,
# #                 y = 1.05,
# #                 label = sample), size = 1.25) +
# #   geom_text(data = df2, aes(x = x, y = y, label = text), size = 1.25, vjust = "inward")
# 
# br_gg <- ggplot(gr_dat, aes(reorder(brtest, is.na(Synurbic)))) + 
#   geom_bar(aes(fill = Synurbic), position = position_fill(reverse = TRUE)) +
#   scale_y_continuous(labels = scales::percent) +
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = "coverage", fill = "Anthropogenic Roosting") +
#   theme(legend.position = "none",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
#   geom_text(data = br_cov,
#             aes(x = brtest,
#                 y = 1.05,
#                 label = sample), size = 1.25) +
#   geom_text(data = df2, aes(x = x, y = y, label = text), size = 1.25, vjust = "inward")
# 
# br_gg <- br_gg + theme(axis.text.y=element_text(size=6),
#                        axis.text.x=element_text(size=6),
#                        axis.title=element_text(size=7))
# 
# # multiplot, need to have A-C lables for pub
# fig_1 <- (fam_gg + plot_layout(widths = c(1.75, 1.25), guide = "collect") & 
#             theme(legend.position = "bottom", legend.text = element_text(size = 6), 
#                   legend.title = element_text(size=6), legend.key.size = unit(0.50, 'cm'), plot.tag = element_text(size = 8))) + 
#   (br_gg / iucn_gg & theme(legend.position = "none", plot.tag = element_text(size = 8)))
# 
# 
# png("/Users/brianabetke/Desktop/fig1.png",width=7,height=5,units="in",res=300)
# fig_1 + plot_annotation(tag_levels = 'A')
# dev.off()
