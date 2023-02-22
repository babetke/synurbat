library(tidyverse)

no_NAs <- read_csv("/Volumes/BETKE 2021/synurbat/flat files/grid search without NAs.csv",
                   col_types = cols(interaction.depth = col_factor(), shrinkage = col_factor()))


AUC_gg <- ggplot(no_NAs, aes(x = shrinkage, y = testAUC)) +
  geom_boxplot(aes(fill = interaction.depth), color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Shrinkage", y = "AUC", fill = "Interaction Depth")

spec_gg <- ggplot(no_NAs, aes(x = shrinkage, y = spec)) +
  geom_boxplot(aes(fill = interaction.depth), color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Shrinkage", y = "Specificity", fill = "Interaction Depth")

sen_gg <- ggplot(no_NAs, aes(x = shrinkage, y = sen)) +
  geom_boxplot(aes(fill = interaction.depth), color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Shrinkage", y = "Sensitivity", fill = "Interaction Depth")

library(patchwork)
AUC_gg + spec_gg + sen_gg

