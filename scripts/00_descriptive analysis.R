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




