# Data visualization for BRT analysis
# babetke@utexas.edu

# packages
library(tidyverse)
library(patchwork)
library(gbm)

# read in rds of model results 
noNA <- readRDS("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/noNA_brts.rds")
pseudo <- readRDS("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/pseudo_brts.rds")

#### Variable importance plots
inf_plot <- function(rds_name1, rds_name2, bar_color){

  # pull rinf
  rds_vars1 <- rds_name1[["rinf"]]
  rds_vars2 <- rds_name2[["rinf"]]
  
  # merge by variable name
  
  rds_vars <- merge(rds_name1, rds_name2, by = "var")
  
  # Clean up variable names
  rds_vars$var <- recode(rds_vars$var,
                        "cites" = "Citation Count",
                        "X26.1_GR_Area_km2" = "Geographic Area",
                        "X30.1_AET_Mean_mm" = "Mean Monthly AET",
                        "X26.2_GR_MaxLat_dd" = "Maximum Latitude",
                        "X26.3_GR_MinLat_dd" = "Minimum Latitude",
                        "habitat_breadth_n" = "Habitat Breadth",
                        "litters_per_year_n" = "Litters Per Year",
                        "adult_body_length_mm" = "Adult Body Length",
                        "X28.2_Temp_Mean_01degC" = "Mean Monthly Temperature",
                        "X27.2_HuPopDen_Mean_n.km2" = "Mean Human Population Density",
                        "X28.1_Precip_Mean_mm" = "Mean Monthly Precipitation",
                        "litter_size_n" = "Litter Size",
                        "upper_elevation_m" = "Upper Elevation",
                        "disected_by_mountains" = "Disected by Mountains",
                        "adult_forearm_length_mm" = "Adult Forearm Length",
                        "altitude_breadth_m" = "Altitude Breadth",
                        "X26.4_GR_MidRangeLat_dd" = "Median Latitudinal Range",
                        "foraging_stratum" = "Foraging stratum",
                        "adult_mass_g" = "Adult Mass",
                        "X30.2_PET_Mean_mm" = "Mean Monthly PET",
                        "det_vfish" = "Diet Fish",
                        "X26.5_GR_MaxLong_dd" = "Maximum Longitude",
                        "fam_RHINOLOPHIDAE" = "Rhinolophidae",
                        "det_diet_breadth_n" = "Diet Breadth",
                        "X26.6_GR_MinLong_dd" = "Minimum Longitude",
                        "det_vend" = "Diet Vend",
                        "X27.1_HuPopDen_Min_n.km2" = "Min Human Population Density",
                        "dphy_vertebrate" = "Diet Vertebrate",
                        "lower_elevation_m" = "Lower Elevation Limit",
                        "det_nect" = "Diet Nectar",
                        "X27.4_HuPopDen_Change" = "Human Population Density Change",
                        "X27.3_HuPopDen_5p_n.km2" = "Human Population Density 5th Percentile",
                        "X26.7_GR_MidRangeLong_dd" = "Median Longitudinal Range",
                        "det_fruit" = "Diet Fruit",
                        "fam_PHYLLOSTOMIDAE" = "Phyllostomidae",
                        "det_vect" = "Diet Vect",
                        "trophic_level" = "Trophic Level",
                        "dphy_invertebrate" = "Diet Invertebrate (dypy)",
                        "island_dwelling" = "Island Dwelling",
                        "fam_MOLOSSIDAE" = "Molossidae",
                        "fam_HIPPOSIDERIDAE" = "Hipposideridae",
                        "dphy_plant" = "Diet Plants",
                        "glaciation" = "Glaciation",
                        "fam_VESPERTILIONIDAE" = "Vespertilionidae",
                        "fam_EMBALLONURIDAE" = "Emballonuridae",
                        "fam_PTEROPODIDAE" = "Pteropodidae",
                        "activity_cycle" = "Activity Cycle",
                        "det_seed" = "Diet Seeds",
                        "fam_MORMOOPIDAE" = "Mormoopidae",
                        "fam_NATALIDAE" = "Natalidae",
                        "fam_NYCTERIDAE" = "Nycteridae"
  )
  
  ggplot(rds_vars, aes(x = reorder(var, rel.inf), y = rel.inf)) + 
    #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
    geom_bar(stat = "identity", fill = bar_color) +
    theme(axis.text.x = element_text(size = 10, hjust = 1),
          legend.position = "none") +
    theme_bw() +
    labs(x = " ", y = "Relative Importance") +
    coord_flip()
  
}

png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/variable importance no NAs.png", width=8,height=10,units="in",res=600)
inf_plot(noNA,"grey")
dev.off()

noNA_gg <- inf_plot(noNA,"dark grey")
pseudo_gg <- inf_plot(pseudo, "purple")

noNA_gg + pseudo_gg

na_inf <- pseudo


#### marginal effect plots (separate for now) added rugs but thats not really what I want.
# basic plot function for now
make_pdp <- function(model, type, predictor) {
  
  # return grid
  vals <- plot.gbm(model[["mod"]], i.var = predictor, type = "response", return.grid = TRUE)
  
  # if numeric, run geom line, if not, point?
  if(type == "num"){
    
    fig <- ggplot(vals, aes(!!sym(predictor), y)) + 
      geom_line() +
      geom_rug() +
      labs(x = predictor, y = "Marginal Effect") +
      theme_bw()
    
  } else {
    
    fig <- ggplot(vals, aes(!!sym(predictor), y)) + 
      geom_point(shape=23, fill="blue", size=3) +
      geom_rug() +
      labs(x = predictor, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 10, hjust = 1),
            legend.position = "none")
    
  }

}

hb <- make_pdp(model = noNA, type = "num", predictor = "habitat_breadth_n")
gr <- make_pdp(model = noNA, type = "num", predictor = "X26.1_GR_Area_km2")
ls <- make_pdp(model = noNA, type = "num", predictor = "litter_size_n")
pm <- make_pdp(model = noNA, type = "num", predictor = "X28.1_Precip_Mean_mm")
ml <- make_pdp(model = noNA, type = "num", predictor = "X26.2_GR_MaxLat_dd")
fr <- make_pdp(model = noNA, type = "num", predictor = "det_fruit")

png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/pdp 6 no NAs.png", width=10,height=5,units="in",res=600)
hb + gr + ls + pm + ml + fr + plot_layout(nrow = 2, ncol = 3, byrow = F)
dev.off()

# pseudoabsesnces
hbp <- make_pdp(model = pseudo, type = "num", predictor = "habitat_breadth_n")
grp <- make_pdp(model = pseudo, type = "num", predictor = "X26.1_GR_Area_km2")
lsp <- make_pdp(model = pseudo, type = "num", predictor = "litter_size_n")
mlp <- make_pdp(model = pseudo, type = "num", predictor = "X26.2_GR_MaxLat_dd")
pmp <- make_pdp(model = pseudo, type = "num", predictor = "X28.1_Precip_Mean_mm")
ccp <- make_pdp(model = pseudo, type = "num", predictor = "cites")

png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/pdp 6 pseudo.png", width=10,height=5,units="in",res=600)
hbp + grp + lsp + mlp + pmp + ccp + plot_layout(nrow = 2, ncol = 3, byrow = F)
dev.off()



# try overlapping bars to compare var infl
# pull rinf
rds_vars1 <- noNA[["rinf"]]
rds_vars2 <- pseudo[["rinf"]]

# merge by variable name

rds_vars <- merge(rds_vars1, rds_vars2, by = "var")

rds_vars$var <- recode(rds_vars$var,
                       "cites" = "Citation Count",
                       "X26.1_GR_Area_km2" = "Geographic Area",
                       "X30.1_AET_Mean_mm" = "Mean Monthly AET",
                       "X26.2_GR_MaxLat_dd" = "Maximum Latitude",
                       "X26.3_GR_MinLat_dd" = "Minimum Latitude",
                       "habitat_breadth_n" = "Habitat Breadth",
                       "litters_per_year_n" = "Litters Per Year",
                       "adult_body_length_mm" = "Adult Body Length",
                       "X28.2_Temp_Mean_01degC" = "Mean Monthly Temperature",
                       "X27.2_HuPopDen_Mean_n.km2" = "Mean Human Population Density",
                       "X28.1_Precip_Mean_mm" = "Mean Monthly Precipitation",
                       "litter_size_n" = "Litter Size",
                       "upper_elevation_m" = "Upper Elevation",
                       "disected_by_mountains" = "Disected by Mountains",
                       "adult_forearm_length_mm" = "Adult Forearm Length",
                       "altitude_breadth_m" = "Altitude Breadth",
                       "X26.4_GR_MidRangeLat_dd" = "Median Latitudinal Range",
                       "foraging_stratum" = "Foraging stratum",
                       "adult_mass_g" = "Adult Mass",
                       "X30.2_PET_Mean_mm" = "Mean Monthly PET",
                       "det_vfish" = "Diet Fish",
                       "X26.5_GR_MaxLong_dd" = "Maximum Longitude",
                       "fam_RHINOLOPHIDAE" = "Rhinolophidae",
                       "det_diet_breadth_n" = "Diet Breadth",
                       "X26.6_GR_MinLong_dd" = "Minimum Longitude",
                       "det_vend" = "Diet Vend",
                       "X27.1_HuPopDen_Min_n.km2" = "Min Human Population Density",
                       "dphy_vertebrate" = "Diet Vertebrate",
                       "lower_elevation_m" = "Lower Elevation Limit",
                       "det_nect" = "Diet Nectar",
                       "X27.4_HuPopDen_Change" = "Human Population Density Change",
                       "X27.3_HuPopDen_5p_n.km2" = "Human Population Density 5th Percentile",
                       "X26.7_GR_MidRangeLong_dd" = "Median Longitudinal Range",
                       "det_fruit" = "Diet Fruit",
                       "fam_PHYLLOSTOMIDAE" = "Phyllostomidae",
                       "det_vect" = "Diet Vect",
                       "trophic_level" = "Trophic Level",
                       "dphy_invertebrate" = "Diet Invertebrate (dypy)",
                       "island_dwelling" = "Island Dwelling",
                       "fam_MOLOSSIDAE" = "Molossidae",
                       "fam_HIPPOSIDERIDAE" = "Hipposideridae",
                       "dphy_plant" = "Diet Plants",
                       "glaciation" = "Glaciation",
                       "fam_VESPERTILIONIDAE" = "Vespertilionidae",
                       "fam_EMBALLONURIDAE" = "Emballonuridae",
                       "fam_PTEROPODIDAE" = "Pteropodidae",
                       "activity_cycle" = "Activity Cycle",
                       "det_seed" = "Diet Seeds",
                       "fam_MORMOOPIDAE" = "Mormoopidae",
                       "fam_NATALIDAE" = "Natalidae",
                       "fam_NYCTERIDAE" = "Nycteridae"
)

ggplot(rds_vars) + 
  #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
  geom_bar(aes(x = reorder(var, rel.inf.x), y = rel.inf.y), stat = "identity", fill = "purple",alpha = 0.5) +
  geom_bar(aes(x = reorder(var, rel.inf.x), y = rel.inf.x), stat = "identity", alpha = 0.75) +
  theme(axis.text.x = element_text(size = 10, hjust = 1),
        legend.position = "none") +
  theme_bw() +
  labs(x = " ", y = "Relative Importance") +
  coord_flip()
