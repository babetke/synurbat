# Data visualization for BRT analysis
# babetke@utexas.edu

# clear environment
rm(list=ls()) 
graphics.off()

# packages
library(tidyverse)
library(patchwork)
library(gbm)
library(ggrepel)

# read in rds of model results 
noNA <- readRDS("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/noNA_brts.rds")
pseudo <- readRDS("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/pseudo_brts.rds")

#### Variable importance plots
inf_plot <- function(rds_name1, rds_name2 = NULL){

  if(is.null(rds_name2)){
    # extract variable importance
    rds_vars <- noNA[["rinf"]]
    
    # Create type variable
    rds_vars <- mutate(rds_vars, type = ifelse(startsWith(var, "X"), "Geographic", 
                                               ifelse(startsWith(var, "det"), "Forage", 
                                                      ifelse(var %in% c("Palearctic", "Neotropical", "Afrotropical","Indomalayan","Nearctic",
                                                                        "Oceanian","Australasian","glaciation", "habitat_breadth_n","altitude_breadth_m","disected_by_mountains"),"Geographic",
                                                             ifelse(var %in% c("dphy_invertebrate","dphy_plant","dphy_vertebrate","trophic_level","foraging_stratum","lower_elevation_m","upper_elevation_m"),"Forage",
                                                                    ifelse(startsWith(var, "fam"),"Phylogeny",
                                                                           ifelse(var %in% c("category","population_trend","cites"),"Other", "Life History")))))))
    
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
                           "upper_elevation_m" = "Upper Elevation Limit",
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
                           "fam_MINIOPTERIDAE" = "Miniopteridae",
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
                           "fam_NYCTERIDAE" = "Nycteridae",
                           "category" = "Conservation Status",
                           "population_trend" = "Population Trend"
    )
    
    # plot points
    ggplot(rds_vars) +
      geom_point(aes(x = reorder(var, rel.inf), y = rel.inf, color = type), stat = "identity", alpha = 0.75, size = 3) +
      theme(axis.text.x = element_text(size = 10, hjust = 1)) +
      theme_bw() +
      labs(x = " ", y = "Relative Importance") +
      coord_flip() +
      theme(legend.position = c(0.80, 0.15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
  } else {
  # pull rinf
  rds_vars1 <- rds_name1[["rinf"]]
  rds_vars2 <- rds_name2[["rinf"]]
  
  # merge by variable name
  rds_vars <- merge(rds_vars1, rds_vars2, by = "var")
  
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
                         "upper_elevation_m" = "Upper Elevation Limit",
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
                         "fam_MINIOPTERIDAE" = "Miniopteridae",
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
                         "fam_NYCTERIDAE" = "Nycteridae",
                         "category" = "Conservation Status",
                         "population_trend" = "Population Trend"
  )
  
  # update names
  rds_vars <- rds_vars %>% 
    rename(Initial = rel.inf.x,
           Pseudoabsence = rel.inf.y)
  
  # points
  ggplot(rds_vars) + 
    #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
    geom_point(aes(x = reorder(var, Initial), y = Pseudoabsence, color = "Pseudoabsence"), stat = "identity", alpha = 0.5, size = 3) +
    geom_point(aes(x = reorder(var, Initial), y = Initial, color = "Initial"), stat = "identity", alpha = 0.75, size = 2) +
    theme(axis.text.x = element_text(size = 10, hjust = 1)) +
    theme_bw() +
    labs(x = " ", y = "Relative Importance") +
    coord_flip() +
    scale_color_manual(values = c("Initial" = "black", "Pseudoabsence" = "orange"))
  
  # bars
  # ggplot(rds_vars) + 
  #   #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
  #   geom_bar(aes(x = reorder(var, rel.inf.x), y = rel.inf.y), stat = "identity", fill = "orange",alpha = 0.5) +
  #   geom_bar(aes(x = reorder(var, rel.inf.x), y = rel.inf.x), stat = "identity", alpha = 0.75) +
  #   theme(axis.text.x = element_text(size = 10, hjust = 1),
  #         legend.position = "none") +
  #   theme_bw() +
  #   labs(x = " ", y = "Relative Importance") +
  #   coord_flip()
  }
  
}

# No NA only 
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/variable importance no NAs.png", width=6, height=7.5, units="in", res=300)
inf_plot(noNA)
dev.off()

# both models
rel_gg <- inf_plot(noNA, pseudo)

#### Scatter plot
rds_vars1 <- noNA[["rinf"]]
rds_vars2 <- pseudo[["rinf"]]

# merge by variable name
rds_vars <- merge(rds_vars1, rds_vars2, by = "var")

ranks_gg <- ggplot(rds_vars, aes(x = rel.inf.x, y = rel.inf.y, label = var)) +
  geom_text_repel() +
  geom_jitter() +
  labs(x = "Feature Rank for Initial Model (AUC = 0.93)", y = "Feature rank for Pseudoabsence Model (AUC = 0.95)") +
  theme_bw()

ranks2_gg <- ggplot(rds_vars, aes(x = rel.inf.x, y = rel.inf.y, label = var)) +
  geom_text_repel() +
  geom_jitter() +
  scale_x_reverse(limits = c(12.92, 0)) +
  scale_y_reverse() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Feature Rank for Initial Model (AUC = 0.93)", y = "Feature rank for Pseudoabsence Model (AUC = 0.95)") +
  theme_bw()

# multiplot of ranks and dotplots
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/feature rank comparisons.png", width=15,height=7,units="in",res=600)
rel_gg + ranks2_gg + plot_layout(widths = c(1, 2))
dev.off()

#### Functions for patial dependence plots
## Partial dependence plots for continuous variables 
make_pdp_cont <- function(model, predictor, var_name, pcolor = FALSE) {
  
  # return grid
  vals <- plot.gbm(model[["mod"]], i.var = predictor, type = "response", return.grid = TRUE)
  
  # data for hist
  yrange = range(vals$y, na.rm = TRUE)
  
  # pull histogram values
  hi=hist(model[["testdata"]][[predictor]],breaks=30,plot=F)
  hi=with(hi,data.frame(breaks[1:(length(breaks)-1)],counts))
  names(hi)=c("mids","counts")
  
  if(pcolor == FALSE){
    
    # plot
    ggplot() + 
      geom_line(data = vals, aes(x = !!sym(predictor), y = y)) +
      geom_segment(data=hi,inherit.aes=F,
                   aes(x=mids,xend=mids,
                       y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                   size=2,colour="darkgrey",alpha=0.50) +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    
  }else{
    
    ggplot() + 
      geom_line(data = vals, aes(x = !!sym(predictor), y = y)) +
      geom_segment(data=hi,inherit.aes=F,
                   aes(x=mids,xend=mids,
                       y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                   size=2,colour="Orange",alpha=0.40) +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
  }
  
}

## Function for factor pdp plots
make_pdp_fact <- function(model, predictor, var_name, pcolor = FALSE) {
  
  # return grid
  vals <- plot.gbm(model[["mod"]], i.var = predictor, type = "response", return.grid = TRUE)
  
  # data for hist
  yrange = range(vals$y, na.rm = TRUE)
  
  # pull counts for color
  df_cat <- as.data.frame(table(noNA[["testdata"]][[predictor]]))
  
  # fix y axis point
  df_cat$ymin <- yrange[1]-0.01
  
  if(pcolor == FALSE){ #greys for initial model
    
    ggplot() +
      geom_point(data = vals, size= 2, shape=15, aes(category, y)) +
      geom_point(data = df_cat, aes(Var1, ymin, color = Freq)) +
      scale_color_continuous(high = "#636363", low = "#D9D9D9", guide = "none") +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    
  }else{ #different color for pseudoab variables
    
    ggplot() +
      geom_point(data = vals, size= 2, shape=15, aes(category, y)) +
      geom_point(data = df_cat, aes(Var1, ymin, color = Freq)) +
      scale_color_continuous(high = "#8C2D04", low = "#FEE6CE", guide = "none") +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    
  }
  
}

## Plot partial dependence
# No NA pdps
hb <- make_pdp_cont(noNA,"habitat_breadth_n", "Habitat Breadth", pcolor = FALSE)
gr <- make_pdp_cont(noNA, "X26.1_GR_Area_km2", "Geographic Area (km2)", pcolor = FALSE)
at <- make_pdp_cont(noNA, "X30.1_AET_Mean_mm", "Mean Monthly AET", pcolor = FALSE)
pm <- make_pdp_cont(noNA, "X28.1_Precip_Mean_mm", "Mean Monthly Precipitation (mm)", pcolor = FALSE)
ls <- make_pdp_cont(noNA, "litter_size_n", "Litter Size", pcolor = FALSE)
fr <- make_pdp_cont(noNA, "det_fruit", "Diet Fruit (%)", pcolor = FALSE)
fa <- make_pdp_cont(noNA, "adult_forearm_length_mm", "Adult Forearm Length", pcolor = FALSE)
am <- make_pdp_cont(noNA, "adult_mass_g","Adult Mass (g)", pcolor = FALSE)
bl <- make_pdp_cont(noNA, "adult_body_length_mm", "Adult Body Length", pcolor = FALSE)
ml <- make_pdp_cont(noNA, "X26.3_GR_MinLat_dd", "Minimum Latitude", pcolor = FALSE)
mx <- make_pdp_cont(noNA, "X26.3_GR_MinLat_dd", "Maximum Latitude", pcolor = FALSE)
di <- make_pdp_cont(noNA, "dphy_invertebrate", "Diet Invertebrate (%)", pcolor = FALSE)
hp <- make_pdp_cont(noNA, "X27.2_HuPopDen_Mean_n.km2", "Mean Human Population Density", pcolor = FALSE)
dp <- make_pdp_cont(noNA, "dphy_plant","Diet Plants (%)", pcolor = FALSE)
cs <- make_pdp_fact(noNA, "category", "Conservation Status", pcolor = FALSE)

# Save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/pdp 15 no NAs.png", width=7,height=8,units="in",res=600)
hb + gr + at + pm + ls + fr + fa + am + bl + ml + mx + di + hp + dp + cs + plot_layout(nrow = 5, ncol = 3, byrow = TRUE)
dev.off()

png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/pdp 9 no NAs.png", width=10,height=8,units="in",res=600)
hb + gr + at + pm + ls + fr + fa + am + bl + plot_layout(nrow = 3, ncol = 3, byrow = TRUE)
dev.off()

#  Pseudo model pdps 
phb <- make_pdp_cont(pseudo,"habitat_breadth_n", "Habitat Breadth", pcolor = TRUE)
pgr <- make_pdp_cont(pseudo, "X26.1_GR_Area_km2", "Geographic Area (km2)", pcolor = TRUE)
pat <- make_pdp_cont(pseudo, "X30.1_AET_Mean_mm", "Mean Monthly AET", pcolor = TRUE)
pls <- make_pdp_cont(pseudo, "litter_size_n", "Litter Size", pcolor = TRUE)
ppm <- make_pdp_cont(pseudo, "X28.1_Precip_Mean_mm", "Mean Monthly Precipitation (mm)", pcolor = TRUE)
pcc <- make_pdp_cont(pseudo, "cites","Citation Count", pcolor = TRUE)
pml <- make_pdp_cont(pseudo, "X26.3_GR_MinLat_dd", "Minimum Latitude", pcolor = TRUE)
pmx <- make_pdp_cont(pseudo, "X26.3_GR_MinLat_dd", "Maximum Latitude", pcolor = TRUE)
pam <- make_pdp_cont(pseudo, "adult_mass_g","Adult Mass (g)", pcolor = TRUE)
pfa <- make_pdp_cont(pseudo, "adult_forearm_length_mm", "Adult Forearm Length", pcolor = TRUE)
pfr <- make_pdp_cont(pseudo, "det_fruit", "Diet Fruit (%)", pcolor = TRUE)
php <- make_pdp_cont(pseudo, "X27.2_HuPopDen_Mean_n.km2", "Mean Human Population Density", pcolor = TRUE)
pbl <- make_pdp_cont(pseudo, "adult_body_length_mm", "Adult Body Length", pcolor = TRUE)
pcs <- make_pdp_fact(pseudo, "category", "Conservation Status", pcolor = TRUE)
pdi <- make_pdp_cont(noNA, "dphy_invertebrate", "Diet Invertebrate (%)", pcolor = TRUE)

# Save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/pdp 15 pseudo model.png", width=7,height=8,units="in",res=600)
phb + pgr + pat + pls + ppm + pcc + pml + pmx + pam + pfa + pfr + php + pbl + pcs + pdi + plot_layout(nrow = 5, ncol = 3, byrow = TRUE)
dev.off()



