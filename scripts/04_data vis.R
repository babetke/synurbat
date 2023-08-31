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
                           "X27.2_HuPopDen_Mean_n.km2" = "Mean Human Density",
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
                           "X27.1_HuPopDen_Min_n.km2" = "Min Human Density",
                           "dphy_vertebrate" = "Diet Vertebrate",
                           "lower_elevation_m" = "Lower Elevation Limit",
                           "det_nect" = "Diet Nectar",
                           "X27.4_HuPopDen_Change" = "Human Density Change",
                           "X27.3_HuPopDen_5p_n.km2" = "Human Density 5th Percentile",
                           "X26.7_GR_MidRangeLong_dd" = "Median Longitudinal Range",
                           "det_fruit" = "Diet Fruit",
                           "fam_PHYLLOSTOMIDAE" = "Phyllostomidae",
                           "det_vect" = "Diet Vect",
                           "trophic_level" = "Trophic Level",
                           "dphy_invertebrate" = "Diet Invertebrate",
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
      geom_point(aes(x = reorder(var, rel.inf/100), y = rel.inf/100, color = type), stat = "identity", alpha = 0.75, size = 3) +
      theme(axis.text.x = element_text(size = 10, hjust = 1)) +
      theme_bw() +
      labs(x = " ", y = "Relative Importance") +
      coord_flip() +
      theme(legend.position = c(0.80, 0.15),
            legend.text = element_text(size=6.5),
            legend.title = element_text(size=6.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(axis.text.y=element_text(size=6.5),
            axis.text.x=element_text(size=6.5),
            axis.title=element_text(size=8.5)) +
      scale_y_continuous(labels = scales::percent)
  
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
                         "X27.2_HuPopDen_Mean_n.km2" = "Mean Human Density",
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
                         "X27.1_HuPopDen_Min_n.km2" = "Min Human Density",
                         "dphy_vertebrate" = "Diet Vertebrate",
                         "lower_elevation_m" = "Lower Elevation Limit",
                         "det_nect" = "Diet Nectar",
                         "X27.4_HuPopDen_Change" = "Human Density Change",
                         "X27.3_HuPopDen_5p_n.km2" = "Human Density 5th Percentile",
                         "X26.7_GR_MidRangeLong_dd" = "Median Longitudinal Range",
                         "det_fruit" = "Diet Fruit",
                         "fam_PHYLLOSTOMIDAE" = "Phyllostomidae",
                         "det_vect" = "Diet Vect",
                         "trophic_level" = "Trophic Level",
                         "dphy_invertebrate" = "Diet Invertebrate",
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
  
  # correlation coeff
  rho <- cor.test(rds_vars$rel.inf.x, rds_vars$rel.inf.y, method = "pearson")
  
  # update names
  rds_vars <- rds_vars %>% 
    rename(Initial = rel.inf.x,
           Pseudoabsence = rel.inf.y)
  
  rds_vars <- rds_vars %>% 
    mutate(i_imp = Initial/100, 
           p_imp = Pseudoabsence/100)
  
  # points
  var.inf <- ggplot(rds_vars) + 
    #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
    geom_point(aes(x = reorder(var, i_imp), y = p_imp, color = "Pseudoabsence"), stat = "identity", alpha = 0.5, size = 3) +
    geom_point(aes(x = reorder(var, i_imp), y = i_imp, color = "Initial"), stat = "identity", alpha = 0.75, size = 2) +
    theme(axis.text.x = element_text(size = 10, hjust = 1)) +
    theme_bw() +
    labs(x = " ", y = "Relative Importance") +
    coord_flip() +
    scale_color_manual(values = c("Initial" = "black", "Pseudoabsence" = "orange")) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    scale_y_continuous(labels = scales::percent)
  
  var.inf <- var.inf + theme(legend.position = "bottom",
                             legend.text = element_text(size = 6),
                             legend.title = element_text(size = 6),
                             axis.text = element_text(size = 6),
                             axis.title=element_text(size = 8))
  
  # scatter plot
  ranks_gg <- ggplot(rds_vars, aes(x = i_imp, y = p_imp, label = var)) +
    geom_text_repel(size = 2) +
    geom_point() +
    #xlim(limits = c(0, 12.36)) +
    #geom_jitter() +
    #scale_x_reverse(limits = c(12.48, 0)) +
    #scale_y_reverse() +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Initial Relative Importance", y = "Pseudoabsence Relative Importance") +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.text = element_text(size=6),
          axis.title=element_text(size=7)) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 0.1252)) +
    scale_y_continuous(labels = scales::percent)
  
  # ranks_gg <- ranks_gg + theme(axis.text = element_text(size=5),
  #                              axis.title=element_text(size=6))
  
  #multi <- (var.inf + ranks_gg + plot_layout(widths = c(0.75, 2.25)) & theme(plot.tag = element_text(size = 7))) + plot_annotation(tag_levels = "A")
  
  
  # multi plot
  #multi <- (var.inf + ranks_gg + plot_layout(widths = c(1, 2)) & theme(plot.tag = element_text(size = 8)))
  # left <- var.inf & theme(plot.tag = element_text(size = 8))
  # right <- ranks_gg & theme(plot.tag = element_text(size = 8))
  # multi <- wrap_plots(left + right) + plot_annotation(tag_levels = 'A')
  
  #multi <- multi + plot_annotation(tag_levels = 'A')
  
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
  
  return(list(rho = rho,
              var.inf = var.inf,
              scatter = ranks_gg))
  }
  
}

# No NA only 
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure 2.png", width=4.5, height=5.75, units="in", res=300)
inf_plot(noNA)
dev.off()

# both models, scatter plot, and pearons rho
rel_gg <- inf_plot(noNA, pseudo)

# view cor output
rel_gg[["rho"]]

# Pearson's product-moment correlation
# 
# data:  rds_vars$rel.inf.x and rds_vars$rel.inf.y
# t = 24.337, df = 59, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9236110 0.9720261
# sample estimates:
#       cor 
# 0.9536318 

# or inset plot
library(cowplot)
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure S1.png", width=4.5,height=6,units="in",res=300)
ggdraw(rel_gg[["var.inf"]]) +
  draw_plot(rel_gg[["scatter"]], .55, .17, .42, .32) + # side to side, up and down, width, height
  draw_plot_label(
    c("A", "B"),
    c(0, 0.50), # moves A label around
    c(1, 0.50), # moves B label around
    size = 10
  )
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
  df_cat <- as.data.frame(table(model[["testdata"]][[predictor]]))
  
  # fix y axis point
  df_cat$ymin <- yrange[1]-0.01
  
  if(pcolor == FALSE){ #greys for initial model
    
    ggplot() +
      geom_point(data = vals, size= 2, shape=15, aes(!!sym(predictor), y)) +
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
      geom_point(data = vals, size= 2, shape=15, aes(!!sym(predictor), y)) +
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
pm <- make_pdp_cont(noNA, "X28.1_Precip_Mean_mm", "Mean Monthly Precipitation (mm)", pcolor = FALSE)
at <- make_pdp_cont(noNA, "X30.1_AET_Mean_mm", "Mean Monthly AET", pcolor = FALSE)
ls <- make_pdp_cont(noNA, "litter_size_n", "Litter Size", pcolor = FALSE)
mp <- make_pdp_cont(noNA, "X30.2_PET_Mean_mm", "Mean Monthly PET", pcolor = FALSE)
bl <- make_pdp_cont(noNA, "adult_body_length_mm", "Adult Body Length", pcolor = FALSE)
am <- make_pdp_cont(noNA, "adult_mass_g","Adult Mass (g)", pcolor = FALSE)
dp <- make_pdp_cont(noNA, "dphy_plant","Diet Plants (%)", pcolor = FALSE)
fr <- make_pdp_cont(noNA, "det_fruit", "Diet Fruit (%)", pcolor = FALSE)
mx <- make_pdp_cont(noNA, "X26.2_GR_MaxLat_dd", "Maximum Latitude", pcolor = FALSE)
fa <- make_pdp_cont(noNA, "adult_forearm_length_mm", "Adult Forearm Length", pcolor = FALSE)
cs <- make_pdp_fact(noNA, "category", "Conservation Status", pcolor = FALSE)
vt <- make_pdp_cont(noNA, "dphy_invertebrate","Diet Invertebrate", pcolor = FALSE)
hp <- make_pdp_cont(noNA, "X27.2_HuPopDen_Mean_n.km2", "Mean Human Density", pcolor = FALSE)

# Save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure 3.png", width=7,height=7.5,units="in",res=300)
hb + gr + pm + at + ls + mp + bl + am + dp + fr + mx + fa + cs + vt + hp + plot_layout(nrow = 5, ncol = 3, byrow = TRUE)
dev.off()

png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/pdp 9 no NAs.png", width=10,height=8,units="in",res=300)
hb + gr + pm + at + ls + mp + bl + am + dp + plot_layout(nrow = 3, ncol = 3, byrow = TRUE)
dev.off()

#  Pseudo model pdps 
phb <- make_pdp_cont(pseudo,"habitat_breadth_n", "Habitat Breadth", pcolor = TRUE)
pgr <- make_pdp_cont(pseudo, "X26.1_GR_Area_km2", "Geographic Area (km2)", pcolor = TRUE)
pat <- make_pdp_cont(pseudo, "X30.1_AET_Mean_mm", "Mean Monthly AET", pcolor = TRUE)
ppm <- make_pdp_cont(pseudo, "X28.1_Precip_Mean_mm", "Mean Monthly Precipitation (mm)", pcolor = TRUE)
pls <- make_pdp_cont(pseudo, "litter_size_n", "Litter Size", pcolor = TRUE)
pmx <- make_pdp_cont(pseudo, "X26.2_GR_MaxLat_dd", "Maximum Latitude", pcolor = TRUE)
pcc <- make_pdp_cont(pseudo, "cites","Citation Count", pcolor = TRUE)
pcs <- make_pdp_fact(pseudo, "category","Conservation Status", pcolor = TRUE)
pam <- make_pdp_cont(pseudo, "adult_mass_g","Adult Mass (g)", pcolor = TRUE)
ppt <- make_pdp_cont(pseudo, "X30.2_PET_Mean_mm", "Mean Monthly PET", pcolor = TRUE)
pdp <- make_pdp_cont(pseudo, "dphy_plant","Diet Plants (%)", pcolor = TRUE)
pml <- make_pdp_cont(pseudo, "X26.3_GR_MinLat_dd", "Minimum Latitude", pcolor = TRUE)
pfa <- make_pdp_cont(pseudo, "adult_forearm_length_mm", "Adult Forearm Length", pcolor = TRUE)
pbl <- make_pdp_cont(pseudo, "adult_body_length_mm", "Adult Body Length", pcolor = TRUE)
pfr <- make_pdp_cont(pseudo, "det_fruit", "Diet Fruit (%)", pcolor = TRUE)


# Save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure S2.png", width=7,height=7.5,units="in",res=300)
phb + pgr + pat + ppm + pls + pmx + pcc + pcs + pam + ppt + pdp + pml + pfa + pbl + pfr + plot_layout(nrow = 5, ncol = 3, byrow = TRUE)
dev.off()
