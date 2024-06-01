# BRT analysis - Trait determinants of anthropogenic roosting
# babetke@utexas.edu
# Updated 06/15/2023

# Script for data cleaning
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
library(gbm)
library(rsample)
library(ROCR)
library(caret) 
library(InformationValue)

# read in data 
# data <- read_csv("~/Desktop/Synurbic_Bats/synurbat/flat files/cleaned dataset 30 cutoff.csv") Directory for personal comp
# data <- read_csv("/Volumes/BETKE 2021/synurbat/flat files/cleaned dataset 30 cutoff.csv") # directory for lab comp

# data classes didn't hold in the cleaned file (may just want to remove that and factorize here OR save as RDS and read in here)
# data <- data %>% # Synurbic and variables that are factors according to COMBINE
#   mutate(across(c("trophic_level","foraging_stratum","activity_cycle",
#                   "disected_by_mountains", "glaciation", "biogeographical_realm","fam_EMBALLONURIDAE":"fam_VESPERTILIONIDAE"), 
#                 factor))

# read in rds
data <- readRDS("~/Desktop/Synurbic_Bats/synurbat/flat files/cleaned dataset 30 cutoff.rds")

# lab comp
# data <- readRDS("/Volumes/BETKE 2021/synurbat/flat files/cleaned dataset 30 cutoff.rds")

# look into distribution of continuous variables
# pull numeric vars
num <- select(data, where(is.numeric))

# remove % diet variables
num <- select(num, !starts_with(c("det","dphy")))

# remove ordinal type variables
num <- num %>% select(!c(litters_per_year_n, litter_size_n, island_dwelling, habitat_breadth_n))

# histograms
Hmisc::hist.data.frame(num)

# log plus constant
log_data <- data
log_data$log_cites <- log1p(log_data$cites)
log_data$log_lower_elevation_m <- log1p(log_data$lower_elevation_m)
log_data$log_X26.1_GR_Area_km2 <- log1p(log_data$X26.1_GR_Area_km2)
log_data$log_X27.1_HuPopDen_Min_n.km2 <- log1p(log_data$X27.1_HuPopDen_Min_n.km2)
log_data$log_X27.2_HuPopDen_Mean_n.km2 <- log1p(log_data$X27.2_HuPopDen_Mean_n.km2)
log_data$log_X27.3_HuPopDen_5p_n.km2 <- log1p(log_data$X27.3_HuPopDen_5p_n.km2)

# log no constant
log_data$log_adult_body_length_mm <- log10(log_data$adult_body_length_mm)
log_data$log_adult_forearm_length_mm <- log10(log_data$adult_forearm_length_mm)
log_data$log_adult_mass_g <- log10(log_data$adult_mass_g)

# # remove species col
# data <- data %>% 
#   select(-c(species, det_inv, biogeographical_realm)) 

# remove species col and non-log versions
log_data <- log_data %>%
  select(!c(species, det_inv, biogeographical_realm, cites, lower_elevation_m, X26.1_GR_Area_km2,
            X27.1_HuPopDen_Min_n.km2, X27.2_HuPopDen_Mean_n.km2, X27.3_HuPopDen_5p_n.km2,
            adult_body_length_mm, adult_forearm_length_mm, adult_mass_g))

# how many NAs are there - 234
length(log_data$Synurbic[is.na(log_data$Synurbic)])

# # make Synurbic numeric (gbm will crash)
# data$Synurbic <- as.numeric(as.character(data$Synurbic))
# data$pseudo <- as.numeric(as.character(data$pseudo))

# synrubic numeric for log data
log_data$Synurbic <- as.numeric(as.character(log_data$Synurbic))
log_data$pseudo <- as.numeric(as.character(log_data$pseudo))

# clean out NAs
log_na.data <- log_data[!is.na(log_data$Synurbic),] # 1045 species

# remove 
rm(num)

# # dataset without NAs in synurbic. NAs cannot be in response.
# na.data <- data[!is.na(data$Synurbic),] # 1042 species

# Set up BRT tuning via search grid
# function to make grids?
## hyperparameter grid, maybe allow the number of seeds to change?
makegrid <- function(seed,trees) {
  
  # create grid
  tgrid <- expand.grid(n.trees = trees,
                       interaction.depth = c(2,3,4),
                       shrinkage = c(0.01,0.001,0.0005),
                       n.minobsinnode = 4, 
                       seed = seq(1,seed,by = 1))
  
  ## trees, depth, shrink, min, prop
  tgrid$id <- with(tgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))
  
  ## sort by id then seed
  tgrid <- tgrid[order(tgrid$id,tgrid$seed),]
  
  ## now add rows
  tgrid$row <- 1:nrow(tgrid)
  
  ## factor id
  tgrid$id2 <- factor(as.numeric(factor(tgrid$id)))
  
  return(tgrid)
}

#### Function to assess each hyperparameter combination
# this function takes a search grid and runs them through gbms for a given split of data
# then assess performance of each combination of parameters to get the optimal parameters for
# final gbms.


# I added set and split_prop so could change datasets and proportion of splits
# row <- row in the hgrid (I assume)
# set <- the subset of data I want to use (two full datasets and two trimmed to complete)
# response <- indicate which response I want to use. Either dum_virus or dum_zvirus
# folds <- indicating either 10 or 5 folds depending on the dataset (5 for subset of 544 and 10 for full)
# nsplit <- allow for randomized test and training sets or only data
grid_search <- function(row, data_df, response, folds, nsplit){
  
  # new data
  ndata <- data_df
  
  # correct the response and remove raw response variables
  ndata <- ndata %>% 
    mutate(response = !!sym(response)) %>%
    select(-c("pseudo", "Synurbic"))
  
  # how to indicate no splits
  nsplit = nsplit
  if(nsplit == "yes"){
    
    ## test and train for data
    dataTrain <- ndata
    dataTest <- ndata
    
    # pull response test and train
    yTrain <- dataTrain$response
    yTest <- dataTest$response
    
  } else {
    
    ## use rsample to split (allow the proportion to be changed incase)
    set.seed(hgrid$seed[row])
    split=initial_split(ndata,prop=0.7,strata="response")
    
    ## test and train
    dataTrain=training(split)
    dataTest=testing(split)
    
    ## yTest and yTrain
    yTrain=dataTrain$response
    yTest=dataTest$response
    
  }
  
  ## BRT
  set.seed(1)
  gbmOut=gbm(response ~ . ,data=dataTrain,
             n.trees=hgrid$n.trees[row],
             distribution="bernoulli",
             shrinkage=hgrid$shrinkage[row],
             interaction.depth=hgrid$interaction.depth[row],
             n.minobsinnode=hgrid$n.minobsinnode[row],
             cv.folds=folds,
             class.stratify.cv=TRUE,
             bag.fraction=0.5,
             train.fraction=1,
             n.cores=1,
             verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## known
  result=dataTest$response
  
  ## sensitiviy and specificity
  sen=InformationValue::sensitivity(result,preds)
  spec=InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## print
  print(paste("hpar row ",row," done; test AUC is ",auc_test,sep=""))
  
  ## save outputs
  return(list(best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              wrow=row))
}

# ifelse for running grid search function
gsrun = "no" 
if(gsrun == "yes"){# run grid search 
  
  hgrid <- makegrid(1, c(seq(5000, 25000, 5000)))
  
  # Removing combinations that constantly max out
  hgrid %>% 
    filter(!(n.trees == 5000 & shrinkage < 0.001)) %>%
    filter(!(n.trees == 10000 & shrinkage < 0.001)) -> hgrid
  
  # renumber the rows for matching 
  hgrid$row <- 1:nrow(hgrid)
  
  # run for the two types of data?
  na.pars <- lapply(1:nrow(hgrid),function(x) grid_search(x, data_df = na.data, response="Synurbic", folds = 10, nsplit = "yes"))
  p.pars <- lapply(1:nrow(hgrid),function(x) grid_search(x, data_df = data, response="pseudo", folds = 10, nsplit = "yes"))
  
  ## get results for main text/initial model
  na.results <- data.frame(sapply(na.pars,function(x) x$trainAUC),
                        sapply(na.pars,function(x) x$testAUC),
                        sapply(na.pars,function(x) x$spec),
                        sapply(na.pars,function(x) x$sen),
                        sapply(na.pars,function(x) x$wrow),
                        sapply(na.pars,function(x) x$best))
  names(na.results) <- c("trainAUC","testAUC",
                      "spec","sen","row","best")
  
  # Merge with hgrid
  na.complete <- merge(na.results, hgrid, by = "row")
  
  # write as csv
  write_csv(na.complete, "/Volumes/BETKE 2021/synurbat/flat files/grid search without NAs.csv")
  
  ## get results for pseudo model
  p.results <- data.frame(sapply(p.pars,function(x) x$trainAUC),
                          sapply(p.pars,function(x) x$testAUC),
                          sapply(p.pars,function(x) x$spec),
                          sapply(p.pars,function(x) x$sen),
                          sapply(p.pars,function(x) x$wrow),
                          sapply(p.pars,function(x) x$best))
  names(p.results) <- c("trainAUC","testAUC",
                        "spec","sen","row","best")
  
  # Merge with hgrid
  p.complete <- merge(p.results, hgrid, by = "row")
  
  # write as csv
  write_csv(p.complete, "/Volumes/BETKE 2021/synurbat/flat files/grid search pseudo.csv")
  
} else {# read in grid search results
  na.complete <- read_csv("~/Desktop/Synurbic_Bats/synurbat/flat files/grid search without NAs.csv")
  p.complete <- read_csv("~/Desktop/Synurbic_Bats/synurbat/flat files/grid search pseudo.csv")
  # na.complete <- read_csv("/Volumes/BETKE 2021/synurbat/flat files/grid search without NAs.csv")
}

# Plot of parameters
# sen_gg <- ggplot(na.complete, aes(x = factor(shrinkage), y = sen)) +
#   geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
#   theme_classic() +
#   labs(x = "Shrinkage", y = "Sensitivity", fill = "Interaction Depth") +
#   theme(legend.position = "none")
# 
# spec_gg <- ggplot(na.complete, aes(x = factor(shrinkage), y = spec)) +
#   geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
#   theme_classic() +
#   labs(x = "Shrinkage", y = "Specificity", fill = "Interaction Depth")

auc_gg <- ggplot(na.complete, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = "AUC", fill = "Interaction Depth", title = "Initial Model") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")

tree_gg <- ggplot(na.complete, aes(x = factor(n.trees), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "No.Trees", y = "AUC", fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")
  
pauc_gg <- ggplot(p.complete, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "Learning Rate", y = NULL, fill = "Interaction Depth", title = "Pseudoabsence Model") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")

ptree_gg <- ggplot(p.complete, aes(x = factor(n.trees), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_bw() +
  labs(x = "No.Trees", y = NULL, fill = "Interaction Depth") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10)) +
  scale_fill_brewer(palette="Accent")

library(patchwork) # multiplot and save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure S4.png", width=6, height=6,units="in",res=300)
guide_area () + 
  (auc_gg + pauc_gg) / (tree_gg + ptree_gg) + 
  plot_layout(guides = "collect", heights = c(1, 15)) & theme(legend.position = "top")
dev.off()

# remove
rm(auc_gg, tree_gg, pauc_gg, ptree_gg)

# unload patchwork
detach("package:patchwork", unload = TRUE)

# Sort output to view top model combinations
sort <- na.complete %>% 
  arrange(desc(testAUC))

# Sort output to view top model combinations
psort <- p.complete %>% 
  arrange(desc(testAUC))

# remove sorts
rm(sort, psort)

# Define BRT function 
# take a specified dataset
# label response
# maybe build in some options to change the split section - keep from analysis before but maybe add a if else 

# Basic gbm functions from virus analysis - the innerds of partition function w/out partitions 
get_brt <- function(data_df, response, nt, shr, int.d, nsplit, seed=NULL) {
  
  # rename dataset
  ndata <- data_df
  
  # correct the response and remove raw response variables
  ndata <- ndata %>% 
    mutate(response = !!sym(response)) %>%
    select(-c("pseudo", "Synurbic"))
  
  # how to indicate no splits
  nsplit = nsplit
  if(nsplit == "yes"){
  
  ## test and train for data, same because no splits
  train <- ndata
  test <- ndata
  
  # pull response test and train
  yTrain <- train$response
  yTest <- test$response
  
  } else {
    
    ## use rsample to split (allow the proportion to be changed incase)
    set.seed(seed)
    split <- initial_split(ndata,prop=0.7,strata="response")
    
    ## test and train
    train <- training(split)
    test <- testing(split)
    
    ## yTest and yTrain
    yTrain <- train$response
    yTest <- test$response
    
  }
  
  ## parameters from parameter_df - search grid data
  # nt <- params$n.trees
  # shr <- params$shrinkage
  # int.d <- params$interaction.depth
  
  ## BRT
  set.seed(1)
  gbmOut <- gbm(response ~ . ,data=train,
                n.trees=nt,
                distribution="bernoulli",
                shrinkage=shr,
                interaction.depth=int.d,
                n.minobsinnode=4,
                cv.folds=10,class.stratify.cv=TRUE,
                bag.fraction=0.5,train.fraction=1,
                n.cores=1,
                verbose=F)
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter <- gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds <- predict(gbmOut,test,n.trees=best.iter,type="response")
  
  ## known
  result <- test$response
  
  ## sensitivity and specificity
  sen <- InformationValue::sensitivity(result,preds)
  spec <- InformationValue::specificity(result,preds)
  
  ## AUC on train
  auc_train <- gbm.roc.area(yTrain,predict(gbmOut,train,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test <- gbm.roc.area(yTest,predict(gbmOut,test,n.trees=best.iter,type="response"))
  
  ## relative importance
  bars <- summary(gbmOut,n.trees=best.iter,plotit=F)
  bars$rel.inf <- round(bars$rel.inf,2)
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              preds=preds,
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              rinf=bars,
              traindata=train,
              testdata=test))
  
}

# # Run BRTs
# noNA_gbm <- get_brt(data_df = na.data, response = "Synurbic", nt = 20000, shr = 0.001, int.d = 4, nsplit = "yes")
# pseudo_gbm <- get_brt(data_df = data, response = "pseudo", nt = 25000, shr = 0.001, int.d = 4, nsplit = "yes")
# 
# # Save 
# saveRDS(noNA_gbm,"/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/noNA_brts.rds")
# saveRDS(pseudo_gbm,"/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/pseudo_brts.rds")

# # save lab computer directory
# saveRDS(noNA_gbm,"/Volumes/BETKE 2021/synurbat/flat files/noNA_brts.rds")
# saveRDS(pseudo_gbm,"/Volumes/BETKE 2021/synurbat/flat files/pseudo_brts.rds")

# brts with transformed dataframe
log_gbm <- get_brt(data_df = log_na.data, response = "Synurbic", nt = 20000, shr = 0.001, int.d = 4, nsplit = "yes")
log_pseudo_gbm <- get_brt(data_df = log_data, response = "pseudo", nt = 20000, shr = 0.001, int.d = 4, nsplit = "yes")

# save
saveRDS(log_gbm,"/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/log_brts.rds")
saveRDS(log_pseudo_gbm,"/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/log_pseudo_brts.rds")

