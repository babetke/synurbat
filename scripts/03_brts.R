# BRT analysis - Trait determinants of anthropogenic roosting
# babetke@utexas.edu

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

# data classes didn't hold in the cleaning file (may just want to remove that and factorize here OR save as RDS and read in here)
# data <- data %>% # Synurbic and variables that are factors according to COMBINE
#   mutate(across(c("trophic_level","foraging_stratum","activity_cycle",
#                   "disected_by_mountains", "glaciation", "biogeographical_realm","fam_EMBALLONURIDAE":"fam_VESPERTILIONIDAE"), 
#                 factor))

# read in rds
data <- readRDS("~/Desktop/Synurbic_Bats/synurbat/flat files/cleaned dataset 30 cutoff.rds")

# remove species col
data <- data %>% 
  select(-c(species, det_inv, biogeographical_realm)) 

# how many NAs are there - 253
length(data$Synurbic[is.na(data$Synurbic)])

# make Synurbic numeric (gbm will crash)
data$Synurbic <- as.numeric(as.character(data$Synurbic))
data$pseudo <- as.numeric(as.character(data$pseudo))

# dataset without NAs in synurbic. NAs cannot be in response.
na.data <- data[!is.na(data$Synurbic),]

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
             class.stratify.cv=FALSE,
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

hgrid <- makegrid(1, trees = seq(15000,50000,5000))
rgrid <- makegrid(1, c(15000, seq(20000, 50000, 10000))) #45, trimmed to not have 5000, or 10000 because they max out on

# run for the two types of data?
na.pars <- lapply(1:nrow(hgrid),function(x) hfit(x, data_df = na.data, response="Synurbic", folds = 10, nsplit = "yes"))
f.pars <- lapply(1:nrow(hgrid),function(x) hfit(x, data_df = data, response="pseudo", folds = 10, nsplit = "yes"))

## get results
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

# Plot of parameters
auc_gg <- ggplot(na.complete, aes(x = factor(shrinkage), y = testAUC)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Shrinkage", y = "AUC", fill = "Interaction Depth") +
  theme(legend.position = "none")

sen_gg <- ggplot(na.complete, aes(x = factor(shrinkage), y = sen)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Shrinkage", y = "Sensitivity", fill = "Interaction Depth") +
  theme(legend.position = "none")

spec_gg <- ggplot(na.complete, aes(x = factor(shrinkage), y = spec)) +
  geom_boxplot(aes(fill = factor(interaction.depth)), color = "black", alpha = 0.5) +
  theme_classic() +
  labs(x = "Shrinkage", y = "Specificity", fill = "Interaction Depth")

library(patchwork) # multiplot and save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/parameters no NAs.png", width=15,height=5,units="in",res=600)
auc_gg + sen_gg + spec_gg 
dev.off()

# remove
rm(auc_gg, sen_gg, spec_gg)

# unload patchwork
detach("package:patchwork", unload = TRUE)

# write as csv
write_csv(na.complete, "/Volumes/BETKE 2021/synurbat/flat files/grid search without NAs.csv")

# Define BRT function 
# take a specified dataset
# label response
# maybe build in some options to change the split section - keep from analysis before but maybe add a if else 

# Basic gbm functions from virus analysis - the innerds of partition function w/out partitions 
get_brt <- function(data_df, response, nt, shr, int.d, folds, nsplit, seed=NULL) {
  
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

# Run BRTs
noNA_gbm <- get_brt(data_df = na.data, response = "Synurbic", nt = 10000, shr = 0.001, int.d = 4, nsplit = "yes")
pseudo_gbm <- get_brt(data_df = data, response = "pseudo", nt = 10000, shr = 0.001, int.d = 4, nsplit = "yes")

# Save 
saveRDS(noNA_gbm,"/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/noNA_brts.rds")
saveRDS(pseudo_gbm,"/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/flat files/pseudo_brts.rds")




