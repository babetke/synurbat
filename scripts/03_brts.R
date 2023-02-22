# BRT analysis - Determinants of anthropogenic roosting
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

# read in data (do I want a separate file for data processing?)
data <- read_csv("~/Desktop/Synurbic_Bats/synurbat/flat files/cleaned dataset 30 cutoff.csv")

# data classes didn't hold in the cleaning file (may just want to remove that and factorize here)
data <- data %>% # Synurbic and variables that are factors according to COMBINE
  mutate(across(c("trophic_level","foraging_stratum","activity_cycle",
                  "disected_by_mountains", "glaciation", "biogeographical_realm","fam_EMBALLONURIDAE":"fam_VESPERTILIONIDAE"), 
                factor))

# remove species col
data <- data %>% select(-species) 

# how many NAs are there - 264 as of 02/18/23
length(data$Synurbic[is.na(data$Synurbic)])

# I also need a second dataset without NAs in synurbic. NAs cannot be in response
na.data <- data[!is.na(data$Synurbic),]

# Set up BRT tuning via search grid
# function to make grids?
## hyperparameter grid, maybe allow the number of seeds to change?
makegrid <- function(seed) {
  
  # create grid
  tgrid <- expand.grid(n.trees = seq(5000,15000,5000),
                       interaction.depth = c(2,3,4),
                       shrinkage = c(0.01,0.001,0.0005),
                       n.minobsinnode = 4, 
                       seed = seq(1,seed,by = 1))
  
  ## trees, depth, shrink, min, prop
  tgrid$id=with(tgrid,paste(n.trees,interaction.depth,shrinkage,n.minobsinnode))
  
  ## sort by id then seed
  tgrid=tgrid[order(tgrid$id,tgrid$seed),]
  
  ## now add rows
  tgrid$row=1:nrow(tgrid)
  
  ## factor id
  tgrid$id2=factor(as.numeric(factor(tgrid$id)))
  
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
hfit=function(row, data_df, response, folds, nsplit){
  
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

hgrid <- makegrid(1)

# trim for now
hgrid <- hgrid[1,]

spars <- lapply(1:nrow(hgrid),function(x) hfit(x, data_df = na.data, response="Synurbic", folds = 5, nsplit = "yes"))

## get results
sresults=data.frame(sapply(spars,function(x) x$trainAUC),
                    sapply(spars,function(x) x$testAUC),
                    sapply(spars,function(x) x$spec),
                    sapply(spars,function(x) x$sen),
                    sapply(spars,function(x) x$wrow),
                    sapply(spars,function(x) x$best))
names(sresults)=c("trainAUC","testAUC",
                  "spec","sen","row","best")

## combine and save
vsearch_complete = merge(sresults,hgrid,by="row")

# define BRT function 
# take a specified dataset
# label response
# maybe build in some options to change the split section - keep from analysis before but maybe add a if else 

# Basic gbm functions from virus analysis - the innerds of partition function w/out partitions 
set.seed(seed)
testbrt <- function(data_df, response, nt, shr, int.d) {
  
  # rename dataset
  ndata <- data_df
  
  # correct the response and remove raw response variables
  ndata <- ndata %>% 
    mutate(response = !!sym(response)) %>%
    select(-c("pseudo", "Synurbic"))

  ## test and train for data, same because no splits
  train <- ndata
  test <- ndata
  
  # pull response test and train
  yTrain <- train$response
  yTest <- test$response
  
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
  
  ## performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter <- gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds <- predict(gbmOut,test,n.trees=best.iter,type="response")
  
  ## known
  result <- test$response
  
  ## sensitiviy and specificity
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
              trainAUC=auc_train,
              testAUC=auc_test,
              spec=spec,
              sen=sen,
              roc=perf,
              rinf=bars,
              #predict=pred_data,
              traindata=train,
              testdata=test,
              seed=seed))
}

# change name of response
ndata <- data_df

# correct the response and remove raw response variables
ndata <- ndata %>% 
  mutate(response = !!sym(response)) %>%
  select(-c("virus", "zvirus", "dum_zvirus", "dum_virus"))

## test and train for data, same because no splits
train <- ndata
test <- ndata

# pull response test and train
yTrain <- train$response
yTest <- test$response

## parameters from parameter_df - search grid data
# nt <- params$n.trees
# shr <- params$shrinkage
# int.d <- params$interaction.depth

nt <- nt
shr <- shr
int.d <- int.d

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

## performance
par(mfrow=c(1,1),mar=c(4,4,1,1))
best.iter <- gbm.perf(gbmOut,method="cv")

## predict with test data
preds <- predict(gbmOut,test,n.trees=best.iter,type="response")

## known
result <- test$response

## sensitiviy and specificity
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
            trainAUC=auc_train,
            testAUC=auc_test,
            spec=spec,
            sen=sen,
            roc=perf,
            rinf=bars,
            #predict=pred_data,
            traindata=train,
            testdata=test,
            seed=seed))
