## phylogenetic analyses of bat roosting data
## danbeck@ou.edu
## last updated 03/12/2023

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(plyr)
library(ggplot2)
library(ape)
library(caper)
library(treeio)
library(ggtree)
library(car)
library(phylofactor)
library(phytools)

## load in data
setwd("~/Desktop/synurbat/flat files")
data=read.csv("Bat References Spreadsheet.csv")

## load in Upham phylogeny
setwd("~/Desktop/bathaus/phylos")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa=taxa[taxa$ord=="CHIROPTERA",]
taxa$tip=taxa$Species_Name

## trim phylo to bats
tree=keep.tip(tree,taxa$tiplabel)

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## merge data and taxa by species
data=merge(data,taxa[c("species","fam","gen","clade")],by="species")
rm(taxa)

## make label
data$label=data$species

## merge
cdata=comparative.data(phy=tree,data=data,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## species
cdata$data$Species=cdata$data$species

## trim to tree with data
set=cdata[!is.na(cdata$data$Synurbic),]

## phylogenetic signal in response
## D of 0 = Brownian model, D of 1 = random (no phylogenetic signal)
set.seed(1)
mod=phylo.d(set,binvar=Synurbic,permut=1000); mod

## extract states
state=setNames(set$data$Synurbic,rownames(set$data))

## fit models
ER.model=fitMk(set$phy,state,model="ER") 
ARD.model=fitMk(set$phy,state,model="ARD") 
Irr1.model=fitMk(set$phy,state,model=matrix(c(0,1,0,0),2,2,byrow=TRUE)) 
Irr2.model=fitMk(set$phy,state,model=matrix(c(0,0,1,0),2,2,byrow=TRUE))

## model comparison
state.aov=anova(ER.model,ARD.model,Irr1.model,Irr2.model)

## simmap
state.simmap=simmap(state.aov,nsim=100)

## set colors and plot
cols=setNames(viridisLite::viridis(n=2),levels(state)) 
plot(summary(state.simmap),ftype="i", fsize=0.7,colors=cols,cex=c(0.6,0.3))

## get density
state.density=density(state.simmap)

## plot
par(mfrow=c(1,2),mar=c(4.5,4.5,0.5,0.5))
COLS<-setNames(cols,state.density$trans) 
plot(state.density,ylim=c(0,0.6), transition=names(COLS)[1],colors=COLS[1],main="") 
mtext("a)transitionstopiscivory",line=1, adj=0,cex=0.8) 

## other
plot(state.density,ylim=c(0,0.6), transition=names(COLS)[2],colors=COLS[2], main="") 
mtext("b)transitionstonon-piscivory", line=1,adj=0,cex=0.8)

## density map
state.densityMap=densityMap(state.simmap, plot=FALSE,res=100) 

## plot
state.densityMap=setMap(state.densityMap, viridisLite::viridis(n=10)) 
par(mfrow=c(1,1),mar=c(0.5,0.5,0.5,0.5))
plot(state.densityMap,lwd=1,outline=F)

## taxonomy
set$data$taxonomy=paste(set$data$fam,set$data$gen,set$data$Species,sep='; ')

## set taxonomy
taxonomy=data.frame(set$data$taxonomy)
names(taxonomy)="taxonomy"
taxonomy$Species=rownames(set$data)
taxonomy=taxonomy[c("Species","taxonomy")]
taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf){
  
  ## get formula
  #chars=as.character(pf$frmla.phylo)[-1]
  chars=strsplit(as.character(pf$frmla.phylo)," ~ ")[[1]]
  
  ## response
  resp=chars[1]
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## binary
set.seed(1)
bpf=gpf(Data=set$data,tree=set$phy,
        frmla.phylo=Synurbic~phylo,
        family=binomial,algorithm='phylo',nfactors=3)

## summarize
bpf_results=pfsum(bpf)$results

## save
set$data$status=factor(set$data$Synurbic)

## save tree
set$data$label=set$data$species
dtree=treeio::full_join(as.treedata(set$phy),set$data,by="label")

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## set x max
plus=1
pplus=plus+1

## ggtree
gg=ggtree(dtree,size=0.2,colour="grey30")+
  geom_tippoint(aes(x=x+0.5,colour=status),shape=15,size=1)+
  scale_colour_manual(values=c("grey80","black"))+
  theme(legend.position = "None")

## reset pplus
pplus=1.5

## add clades
for(i in 1:nrow(bpf_results)){
  
  gg=gg+
    geom_hilight(node=bpf_results$node[i],
                 alpha=0.25,
                 fill=ifelse(bpf_results$clade>
                               bpf_results$other,pcols[2],pcols[1])[i])
    # geom_cladelabel(node=bpf_results$node[i],
    #                 label=bpf_results$taxa[i],
    #                 offset=pplus,
    #                 hjust=0.5,
    #                 offset.text=pplus*1.25,
    #                 parse=T,
    #                 angle=90)
}
bgg=gg
