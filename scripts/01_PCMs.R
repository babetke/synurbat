## 01_phylogenetic analyses of bat roosting data
## danbeck@ou.edu
## last updated 03/21/2023

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

## load in roosting data
setwd("~/Desktop/synurbat/flat files")
data=read.csv("Bat References Spreadsheet.csv")

## load in taxonomy
setwd("~/Desktop/bathaus/phylos")
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa=taxa[taxa$ord=="CHIROPTERA",]
taxa$tip=taxa$Species_Name
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## merge data and taxa by species
data=merge(data,taxa[c("species","fam","gen","clade","extinct.","tiplabel")],by="species")
rm(taxa)

## remove extinct
data=data[!data$extinct.==1,]
data$extinct.=NULL

## load in Upham phylogeny
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## trim phylo to bats
tree=keep.tip(tree,data$tiplabel)

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))

## make label
data$label=data$species

## merge
cdata=comparative.data(phy=tree,data=data,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## species
cdata$data$Species=cdata$data$species

## make Synurbic with pseudoabsences
cdata$data$Synurbic_pseudo=ifelse(is.na(cdata$data$Synurbic),0,cdata$data$Synurbic)

## taxonomy
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')

## log1p citations
hist(cdata$data$cites)
hist(log1p(cdata$data$cites))
hist(sqrt(cdata$data$cites))
cdata$data$logcites=log1p(cdata$data$cites)

## save status and label
cdata$data$status=factor(cdata$data$Synurbic)
cdata$data$label=cdata$data$species

## trim to tree with data
set=cdata[!is.na(cdata$data$Synurbic),]

## phylogenetic signal in response for true and pseudo
## D of 0 = Brownian model, D of 1 = random (no phylogenetic signal)
set.seed(1)
mod1=phylo.d(set,binvar=Synurbic,permut=1000); mod1
set.seed(1)
mod2=phylo.d(cdata,binvar=Synurbic_pseudo,permut=1000); mod2

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>2,cs[3],cs[1])
  
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
  chars=as.character(pf$frmla.phylo)[-1]
  #chars=strsplit(as.character(pf$frmla.phylo)," ~ ")[[1]]
  
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

## binary, adjust for log1p citations
set.seed(1)
bpf=gpf(Data=set$data,tree=set$phy,
        frmla.phylo=Synurbic~phylo+cites,
        family=binomial,algorithm='phylo',nfactors=3,min.group.size=10)

## summarize
bpf_results=pfsum(bpf)$results

## repeat for pseudo
set.seed(1)
bpf2=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=Synurbic_pseudo~phylo+cites,
        family=binomial,algorithm='phylo',nfactors=3,min.group.size=10)

## summarize
bpf2_results=pfsum(bpf2)$results

## save tree
dtree=treeio::full_join(as.treedata(set$phy),set$data,by="label")
dtree2=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

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
gg=ggtree(dtree,size=0.2,colour="grey30",layout="fan")+
  geom_tippoint(aes(x=x+0.5,colour=status),shape=16,size=0.5)+
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

## repeat for pseudo
gg=ggtree(dtree2,size=0.2,colour="grey30",layout="fan")+
  geom_tippoint(aes(x=x+0.5,colour=status),shape=16,size=0.5)+
  scale_colour_manual(values=c("grey80","black"))+
  theme(legend.position = "None")

## add clades
for(i in 1:nrow(bpf2_results)){
  
  gg=gg+
    geom_hilight(node=bpf2_results$node[i],
                 alpha=0.25,
                 fill=ifelse(bpf2_results$clade>
                               bpf2_results$other,pcols[2],pcols[1])[i])
  # geom_cladelabel(node=bpf_results$node[i],
  #                 label=bpf_results$taxa[i],
  #                 offset=pplus,
  #                 hjust=0.5,
  #                 offset.text=pplus*1.25,
  #                 parse=T,
  #                 angle=90)
}
bgg2=gg

## extract states
state=setNames(set$data$status,rownames(set$data))

## set prior on root node as non-anthropogenic
rootn=c(1,0)
names(rootn)=levels(set$data$status)

## fit models
ER.model=fitMk(set$phy,state,model="ER",pi=as.matrix(rootn)) 
ARD.model=fitMk(set$phy,state,model="ARD",pi=as.matrix(rootn)) 
SYM.model=fitMk(set$phy,state,model="SYM",pi=as.matrix(rootn)) 
Irr.model=fitMk(set$phy,state,model=matrix(c(0,1,0,0),2,2,byrow=TRUE),
                 pi=as.matrix(rootn)) 

## model comparison
state.aov=anova(ER.model,ARD.model,SYM.model,Irr.model)
round(state.aov$weight,2)

## delta
state.aov=state.aov[order(state.aov$AIC,decreasing=F),]
state.aov$delta=state.aov$AIC-state.aov$AIC[1]

## plot
plot(ARD.model,width=T,color=T,offset=0.1,tol=0)

## simmap
state.simmap=simmap(state.aov,nsim=1000)

## summarize
simsum=summary(state.simmap)

## set colors and plot
cols=setNames(viridisLite::viridis(n=2),levels(set$data$synroot)) 
plot(simsum,ftype="i", fsize=0.7,colors=cols,cex=0.5)

## get density
state.density=density(state.simmap)

## plot
par(mfrow=c(1,2),mar=c(4.5,4.5,0.5,0.5))
COLS<-setNames(cols,state.density$trans) 
plot(state.density,transition=names(COLS)[1],colors=COLS[1],main="") 
mtext("a) transitions to anthropogenic roosting",line=-1) 

## other
plot(state.density,transition=names(COLS)[2],colors=COLS[2], main="") 
mtext("b) transitions to natural roosting", line=-1,adj=0)

## density map
state.densityMap=densityMap(state.simmap, plot=FALSE,res=20) 

## plot
state.densityMap=setMap(state.densityMap, viridisLite::viridis(n=10)) 
par(mfrow=c(1,1),mar=c(0.5,0.5,0.5,0.5))
plot(state.densityMap,lwd=1,outline=F,ftype="off",type="fan")

## add nodes
nodelabels(pie=simsum$ace,piecol=cols,cex=0.3)
