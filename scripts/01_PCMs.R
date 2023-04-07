## 01_phylogenetic analyses of bat roosting data
## danbeck@ou.edu
## last updated 04/06/2023

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
library(igraph)
library(ggraph)
library(tibble)

## load in roosting data
setwd("~/Desktop/synurbat/flat files")
data=readRDS("synurbic and traits only.rds")

## load in Upham phylogeny
setwd("~/Desktop/synurbat/phylos")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

## trim phylo to bats
tree=keep.tip(tree,data$tip)

## make label
data$label=data$tip

## merge
cdata=comparative.data(phy=tree,data=data,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## species
cdata$data$Species=cdata$data$tip

## make Synurbic with pseudoabsences
cdata$data$Synurbic=as.numeric(as.character(cdata$data$Synurbic))
cdata$data$Synurbic_pseudo=ifelse(is.na(cdata$data$Synurbic),0,cdata$data$Synurbic)

## taxonomy
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')

## save status and label
cdata$data$status=factor(cdata$data$Synurbic)
cdata$data$label=cdata$data$tip

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
        frmla.phylo=Synurbic~phylo+log1p(cites),
        family=binomial,algorithm='phylo',nfactors=3,min.group.size=5)

## summarize
bpf_results=pfsum(bpf)$results

## repeat for pseudo
set.seed(1)
bpf2=gpf(Data=cdata$data,tree=cdata$phy,
        frmla.phylo=Synurbic_pseudo~phylo+log1p(cites),
        family=binomial,algorithm='phylo',nfactors=3,min.group.size=5)

## summarize
bpf2_results=pfsum(bpf2)$results

## extract states for known status
state=setNames(set$data$status,rownames(set$data))

## set prior on root node as non-anthropogenic
rootn=c(1,0)
names(rootn)=levels(set$data$status)

## fit models
ER.model=fitMk(set$phy,state,model="ER",pi=as.matrix(rootn)) 
ARD.model=fitMk(set$phy,state,model="ARD",pi=as.matrix(rootn)) 
Irr.model=fitMk(set$phy,state,model=matrix(c(0,1,0,0),2,2,byrow=TRUE),
                pi=as.matrix(rootn)) 

## model comparison
state.aov=anova(ER.model,ARD.model,Irr.model)

## delta
state.aov=state.aov[order(state.aov$AIC,decreasing=F),]
state.aov$delta=state.aov$AIC-state.aov$AIC[1]
round(state.aov,2)

## plot
plot(ARD.model,width=T,color=T,offset=0.1,tol=0)

## plot in ggplot
qm=as.Qmatrix(ARD.model)
qm=matrix(as.numeric(qm),
          ncol=length(ARD.model$states),
          nrow=length(ARD.model$states))
rownames(qm)=ARD.model$states
colnames(qm)=ARD.model$states
diag(qm)=NA

## make edges
edges=data.frame(from=rev(ARD.model$states),
                 to=(ARD.model$states),
                 rate=qm[!is.na(qm)])

## convert to network with meta
g=graph_from_data_frame(edges,directed=T)
V(g)$name=c("anthropogenic","natural")
E(g)$rates=paste0(round(edges$rate,2))
E(g)$type=c("natural","anthropogenic")

## plot
asr_plot=ggraph(g, layout="linear")+
  geom_edge_arc(aes(edge_width=rate,
                    label=rates,
                    colour=type),
                label_dodge=unit(-5,'mm'),
                strength=0.65,
                label_size=2.5,
                angle_calc='along',
                arrow=arrow(length=unit(3,'mm')),
                start_cap=square(15,'mm'),
                end_cap=square(15,'mm'))+
  #geom_node_point(size=30,shape=16)+
  geom_node_text(aes(label=name,
                     colour=name),
                 size=2.5,fontface="bold",
                 nudge_x = c(0.06,-0.02))+
  scale_colour_manual(values=rev(c("palegreen3","wheat4")))+
  scale_edge_width_continuous(range=c(0.25,1.5))+
  scale_edge_color_manual(values=rev(c("palegreen3","wheat4")))+
  theme_void()+
  theme(legend.position = "none")

## simmap
state.simmap=simmap(state.aov,nsim=10)
#state.simmap=simmap(state.aov,nsim=1000)

## lineage-through-time
lobj=ltt(state.simmap,plot=F)
tmp=data.frame(time=lobj[[1]]$times,
               lobj[[1]]$ltt)
names(tmp)=c("time","natural","anthropogenic","total")
tmp$total=NULL

## wide to long
tmp=tidyr::gather(tmp,lin,number,natural:anthropogenic)

## visualize
ggplot(tmp,aes(time,number,colour=lin,group=lin))+
  geom_line()+
  scale_y_continuous(trans="log10")+
  theme_bw()

## summarize simmap
simsum=summary(state.simmap,plot=F)

## get states
top_prob=data.frame(simsum$ace)
names(top_prob)=c("neg","pos")

## make tibble
asrs=tibble(node=as.numeric(rownames(top_prob)),
            asr=top_prob$pos)
asrs=asrs[!is.na(asrs$node),]

## save tree
dtree=treeio::full_join(as.treedata(set$phy),set$data,by="label")

## join with asrs
dtree=treeio::full_join(dtree,asrs,by="node")

## circular tree for known data, with colors
circ=ggtree(dtree, layout="circular",aes(colour=asr),size=0.25)+
  scale_colour_gradient(low="palegreen3",high="wheat4")+
  guides(colour="none")

## add raw data into heatmap
tdat=as.data.frame(set$data$Synurbic)
rownames(tdat)=set$phy$tip.label
asr_tree=gheatmap(circ,tdat,offset=0.1,width=0.05,colnames=F,
                  low="palegreen3",high="wheat4")+
  theme(legend.position = "none")

## add phylofactor
asr_tree=asr_tree+
  geom_hilight(node=bpf_results$node[1],
               alpha=0.5,fill="grey90")+
  geom_hilight(node=bpf_results$node[2],
               alpha=0.5,fill="grey90")

## label
asr_tree=asr_tree+
  geom_cladelabel(node=bpf_results$node[1],
                  label="Pteropodidae", 
                  offset=5, offset.text = 5,
                  angle=-75, fontsize = 2)+
  geom_cladelabel(node=bpf_results$node[2],
                  label="Noctilionoidea", 
                  offset=5, offset.text = 5,
                  angle=45, fontsize = 2)

## combine plots
library(patchwork)
setwd("~/Desktop")
png("test.png",width=4,height=5,units="in",res=300)
asr_tree+asr_plot+plot_layout(ncol=1, heights=c(3,1))
dev.off()

## repeat for pseudoabsence dataset
cdata$data$status=as.character(cdata$data$Synurbic_pseudo)

## extract states for known status
state=setNames(cdata$data$status,rownames(cdata$data))

## set prior on root node as non-anthropogenic
rootn=c(1,0)
names(rootn)=levels(cdata$data$status)

## fit models
ER.model=fitMk(cdata$phy,state,model="ER",pi=as.matrix(rootn)) 
ARD.model=fitMk(cdata$phy,state,model="ARD",pi=as.matrix(rootn)) 
Irr.model=fitMk(cdata$phy,state,model=matrix(c(0,1,0,0),2,2,byrow=TRUE),
                pi=as.matrix(rootn)) 

## model comparison
state.aov=anova(ER.model,ARD.model,Irr.model)

## delta
state.aov=state.aov[order(state.aov$AIC,decreasing=F),]
state.aov$delta=state.aov$AIC-state.aov$AIC[1]
round(state.aov,2)

## plot
plot(ARD.model,width=T,color=T,offset=0.1,tol=0)

## plot in ggplot
qm=as.Qmatrix(ARD.model)
qm=matrix(as.numeric(qm),
          ncol=length(ARD.model$states),
          nrow=length(ARD.model$states))
rownames(qm)=ARD.model$states
colnames(qm)=ARD.model$states
diag(qm)=NA

## make edges
edges=data.frame(from=rev(ARD.model$states),
                 to=(ARD.model$states),
                 rate=qm[!is.na(qm)])

## convert to network with meta
g=graph_from_data_frame(edges,directed=T)
V(g)$name=c("anthropogenic","natural")
E(g)$rates=paste0(round(edges$rate,2))
E(g)$type=c("natural","anthropogenic")

## plot
asr_plot=ggraph(g, layout="linear")+
  geom_edge_arc(aes(edge_width=rate,
                    label=rates,
                    colour=type),
                label_dodge=unit(-5,'mm'),
                strength=0.65,
                label_size=2.5,
                angle_calc='along',
                arrow=arrow(length=unit(3,'mm')),
                start_cap=square(15,'mm'),
                end_cap=square(15,'mm'))+
  #geom_node_point(size=30,shape=16)+
  geom_node_text(aes(label=name,
                     colour=name),
                 size=2.5,fontface="bold",
                 nudge_x = c(0.06,-0.02))+
  scale_colour_manual(values=rev(c("palegreen3","wheat4")))+
  scale_edge_width_continuous(range=c(0.25,1.5))+
  scale_edge_color_manual(values=rev(c("palegreen3","wheat4")))+
  theme_void()+
  theme(legend.position = "none")

## simmap
state.simmap=simmap(state.aov,nsim=100)

## summarize
simsum=summary(state.simmap,plot=F)

## get states
top_prob=data.frame(simsum$ace)
names(top_prob)=c("neg","pos")

## make tibble
asrs=tibble(node=as.numeric(rownames(top_prob)),
            asr=top_prob$pos)
asrs=asrs[!is.na(asrs$node),]

## save tree
dtree=treeio::full_join(as.treedata(cdata$phy),set$data,by="label")

## join with asrs
dtree=treeio::full_join(dtree,asrs,by="node")

## circular tree for known data, with colors
circ=ggtree(dtree, layout="circular",aes(colour=asr),size=0.25)+
  scale_colour_gradient(low="palegreen3",high="wheat4")+
  guides(colour="none")

## add raw data into heatmap
tdat=as.data.frame(cdata$data$Synurbic_pseudo)
rownames(tdat)=cdata$phy$tip.label
asr_tree=gheatmap(circ,tdat,offset=0.1,width=0.05,colnames=F,
                  low="palegreen3",high="wheat4")+
  theme(legend.position = "none")

## add phylofactor
asr_tree=asr_tree+
  geom_hilight(node=bpf2_results$node[1],
               alpha=0.5,fill="grey90")+
  geom_hilight(node=bpf2_results$node[2],
               alpha=0.5,fill="grey90")

## label
asr_tree=asr_tree+
  geom_cladelabel(node=bpf2_results$node[1],
                  label="Pteropodidae", 
                  offset=5, offset.text = 5,
                  angle=-75, fontsize = 2)+
  geom_cladelabel(node=bpf2_results$node[2],
                  label="sub-Phyllostomidae", 
                  offset=5, offset.text = 5,
                  angle=45, fontsize = 2)

## combine plots
library(patchwork)
setwd("~/Desktop")
png("test2.png",width=4,height=5,units="in",res=300)
asr_tree+asr_plot+plot_layout(ncol=1, heights=c(3,1))
dev.off()
