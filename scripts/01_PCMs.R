## 01_phylogenetic analyses of bat roosting data
## danbeck@ou.edu
## last updated 03/24/2024

## clean environment & plots
rm(list=ls()) 
graphics.off()
gc()

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
library(phyloregion)
library(diversitree)
library(tidyr)
library(ggdist)
library(MuMIn)

## save theme
th=theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=11))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## load in roosting data
setwd("~/Desktop/synurbat/flat files")
setwd("/Users/danielbecker/Desktop/GitHub/synurbat/flat files")
data=readRDS("synurbic and traits only.rds")

## load in Upham phylogeny
setwd("~/Desktop/synurbat/phylos")
setwd("/Users/danielbecker/Desktop/GitHub/synurbat/phylos")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

## trim phylo to bats
tree=keep.tip(tree,data$tip)

## make label
data$label=data$tip

## state colors
scols=c("#8470FF","#9DD866")

## merge
cdata=comparative.data(phy=tree,data=data,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## species
cdata$data$Species=cdata$data$tip

## make Synurbic with pseudoabsences
cdata$data$Synurbic=as.numeric(as.character(cdata$data$Synurbic))
cdata$data$Synurbic_pseudo=ifelse(is.na(cdata$data$Synurbic),0,cdata$data$Synurbic)

## tally
tab=table(cdata$data$Synurbic)
round(tab["1"]/nrow(cdata$data),3)
round(tab["0"]/nrow(cdata$data),2)
tab

## taxonomy
cdata$data$taxonomy=paste(cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')

## save status and label
cdata$data$status=factor(cdata$data$Synurbic)
cdata$data$label=cdata$data$tip

## factor
cdata$data$sgroup=factor(cdata$data$Synurbic)

## pgls models for geographic range and size
mod=pgls(log10(X26.1_GR_Area_km2)~log10(adult_mass_g),data=cdata,lambda="ML"); summary(mod)
mod=pgls(log10(X26.1_GR_Area_km2)~adult_forearm_length_mm,data=cdata,lambda="ML"); summary(mod)

## pgls models for forearm
mod=pgls(adult_forearm_length_mm~det_fruit,data=cdata,lambda="ML"); summary(mod)
mod=pgls(adult_forearm_length_mm~det_inv,data=cdata,lambda="ML"); summary(mod)

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

## BiSSE
bstate=setNames(as.numeric(as.character(set$data$Synurbic)),rownames(set$data))
uphy=force.ultrametric(set$phy)
p=starting.point.bisse(uphy)
bmodel=make.bisse(tree=uphy,
                  states=bstate)
fit1=find.mle(bmodel,p)

## constrain all (equal rates model)
cmodel1=constrain(bmodel,lambda1~lambda0,mu1~mu0,q01~q10)

## constrain lambda and mu but not q (all rates different model)
cmodel2=constrain(bmodel,lambda1~lambda0,mu1~mu0)

## constrain lambda and mu, 0 to 1 only (irreversible model)
cmodel3=constrain(bmodel,lambda1~lambda0,mu1~mu0,q10~0)

## constrain  mu, 0 to 1 only (irreversible model)
cmodel4=constrain(bmodel,lambda1~lambda0,q10~0)

## constrain  lambda, 0 to 1 only (irreversible model)
cmodel5=constrain(bmodel,mu1~mu0,q10~0)

## constrain q only (irreversible)
cmodel6=constrain(bmodel,q10~0)

## constrain lambda only
cmodel7=constrain(bmodel,lambda1~lambda0)

## constrain mu only
cmodel8=constrain(bmodel,mu1~mu0)

## fit
cfit1=find.mle(cmodel1,p[argnames(cmodel1)])
cfit2=find.mle(cmodel2,p[argnames(cmodel2)])
cfit3=find.mle(cmodel3,p[argnames(cmodel3)])
cfit4=find.mle(cmodel4,p[argnames(cmodel4)])
cfit5=find.mle(cmodel5,p[argnames(cmodel5)])
cfit6=find.mle(cmodel6,p[argnames(cmodel6)])
cfit7=find.mle(cmodel7,p[argnames(cmodel7)])
cfit8=find.mle(cmodel8,p[argnames(cmodel8)])

## save in list
blist=list(fit1,cfit1,cfit2,cfit3,cfit4,cfit5,cfit6,cfit7,cfit8)
#anova(fit1,cfit1,cfit2,cfit3,cfit4,cfit5,cfit6,cfit7,cfit8)

## model comparison table
bcomp=data.frame(AIC=sapply(blist,AIC),
                 k=sapply(blist,function(x) length(coef(x))),
                 model=sapply(blist,function(x) paste(names(x$par),collapse=", ")))
bcomp=bcomp[order(bcomp$AIC),]
bcomp$delta=round(bcomp$AIC-bcomp$AIC[1],2)
bcomp$wi=round(Weights(bcomp$AIC),2)

## rearrange
bcomp=bcomp[c("model","k","delta","wi")]

## compare top models
anova(blist[[as.numeric(rownames(bcomp[1,]))]],
      blist[[as.numeric(rownames(bcomp[2,]))]])

## set top via parsimony
tmodel=attr(blist[[as.numeric(rownames(bcomp[1,]))]],"func")
top=find.mle(tmodel,p[argnames(tmodel)])

## MCMC framework, exponential prior as character-independent rate
prior=make.prior.exponential(1/(2*(p[1]-p[3])))

## assign w
set.seed(1)
tmp=mcmc(tmodel,top$par,nsteps=10,prior=prior,lower=0,w=rep(1,length(top$par)),print.every=0) 
w=diff(sapply(tmp[names(top$par)],range))

## run the chain to get 95% CIs
set.seed(1)
samples=mcmc(tmodel,top$par,nsteps=2000,w=w,lower=0,prior=prior,print.every=0)
samples$chain=1:nrow(samples)

## save
samples_raw=samples

## cut first 20% as burn-in
end=0.2*nrow(samples)
keep=(end+1):nrow(samples)
samples=samples[keep,]

## posterior package to summarize
library(posterior)
sdraw=data.frame(summarise_draws(samples,mean=mean,~ quantile(.x, probs = c(.025, .975)),rhat=rhat))

## minimize 
sdraw=sdraw[sdraw$variable%in%c("lambda0","lambda1","mu0","mu1","q01","q10"),]
rownames(sdraw)=sdraw$variable

## save q matrix for stochastic character mapping
qm=matrix(c(-sdraw["q01","mean"],sdraw["q10","mean"],sdraw["q01","mean"],-sdraw["q10","mean"]),2)
rownames(qm)=c("0","1")
colnames(qm)=rownames(qm)

## simplify
tmp=samples[names(top$par)]
tmp$sample=as.numeric(rownames(tmp))
tmp=gather(tmp,par,est,head(names(top$par),1):tail(names(top$par),1))

## check chains
ggplot(tmp,aes(sample,est))+
  facet_wrap(~par)+
  th+
  geom_line()

## density
ggplot(tmp,aes(est))+
  facet_wrap(~par)+
  th+
  geom_density(fill="grey")

## fix par
tmp$par2=revalue(tmp$par,
                 c("lambda0"="lambda[0]",
                   "lambda1"="lambda[1]",
                   "mu0"="mu[0]",
                   "mu1"="mu[1]",
                   "q01"="italic(q)[0][1]",
                   "q10"="italic(q)[1][0]"))

## order parameters
tmp$par=factor(tmp$par,levels=c("lambda0","lambda1","mu0","mu1","q10","q01"))

## assign similar levels for par2
tmp$par2=factor(tmp$par2,levels=unique(tmp[order(tmp$par),"par2"]))

## assign type
tmp$type=ifelse(tmp$par%in%c("lambda0","mu0","q10"),"natural","anthropogenic")
tmp$type=factor(tmp$type,levels=rev(unique(tmp$type)))

## parameter type
tmp$ptype=revalue(tmp$par,
                  c("mu0"="extinction",
                    "mu1"="extinction",
                    "lambda0"="speciation",
                    "lambda1"="speciation",
                    "q10"="transition",
                    "q01"="transition"))

## factor
tmp$ptype=factor(tmp$ptype,levels=c("speciation","extinction","transition"))

## ggdist
bisse_plot=ggplot(tmp,aes(par2,est,colour=type,fill=type))+
  coord_flip()+
  facet_wrap(~ptype,ncol=1,scales="free_y",strip.position="right",shrink=T)+
  stat_halfeye(slab_alpha=0.25,point_size=2,point_interval="mean_qi")+
  th+
  labs(y="Posterior Samples",x="BiSSE Parameters")+
  scale_x_discrete(labels=scales::label_parse())+
  theme(axis.text.y=element_text(size=12))+
  guides(colour="none",fill="none")+
  scale_colour_manual(values=rev(scols))+
  scale_fill_manual(values=rev(scols))

## set prior on root node as non-anthropogenic
rootn=c(1,0)
names(rootn)=levels(factor(set$data$Synurbic))

## clean large files
rm(mod1,fit1,cfit1,cfit2,cfit3,cfit4,cfit5,cfit6,cfit7,cfit8)

## test
set.seed(1)
mk=make.simmap(uphy,bstate,Q=qm,nsim=1000,pi=as.matrix(rootn))

## summarize simmap
simsum=summary(mk,plot=F)

## lineage-through-time
lobj=ltt(mk,plot=F)

## clean
rm(mk)

## add sim identity and save as data
lset=lapply(1:length(lobj),function(x){
  
  ## get from list
  tt=lobj[[x]]
  
  ## times
  times=tt$times
  
  ## data
  lset=data.frame(tt$ltt)
  lset$total=NULL
  
  ## data frame
  lset=data.frame(times=times,
             lset)
  names(lset)=c("times","natural","anthropogenic")
  
  ## add sim
  lset$sim=x
  
  ## return
  return(lset)
  ## resave as list
  #return(list(times=times,ltt=ltt))

})
  
## save as dataset
lset=do.call(rbind.data.frame,lset)

## wide to long
lset=tidyr::gather(lset,lin,number,natural:anthropogenic)
lset$lin=factor(lset$lin,levels=c("natural","anthropogenic"))

## get means
H=max(lset$times)
TIMES<-seq(0,max(H),length.out=10000)
LINEAGES<-matrix(0,length(TIMES),length(levels(set$data$sgroup))+1)

## iterate over times and ltts
for(i in 1:length(TIMES)){
  for(j in 1:length(lobj)){
    ii<-which(lobj[[j]]$times<=TIMES[i])
    ADD<-if(length(ii)==0) rep(0,length(levels(set$data$sgroup))) else 
      lobj[[j]]$ltt[max(ii),]/length(lobj)
    LINEAGES[i,]<-LINEAGES[i,]+ADD
  }
}

## save means
mset=data.frame(TIMES,LINEAGES)
mset$X3=NULL
names(mset)=c("times","natural","anthropogenic")
mset=tidyr::gather(mset,lin,number,natural:anthropogenic)
mset$lin=factor(mset$lin,levels=c("natural","anthropogenic"))

## visualize
ltt_plot=ggplot()+
  geom_line(data=lset,aes(times,number,colour=lin,group=sim),alpha=0.1,linewidth=0.05)+
  geom_line(data=mset,aes(times,number,colour=lin),linewidth=1.25)+
  scale_y_continuous(trans = scales::pseudo_log_trans(2, 1000), breaks=c(1,10,100,1000), limits=c(0,1000), 
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  th+
  scale_colour_manual(values=scols)+
  guides(colour="none")+
  labs(x="Relative Time",y="Number of Lineages")

## extract statistics
mean(round(sapply(lobj,function(x) x$gamma),2))
var(round(sapply(lobj,function(x) x$gamma),2))
mean(round(sapply(lobj,function(x) x$p),3))

## get states from simsum
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
circ=ggtree(dtree, layout="rectangular",
            ladderize=T,right=T,
            aes(colour=asr),size=0.25)+
  scale_colour_gradient(low=scols[1],high=scols[2])+
  guides(colour="none")

## add raw data into heatmap
tdat=as.data.frame(set$data$Synurbic)
rownames(tdat)=set$phy$tip.label
asr_tree=gheatmap(circ,tdat,offset=0.1,width=0.025,colnames=F,
                  low=scols[1],high=scols[2],color=NA)+
  theme(legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.margin=margin(t=-2.5,b=-2.5),
        legend.title=element_text(size=10))+
  guides(colour=guide_colorbar(barwidth=8,barheight=0.75,
                               title=expression(paste(italic(P),"(Anthropogenic Roosting)")),
                               title.vjust=1),
         fill="none")

## add phylofactor
asr_tree=asr_tree+
  geom_hilight(node=bpf_results$node[1],
               fill=NA,colour="black",linewith=2)+
  geom_hilight(node=bpf_results$node[2],
               fill=NA,colour="black",linewith=2)

## label
asr_tree=asr_tree+
  geom_cladelabel(node=bpf_results$node[1],
                  label="Pteropodidae",hjust=0.55,
                  offset=2.5, offset.text = 2.5,
                  angle=90, fontsize = 3)+
  geom_cladelabel(node=bpf_results$node[2],
                  label="Stenodermatinae",hjust=0.5,
                  offset=2.5, offset.text = 2.5,
                  angle=90, fontsize = 3)

## combine bisse and ltt
library(ggpubr)
#bt=egg::ggarrange(bisse_plot,ltt_plot,ncol=1,heights=c(1.5,1.1))
bt=ggpubr::ggarrange(bisse_plot,ltt_plot,
                     ncol=1,heights=c(1.5,1.1),
                     labels=c("B","C"),
                     font.label=list(size=11,face="plain"),
                     hjust=-0.7,
                     vjust=c(2,0.1))

## combine all plots
library(patchwork)
setwd("/Users/danielbecker/Desktop/synurbat/figures")
setwd("/Users/danielbecker/Desktop/GitHub/synurbat/figures")
png("Figure 4.png",width=7,height=5,units="in",res=300)
#(asr_tree|bt)+plot_layout(widths=c(1.5,1))
ggpubr::ggarrange(asr_tree,bt,
                  ncol=2,widths=c(1.5,1),
                  labels=c("A"),
                  font.label=list(size=11,face="plain"),
                  vjust=2)
dev.off()

## repeat BiSSE for pseudoabsence dataset
bstate=setNames(as.numeric(as.character(cdata$data$Synurbic_pseudo)),rownames(cdata$data))
uphy=as.ultrametric(cdata$phy)
p=starting.point.bisse(uphy)
bmodel=make.bisse(tree=uphy,
                  states=bstate)
fit1=find.mle(bmodel,p)

## constrain all (equal rates model)
cmodel1=constrain(bmodel,lambda1~lambda0,mu1~mu0,q01~q10)

## constrain lambda and mu but not q (all rates different model)
cmodel2=constrain(bmodel,lambda1~lambda0,mu1~mu0)

## constrain lambda and mu, 0 to 1 only (irreversible model)
cmodel3=constrain(bmodel,lambda1~lambda0,mu1~mu0,q10~0)

## constrain  mu, 0 to 1 only (irreversible model)
cmodel4=constrain(bmodel,lambda1~lambda0,q10~0)

## constrain  lambda, 0 to 1 only (irreversible model)
cmodel5=constrain(bmodel,mu1~mu0,q10~0)

## constrain q only (irreversible)
cmodel6=constrain(bmodel,q10~0)

## constrain lambda only
cmodel7=constrain(bmodel,lambda1~lambda0)

## constrain mu only
cmodel8=constrain(bmodel,mu1~mu0)

## fit
cfit1=find.mle(cmodel1,p[argnames(cmodel1)])
cfit2=find.mle(cmodel2,p[argnames(cmodel2)])
cfit3=find.mle(cmodel3,p[argnames(cmodel3)])
cfit4=find.mle(cmodel4,p[argnames(cmodel4)])
cfit5=find.mle(cmodel5,p[argnames(cmodel5)])
cfit6=find.mle(cmodel6,p[argnames(cmodel6)])
cfit7=find.mle(cmodel7,p[argnames(cmodel7)])
cfit8=find.mle(cmodel8,p[argnames(cmodel8)])

## save in list
blist=list(fit1,cfit1,cfit2,cfit3,cfit4,cfit5,cfit6,cfit7,cfit8)

## model comparison table
bcomp2=data.frame(AIC=sapply(blist,AIC),
                 k=sapply(blist,function(x) length(coef(x))),
                 model=sapply(blist,function(x) paste(names(x$par),collapse=", ")))
bcomp2=bcomp2[order(bcomp2$AIC),]
bcomp2$delta=round(bcomp2$AIC-bcomp2$AIC[1],2)
bcomp2$wi=round(Weights(bcomp2$AIC),2)

## rearrange
bcomp2=bcomp2[c("model","k","delta","wi")]

## compare top models
anova(blist[[as.numeric(rownames(bcomp2[1,]))]],
      blist[[as.numeric(rownames(bcomp2[2,]))]])

## set top via parsimony
tmodel=attr(blist[[as.numeric(rownames(bcomp2[1,]))]],"func")
top=find.mle(tmodel,p[argnames(tmodel)])

## MCMC framework, exponential prior as character-independent rate
prior=make.prior.exponential(1/(2*(p[1]-p[3])))

## assign w
set.seed(1)
tmp=mcmc(tmodel,top$par,nsteps=100,prior=prior,lower=0,w=rep(1,length(top$par)),print.every=0) 
w=diff(sapply(tmp[names(top$par)],range))

## run the chain to get 95% CIs
set.seed(1)
samples=mcmc(tmodel,top$par,nsteps=2000,w=w,lower=0,prior=prior,print.every=0)
samples$chain=1:nrow(samples)

## save
samples_raw=samples

## cut first 20% as burn-in
end=0.2*nrow(samples)
keep=(end+1):nrow(samples)
samples=samples[keep,]

## posterior package to summarize
library(posterior)
sdraw=data.frame(summarise_draws(samples,mean=mean,~ quantile(.x, probs = c(.025, .975))))

## minimize 
sdraw=sdraw[sdraw$variable%in%c("lambda0","lambda1","mu0","mu1","q01","q10"),]
rownames(sdraw)=sdraw$variable

## simplify
tmp=samples[names(top$par)]
tmp=gather(tmp,par,est,head(names(top$par),1):tail(names(top$par),1))

## fix par
tmp$par2=revalue(tmp$par,
                 c("lambda0"="lambda[0]",
                   "lambda1"="lambda[1]",
                   "mu0"="mu[0]",
                   "mu1"="mu[1]",
                   "q01"="italic(q)[0][1]",
                   "q10"="italic(q)[1][0]"))

## order parameters
tmp$par=factor(tmp$par,levels=c("lambda0","lambda1","mu0","mu1","q10","q01"))

## assign similar levels for par2
tmp$par2=factor(tmp$par2,levels=unique(tmp[order(tmp$par),"par2"]))

## assign type
tmp$type=ifelse(tmp$par%in%c("lambda0","mu0","q10"),"natural","anthropogenic")
tmp$type=factor(tmp$type,levels=rev(unique(tmp$type)))

## parameter type
tmp$ptype=revalue(tmp$par,
                  c("mu0"="extinction",
                    "mu1"="extinction",
                    "lambda0"="speciation",
                    "lambda1"="speciation",
                    "q10"="transition",
                    "q01"="transition"))

## factor
tmp$ptype=factor(tmp$ptype,levels=c("speciation","extinction","transition"))

## ggdist
#setwd("/Users/danielbecker/Desktop/synurbat/figures")
png("Figure S3.png",width=4,height=4,units="in",res=300)
ggplot(tmp,aes(par2,est,colour=type,fill=type))+
  coord_flip()+
  facet_wrap(~ptype,ncol=1,scales="free_y",strip.position="right",shrink=T)+
  stat_halfeye(slab_alpha=0.25,point_size=2,point_interval="mean_qi")+
  th+
  labs(y="Posterior Samples",x="BiSSE Parameters")+
  scale_x_discrete(labels=scales::label_parse())+
  theme(axis.text.y=element_text(size=12),
        strip.text=element_blank(),
        strip.background=element_blank())+
  theme(legend.position="top")+
  #guides(colour="none",fill="none")+
  scale_colour_manual(values=rev(scols))+
  scale_fill_manual(values=rev(scols))
dev.off()