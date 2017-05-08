setwd("/Users/koenvandenberge/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/fpr/")
library(BiocParallel) ; library(Biobase)
load("../esetUsoskin.RData")


#### low variability setting: one picking session, one celltype
## only get cells from picking session Cold
eset=eset[,pData(eset)[,2]=="Cold"]
## celltype TH
eset=eset[,pData(eset)[,7]=="TH"]
cellType=factor(pData(eset)[,7], levels="TH")
## remove one sample from different lib
eset=eset[,!pData(eset)[,4]=="L281"]
condition=factor(rep(c("A","B"),each=45))
pData(eset)$condition=factor(condition)
#exprs(eset)=exprs(eset)[rowSums(exprs(eset))>0,]
eset=eset[rowSums(exprs(eset)>0)>=10,]
exprs(eset) <- apply(exprs(eset),2,function(x) {storage.mode(x) <- 'integer'; x})
subsets <- read.table("subsetMatrixUsoskinFPR.txt")
library(DESeq2) ; library(edgeR) ; library(limma) ; library(scde) ; library(MAST) ; library(genefilter) ; library(metagenomeSeq)
source("runScriptsUsoskin_pickingSession.R")

algos <-  list("zingeR_edgeR"=runEdgeREMLibSize, "edgeR"=runEdgeR, "limma-voom"=runVoom, "DESeq2"=runDESeq2, "MAST"=runMAST, "scde"=runScde, "metagenomeSeq"=runMetagenomeSeq, "zingeR_DESeq2"=runDESeq2Zero, "DESeq2_poscounts"=runDESeq2_poscounts, "MAST_count"=runMAST_count)

namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

nreps <- 3

#register(MulticoreParam(workers=1,verbose=TRUE, type="FORK"))

res <- lapply(1:nreps, function(i) {
#res <- bplapply(1:nreps, function(i) {   		    
    cat(i," ")
    eLoop <- eset[,as.numeric(subsets[i,])]
    condition=factor(rep(c("A","B"),each=45))
    pickingSession=factor(rep(rep(c("Cold","rt1","rt2"),each=15),2))
    pData(eLoop)$condition=condition
    pData(eLoop)$pickingSession=pickingSession  
    resFPR <- lapply(namesAlgos, function(n) algos[[n]](eLoop))
    return(resFPR)
})
#save(res,file="/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/FPR/res30Iter.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/FPR/res30Iter.rda") #precomputed
hlp=lapply(res,function(replication){
	   lapply(replication,function(method){
		      pvalHlp=method$pvals
		      pvalHlp[is.na(pvalHlp)]=1 #independent filtering
		      mean(pvalHlp<=0.01)
})
})

fprHat=Reduce(hlp,f=cbind)
#fprHat[1,]=unlist(lapply(res,function(x){ 
#			     pvalHlp=x$edgeRLibSize[[1]]
#			     pvalHlp[is.na(pvalHlp)]=1
#			     mean(pvalHlp<=0.01)
#}))
fprHat=matrix(unlist(fprHat),nrow=length(algos),ncol=nreps,byrow=FALSE)
rownames(fprHat)=namesAlgos

library(ggplot2)
#boxplotData=data.frame(fpr=unlist(c(hlp)),method=rep(namesAlgos,each=30))
boxplotData=data.frame(fpr=c(t(fprHat)),method=rep(namesAlgos,each=nreps))
p=ggplot(boxplotData,aes(x=reorder(method,fpr,median),y=fpr))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") +
    scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))

## without the hurdles
boxplotDataNoHurdle = boxplotData[!boxplotData$method%in%c("limma-voom_hurdle","edgeR_hurdle"),]
p=ggplot(boxplotDataNoHurdle,aes(x=reorder(method,fpr,median),y=fpr))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + xlab("") + theme(axis.text.y=element_text(size=12)) +
    scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))

### without metagenomeSeq
boxplotDataSub=boxplotDataNoHurdle[!boxplotDataNoHurdle$method=="metagenomeSeq",]

p=ggplot(boxplotDataSub,aes(x=reorder(method,fpr,median),y=fpr))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + xlab("") + theme(axis.text.y=element_text(size=12)) +
    scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))



### aggregate p-values and make histograms
pvalMethods=lapply(res,function(replication){
       do.call(cbind,lapply(replication,function(method){
              pvalHlp=method$pvals
                return(pvalHlp)
}))
})

pDataset=matrix(NA,nrow=nreps*nrow(eset),ncol=length(namesAlgos))
colnames(pDataset)=colnames(pvalMethods[[1]])
for(i in 1:length(namesAlgos)) pDataset[,i] = do.call(cbind,lapply(pvalMethods, function(x) x[,i]))


par(mfrow=c(4,3))
for(i in 1:ncol(pDataset)) hist(pDataset[,i], main=colnames(pDataset)[i])


#### medium variability setting: one celltype, three picking sessions (15 samples each) in every condition
################ NP celltype
setwd("/Users/koenvandenberge/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/fpr/")
library(BiocParallel) ; library(Biobase)
load("../esetUsoskin.RData")
eset=eset[,pData(eset)[,7]=="NP"]
#condition=factor(rep(c("A","B"),each=45))
#pData(eset)$condition=factor(condition)
eset=eset[rowSums(exprs(eset)>0)>=10,]
exprs(eset) <- apply(exprs(eset),2,function(x) {storage.mode(x) <- 'integer'; x})

subsets <- read.table("subsetMatrixUsoskinFPR_cellTypeNP.txt")
library(DESeq2) ; library(edgeR) ; library(limma) ; library(scde) ; library(MAST) ; library(genefilter) ; library(metagenomeSeq)
source("runScriptsUsoskin_pickingSession.R")

algos <-  list("zingeR_edgeR"=runEdgeREMLibSize, "edgeR"=runEdgeR, "limma-voom"=runVoom, "DESeq2"=runDESeq2, "MAST"=runMAST, "scde"=runScde, 
   "metagenomeSeq"=runMetagenomeSeq, "zingeR_DESeq2"=runDESeq2Zero, "DESeq2_poscounts"=runDESeq2_poscounts, "MAST_count"=runMAST_count)

namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

nreps <- 3

res <- lapply(1:nreps, function(i) {
#res <- bplapply(1:nreps, function(i) {             
    cat(i," ")
    eLoop <- eset[,as.numeric(subsets[i,])]
    condition=factor(rep(c("A","B"),each=45))
    pickingSession=factor(rep(rep(c("Cold","rt1","rt2"),each=15),2))
    pData(eLoop)$condition=condition
    pData(eLoop)$pickingSession=pickingSession  
    resFPR <- lapply(namesAlgos, function(n) algos[[n]](eLoop))
    return(resFPR)
})
#save(res,file="/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/FPR/res3IterNP_pickingSession.rda")

hlp=lapply(res,function(replication){
       lapply(replication,function(method){
              pvalHlp=method$pvals
              pvalHlp[is.na(pvalHlp)]=1 #independent filtering
              mean(pvalHlp<=0.01)
})
})

fprHat=Reduce(hlp,f=cbind)
#fprHat[1,]=unlist(lapply(res,function(x){ 
#                pvalHlp=x$edgeRLibSize[[1]]
#                pvalHlp[is.na(pvalHlp)]=1
#                mean(pvalHlp<=0.01)
#}))
fprHat=matrix(unlist(fprHat),nrow=length(algos),ncol=nreps,byrow=FALSE)
rownames(fprHat)=namesAlgos

library(ggplot2)
boxplotData=data.frame(fpr=c(t(fprHat)),method=rep(namesAlgos,each=nreps))
p=ggplot(boxplotData,aes(x=reorder(method,fpr,median),y=fpr))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") +
    scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))

### without metagenomeSeq
boxplotDataSub=boxplotData[!boxplotData$method=="metagenomeSeq",]

p=ggplot(boxplotDataSub,aes(x=reorder(method,fpr,median),y=fpr))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() + geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + xlab("") + theme(axis.text.y=element_text(size=12)) +scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))


### p-value distributions
pvalMethods=lapply(res,function(replication){
       do.call(cbind,lapply(replication,function(method){
              pvalHlp=method$pvals
                return(pvalHlp)
}))
})

pDataset=matrix(NA,nrow=3*nrow(eset),ncol=length(algos))
colnames(pDataset)=colnames(pvalMethods[[1]])
for(i in 1:length(algos)) pDataset[,i] = do.call(cbind,lapply(pvalMethods, function(x) x[,i]))

par(mfrow=c(4,3))
for(i in 1:ncol(pDataset)) hist(pDataset[,i], main=colnames(pDataset)[i])


#### high variability setting: all celltypes, three picking sessions (15 samples each) in every condition.
setwd("/Users/koenvandenberge/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/fpr/")
library(BiocParallel) ; library(Biobase)
load("../esetUsoskin.RData")
#eset=eset[,pData(eset)[,7]=="NP"]
#condition=factor(rep(c("A","B"),each=45))
#pData(eset)$condition=factor(condition)
#at least 20 positive counts
eset=eset[rowSums(exprs(eset)>0)>=20,]
exprs(eset) <- apply(exprs(eset),2,function(x) {storage.mode(x) <- 'integer'; x})

subsets <- read.table("subsetMatrixUsoskinFPR_randomCellTypes.txt")
library(DESeq2) ; library(edgeR) ; library(limma) ; library(scde) ; library(MAST) ; library(genefilter) ; library(metagenomeSeq)
source("runScriptsUsoskin_pickingSession.R")


algos <-  list("zingeR_edgeR"=runEdgeREMLibSize, "edgeR"=runEdgeR, "limma-voom"=runVoom, "DESeq2"=runDESeq2, "MAST"=runMAST, "scde"=runScde, "metagenomeSeq"=runMetagenomeSeq, "zingeR_DESeq2"=runDESeq2Zero, "DESeq2_poscounts"=runDESeq2_poscounts, "MAST_count"=runMAST_count)

namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

nreps <- 30

res <- lapply(1:nreps, function(i) {
#res <- bplapply(1:nreps, function(i) {             
    cat(i," ")
    eLoop <- eset[,as.numeric(subsets[i,])]
    condition=factor(rep(c("A","B"),each=45))
    pickingSession=factor(rep(rep(c("Cold","rt1","rt2"),each=15),2))
    pData(eLoop)$condition=condition
    pData(eLoop)$pickingSession=pickingSession  
    resFPR <- lapply(namesAlgos, function(n) algos[[n]](eLoop))
    return(resFPR)
})
#save(res,file="/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/FPR/res30Iter_highVariability_pickingSession.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/FPR/res30Iter_highVariability_pickingSession.rda")

### FPR evaluation boxplots
hlp=lapply(res,function(replication){
     lapply(replication,function(method){
          pvalHlp=method$pvals
          pvalHlp[is.na(pvalHlp)]=1 #independent filtering
          mean(pvalHlp<=0.01)
})
})

fprHat=Reduce(hlp,f=cbind)
fprHat=matrix(unlist(fprHat),nrow=length(algos),ncol=nreps,byrow=FALSE)
rownames(fprHat)=namesAlgos

library(ggplot2)
#boxplotData=data.frame(fpr=unlist(c(hlp)),method=rep(namesAlgos,each=30))
boxplotData=data.frame(fpr=c(t(fprHat)),method=rep(namesAlgos,each=nreps))
p=ggplot(boxplotData,aes(x=reorder(method,fpr,median),y=fpr))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))

### without metagenomeSeq
boxplotDataSub=boxplotData[!boxplotData$method=="metagenomeSeq",]

p=ggplot(boxplotDataSub,aes(x=reorder(method,fpr,median),y=fpr))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + xlab("") + theme(axis.text.y=element_text(size=12)) +
    scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))


### aggregate p-values and make histograms
pvalMethods=lapply(res,function(replication){
       do.call(cbind,lapply(replication,function(method){
              pvalHlp=method$pvals
                return(pvalHlp)
}))
})

pDataset=matrix(NA,nrow=nreps*nrow(eset),ncol=length(namesAlgos))
colnames(pDataset)=colnames(pvalMethods[[1]])
for(i in 1:length(namesAlgos)) pDataset[,i] = do.call(cbind,lapply(pvalMethods, function(x) x[,i]))


par(mfrow=c(4,3))
for(i in 1:ncol(pDataset)) hist(pDataset[,i], main=colnames(pDataset)[i])

## p-value histogram for the first iteration
par(mfrow=c(5,2))
o=c(2,1,3,6,5,10,7,4,9,8)
sapply(1:length(res[[1]]), function(i) hist(res[[1]][[o[i]]]$pvals, main=names(res[[1]])[o[i]], xlab=""))




