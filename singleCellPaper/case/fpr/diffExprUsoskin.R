setwd("/Users/koenvandenberge/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/fpr/")
library(BiocParallel) ; library(Biobase)
load("../esetUsoskin.RData")
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
source("runScriptsUsoskin.R")

algos <-  list("zingeR_edgeR"=runEdgeREMLibSize, "edgeR"=runEdgeR, "limma-voom"=runVoom, "DESeq2"=runDESeq2, "MAST"=runMAST, "scde"=runScde, "limma-voom_hurdle"=runLimmaHurdle, "edgeR_hurdle"=runEdgeRHurdle, "metagenomeSeq"=runMetagenomeSeq, "zingeR_DESeq2"=runDESeq2Zero, "DESeq2_poscounts"=runDESeq2_poscounts, "MAST_count"=runMAST_count)

namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

nreps <- 30

#register(MulticoreParam(workers=1,verbose=TRUE, type="FORK"))

res <- lapply(1:nreps, function(i) {
#res <- bplapply(1:nreps, function(i) {   		    
    cat(i," ")
    eLoop <- eset[,as.numeric(subsets[i,])]
    pData(eLoop)$condition=factor(condition)    
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
boxplotData=data.frame(fpr=c(t(fprHat)),method=rep(namesAlgos,each=30))
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





