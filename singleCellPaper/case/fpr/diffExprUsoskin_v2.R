setwd("/Users/koenvandenberge/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/fpr/")
library(BiocParallel) ; library(Biobase)
load("../esetUsoskin.RData")


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


algos <-  list("zingeR_edgeR"=runEdgeREMLibSize, "edgeR"=runEdgeR, "limma-voom"=runVoom, "DESeq2"=runDESeq2, "MAST"=runMAST, "scde"=runScde, "metagenomeSeq"=runMetagenomeSeq, "zingeR_DESeq2"=runDESeq2Zero, "DESeq2_poscounts"=runDESeq2_poscounts)

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
#save(res,file="/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/FPR/res30Iter_highVariability_pickingSession_final.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/usoskin/FPR/res30Iter_highVariability_pickingSession_final.rda")

### FPR evaluation boxplots
hlp=lapply(res,function(replication){
     lapply(replication,function(method){
          pvalHlp=method$pvals
          pvalHlp[is.na(pvalHlp)]=1 #independent filtering
          mean(pvalHlp<=0.05)
})
})

fprHat=Reduce(hlp,f=cbind)
fprHat=matrix(unlist(fprHat),nrow=length(algos),ncol=nreps,byrow=FALSE)
rownames(fprHat)=namesAlgos

library(ggplot2)
#boxplotData=data.frame(fpr=unlist(c(hlp)),method=rep(namesAlgos,each=30))
boxplotData=data.frame(fpr=c(t(fprHat)),method=rep(namesAlgos,each=nreps))
p=ggplot(boxplotData,aes(x=reorder(method,fpr,median),y=fpr))
png("~/Dropbox/phdKoen/singleCell/figures/supplementary/fprCallsUsoskin_highVar_pickingSession.png", width=7,height=5, units="in", res=300)
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))
dev.off()

### without metagenomeSeq
boxplotDataSub=boxplotData[!boxplotData$method=="metagenomeSeq",]

p=ggplot(boxplotDataSub,aes(x=reorder(method,fpr,median),y=fpr))
png("~/Dropbox/phdKoen/singleCell/figures/supplementary/fprCallsUsoskin_highVar_pickingSession_noMetagenomeSeq.png", width=7,height=5, units="in", res=300)
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
    geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + xlab("") + theme(axis.text.y=element_text(size=12)) +
    scale_colour_discrete(guide="none")  + ylab("False Positive Rate") + geom_hline(aes(yintercept=0.01,colour="red"))
dev.off()

## combined
png("~/Dropbox/phdKoen/singleCell/figures/supplementary/fprCallsUsoskin_highVar_pickingSession_combined.png", width=11,height=9, units="in", res=300)
boxplotData$method=factor(boxplotData$method,levels=c("scde", "zingeR_DESeq2", "DESeq2_poscounts", "DESeq2", "edgeR", "zingeR_edgeR", "limma-voom", "MAST", "metagenomeSeq"))
hlp=boxplotData
hlp$fpr[hlp$method=="metagenomeSeq"]=NA
#boxplot(fpr~reorder(method,fpr,median),data=hlp)
par(mar=c(8,4.2,4,2.5))
boxplot(fpr~method,data=hlp, ylab="False positive rate", outline=FALSE, ylim=c(0,0.22), bty="l", cex.lab=1.5, xaxt="n")
#axis(1,at=1:9,labels=levels(boxplotData$method), las=2.5, cex.axis=1.33)
axis(1,at=1:9,labels=rep(" ",9))
text(cex=1.33,x=c(1:9)-0.15,y=-0.04,levels(boxplotData$method),xpd=TRUE,srt=45)
stripchart(fpr~method,data=hlp,vertical=TRUE,method="jitter",add=TRUE, pch=16, col="dimgray")
lines(x=c(0,8.5),y=c(0.05,0.05), col=2,lwd=2)

hlp2=boxplotData
hlp2$fpr[hlp2$method!="metagenomeSeq"]=NA
par(new="T")
boxplot(fpr~method,data=hlp2, yaxt="n", bty="l", xaxt="n")
stripchart(fpr~method,data=hlp2,vertical=TRUE,method="jitter",add=TRUE, pch=16, col="dimgray")
axis(4,at=seq(0,1,by=0.05), cex=1.33)
abline(v=8.5, lty=2, lwd=2)
lines(x=c(0,8.5),y=c(0.05,0.05), col=2,lwd=2)
dev.off()



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


par(mfrow=c(3,3))
for(i in 1:ncol(pDataset)) hist(pDataset[,i], main=colnames(pDataset)[i], xlim=c(0,1), breaks=seq(0,1,by=0.05))

## p-value histogram for the first iteration
png("~/Dropbox/phdKoen/singleCell/figures/supplementary/pValueHistogram_highVar_pickingSession_1iter.png", width=7,height=5, units="in", res=300)
par(mfrow=c(3,3))
o=c(2,1,3,4,9,8,5,6,7)
sapply(1:length(res[[1]]), function(i) hist(res[[1]][[o[i]]]$pvals, main=names(res[[1]])[o[i]], xlab="", breaks=seq(0,1,by=0.05)))
dev.off()


######## composite case study plot
library(Biobase)
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/esetUsoskin.RData")
usoskin=exprs(eset)
usoskin=usoskin[!rowSums(usoskin)==0,]
#cellTypeUsoskin=factor(pData(eset)[,7], levels=c("NF", "NP", "PEP", "TH"))
pickingSession=factor(pData(eset)[,2])
cellTypeAll=factor(pData(eset)[,"Level 3"],exclude=TRUE)
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v6.R") ## for weight functions
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

png("~/Dropbox/phdKoen/singleCell/figures/caseUsoskin_composite.png", width=12,height=11, units="in", res=350)
layout(matrix(c(1,2,3,3),nrow=2,ncol=2,byrow=TRUE))
## P(zero) ~ library size
expit <- function(x) 1/(1+exp(-x))
## not including batch effect
par(mar=c(5,4.5,4,1))
logLib=log(colSums(usoskin))
pZero=colMeans(usoskin==0)
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zeros", main= "", cex.lab=1.5, cex.main=1.5, bty="l", pch=1, cex.axis=1.33, col=as.numeric(pickingSession))
m <- glm(pZero~logLib,family="binomial")
a <- coef(m)[1]
b <- coef(m)[2]
lines(x=sort(logLib),y=expit(a+b*sort(logLib)),lwd=3,col="steelblue")
#legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,pch=1,cex=1.5)

## including main batch effect
#plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m2 <- glm(pZero~logLib+pickingSession,family="binomial")
a <- coef(m2)[1]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="Cold"]),y=expit(a+b*sort(logLib[pickingSession=="Cold"])),lwd=2,col=1)
a <- coef(m2)[1]+coef(m2)[3]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="RT-1"]),y=expit(a+b*sort(logLib[pickingSession=="RT-1"])),col=2,lwd=2)
a <- coef(m2)[1]+coef(m2)[4]
lines(x=sort(logLib[pickingSession=="RT-2"]),y=expit(a+b*sort(logLib[pickingSession=="RT-2"])),col=3,lwd=2)
legend("bottomleft",c("Global","Cold","RT-1","RT-2"),col=c("steelblue","black","red","green"),lty=1,cex=1.5, lwd=2, bty="n")
mtext("a" ,at=5.1, font=2, cex=4/3)

# Posterior probabilities
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/zingeR_edgeR/weightsConvergedZingeREdgeRUsoskin.rda")
par(mar=c(5,4.5,4,1))
hist(w[usoskin==0],xlab="Posterior probability",main="",cex.lab=1.5,cex.axis=1.33, ylim=c(0,3.5e6), yaxt="n", ylab="Frequency (x 1e5)")
axis(2,at=seq(0,3.5e6,by=5e5), labels=seq(0,3.5e6,by=5e5)/1e5, cex.axis=1.33)
#load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/weightsUsoskinNoBatchLibSizeDispFast300Iter.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/zingeR_edgeR_noBatch/weightsDESeq2ZeroUsoskin45.RData")
hist(w[usoskin==0],add=TRUE,col=rgb(0.1,0.8,0.1,.2))
legend("topleft",c("picking session + eff. lib. size","eff. lib. size"),fill=c(0,rgb(0.1,0.8,0.1,.2)), bty="n",cex=1.5)
mtext("b" ,at=-.13, font=2, cex=4/3)

# FPR calls
boxplotData$method=factor(boxplotData$method,levels=c("scde", "zingeR_DESeq2", "DESeq2_poscounts", "DESeq2", "edgeR", "zingeR_edgeR", "limma-voom", "MAST", "metagenomeSeq"))
hlp=boxplotData
hlp$fpr[hlp$method=="metagenomeSeq"]=NA
#boxplot(fpr~reorder(method,fpr,median),data=hlp)
par(mar=c(8,4.2,4,2.5))
boxplot(fpr~method,data=hlp, ylab="False positive rate", outline=FALSE, ylim=c(0,0.225), bty="l", cex.lab=1.5, xaxt="n", cex.axis=1.33)
#axis(1,at=1:9,labels=levels(boxplotData$method), las=2.5, cex.axis=1.33)
axis(1,at=1:9,labels=rep(" ",9))
text(cex=1.33,x=c(1:9)-0.15,y=-0.05,levels(boxplotData$method),xpd=TRUE,srt=45)
stripchart(fpr~method,data=hlp,vertical=TRUE,method="jitter",add=TRUE, pch=16, col="dimgray")
lines(x=c(0,8.5),y=c(0.05,0.05), col=2,lwd=2)

hlp2=boxplotData
hlp2$fpr[hlp2$method!="metagenomeSeq"]=NA
par(new="T")
boxplot(fpr~method,data=hlp2, yaxt="n", bty="l", xaxt="n")
stripchart(fpr~method,data=hlp2,vertical=TRUE,method="jitter",add=TRUE, pch=16, col="dimgray")
axis(4,at=seq(0,1,by=0.05), cex.axis=1.33)
abline(v=8.5, lty=2, lwd=2)
lines(x=c(0,8.5),y=c(0.05,0.05), col=2,lwd=2)
mtext("c" ,at=-.12, font=2, cex=4/3)
dev.off()
