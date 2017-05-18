library(Biobase)
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/esetUsoskin.RData")
usoskin=exprs(eset)
usoskin=usoskin[!rowSums(usoskin)==0,]
#cellTypeUsoskin=factor(pData(eset)[,7], levels=c("NF", "NP", "PEP", "TH"))
pickingSession=factor(pData(eset)[,2])
cellTypeAll=factor(pData(eset)[,"Level 3"],exclude=TRUE)
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v6.R") ## for weight functions
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")


## P(zero) ~ library size
expit <- function(x) 1/(1+exp(-x))
## not including batch effect
par(mar=c(5,4.5,4,1))
logLib=log(colSums(usoskin))
pZero=colMeans(usoskin==0)
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m <- glm(pZero~logLib,family="binomial")
a <- coef(m)[1]
b <- coef(m)[2]
lines(x=sort(logLib),y=expit(a+b*sort(logLib)),lwd=2,col="steelblue")
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,pch=1,cex=1.5)

## including main batch effect
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m2 <- glm(pZero~logLib+pickingSession,family="binomial")
a <- coef(m2)[1]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="Cold"]),y=expit(a+b*sort(logLib[pickingSession=="Cold"])),lwd=2,col=1)
a <- coef(m2)[1]+coef(m2)[3]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="RT-1"]),y=expit(a+b*sort(logLib[pickingSession=="RT-1"])),col=2,lwd=2)
a <- coef(m2)[1]+coef(m2)[4]
lines(x=sort(logLib[pickingSession=="RT-2"]),y=expit(a+b*sort(logLib[pickingSession=="RT-2"])),col=3,lwd=2)
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,lty=1,cex=1.5)


## including batch effect as interaction with logLib
plot(x=logLib,y=pZero, xlab="log library size", ylab="fraction of zeroes", main= "Usoskin et al. 2015", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m2 <- glm(pZero~logLib*pickingSession,family="binomial")
a <- coef(m2)[1]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="Cold"]),y=expit(a+b*sort(logLib[pickingSession=="Cold"])),lwd=2,col=1)
a <- coef(m2)[1]+coef(m2)[3]
b <- coef(m2)[2]+coef(m2)[5]
lines(x=sort(logLib[pickingSession=="RT-1"]),y=expit(a+b*sort(logLib[pickingSession=="RT-1"])),col=2,lwd=2)
a <- coef(m2)[1]+coef(m2)[4]
b <- coef(m2)[2]+coef(m2)[6]
lines(x=sort(logLib[pickingSession=="RT-2"]),y=expit(a+b*sort(logLib[pickingSession=="RT-2"])),col=3,lwd=2)


### ANALYSES
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}
zeroWeightsLibSizeDispFastUsoskin <- function(counts, design, designFormula, initialWeightAt0=TRUE, maxit=100, plot=FALSE, plotW=FALSE, designZI=NULL, llTol=1e-4, normalization="TMM"){
    require(edgeR) ; require(DESeq2)
    if(plot | plotW) par(mfrow=c(1,plot+plotW))    
    counts <- DGEList(counts)
    #counts <- edgeR::calcNormFactors(counts)
    if(normalization=="TMM"){
    	counts = edgeR::calcNormFactors(counts)
    } else if(normalization=="DESeq2_default"){
    	dse = DESeqDataSetFromMatrix(counts$counts, colData=data.frame(cellType=cellTypeAll, pickingSession=pickingSession), design=designFormula)
    	 dse = DESeq2::estimateSizeFactors(dse)    	 
    	 counts$samples$norm.factors = 1/dse$sizeFactor
    } else if(normalization=="DESeq2_poscounts"){
		dse = DESeqDataSetFromMatrix(counts$counts, colData=data.frame(cellType=cellTypeAll, pickingSession=pickingSession), design=designFormula)
    	 dse = DESeq2::estimateSizeFactors(dse, type="poscounts")
    	 counts$samples$norm.factors = 1/dse$sizeFactor
    }
    effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
    logEffLibSize <- log(effLibSize)
    zeroId <- counts$counts==0
    w <- matrix(1,nrow=nrow(counts),ncol=ncol(counts), dimnames=list(c(1:nrow(counts)), NULL))
      ## starting values based on P(zero) in the library
    for(k in 1:ncol(w)) w[counts$counts[,k]==0,k] <- 1-mean(counts$counts[,k]==0)
    
    llOld <- matrix(-1e4,nrow=nrow(counts),ncol=ncol(counts))
    likCOld <- matrix(0,nrow=nrow(counts),ncol=ncol(counts))
    converged=FALSE
    j=0

    for(i in 1:maxit){
	j=j+1
    zeroId <- counts$counts==0	
	counts$weights <- w
	
	### M-step counts
	#only estimate dispersions every 5 iterations
	#if(i==1 | is.wholenumber(i/10)){
	if(i==1 | converged){
	#counts <- estimateGLMCommonDisp(counts, design, interval=c(0,10))
	#counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0, min.row.sum=1)
	counts = estimateDisp(counts, design, prior.df=0, min.row.sum=1)
	}
	if(plot) plotBCV(counts)	
	fit <- glmFit(counts, design)
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)
	
	### M-step mixture parameter: model zero probability
	successes <- colSums(1-w) #P(zero)
	failures <- colSums(w) #1-P(zero)
	if(is.null(designZI)){
	zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")} else{
	zeroFit <- glm(cbind(successes,failures) ~-1+designZI, family="binomial")}
	pi0Hat <- predict(zeroFit,type="response") 
	
	## E-step: Given estimated parameters, calculate expected value of weights
	pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
	w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)

	## data log-likelihood
	if(i>1) llOld=ll
	ll <- log(pi0HatMat*zeroId + (1-pi0HatMat)*likC + 1e-6)

	delta <- (rowSums(ll)-rowSums(llOld))/(rowSums(llOld)+llTol)
	if(mean(abs(delta) < llTol)>.999){ #if 99.9% has converged
	    if(j==1 & mean(abs(delta) < llTol)>.999){ #final convergence?
	    	cat(paste0("converged. \n")) ; return(w)}
	    j=0 
	    converged=TRUE} else {converged=FALSE}
	cat(paste0("iteration: ",i,". mean conv.: ",mean(abs(delta) < llTol),"\n"))

	#plot weights and make new pdf doc every 10 iterations to check convergence
	if(plotW){
		if(i==1) pdf(paste0("weightsHistogramDESeq2Zero",i,".pdf"))
		if(is.wholenumber(i/10)){
			dev.off()
			pdf(paste0("weightsHistogramDESeq2Zero",(i/10)+1,".pdf"))
			} 
		hist(w[zeroId],main=paste0("iteration: ",i,". mean conv.: ",mean(abs(delta) < llTol)))
	}

	if(is.wholenumber(i/10)){
		save(w, file=paste0("weightsDESeq2ZeroUsoskin",(i/10)+1,".RData"))
	}

    }
    dev.off()
    return(w)
}


### including the batch effect
dUsoskin=DGEList(usoskin)
dUsoskin=calcNormFactors(dUsoskin)
effLogLib=log(dUsoskin$samples$lib.size*dUsoskin$samples$norm.factors)
setwd("~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/zingeR_edgeR/")
w = zeroWeightsLibSizeDispFastUsoskin(counts=usoskin, designZI=model.matrix(~effLogLib+pickingSession), design=model.matrix(~cellTypeAll+pickingSession), maxit=500, plotW=TRUE, normalization="TMM")
#save(w,file="weightsConvergedZingeREdgeRUsoskin.rda")
#convergence after 476 iterations.
load("weightsConvergedZingeREdgeRUsoskin.rda")
dUsoskin$weights=w
design=model.matrix(~cellTypeAll+pickingSession)
dUsoskin=estimateDisp(dUsoskin,design)
plotBCV(dUsoskin)
fitBatch = glmFit(dUsoskin,design)
## testing every celltype against the average of the other 10 celltypes
L <- matrix(0,nrow=ncol(fitBatch$coefficients),ncol=11)
rownames(L) <- colnames(fitBatch$coefficients)
colnames(L) <- c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH")
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
lrtList=list()
fitBatch$df.residual <- rowSums(fitBatch$weights)-ncol(design)
for(i in 1:ncol(L)) lrtList[[i]] <- glmLRTOld(fitBatch,contrast=L[,i],test="F")
#### independent filtering on the p-values
#padjListBatch=lapply(lrtList, function(x) p.adjust(x$table$PValue,"BH"))
baseMean = unname(rowMeans(sweep(dUsoskin$counts,2,dUsoskin$samples$norm.factors,FUN="*")))
library(genefilter)
padjListBatch <- lapply(lrtList, function(x) pvalueAdjustment_kvdb(baseMean=baseMean, pValue=x$table$PValue))
unlist(lapply(padjListBatch,function(x) sum(x$padj<.05,na.rm=TRUE)))




### not including batch in zero model
dUsoskin=DGEList(usoskin)
dUsoskin=calcNormFactors(dUsoskin)
effLogLib=log(dUsoskin$samples$lib.size*dUsoskin$samples$norm.factors)
setwd("~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/zingeR_edgeR_noBatch/")
#pdf("~/Dropbox/phdKoen/singleCell/usoskinweightsAllCellTypesNoBatchInZeroModel_v2.pdf")
weightsUsoskinNoBatch = zeroWeightsLibSizeDispFastUsoskin(counts=usoskin, designZI=model.matrix(~effLogLib), design=model.matrix(~cellTypeAll+pickingSession), maxit=500, plotW=TRUE, normalization="TMM")
#dev.off()
#save(weightsUsoskinNoBatch,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/weightsUsoskinNoBatchLibSizeDispFast300Iter.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/weightsUsoskinNoBatchLibSizeDispFast300Iter.rda")


######### PLOT
## P(zero) ~ library size
png("~/usoskin.png", width=10,height=5, units="in", res=330)
par(mfrow=c(1,3))
expit <- function(x) 1/(1+exp(-x))
## not including batch effect
par(mar=c(5,4.5,4,1))
logLib=log(colSums(usoskin))
pZero=colMeans(usoskin==0)
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m <- glm(pZero~logLib,family="binomial")
a <- coef(m)[1]
b <- coef(m)[2]
lines(x=sort(logLib),y=expit(a+b*sort(logLib)),lwd=2,col="steelblue")
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,pch=1,cex=1.5)

## including main batch effect
plot(x=logLib,y=pZero, xlab="Log library size", ylab="Fraction of zero's", main= "", cex.lab=2, cex.main=2, bty="l", pch=1, cex.axis=1.5, col=as.numeric(pickingSession))
m2 <- glm(pZero~logLib+pickingSession,family="binomial")
a <- coef(m2)[1]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="Cold"]),y=expit(a+b*sort(logLib[pickingSession=="Cold"])),lwd=2,col=1)
a <- coef(m2)[1]+coef(m2)[3]
b <- coef(m2)[2]
lines(x=sort(logLib[pickingSession=="RT-1"]),y=expit(a+b*sort(logLib[pickingSession=="RT-1"])),col=2,lwd=2)
a <- coef(m2)[1]+coef(m2)[4]
lines(x=sort(logLib[pickingSession=="RT-2"]),y=expit(a+b*sort(logLib[pickingSession=="RT-2"])),col=3,lwd=2)
legend("bottomleft",c("Cold","RT-1","RT-2"),col=1:3,lty=1,cex=1.5)

load("~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/zingeR_edgeR/weightsConvergedZingeREdgeRUsoskin.rda")
par(mar=c(5,4.5,4,1))
hist(w[usoskin==0],xlab="Posterior probability",main="",cex.lab=2,cex.axis=2, ylim=c(0,3e6))
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/weightsUsoskinNoBatchLibSizeDispFast300Iter.rda")
hist(weightsUsoskinNoBatch[usoskin==0],add=TRUE,col=rgb(0.1,0.8,0.1,.2))
legend("topleft",c("picking session + effective lib. size","effective lib. size"),fill=c(0,rgb(0.1,0.8,0.1,.2)), bty="n",cex=1.25)
dev.off()

### edgeR analysis
design=model.matrix(~cellTypeAll+pickingSession)
dNoWeights=DGEList(usoskin)
dNoWeights=edgeR::calcNormFactors(dNoWeights)
dNoWeights=estimateDisp(dNoWeights,design)
#save(dNoWeights,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/dUsoskinNoWeights.rda")
fitNoWeights=glmFit(dNoWeights,design)
L <- matrix(0,nrow=ncol(fitNoWeights$coefficients),ncol=11)
rownames(L) <- colnames(fitNoWeights$coefficients)
colnames(L) <- c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH")
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
lrtListNoWeights=list()
for(i in 1:ncol(L)) lrtListNoWeights[[i]] <- glmLRTOld(fitNoWeights,contrast=L[,i],test="F", ZI=FALSE)
padjListNoWeights=lapply(lrtListNoWeights, function(x) p.adjust(x$table$PValue,"BH"))
deGenesNoWeights=unlist(lapply(padjListNoWeights,function(x) sum(x<.05)))

### DESeq2 analysis
library(DESeq2)
colData <- data.frame(cellType=cellTypeAll,pickingSession=pickingSession)
dse <- DESeqDataSetFromMatrix(countData = usoskin, colData = colData, design = ~ cellType+pickingSession)
dse <- DESeqDataSetFromMatrix(countData = usoskin, colData = colData, design = ~ cellType+pickingSession)
dse = estimateSizeFactors(dse)
dse = estimateDispersions(dse)
dse = nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE)
#save(dse,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/dseUsoskin.rda")
resultsNames(dse) #for building contrasts, see ?results
# L2=matrix(0,nrow=length(resultsNames(dse)),ncol=11)
# rownames(L2)=resultsNames(dse)
# L2[2:12,1] = c(1,rep(-1/10,10))
# L2[2:12,2] = c(-1/10,1,rep(-1/10,9))
# L2[2:12,3] = c(rep(-1/10,2),1,rep(-1/10,8))
# L2[2:12,4] = c(rep(-1/10,3),1,rep(-1/10,7))
# L2[2:12,5] = c(rep(-1/10,4),1,rep(-1/10,6))
# L2[2:12,6] = c(rep(-1/10,5),1,rep(-1/10,5))
# L2[2:12,7] = c(rep(-1/10,6),1,rep(-1/10,4))
# L2[2:12,8] = c(rep(-1/10,7),1,rep(-1/10,3))
# L2[2:12,9] = c(rep(-1/10,8),1,rep(-1/10,2))
# L2[2:12,10] = c(rep(-1/10,9),1,rep(-1/10,1))
# L2[2:12,11] = c(rep(-1/10,10),1)
L=matrix(0,nrow=length(resultsNames(dse)),ncol=11)
rownames(L)=resultsNames(dse)
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
resList=list()
for(i in 1:ncol(L)) resList[[i]] = results(dse,contrast=L[,i])
#save(resList,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/resListDESeq2Usoskin.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/resListDESeq2Usoskin.rda")
cbind(deGenesBatch,deGenesNoWeights,deGenesDESeq2)

### DESeq2 poscounts analysis
library(DESeq2)
colData <- data.frame(cellType=cellTypeAll,pickingSession=pickingSession)
dse <- DESeqDataSetFromMatrix(countData = usoskin, colData = colData, design = ~ cellType+pickingSession)
dse = estimateSizeFactors(dse, type="poscounts")
dse = estimateDispersions(dse)
dse = nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE)
#save(dse,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/dsePoscountsUsoskin.rda")
L=matrix(0,nrow=length(resultsNames(dse)),ncol=11)
rownames(L)=resultsNames(dse)
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
resList=list()
for(i in 1:ncol(L)) resList[[i]] = results(dse,contrast=L[,i])
#save(resList,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/resListPoscountsUsoskin.rda")


### MAST analysis on cpm
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/case/esetUsoskinCpm.RData")
library(MAST)
tpm=exprs(esetCpm)
sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(cellType=cellTypeAll, pickingSession=pickingSession))
ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
CD <- cData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
cData(sca) <- CD  
## differential expression 
fit <- zlm.SingleCellAssay(~cngeneson+cellType+pickingSession,sca=sca,method="bayesglm",ebayes=TRUE)
#how many genes have all-zero counts in at least one condition?
mean(apply(fit@coefC[,c(1,3:12)],1,function(row) any(is.na(row)))) 
L <- matrix(0,nrow=ncol(fit@coefC),ncol=11)
rownames(L) <- colnames(fit@coefC)
colnames(L) <- c("NF1","NF2","NF3","NF4","NF5","NP1","NP2","NP3","PEP1","PEP2","TH")
L[3:12,1] <- -1/10 #NF1 vs. others
L[3:12,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[3:12,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[3:12,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[3:12,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[3:12,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[3:12,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[3:12,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[3:12,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[3:12,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[3:12,11] <- c(rep(-1/10,9),1) #TH vs. others
lrList=list()
for(i in 1:ncol(L)) lrList[[i]] = lrTest(fit,hypothesis=matrix(L[,i],ncol=1))
#save(lrList,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/lrListMASTUsoskin.rda")
#hurdle model results
padjListHurdle <- lapply(lrList,function(x) p.adjust(x[,'hurdle','Pr(>Chisq)'],"BH"))
unlist(lapply(padjListHurdle,function(x) sum(x<=0.05))) #nr DE genes at 5%
# count component results
padjListCont <- lapply(lrList,function(x) p.adjust(x[,'cont','Pr(>Chisq)'],"BH"))
unlist(lapply(padjListCont,function(x) sum(x<=0.05))) #nr DE genes at 5%


###### zingeR-DESeq2
library(DESeq2)
dse = DESeqDataSetFromMatrix(usoskin, colData=data.frame(cellType=cellTypeAll, pickingSession=pickingSession), design=~cellType+pickingSession)
dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
effLogLib=colSums(usoskin)*(1/dse$sizeFactor)


setwd("~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/DESeq2Zero_v2/")
weightsUsoskinBatch = zeroWeightsLibSizeDispFastUsoskin(counts=usoskin, designZI=model.matrix(~effLogLib+pickingSession), design=model.matrix(~cellTypeAll+pickingSession), designFormula=~cellType+pickingSession , maxit=500, plotW=TRUE, normalization="DESeq2_poscounts")
#save(weightsUsoskinBatch, file="~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/DESeq2Zero_v2/weightsZingeRDESeq2Usoskin500iter.RData")
## analysis
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/usoskin/DESeq2Zero_v2/weightsZingeRDESeq2Usoskin500iter.RData") 
hist(weightsUsoskinBatch[usoskin==0])
library(DESeq2) ; library(genefilter)
dse = DESeqDataSetFromMatrix(usoskin, colData=data.frame(cellType=cellTypeAll, pickingSession=pickingSession), design=~cellType+pickingSession)
dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
dimnames(weightsUsoskinBatch) = NULL
assays(dse)[["weights"]] = weightsUsoskinBatch
#dse <- DESeq(dse, betaPrior=TRUE, modelMatrixType="standard")
dse = estimateDispersions(dse)
dse = nbinomWaldTest(dse, betaPrior=TRUE, modelMatrixType="standard")
#save(dse,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/dseZingeRDESeq2500IterUsoskin.rda")
L=matrix(0,nrow=length(resultsNames(dse)),ncol=11)
rownames(L)=resultsNames(dse)
L[2:11,1] <- -1/10 #NF1 vs. others
L[2:11,2] <- c(1,rep(-1/10,9)) #NF2 vs. others
L[2:11,3] <- c(-1/10,1,rep(-1/10,8)) #NF3 vs. others
L[2:11,4] <- c(rep(-1/10,2),1,rep(-1/10,7)) #NF4 vs. others
L[2:11,5] <- c(rep(-1/10,3),1,rep(-1/10,6)) #NF5 vs. others
L[2:11,6] <- c(rep(-1/10,4),1,rep(-1/10,5)) #NP1 vs. others
L[2:11,7] <- c(rep(-1/10,5),1,rep(-1/10,4)) #NP2 vs. others
L[2:11,8] <- c(rep(-1/10,6),1,rep(-1/10,3)) #NP3 vs. others
L[2:11,9] <- c(rep(-1/10,7),1,rep(-1/10,2)) #PEP1 vs. others
L[2:11,10] <- c(rep(-1/10,8),1,rep(-1/10,1)) #PEP2 vs. others
L[2:11,11] <- c(rep(-1/10,9),1) #TH vs. others
resList=list()
for(i in 1:ncol(L)) resList[[i]] = results(dse,contrast=L[,i])
#save(resList,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/resListZingeRDESeq2Usoskin.rda")
unlist(lapply(resList, function(x) sum(x$padj<=0.05,na.rm=TRUE)))












###### JUNK
### slow EM-algorithm
pdf("~/Dropbox/PhD/Research/zeroInflation/singleCell/
   ")
weightsUsoskinBatch = zeroWeightsLibSizeDispFast(counts=usoskin, designZI=model.matrix(~effLogLib+pickingSession), design=model.matrix(~cellTypeAll+pickingSession), maxit=300, plotW=TRUE)
dev.off()



### next 5x10 iterations
w=weightsUsoskinBatch10Iterations
#first I did from 1 to 5 (so 1 is iteration 10-20, 2 is 21-30,...)
#then from 6 to 10
#then from 11 to 15
for(j in 11:15){
    pdf(paste0("~/Dropbox/phdKoen/singleCell/usoskinweightsBatchInMainAndZeroModelAllCellTypes",j,".pdf"))
counts=DGEList(usoskin)
designZI=model.matrix(~effLogLib+pickingSession)
design=model.matrix(~cellTypeAll+pickingSession)
maxit=10
plot=TRUE
plotW=TRUE
    if(plot | plotW) par(mfrow=c(1,plot+plotW))    
    zeroId <- counts$counts==0
    #w <- weightsUsoskinBatch10Iterations
    wFinal <- matrix(NA,nrow=nrow(counts),ncol=ncol(counts))
    active <- rowSums(counts$counts==0)>0 #work with genes with at least 1 zero
    #wFinal[!active,]=w[!active,]
    llOld <- matrix(-1e4,nrow=nrow(counts),ncol=ncol(counts))
    likCOld <- matrix(0,nrow=nrow(counts),ncol=ncol(counts))
    maxit=10
    for(i in 1:maxit){
        zeroId <- counts$counts==0	
	counts$weights <- w
	
	### M-step counts
	counts <- estimateGLMCommonDisp(counts, design, interval=c(0,10))
	counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0, min.row.sum=1)
	if(plot) plotBCV(counts)
	fit <- glmFit(counts, design)
	if(i>1) likCOld <- likC[!converged,]	
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)
	
	### M-step mixture parameter: model zero probability
	successes <- colSums(1-w) #P(zero)
	failures <- colSums(w) #1-P(zero)
	if(is.null(designZI)){
	zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")} else{
	zeroFit <- glm(cbind(successes,failures) ~-1+designZI, family="binomial")}
	pi0Hat <- predict(zeroFit,type="response") 
	
	## E-step: Given estimated parameters, calculate expected value of weights
	pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
	wOld <- w
	w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)
	rownames(w) <- rownames(wOld)
	if(plotW) hist(w[zeroId])

	## expected complete data log-likelihood
	if(i>1) llOld <- ll[!converged,]
	ll <- w*log(pi0HatMat) + (1-w)*log(1-pi0HatMat) + (1-w)*log(likC)

	converged=FALSE
	if(any(converged)){
	    wFinal[as.numeric(rownames(w)[converged]),] = w[converged,]
	    w <- matrix(w[!converged,],ncol=ncol(counts), dimnames=list(c(rownames(w)[!converged]),NULL))
	    counts <- counts[!converged,]
	}
	#cat(paste0("mean diff in L: ",round(mean(rowSums(exp(ll)-exp(llOld))),2),". active features: ",nrow(counts),"\n"))
	cat(paste0("iteration ",i))
	#metagenomeSeq for example seems to report mean(ll) instead of difference.
	if(all(converged)) break
	if(i==maxit){
	    wFinal[apply(wFinal,1,function(row) any(is.na(row))),] = w
	    break
	}
    }
    dev.off()
    save(w,file=paste0("/Users/koenvandenberge/PhD_Data/singleCell/usoskinCase/weightsUsoskinBatch",j,".rda"))
}



### ignoring the batch effect
pdf("~/Dropbox/phdKoen/singleCell/usoskinweights100IterationsBatchInMainNotInZeroModelAllCellTypes.pdf")
weightsUsoskinNoBatch50Iterations = zeroWeightsLibSizeFast(counts=usoskin, design=model.matrix(~cellTypeAll+pickingSession), maxit=100, plot=TRUE, plotW=TRUE)
dev.off()



