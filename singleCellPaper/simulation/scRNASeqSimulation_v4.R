source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v5.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")
library(DESeq2) ; library(edgeR) ; library(limma) ; library(scde) ; library(MAST) ; library(mgcv) ; library(MultiAssayExperiment) ; library(SummarizedExperiment) ; library(scales) ; library(iCOBRA)

##################################################################
########################### ISLAM ################################
##################################################################
library(GEOquery)
data = read.delim("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/expressionTabAdapted_kvdb.txt")
#seriesMatrix=getGEO(filename="~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/GSE29087_series_matrix.txt")
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/seriesMatrix.rda")
countData = data[,8:ncol(data)]
rownames(countData)=data[,1]
#(pData(seriesMatrix)$title)[93:96] #negative controls added by the authors: can be discarded.
countData=countData[,1:92]
well = factor(substr(colnames(countData),1,1))
fibroID <- grep(x=pData(seriesMatrix[[1]])$title,pattern="fibroblast")
stemCellID <- grep(x=pData(seriesMatrix[[1]])$title,pattern="stem cell")
colnames(countData)[fibroID] = paste("fibro",1:length(fibroID),sep="_")
colnames(countData)[stemCellID] = paste("stemCell",1:length(stemCellID),sep="_")
cellType=vector(length=ncol(countData))
cellType[fibroID] = "fibro"
cellType[stemCellID] <- "stemCell"
islam = as.matrix(countData) 
islam = islam[!rowSums(islam>0)<5,]

## weight distribution on real data
design=model.matrix(~cellType)
weightsIslamIter50=zeroWeightsLibSize(counts=islam, niter=50, design=design, plotW=FALSE)
weightsIslamFast=zeroWeightsLibSizeDispFast(counts=islam, maxit=200, design=design)
par(mfrow=c(1,2))
hist(weightsIslamIter50[islam==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))
hist(weightsIslamFast[islam==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))

#axis(2,at=c(0,2e4,4e4,6e4),labels=c("0","2e4","4e4","6e4"))

## get gene-wise parameters acc to zero-truncated NB distribution
#paramsIslamAllDesignAveLogCPM=getDatasetZTNB(counts=islam, design=model.matrix(~cellType))
#save(paramsIslamAllDesignAveLogCPM,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsIslamAllDesignAveLogCPM.rda") #with lambda as mean(count/libSize)
#load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/paramsIslamAllDesignAveLogCPM.rda")
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsIslamAllDesignAveLogCPM.rda")
#save(paramsIslamAllDesignAveLogCPM,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsIslamAllDesignAveLogCPM_lambda.rda") #with ZTNB lambda
#load("~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsIslamAllDesignAveLogCPM_lambda.rda") 

####### 40 samples / condition
set.seed(12)
nSamp <- 80
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 15e3
DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(11)
libSizes=sample(colSums(islam),nSamp,replace=TRUE)
dataIslamAllAveLogCPM <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = islam, nTags = nTags, group = grp, verbose = TRUE, params=paramsIslamAllDesignAveLogCPM, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM")
means = sweep(dataIslamAllAveLogCPM$Lambda,2,colSums(dataIslamAllAveLogCPM$counts),"*")
simFC = log2(dataIslamAllAveLogCPM$foldDiff)

### compare simulated and empirical data
#pdf("~/Dropbox/phdKoen/singleCell/figures/supplementary/simulatedDataIslamNoNoise.pdf")
#empirical BCV
condition = cellType
design = model.matrix(~condition)
dIslam=DGEList(islam)
dIslam=calcNormFactors(dIslam)
dIslam=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dIslam,design, interval=c(0,10)),design,prior.df=0)
#simulated BCV
dataIslam=dataIslamAllAveLogCPM
design=model.matrix(~grp)
dSim=DGEList(dataIslam$counts)
dSim=calcNormFactors(dSim)
dSim=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dSim,design,interval=c(0,10)),design,prior.df=0)

par(mfrow=c(2,2), mar=c(5,4,1,1))
plotBCV(dIslam,main="Empirical", ylim=c(0,12), ylab="BCV")
plotBCV(dSim, main="Simulated", ylim=c(0,12), ylab="BCV",xlim=c(1,15))
## marginal density
#plot(density(log(islam+1)))
#lines(density(log(dataIslam$counts+1)),col=2)
#plot(density(log(dataIslam$counts+1)))
#legend("topright",c("empirical","simulated"),lty=1,col=c("black","red"))
## zero ~ libSize
plot(x=log(colSums(islam)),y=colMeans(islam==0), xlab="Log library size", ylab="Fraction of zeros", xlim=c(9.5,15.5), ylim=c(0.25,0.95))
points(x=log(colSums(dataIslamAllAveLogCPM$counts)),y=colMeans(dataIslamAllAveLogCPM$counts==0), xlab="Log library size", ylab="Fraction of zeros",col=2)
legend("topright",c("real","simulated"),pch=1,col=1:2, bty="n", cex=.6)
## zero ~ cpm
#par(mar=c(7,4,1,1))
plot(x=aveLogCPM(islam),y=rowMeans(islam==0), xlab="Average log CPM", ylab="Fraction of zeros", xlim=c(0.5,16),col=alpha(1,1/2),pch=19,cex=.2)
points(x=aveLogCPM(dataIslamAllAveLogCPM$counts),y=rowMeans(dataIslamAllAveLogCPM$counts==0), xlab="aCPM", ylab="Fraction of zeros",col=alpha(2,1/2),pch=19,cex=.2)
legend("topright",c("real","simulated"),pch=19,col=1:2, bty="n", cex=.6)
#dev.off()

### performance characteristics
selectedMethods <- c("edgeREMLibSizeDispFastOldFFilteredEdgeR", "MAST", "limma_voomZeroFiltered",  "DESeq2", "DESeq2_poscounts", "DESeq2Zero_adjustedDf_posCountsNormZeroWeights", "DESeq2Zero_wald_posCountsNormZeroWeights", "edgeR_robust", "edgeR", "edgeROldF", "limma_voom", "limma_voomFiltered", "metagenomeSeq","edgeRFiltered", "NODES", "MAST_count")

group=grp
pvalsIslam <- pval(dataIslamAllAveLogCPM, method=selectedMethods, mc.cores=2, niter=200)
#save(pvalsIslam,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamPaper.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamPaper.rda")
for(k in 1:length(pvalsIslam$padj)) pvalsIslam$padj[[k]][is.na(pvalsIslam$padj[[k]])]=1

## scde
#scdePIslam=scde.pfun(dataIslamAllAveLogCPM$counts,group=grp,mc.cores=1) #gives trouble in parallellization
#save(scdePIslam,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslam.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslam.rda")


resDeseq2Wald= DESeq2Zero_wald_posCountsNormZeroWeights.pfun(counts=dataIslamAllAveLogCPM$counts, group=grp, niter=200)
resDeseq2Wald[is.na(resDeseq2Wald)]=1

#FDR-TPR
truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslamAllAveLogCPM))
truthIslam[dataIslamAllAveLogCPM$indDE,"status"]=1
cobraIslam <- COBRAData(pval =data.frame(
					zingeR_edgeR=pvalsIslam$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					 DESeq2=pvalsIslam$pval$DESeq2,
					 DESeq2_poscounts=pvalsIslam$pval$DESeq2_poscounts,
					 zingeR_DESeq2=pvalsIslam$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					 zingeR_DESeq2_wald=resDeseq2Wald[,"pval"],
					 limma_voom=pvalsIslam$pval$limma_voom,
					 zingeR_limma_voom=pvalsIslam$pval$limma_voomZero,
					 edgeR=pvalsIslam$pval$edgeR,
					 #edgeRFiltered=pvalsIslam$pval$edgeRFiltered,
					 MAST=pvalsIslam$pval$MAST,
					 MAST_count=pvalsIslam$pval$MAST_count,
					 metagenomeSeq=pvalsIslam$pval$metagenomeSeq,
					 scde=scdePIslam[,"pval"],
					 NODES=pvalsIslam$pval$NODES,
					 row.names = rownames(dataIslamAllAveLogCPM)), 
		   padj = data.frame(
		   			 zingeR_edgeR=pvalsIslam$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
				     DESeq2=pvalsIslam$padj$DESeq2,
				     DESeq2_poscounts=pvalsIslam$padj$DESeq2_poscounts,
					 zingeR_DESeq2=pvalsIslam$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					 zingeR_DESeq2_wald=resDeseq2Wald[,"padj"],
				     limma_voom=pvalsIslam$padj$limma_voom, 
					 zingeR_limma_voom=pvalsIslam$padj$limma_voomZero,	
				     edgeR=pvalsIslam$padj$edgeR,
					 #edgeRFiltered=pvalsIslam$padj$edgeRFiltered,
					MAST=pvalsIslam$padj$MAST, 
					MAST_count=pvalsIslam$padj$MAST_count,
					metagenomeSeq=pvalsIslam$padj$metagenomeSeq,
					scde=scdePIslam[,"padj"],
					NODES=pvalsIslam$padj$NODES,
				     	row.names = rownames(dataIslamAllAveLogCPM)),
                   truth = truthIslam)
cobraperf <- calculate_performance(cobraIslam, binary_truth = "status")
colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", MAST_count="purple", metagenomeSeq="green", scde="grey", NODES="black", zingeR_DESeq2_wald="steelblue")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslam.rda")
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamLimma.rda")
plot_fdrtprcurve(cobraplot, pointsize=2)


#high FC FDR-TPR
fcSim=(3 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(11)
libSizes=sample(colSums(islam),nSamp,replace=TRUE)
dataIslamAllAveLogCPMHighFC <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = islam, nTags = nTags, group = grp, verbose = TRUE, params=paramsIslamAllDesignAveLogCPM, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM")
pvalsIslamHighFC <- pval(dataIslamAllAveLogCPMHighFC, method=selectedMethods, mc.cores=2, niter=200)
#save(pvalsIslamHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamHighFCPaper.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamHighFCPaper.rda")
for(k in 1:length(pvalsIslamHighFC$padj)) pvalsIslamHighFC$padj[[k]][is.na(pvalsIslamHighFC$padj[[k]])]=1
#scdePHighFC=scde.pfun(dataIslamAllAveLogCPMHighFC$counts,group=grp)
#save(scdePHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslamHighFC.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePIslamHighFC.rda")
truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslamAllAveLogCPMHighFC))
truthIslam[dataIslamAllAveLogCPMHighFC$indDE,"status"]=1
cobraIslam <- COBRAData(pval =data.frame(
					zingeR_edgeR=pvalsIslamHighFC$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
					 DESeq2=pvalsIslamHighFC$pval$DESeq2,
					 DESeq2_poscounts=pvalsIslamHighFC$pval$DESeq2_poscounts,
					 zingeR_DESeq2=pvalsIslamHighFC$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					 limma_voom=pvalsIslamHighFC$pval$limma_voom,
					 zingeR_limma_voom=pvalsIslamHighFC$pval$limma_voomZero,
					 edgeR=pvalsIslamHighFC$pval$edgeR,
					 #edgeRFiltered=pvalsIslamHighFC$pval$edgeRFiltered,
					 #edgeREMLRT=hlpEdgeRZeroLRT[,"pval"],
					 MAST=pvalsIslamHighFC$pval$MAST,
					 MAST_count=pvalsIslamHighFC$pval$MAST_count,
					 metagenomeSeq=pvalsIslamHighFC$pval$metagenomeSeq,
					 scde=scdePHighFC[,"pval"],
					 NODES=pvalsIslamHighFC$pval$NODES,
					 row.names = rownames(dataIslamAllAveLogCPMHighFC)), 
		   padj = data.frame(
		   			zingeR_edgeR=pvalsIslamHighFC$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
				     DESeq2=pvalsIslamHighFC$padj$DESeq2,
				     DESeq2_poscounts=pvalsIslamHighFC$padj$DESeq2_poscounts,
					 zingeR_DESeq2=pvalsIslamHighFC$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
				     limma_voom=pvalsIslamHighFC$padj$limma_voom, 
					 zingeR_limma_voom=pvalsIslamHighFC$padj$limma_voomZero,	
				     edgeR=pvalsIslamHighFC$padj$edgeR,
					 #edgeRFiltered=pvalsIslamHighFC$padj$edgeRFiltered,
					 #edgeREMLRT=hlpEdgeRZeroLRT[,"padj"],
					MAST=pvalsIslamHighFC$padj$MAST,
					MAST_count=pvalsIslamHighFC$padj$MAST_count,
					 metagenomeSeq=pvalsIslamHighFC$padj$metagenomeSeq,
					scde=scdePHighFC[,"padj"],
					NODES=pvalsIslamHighFC$padj$NODES,
				     	row.names = rownames(dataIslamAllAveLogCPMHighFC)),
                   truth = truthIslam)
cobraperf <- calculate_performance(cobraIslam, binary_truth = "status")
colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", MAST_count="purple", metagenomeSeq="green", scde="grey", NODES="black")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamHighFC.rda")
plot_fdrtprcurve(cobraplot, pointsize=2)


### compare DESeq2 variants
selectedMethodsDESeq2 <- c(#"DESeq2",
			 "DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noImputation",
			 "DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noFiltering",
			 "DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noShrink",
		     "DESeq2_noImputation", #no imputation for outliers
		     "DESeq2_noFiltering", #no independent filtering step
		     "DESeq2_noShrink", #no shrinkage to mean zero prior of coefs
		     "DESeq2_poscounts_noShrink",
		     "DESeq2_poscounts_noFiltering",
		     "DESeq2_poscounts_noImputation"
		     )

pvalsIslamDESeq2 <- pval(dataIslamAllAveLogCPM, method=selectedMethodsDESeq2, mc.cores=2, niter=200)
#save(pvalsIslamDESeq2,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamDESeq2.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsIslamDESeq2.rda")
for(k in 1:length(pvalsIslamDESeq2$padj)) pvalsIslamDESeq2$padj[[k]][is.na(pvalsIslamDESeq2$padj[[k]])]=1


# FDR-TPR DESeq2 
library(iCOBRA)
truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslamAllAveLogCPM))
truthIslam[dataIslamAllAveLogCPM$indDE,"status"]=1
cobraDESeq2 <- COBRAData(pval =data.frame(
					DESeq2=pvalsIslam$pval$DESeq2,
					DESeq2_poscounts=pvalsIslam$pval$DESeq2_poscounts,
					DESeq2_poscounts_noShrink=pvalsIslamDESeq2$pval$DESeq2_poscounts_noShrink,
					DESeq2_poscounts_noFiltering=pvalsIslamDESeq2$pval$DESeq2_poscounts_noFiltering,
					DESeq2_poscounts_noImputation=pvalsIslamDESeq2$pval$DESeq2_poscounts_noShrink,
					DESeq2_noImputation=pvalsIslamDESeq2$pval$DESeq2_noImputation,
				    DESeq2_noFiltering=pvalsIslamDESeq2$pval$DESeq2_noFiltering,
					DESeq2_noShrink=pvalsIslamDESeq2$pval$DESeq2_noShrink,
					#edgeR=pvalsIslam$pval$edgeR,
					#zingeR_edgeR=pvalsIslam$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					zingeR_DESeq2=pvalsIslam$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					zingeR_DESeq2_noShrink=pvalsIslamDESeq2$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noShrink,
					#zingeR_DESeq2_noImputation=pvalsIslamDESeq2$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noImputation,
					zingeR_DESeq2_noFiltering=pvalsIslamDESeq2$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noFiltering,
					row.names = rownames(dataIslamAllAveLogCPM)), 
		   padj = data.frame(
		   			DESeq2=pvalsIslam$padj$DESeq2,
		   			DESeq2_poscounts=pvalsIslam$padj$DESeq2_poscounts,
		   			DESeq2_poscounts_noShrink=pvalsIslamDESeq2$padj$DESeq2_poscounts_noShrink,
					DESeq2_poscounts_noFiltering=pvalsIslamDESeq2$padj$DESeq2_poscounts_noFiltering,
					DESeq2_poscounts_noImputation=pvalsIslamDESeq2$padj$DESeq2_poscounts_noShrink,
					DESeq2_noImputation=pvalsIslamDESeq2$padj$DESeq2_noImputation,
				    DESeq2_noFiltering=pvalsIslamDESeq2$padj$DESeq2_noFiltering,
					DESeq2_noShrink=pvalsIslamDESeq2$padj$DESeq2_noShrink,
					#edgeR=pvalsIslam$padj$edgeR,
					#zingeR_edgeR=pvalsIslam$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					zingeR_DESeq2=pvalsIslam$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					zingeR_DESeq2_noShrink=pvalsIslamDESeq2$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noShrink,
					#zingeR_DESeq2_noImputation=pvalsIslamDESeq2$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noImputation,
					zingeR_DESeq2_noFiltering=pvalsIslamDESeq2$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noFiltering,
				   row.names = rownames(dataIslamAllAveLogCPM)),
                   truth = truthIslam)
cobraperfDESeq2 <- calculate_performance(cobraDESeq2, binary_truth = "status")
cobraplotDESeq2 <- prepare_data_for_plot(cobraperfDESeq2)
#save(cobraplotDESeq2,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotDESeq2Islam.rda")
plot_fdrtprcurve(cobraplotDESeq2,pointsize=2)

#### compare Gaussian models versus limma
#selectedMethods <- c("regularHurdleTPM","positiveHurdleTPM","MAST","MAST_count","limma_voomHurdle", "metagenomeSeq")
selectedMethods <- c("regularHurdleTPM","positiveHurdleTPM","positiveHurdleTPMCDR","positiveHurdleLimmaTPM","MAST_count","limma_voomHurdle", "limma_voomHurdleHeteroscedastic", "edgeRHurdle")
pvalGaussian <- pval(dataIslamAllAveLogCPM,method=selectedMethods,mc.cores=3,niter=200)
for(k in 1:length(pvalGaussian$padj)) pvalGaussian$padj[[k]][is.na(pvalGaussian$padj[[k]])]=1
#hlpLimmaZeroHomo=limma_voomZeroFilteredHomo.pfun(counts=dataIslamAllAveLogCPM$counts, group=grp, design=model.matrix(~grp), niter=200)
#hlpLimmaZeroHomo[is.na(hlpLimmaZeroHomo[,"padj"]),"padj"]=1


#Gaussian models FDR-TPR
truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslamAllAveLogCPM))
truthIslam[dataIslamAllAveLogCPM$indDE,"status"]=1
cobraIslam <- COBRAData(pval =data.frame(
					regularHurdle_CPM=pvalGaussian$pval$regularHurdleTPM,
					 #positiveHurdleTPM=pvalGaussian$pval$positiveHurdleTPM,
					 #positiveHurdleTPMCDR=pvalGaussian$pval$positiveHurdleTPMCDR,
					 positiveHurdle_CPM=pvalGaussian$pval$positiveHurdleLimmaTPM,
					 MAST=pvalsIslam$pval$MAST,
					 MAST_count=pvalGaussian$pval$MAST_count,
					 #limma_voomPositiveHurdle=pvalGaussian$pval$limma_voomHurdle,
					 limma_voom_hurdle=pvalGaussian$pval$limma_voomHurdleHeteroscedastic,
					 #limma_voomZeroHomo=hlpLimmaZeroHomo[,"pval"],
					 metagenomeSeq=pvalsIslam$pval$metagenomeSeq,
					 zingeR=pvalsIslam$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					 edgeRHurdle=pvalGaussian$pval$edgeRHurdle,
					 row.names = rownames(dataIslamAllAveLogCPM)), 
		   padj = data.frame(
		   			regularHurdle_CPM=pvalGaussian$padj$regularHurdleTPM,
					 #positiveHurdleTPM=pvalGaussian$padj$positiveHurdleTPM,
					 #positiveHurdleTPMCDR=pvalGaussian$padj$positiveHurdleTPMCDR,
					 positiveHurdle_CPM=pvalGaussian$padj$positiveHurdleLimmaTPM,					 
					 MAST=pvalsIslam$padj$MAST,
					 MAST_count=pvalGaussian$padj$MAST_count,
					 #limma_voomPositiveHurdle=pvalGaussian$padj$limma_voomHurdle,
					 limma_voom_hurdle=pvalGaussian$padj$limma_voomHurdleHeteroscedastic,
					 #limma_voomZeroHomo=hlpLimmaZeroHomo[,"padj"],
					 metagenomeSeq=pvalsIslam$padj$metagenomeSeq,
					 zingeR=pvalsIslam$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					 edgeRHurdle=pvalGaussian$padj$edgeRHurdle,
				     row.names = rownames(dataIslamAllAveLogCPM)),
                   truth = truthIslam)
cobraperf <- calculate_performance(cobraIslam, binary_truth = "status")
cobraplotGaussian <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplotGaussian)



##################################################################
########################### TRAPNELL #############################
##################################################################
trapnellAssay72 <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE52529-GPL11154.rds")
trapnellAssay <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE52529-GPL16791.rds")
# get the 48h timepoint from large assay set
trapnellAssay48 <- trapnellAssay[,pData(trapnellAssay)[,"characteristics_ch1.1"] == "hour post serum-switch: 48"]
rm(trapnellAssay)
#remove wells containing debris
trapnellAssay72 <- trapnellAssay72[,!pData(trapnellAssay72)[,"characteristics_ch1.2"]=="debris: TRUE"]
trapnellAssay48 <- trapnellAssay48[,!pData(trapnellAssay48)[,"characteristics_ch1.2"]=="debris: TRUE"]
#remove wells that did not contain one cell
trapnellAssay72 <- trapnellAssay72[,!pData(trapnellAssay72)[,"characteristics_ch1.4"]!="cells in well: 1"]
trapnellAssay48 <- trapnellAssay48[,!pData(trapnellAssay48)[,"characteristics_ch1.4"]!="cells in well: 1"]

countsTrapnell48 <- round(assay(experiments(trapnellAssay48)$gene,"count"))
countsTrapnell72 <- round(assay(experiments(trapnellAssay72)$gene,"count"))
countsTrapnell <- cbind(countsTrapnell48,countsTrapnell72)
countsTrapnell <- countsTrapnell[rowSums(countsTrapnell>0)>9,] #expression in at least 10 out of 149 samples. Remains 24,576 genes and 149 samples.
timePoint=factor(c(rep(48,85),rep(72,64)))

## weight distribution on real data
design=model.matrix(~timePoint)
weightsTrapnellIter50=zeroWeightsLibSize(counts=countsTrapnell, niter=50, design=design, plotW=FALSE)
weightsTrapnellFast=zeroWeightsLibSizeDispFast(counts=countsTrapnell, maxit=200, design=design)
par(mfrow=c(1,2))
hist(weightsTrapnellIter50[countsTrapnell==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))
hist(weightsTrapnellFast[countsTrapnell==0],xlab="Posterior probability",main="",yaxt="n",xlim=c(0,1), breaks=seq(0,1,by=0.05))

#paramsTrapnellAllDesignAveLogCPM=getDatasetZTNB(counts=countsTrapnell, design=model.matrix(~timePoint))
#save(paramsTrapnellAllDesignAveLogCPM,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsTrapnellAllDesignAveLogCPM.rda")
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/paramsTrapnellAllDesignAveLogCPM.rda")


## simulate for 80 samples
set.seed(12)
nSamp <- 160
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 15e3
DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
fcSim=(2 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(11)
libSizes=sample(colSums(countsTrapnell),nSamp,replace=TRUE)
dataTrapnell <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = countsTrapnell, nTags = nTags, group = grp, verbose = TRUE, params=paramsTrapnellAllDesignAveLogCPM, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM")

plot(x=log(colSums(countsTrapnell)),y=colMeans(countsTrapnell==0), xlab="Log library size", ylab="Fraction of zeros")
points(x=log(colSums(dataTrapnell$counts)),y=colMeans(dataTrapnell$counts==0), xlab="Log library size", ylab="Fraction of zero's",col=2)
## zero ~ cpm
#par(mar=c(7,4,1,1))
plot(x=aveLogCPM(countsTrapnell),y=rowMeans(countsTrapnell==0), xlab="Average log CPM", ylab="Fraction of zeros", xlim=c(0.5,16),col=alpha(1,1/2),pch=19,cex=.2)
points(x=aveLogCPM(dataTrapnell$counts),y=rowMeans(dataTrapnell$counts==0), xlab="Average Log CPM", ylab="Fraction of zeroes",col=alpha(2,.8),pch=19,cex=.2)

selectedMethods <- c("edgeREMLibSizeDispFastOldFFilteredEdgeR", "MAST", "limma_voomZeroFiltered",  "DESeq2", "DESeq2_poscounts", "DESeq2Zero_adjustedDf_posCountsNormZeroWeights", "DESeq2Zero_wald_posCountsNormZeroWeights", "edgeR_robust", "edgeR", "edgeROldF", "limma_voom", "limma_voomFiltered", "metagenomeSeq","edgeRFiltered", "NODES", "MAST_count")
group=grp
pvalsTrapnell = pval(dataTrapnell,method=selectedMethods,niter=200, mc.cores=1)
#save(pvalsTrapnell,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellPaper.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellPaper.rda")
for(k in 1:length(pvalsTrapnell$padj)) pvalsTrapnell$padj[[k]][is.na(pvalsTrapnell$padj[[k]])]=1
#scdePTrapnell=scde.pfun(dataTrapnell$counts,group=grp,mc.cores=1) #gives trouble in parallellization
#save(scdePTrapnell,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePTrapnell.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePTrapnell.rda")

resDeseq2Wald= DESeq2Zero_wald_posCountsNormZeroWeights.pfun(counts=dataTrapnell$counts, group=grp, niter=200)
resDeseq2Wald[is.na(resDeseq2Wald)]=1

#FDR-TPR
library(iCOBRA)
truthTrapnell=data.frame(status=rep(0,nTags), row.names=rownames(dataTrapnell))
truthTrapnell[dataTrapnell$indDE,"status"]=1
cobraTrapnell <- COBRAData(pval =data.frame(#edgeRRobust=pvalsTrapnell$pval$edgeR_robust, 
					    limma_voom=pvalsTrapnell$pval$limma_voom,
					    zingeR_limma_voom=pvalsTrapnell$pval$limma_voomZero, 
					    edgeR=pvalsTrapnell$pval$edgeR, 
					    #edgeRFiltered=pvalsTrapnell$pval$edgeRFiltered,
					    zingeR_edgeR=pvalsTrapnell$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
					    DESeq2=pvalsTrapnell$pval$DESeq2,
					    DESeq2_poscounts=pvalsTrapnell$pval$DESeq2_poscounts,
					    zingeR_DESeq2=pvalsTrapnell$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					    zingeR_DESeq2_wald=resDeseq2Wald[,"pval"],
					    MAST=pvalsTrapnell$pval$MAST,
					    MAST_count=pvalsTrapnell$pval$MAST_count,
					    metagenomeSeq=pvalsTrapnell$pval$metagenomeSeq,
					    scde=scdePTrapnell[,"pval"],
					    NODES=pvalsTrapnell$pval$NODES,
					    row.names = rownames(dataTrapnell)), 
		   padj = data.frame(#edgeRRobust=pvalsTrapnell$padj$edgeR_robust, 
				     limma_voom=pvalsTrapnell$padj$limma_voom,
				     zingeR_limma_voom=pvalsTrapnell$padj$limma_voomZero,
				     edgeR=pvalsTrapnell$padj$edgeR,
				     #edgeRFiltered=pvalsTrapnell$padj$edgeRFiltered,
				     zingeR_edgeR=pvalsTrapnell$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
				     DESeq2=pvalsTrapnell$padj$DESeq2,
				     DESeq2_poscounts=pvalsTrapnell$padj$DESeq2_poscounts, 
					 zingeR_DESeq2=pvalsTrapnell$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					 zingeR_DESeq2_wald=resDeseq2Wald[,"padj"],
				     MAST=pvalsTrapnell$padj$MAST,
				     MAST_count=pvalsTrapnell$padj$MAST_count,
				     metagenomeSeq=pvalsTrapnell$padj$metagenomeSeq,
				     scde=scdePTrapnell[,"padj"],
				     NODES=pvalsTrapnell$padj$NODES,
				     row.names = rownames(dataTrapnell)),
                   truth = truthTrapnell)
cobraperf <- calculate_performance(cobraTrapnell, binary_truth = "status")
colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", MAST_count="purple", metagenomeSeq="green", scde="grey", NODES="black", zingeR_DESeq2_wald="steelblue")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnell.rda")
#save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellLimma.rda")
plot_fdrtprcurve(cobraplot, pointsize=2)

#high FC FDR-TPR
set.seed(12)
nSamp <- 160
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 15e3
DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
fcSim=(3 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(11)
libSizes=sample(colSums(countsTrapnell),nSamp,replace=TRUE)
dataTrapnellHighFC <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = countsTrapnell, nTags = nTags, group = grp, verbose = TRUE, params=paramsTrapnellAllDesignAveLogCPM, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM")
pvalsTrapnellHighFC = pval(dataTrapnellHighFC,method=selectedMethods,niter=200, mc.cores=2)
#save(pvalsTrapnellHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellHighFCPaper.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellHighFCPaper.rda")
for(k in 1:length(pvalsTrapnellHighFC$padj)) pvalsTrapnellHighFC$padj[[k]][is.na(pvalsTrapnellHighFC$padj[[k]])]=1
#scdePHighFC=scde.pfun(dataTrapnellHighFC$counts,group=grp)
#save(scdePHighFC,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePHighFCTrapnell.rda")
load("/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/scdePHighFCTrapnell.rda")

truthTrapnell=data.frame(status=rep(0,nTags), row.names=rownames(dataTrapnellHighFC))
truthTrapnell[dataTrapnellHighFC$indDE,"status"]=1
cobraTrapnell <- COBRAData(pval =data.frame(#edgeRRobust=pvalsTrapnell$pval$edgeR_robust, 
					    limma_voom=pvalsTrapnellHighFC$pval$limma_voom,
					    zingeR_limma_voom=pvalsTrapnellHighFC$pval$limma_voomZero, 
					    edgeR=pvalsTrapnellHighFC$pval$edgeR, 
					    #edgeRFiltered=pvalsTrapnellHighFC$pval$edgeRFiltered,
					    zingeR_edgeR=pvalsTrapnellHighFC$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
					    DESeq2=pvalsTrapnellHighFC$pval$DESeq2,
					    DESeq2_poscounts=pvalsTrapnellHighFC$pval$DESeq2_poscounts,
					    zingeR_DESeq2=pvalsTrapnellHighFC$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights, 
					    MAST=pvalsTrapnellHighFC$pval$MAST,
					    MAST_count=pvalsTrapnellHighFC$pval$MAST_count,
					    metagenomeSeq=pvalsTrapnellHighFC$pval$metagenomeSeq,
					    scde=scdePHighFC[,"pval"],
					    NODES=pvalsTrapnellHighFC$pval$NODES,
					    row.names = rownames(dataTrapnellHighFC)), 
		   padj = data.frame(#edgeRRobust=pvalsTrapnell$padj$edgeR_robust, 
				     limma_voom=pvalsTrapnellHighFC$padj$limma_voom,
				     zingeR_limma_voom=pvalsTrapnellHighFC$padj$limma_voomZero,
				     edgeR=pvalsTrapnellHighFC$padj$edgeR,
				     #edgeRFiltered=pvalsTrapnellHighFC$padj$edgeRFiltered,
				     zingeR_edgeR=pvalsTrapnellHighFC$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
				     DESeq2=pvalsTrapnellHighFC$padj$DESeq2, 
				     DESeq2_poscounts=pvalsTrapnellHighFC$padj$DESeq2_poscounts,
				     zingeR_DESeq2=pvalsTrapnellHighFC$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights, 
				     MAST=pvalsTrapnellHighFC$padj$MAST,
				     MAST_count=pvalsTrapnellHighFC$padj$MAST_count,
				     metagenomeSeq=pvalsTrapnellHighFC$padj$metagenomeSeq,
				     scde=scdePHighFC[,"padj"],
				     NODES=pvalsTrapnellHighFC$padj$NODES,
				     row.names = rownames(dataTrapnellHighFC)),
                   truth = truthTrapnell)
cobraperf <- calculate_performance(cobraTrapnell, binary_truth = "status")
colors=c(limma_voom="blue", zingeR_limma_voom="steelblue", edgeR="red", zingeR_edgeR="salmon", edgeRFiltered="pink", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Zero_phyloNorm="forestgreen", MAST="darkturquoise", MAST_count="purple", metagenomeSeq="green", scde="grey", NODES="black")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
save(cobraplot,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellHighFC.rda")
plot_fdrtprcurve(cobraplot, pointsize=2)


#### empirical vs simulated
#pdf("~/Dropbox/phdKoen/singleCell/figures/supplementary/simulatedDataTrapnellNoNoise.pdf")
#empirical BCV
design=model.matrix(~timePoint)
dEmp=DGEList(countsTrapnell)
dEmp=calcNormFactors(dEmp)
dEmp=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dEmp,design,interval=c(0,10)),design,prior.df=0)
#plotBCVNoLeg(dEmp,col.tagwise=ifelse(rowSums(countsTrapnell>0)>16,1,rowSums(countsTrapnell>0))) #colored
#simulated BCV
design=model.matrix(~grp)
dSim=DGEList(dataTrapnell$counts)
dSim=calcNormFactors(dSim)
dSim=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dSim,design,interval=c(0,10)),design,prior.df=0)

par(mfrow=c(2,2), mar=c(5,4,1,1))
plotBCV(dEmp,main="Empirical", ylab="BCV", ylim=c(0,10))
plotBCV(dSim, main="Simulated", ylab="BCV", ylim=c(0,10))
## density
#plot(density(log(countsTrapnell+1)))
#lines(density(log(dataTrapnell$counts+1)),col=2)
#legend("topright",c("empirical","simulated"),lty=1,col=c("black","red"))
## zero ~ libSize
plot(x=log(colSums(countsTrapnell)),y=colMeans(countsTrapnell==0), xlab="Log library size", ylab="Fraction of zero's")
points(x=log(colSums(dataTrapnell$counts)),y=colMeans(dataTrapnell$counts==0), xlab="Log library size", ylab="Fraction of zero's",col=2)
legend("topright",c("real","simulated"),pch=1,col=1:2, bty="n", cex=.6)
## zero ~ cpm
plot(x=aveLogCPM(countsTrapnell),y=rowMeans(countsTrapnell==0), xlab="average log CPM", ylab="Fraction of zero's",col=alpha(1,1/2),pch=19,cex=.2)
points(x=aveLogCPM(dataTrapnell$counts),y=rowMeans(dataTrapnell$counts==0), col=alpha(2,1/2),pch=19,cex=.2)
legend("topright",c("real","simulated"),pch=19,col=1:2, bty="n", cex=.6)
#dev.off()


### DESeq2 comparisons
selectedMethodsDESeq2 <- c(#"DESeq2",
			 "DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noImputation",
			 "DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noFiltering",
			 "DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noShrink",
		     "DESeq2_noImputation", #no imputation for outliers
		     "DESeq2_noFiltering", #no independent filtering step
		     "DESeq2_noShrink", #no shrinkage to mean zero prior of coefs
		     "DESeq2_poscounts_noShrink",
		     "DESeq2_poscounts_noFiltering",
		     "DESeq2_poscounts_noImputation"
		     )

pvalsTrapnellDESeq2 <- pval(dataTrapnell, method=selectedMethodsDESeq2, mc.cores=2, niter=200)
#save(pvalsTrapnellDESeq2,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsTrapnellDESeq2.rda")
for(k in 1:length(pvalsTrapnellDESeq2$padj)) pvalsTrapnellDESeq2$padj[[k]][is.na(pvalsTrapnellDESeq2$padj[[k]])]=1


# FDR-TPR DESeq2 
library(iCOBRA)
truthTrapnell=data.frame(status=rep(0,nTags), row.names=rownames(dataTrapnell))
truthTrapnell[dataTrapnell$indDE,"status"]=1
cobraDESeq2 <- COBRAData(pval =data.frame(
					DESeq2=pvalsTrapnell$pval$DESeq2,
					DESeq2_noImputation=pvalsTrapnellDESeq2$pval$DESeq2_noImputation,
				    DESeq2_noFiltering=pvalsTrapnellDESeq2$pval$DESeq2_noFiltering,
					DESeq2_noShrink=pvalsTrapnellDESeq2$pval$DESeq2_noShrink,
					DESeq2_poscounts_noShrink=pvalsTrapnellDESeq2$pval$DESeq2_poscounts_noShrink,
					DESeq2_poscounts_noFiltering=pvalsTrapnellDESeq2$pval$DESeq2_poscounts_noFiltering,
					DESeq2_poscounts_noImputation=pvalsTrapnellDESeq2$pval$DESeq2_poscounts_noImputation,
				    #edgeR=pvalsTrapnell$pval$edgeR,
					#zingeR=pvalsTrapnell$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					zingeR_DESeq2=pvalsTrapnell$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					zingeR_DESeq2_noImputation=pvalsTrapnellDESeq2$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noImputation,
					zingeR_DESeq2_noShrink=pvalsTrapnellDESeq2$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noShrink,
					zingeR_DESeq2_noFiltering=pvalsTrapnellDESeq2$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noFiltering,
					 row.names = rownames(dataTrapnell)), 
		   padj = data.frame(DESeq2=pvalsTrapnell$padj$DESeq2,
					DESeq2_noImputation=pvalsTrapnellDESeq2$padj$DESeq2_noImputation,
				    DESeq2_noFiltering=pvalsTrapnellDESeq2$padj$DESeq2_noFiltering,
					DESeq2_noShrink=pvalsTrapnellDESeq2$padj$DESeq2_noShrink,
					DESeq2_poscounts_noShrink=pvalsTrapnellDESeq2$padj$DESeq2_poscounts_noShrink,
					DESeq2_poscounts_noFiltering=pvalsTrapnellDESeq2$padj$DESeq2_poscounts_noFiltering,
					DESeq2_poscounts_noImputation=pvalsTrapnellDESeq2$padj$DESeq2_poscounts_noImputation,
				    #edgeR=pvalsTrapnell$padj$edgeR,
					#zingeR=pvalsTrapnell$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					zingeR_DESeq2=pvalsTrapnell$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					zingeR_DESeq2_noImputation=pvalsTrapnellDESeq2$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noImputation,
					zingeR_DESeq2_noShrink=pvalsTrapnellDESeq2$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noShrink,
					zingeR_DESeq2_noFiltering=pvalsTrapnellDESeq2$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noFiltering,
				     	row.names = rownames(dataTrapnell)),
                   truth = truthTrapnell)
cobraperfDESeq2 <- calculate_performance(cobraDESeq2, binary_truth = "status")
cobraplotDESeq2 <- prepare_data_for_plot(cobraperfDESeq2)
#save(cobraplotDESeq2,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotDESeq2Trapnell.rda")
plot_fdrtprcurve(cobraplotDESeq2,pointsize=2)



#### compare Gaussian models versus limma
selectedMethods <- c("regularHurdleTPM","positiveHurdleTPM","positiveHurdleTPMCDR","positiveHurdleLimmaTPM","MAST_count","limma_voomHurdle", "limma_voomHurdleHeteroscedastic", "edgeRHurdle")
pvalGaussian <- pval(dataTrapnell,method=selectedMethods,mc.cores=1,niter=1)
#save(pvalGaussian,file="/Users/koenvandenberge/Dropbox/PhD/Research/zeroInflation/singleCell/pvalGaussianTrapnell.rda")

#Gaussian models FDR-TPR
truthTrapnell=data.frame(status=rep(0,nTags), row.names=rownames(dataTrapnell))
truthTrapnell[dataTrapnell$indDE,"status"]=1
cobraTrapnellGaussian <- COBRAData(pval =data.frame(
					regularHurdle_CPM=pvalGaussian$pval$regularHurdleTPM,
					 positiveHurdle_CPM=pvalGaussian$pval$positiveHurdleLimmaTPM,
					 MAST=pvalsTrapnell$pval$MAST,
					 MAST_count=pvalGaussian$pval$MAST_count,
					 limma_voom_hurdle=pvalGaussian$pval$limma_voomHurdleHeteroscedastic,
					 metagenomeSeq=pvalsTrapnell$pval$metagenomeSeq,
					zingeR=pvalsTrapnell$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					edgeRHurdle=pvalGaussian$pval$edgeRHurdle,
					 row.names = rownames(dataTrapnell)), 
		   padj = data.frame(
		   			regularHurdle_CPM=pvalGaussian$padj$regularHurdleTPM,
					 positiveHurdle_CPM=pvalGaussian$padj$positiveHurdleLimmaTPM,
					 MAST=pvalsTrapnell$padj$MAST,
					 MAST_count=pvalGaussian$padj$MAST_count,
					 limma_voom_hurdle=pvalGaussian$padj$limma_voomHurdleHeteroscedastic,
					 metagenomeSeq=pvalsTrapnell$padj$metagenomeSeq,
					zingeR=pvalsTrapnell$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					edgeRHurdle=pvalGaussian$padj$edgeRHurdle,
				     row.names = rownames(dataTrapnell)),
                   truth = truthTrapnell)
cobraperf <- calculate_performance(cobraTrapnellGaussian, binary_truth = "status")
cobraplotGaussian <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplotGaussian)



####### combine two-panel plots with same legend
## main plot
library(cowplot)
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslam.rda")
islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnell.rda")
trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
#plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
prow <- plot_grid( islamPlot + theme(legend.position="none"),
           trapnellPlot + theme(legend.position="none"),
           align = 'vh',
           labels = c("A", "B"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p

## including ZI limma-voom plot
library(cowplot)
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamLimma.rda")
islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellLimma.rda")
trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
#plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
prow <- plot_grid( islamPlot + theme(legend.position="none"),
           trapnellPlot + theme(legend.position="none"),
           align = 'vh',
           labels = c("A", "B"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p

## high fold changes
## including ZI limma-voom plot
library(cowplot)
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotIslamHighFC.rda")
islamPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotTrapnellHighFC.rda")
trapnellPlot=plot_fdrtprcurve(cobraplot, pointsize=2)
#plot_grid(islamPlot,trapnellPlot, labels = c("A", "B"))
prow <- plot_grid( islamPlot + theme(legend.position="none"),
           trapnellPlot + theme(legend.position="none"),
           align = 'vh',
           labels = c("A", "B"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p

## DESeq2 variants
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotDESeq2Islam.rda")
islamPlot=plot_fdrtprcurve(cobraplotDESeq2, pointsize=2)
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/cobraplotDESeq2Trapnell.rda")
trapnellPlot=plot_fdrtprcurve(cobraplotDESeq2, pointsize=2)
prow <- plot_grid( islamPlot + theme(legend.position="none"),
           trapnellPlot + theme(legend.position="none"),
           align = 'vh',
           labels = c("A", "B"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(islamPlot + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p










### OLD
## weight distribution on real data
timePoint=factor(c(rep(48,85),rep(72,64)))
design=model.matrix(~timePoint)
weightsTrapnell4872Iter50=zeroWeightsLibSizeFast(counts=countsTrapnell, maxit=50, design=design, plotW=TRUE)
hist(weightsTrapnell4872Iter50[countsTrapnell==0],xlab="weight",main="",yaxt="n")
axis(2,at=c(1e5,5e5,1e6,1.5e6),labels=c("1e5","5e5","1e6","1.5e6"))

#paramsTrapnell4872lAllAveLogCPM=getDatasetZTNB(countsTrapnell,drop.low.lambda=FALSE,drop.extreme.dispersion=FALSE, cpm="AveLogCPM")
#save(paramsTrapnell4872lAllAveLogCPM,file="~/PhD_Data/singleCell/evaluateSimulations/newSimulation/paramsTrapnell4872AllAveLogCPM.rda")
#load("~/PhD_Data/singleCell/evaluateSimulations/newSimulation/paramsTrapnell4872AllAveLogCPM.rda")

#load("/Users/koenvandenberge/PhD_Data/singleCell/evaluateSimulations/newSimulation/40samplesGroup/dataTrapnellAllAveLogCPM/paramsTrapnellAllAveLogCPM.rda")


#set.seed(10)
#nSamp <- 150
#grp <- as.factor(rep(0:1, each = nSamp/2))
#nTags = 15e3
#DEind = sample(1:nTags,floor(nTags*.1),replace=FALSE) #5% differentially expressed
#fcSim=(2 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
#set.seed(12)
#libSizes=sample(colSums(countsTrapnell),nSamp,replace=TRUE)
#dataTrapnell <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = countsTrapnell, nTags = nTags, group = grp, verbose = TRUE, params=paramsTrapnellAllAveLogCPM, lib.size=libSizes, randomZero=0, noiseCell=0, noiseGene=0, cpm="AveLogCPM")


#### empirical vs simulated
#pdf("~/Dropbox/phdKoen/singleCell/figures/supplementary/simulatedDataTrapnellNoNoise.pdf")
par(mfrow=c(2,2), mar=c(5,4,1,1))
#empirical BCV
#design=matrix(rep(1,ncol(countsTrapnell),ncol=1)) #only estimate intercept
design=model.matrix(~timePoint)
dEmp=DGEList(countsTrapnell)
dEmp=calcNormFactors(dEmp)
dEmp=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dEmp,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dEmp,main="Empirical", ylab="BCV", ylim=c(0,10))
#plotBCVNoLeg(dEmp,col.tagwise=ifelse(rowSums(countsTrapnell>0)>16,1,rowSums(countsTrapnell>0))) #colored
#simulated BCV
design=model.matrix(~grp)
dSim=DGEList(dataTrapnell$counts)
dSim=calcNormFactors(dSim)
dSim=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dSim,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dSim, main="Simulated", ylab="BCV", ylim=c(0,10))
## density
#plot(density(log(countsTrapnell+1)))
#lines(density(log(dataTrapnell$counts+1)),col=2)
#legend("topright",c("empirical","simulated"),lty=1,col=c("black","red"))
## zero ~ libSize
plot(x=log(colSums(countsTrapnell)),y=colMeans(countsTrapnell==0), xlab="Log library size", ylab="Fraction of zeroes")
points(x=log(colSums(dataTrapnell$counts)),y=colMeans(dataTrapnell$counts==0), xlab="Log library size", ylab="Fraction of zeroes",col=2)
## zero ~ cpm
plot(x=aveLogCPM(countsTrapnell),y=rowMeans(countsTrapnell==0), xlab="average log CPM", ylab="Fraction of zeroes")
points(x=aveLogCPM(dataTrapnell$counts),y=rowMeans(dataTrapnell$counts==0), col=alpha(2,.25))
#dev.off()


## FDR-TPR
library(iCOBRA)
selectedMethods <- c("edgeREMLibSizeDispFastOldFFilteredEdgeR", "MAST", "limma_voomZeroFiltered",  "edgeR", "DESeq2", "limma_voom", "limma_voomFiltered", "metagenomeSeq", "edgeR_robust", "edgeRFiltered")
pvalsTrapnell <- pval(dataTrapnell, method=selectedMethods, mc.cores=2, niter=100)
#scdeP=scde.pfun(dataIslam$counts,group=grp,mc.cores=1) #gives trouble in parallellization

#FDR-TPR
library(iCOBRA)
truthTrapnell=data.frame(status=rep(0,nTags), row.names=rownames(dataTrapnell))
truthTrapnell[dataTrapnell$indDE,"status"]=1
cobraTrapnell <- COBRAData(pval =data.frame(#edgeRRobust=pvalsTrapnell$pval$edgeR_robust, 
					    limma_voom=pvalsTrapnell$pval$limma_voom, 
					    edgeR=pvalsTrapnell$pval$edgeR, 
					    edgeRFiltered=pvalsTrapnell$pval$edgeRFiltered,
					    edgeREMLibSize=pvalsTrapnell$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
					    DESeq2=pvalsTrapnell$pval$DESeq2, 
					    MAST=pvalsTrapnell$pval$MAST,
					    metagenomeSeq=pvalsTrapnell$pval$metagenomeSeq,
					    row.names = rownames(dataTrapnell)), 
		   padj = data.frame(#edgeRRobust=pvalsTrapnell$padj$edgeR_robust, 
				     limma_voom=pvalsTrapnell$padj$limma_voom, 
				     edgeR=pvalsTrapnell$padj$edgeR,
				     edgeRFiltered=pvalsTrapnell$padj$edgeRFiltered,
				     edgeREMLibSize=pvalsTrapnell$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR, 
				     DESeq2=pvalsTrapnell$padj$DESeq2, 
				     MAST=pvalsTrapnell$padj$MAST,
				     metagenomeSeq=pvalsTrapnell$padj$metagenomeSeq,				     
				     row.names = rownames(dataTrapnell)),
                   truth = truthTrapnell)
cobraperf <- calculate_performance(cobraTrapnell, binary_truth = "status")
cobraplot <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))



