source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

#########################################################
############# RNA-Seq data simulation ###################
#########################################################
#### no zero inflation simulation
library(Biobase)
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/bottomly_eset.RData")
bottomly=exprs(bottomly.eset)
bottomly=bottomly[!rowSums(bottomly)==0,]
nSamp <- 10
nTags <- 20e3
grp <- as.factor(rep(0:1, each = nSamp/2))
libSize = sample(round(seq(8e6,10e6,length.out=nSamp)))
DEind = sample(1:nTags,floor(nTags/20),replace=FALSE) #5% differentially expressed
fcSim <- (2 + rexp(length(DEind), rate = 1)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(1)
dataNoZI <- NBsim(foldDiff = fcSim, ind=DEind, dataset = bottomly, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = FALSE, lib.size=libSize, drop.low.lambda=TRUE)
selectedMethods <- c("edgeR_robust", "limma_voom","edgeR", "DESeq2", "edgeREMLibSizeFastOldF")
pvalsNoZI <- pval(dataNoZI, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
#rnaSeqPerformanceNoZeroes.pdf
library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataNoZI$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(edgeRRobust=pvalsNoZI$pval$edgeR_robust, limma_voom=pvalsNoZI$pval$limma_voom, edgeR=pvalsNoZI$pval$edgeR,  edgeREMLibSizeFast=pvalsNoZI$pval$edgeREMLibSizeFast, DESeq2=pvalsNoZI$pval$DESeq2, row.names = rownames(dataNoZI)), 
		   padj = data.frame(edgeRRobust=pvalsNoZI$padj$edgeR_robust, limma_voom=pvalsNoZI$padj$limma_voom, edgeR=pvalsNoZI$padj$edgeR, edgeREMLibSizeFast=pvalsNoZI$padj$edgeREMLibSizeFast, DESeq2=pvalsNoZI$padj$DESeq2, row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
cobraplot <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))

# add zeroes: here we see an obvious gain
dataZeroes = dataNoZI
propZeroes=0.05
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
samp=sample(1:length(zeroId),floor(length(zeroId)*propZeroes))
zeroId[samp]=0
dataZeroes$counts = dataZeroes$counts*zeroId
pvalsZeroes = pval(dataZeroes, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
## edgeR with ground truth
pvalEdgeRGroundTruth <- edgeRWeightedOldF.pfun(counts=dataZeroes$counts, group=grp, weights=zeroId)

##performance curves
#rnaSeqPerformanceWithZeroes.pdf
truthZero=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truthZero[dataZeroes$indDE,"status"]=1
cobraZero <- COBRAData(pval =data.frame(edgeRRobust=pvalsZeroes$pval$edgeR_robust, limma_voom=pvalsZeroes$pval$limma_voom, DESeq2=pvalsZeroes$pval$DESeq2, edgeR=pvalsZeroes$pval$edgeR, edgeREMLibSizeFast=pvalsZeroes$pval$edgeREMLibSizeFast, edgeRTruth = pvalEdgeRGroundTruth[,1], row.names = rownames(dataZeroes)), 
		   padj = data.frame(edgeRRobust=pvalsZeroes$padj$edgeR_robust, limma_voom=pvalsZeroes$padj$limma_voom, DESeq2=pvalsZeroes$padj$DESeq2, edgeR=pvalsZeroes$padj$edgeR, edgeREMLibSizeFast=pvalsZeroes$padj$edgeREMLibSizeFast, edgeRTruth = pvalEdgeRGroundTruth[,2], row.names = rownames(dataZeroes)),
                   truth = truthZero)
cobraperf <- calculate_performance(cobraZero, binary_truth = "status")
cobraplotZero <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplotZero)
plot_roc(cobraplotZero,xaxisrange=c(0,0.1), yaxisrange=c(0,0.6))

## plot a histogram of the weights for the introduced zeroes
weights=zeroWeightsLibSizeFast(counts=dataZeroes$counts, design=model.matrix(~grp), maxit=30, plotW=TRUE)
hist(weights[samp[dataNoZI$counts[samp]!=0]], xlab="weight", xlim=c(0,1), main="")



## plot for paper
p1=plot_fdrtprcurve(cobraplot)
p2=plot_fdrtprcurve(cobraplotZero)
library(scater) 
multiplot(p1,p2,cols=2)






























########### OLD CODE

#### no zero inflation simulation
#library(tweeDEseqCountData)
#data(pickrell)
#pickrell <- as.matrix(exprs(pickrell.eset))
nSamp <- 10
nTags <- 20e3
grp <- as.factor(rep(0:1, each = nSamp/2))
libSize = sample(round(seq(5e6,8e6,length.out=nSamp)))
DEind = sample(1:nTags,floor(nTags/20),replace=FALSE) #5% differentially expressed
fcSim <- (2 + rexp(length(DEind), rate = 1)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(1)
dataNoZI <- NBsim(foldDiff = fcSim, ind=DEind, dataset = pickrell, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE)
selectedMethods <- c("edgeR_robust", "limma_voom","edgeR", "DESeq2", "edgeREMLibSizeOldF")
pvalsNoZI <- pval(dataNoZI, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
hlp=edgeREMLibSizeFastOldF.pfun(counts=dataNoZI$counts, group=grp, niter=1e3)
#### performance curves using iCOBRA
#rnaSeqPerformanceNoZeroes.pdf
library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataNoZI$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(edgeRRobust=pvalsNoZI$pval$edgeR_robust, limma_voom=pvalsNoZI$pval$limma_voom, edgeR=pvalsNoZI$pval$edgeR,  edgeREMLibSize=pvalsNoZI$pval$edgeREMLibSize, DESeq2=pvalsNoZI$pval$DESeq2, edgeREMFast=hlp[,1], row.names = rownames(dataNoZI)), 
		   padj = data.frame(edgeRRobust=pvalsNoZI$padj$edgeR_robust, limma_voom=pvalsNoZI$padj$limma_voom, edgeR=pvalsNoZI$padj$edgeR, edgeREMLibSize=pvalsNoZI$padj$edgeREMLibSize, DESeq2=pvalsNoZI$padj$DESeq2,edgeREMFast=hlp[,2], row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
cobraplot <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))


# add zeroes: here we see an obvious gain
dataZeroes = dataNoZI
propZeroes=0.1
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
zeroId[sample(1:length(zeroId),floor(length(zeroId)*propZeroes))]=0
dataZeroes$counts = dataZeroes$counts*zeroId
genesWithAddedZero <- which(rowSums(zeroId)<nSamp)
pvalsZeroes = pval(dataZeroes, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
hlp=edgeREMLibSizeFastOldF.pfun(counts=dataZeroes$counts, group=grp, niter=1e3)
##performance curves
#rnaSeqPerformanceWithZeroes.pdf
truthZero=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truthZero[dataZeroes$indDE,"status"]=1
cobraZero <- COBRAData(pval =data.frame(edgeRRobust=pvalsZeroes$pval$edgeR_robust, limma_voom=pvalsZeroes$pval$limma_voom, edgeR=pvalsZeroes$pval$edgeR, edgeREMLibSize=pvalsZeroes$pval$edgeREMLibSize, DESeq2=pvalsZeroes$pval$DESeq2, edgeREMFast=hlp[,1], row.names = rownames(dataZeroes)), 
		   padj = data.frame(edgeRRobust=pvalsZeroes$padj$edgeR_robust, limma_voom=pvalsZeroes$padj$limma_voom, edgeR=pvalsZeroes$padj$edgeR, edgeREMLibSize=pvalsZeroes$padj$edgeREMLibSize, DESeq2=pvalsZeroes$padj$DESeq2 , edgeREMFast=hlp[,2] , row.names = rownames(dataZeroes)),
                   truth = truthZero)
cobraperf <- calculate_performance(cobraZero, binary_truth = "status")
cobraplotZero <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplotZero)
plot_roc(cobraplotZero,xaxisrange=c(0,0.1), yaxisrange=c(0,0.6))

#performanceRNAseq.pdf
p1=plot_fdrtprcurve(cobraplot)
p2=plot_fdrtprcurve(cobraplotZero)
library(scater) 
multiplot(p1,p2)

