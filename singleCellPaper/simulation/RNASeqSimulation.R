source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

#########################################################
############# RNA-Seq data simulation ###################
#########################################################
#### no zero inflation simulation
library(tweeDEseqCountData)
data(pickrell)
pickrell <- as.matrix(exprs(pickrell.eset))
nSamp <- 10
nTags <- 20e3
grp <- as.factor(rep(0:1, each = nSamp/2))
libSize = sample(round(seq(5e6,8e6,length.out=nSamp)))
DEind = sample(1:nTags,floor(nTags/20),replace=FALSE) #5% differentially expressed
fcSim <- (2 + rexp(length(DEind), rate = 1)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(1)
dataNoZI <- NBsim(foldDiff = fcSim, ind=DEind, dataset = pickrell, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE)
selectedMethods <- c("edgeR_robust", "limma_voom","edgeR", "edgeREMLibSizeOldF", "DESeq2")
pvalsNoZI <- pval(dataNoZI, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
#### performance curves using iCOBRA
library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataNoZI$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(edgeRRobust=pvalsNoZI$pval$edgeR_robust, limma_voom=pvalsNoZI$pval$limma_voom, edgeR=pvalsNoZI$pval$edgeR,  edgeREMLibSize=pvalsNoZI$pval$edgeREMLibSize, DESeq2=pvalsNoZI$pval$DESeq2, row.names = rownames(dataNoZI)), 
		   padj = data.frame(edgeRRobust=pvalsNoZI$padj$edgeR_robust, limma_voom=pvalsNoZI$padj$limma_voom, edgeR=pvalsNoZI$padj$edgeR, edgeREMLibSize=pvalsNoZI$padj$edgeREMLibSize, DESeq2=pvalsNoZI$padj$DESeq2, row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
cobraplot <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))


# add zeroes: here we see an obvious gain
dataZeroes = dataNoZI
propZeroes=0.05
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
zeroId[sample(1:length(zeroId),floor(length(zeroId)*propZeroes))]=0
dataZeroes$counts = dataZeroes$counts*zeroId
genesWithAddedZero <- which(rowSums(zeroId)<nSamp)
pvalsZeroes = pval(dataZeroes, method=selectedMethods, count.type="counts", mc.cores=2, niter=25)
##performance curves
truthZero=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truthZero[dataZeroes$indDE,"status"]=1
cobraZero <- COBRAData(pval =data.frame(edgeRRobust=pvalsZeroes$pval$edgeR_robust, limma_voom=pvalsZeroes$pval$limma_voom, edgeR=pvalsZeroes$pval$edgeR, edgeREMLibSize=pvalsZeroes$pval$edgeREMLibSize, DESeq2=pvalsZeroes$pval$DESeq2, row.names = rownames(dataZeroes)), 
		   padj = data.frame(edgeRRobust=pvalsZeroes$padj$edgeR_robust, limma_voom=pvalsZeroes$padj$limma_voom, edgeR=pvalsZeroes$padj$edgeR, edgeREMLibSize=pvalsZeroes$padj$edgeREMLibSize, DESeq2=pvalsZeroes$padj$DESeq2, row.names = rownames(dataZeroes)),
                   truth = truthZero)
cobraperf <- calculate_performance(cobraZero, binary_truth = "status")
cobraplot <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,0.6))


