source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")
library(DESeq2) ; library(edgeR) ; library(limma) ; library(scde) ; library(MAST) ; library(mgcv) ; library(MultiAssayExperiment) ; library(SummarizedExperiment)

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
islam = islam[!rowSums(islam)==0,]

## get gene-wise parameters acc to zero-truncated NB distribution
#paramsIslamDropDisp=getDatasetNB(islam,drop.low.lambda=FALSE,drop.extreme.dispersion=.1)
#save(paramsIslamDropDisp,file="paramsIslamDropDisp10Percent.RData")
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/paramsIslamDropDisp10Percent.RData")
set.seed(12)
nSamp <- 80
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 20e3
DEind = sample(1:nTags,floor(nTags*.05),replace=FALSE) #5% differentially expressed
fcSim=(3 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
libSizes=round(runif(n=nSamp,min=5e4,max=3e6))
dataIslam <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = islam, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, params=paramsIslamDropDisp, noise=FALSE, lib.size=libSizes, randomZero=0.05)

### compare simulated and empirical data
#empirical BCV
condition = cellType
design = model.matrix(~condition)
dIslam=DGEList(islam)
dIslam=calcNormFactors(dIslam)
dIslam=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dIslam,design, interval=c(0,10)),design,prior.df=0)
plotBCV(dIslam,main="Empirical", ylim=c(0,15))
#simulated BCV
design=model.matrix(~grp)
dSim=DGEList(dataIslam$counts)
dSim=calcNormFactors(dSim)
dSim=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dSim,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dSim, main="Simulated", ylim=c(0,15))
## marginal density
plot(density(log(islam+1)))
lines(density(log(dataIslam$counts+1)),col=2)
legend("topright",c("empirical","simulated"),lty=1,col=c("black","red"))
## zero ~ libSize
plot(x=log(colSums(islam)),y=colMeans(islam==0))
plot(x=log(colSums(dataIslam$counts)),y=colMeans(dataIslam$counts==0))
## zero ~ cpm
plot(x=aveLogCPM(islam),y=rowMeans(islam==0))
plot(x=aveLogCPM(dataIslam$counts),y=rowMeans(dataIslam$counts==0))

### performance characteristics
#selectedMethods <- c("edgeR", "DESeq2", "limma_voom", "edgeR_robust", "scde", "MAST",  "edgeREMLibSizeOldF")
selectedMethods <- c("DESeq2", "edgeREMLibSizeFastOldF")
pvalsIslam <- pval(dataIslam, method=selectedMethods, count.type="counts", mc.cores=2, niter=30)
#FDR-TPR
library(iCOBRA)
truthIslam=data.frame(status=rep(0,nTags), row.names=rownames(dataIslam))
truthIslam[dataIslam$indDE,"status"]=1
cobraIslam <- COBRAData(pval =data.frame(edgeRRobust=pvalsIslam$pval$edgeR_robust, limma_voom=pvalsIslam$pval$limma_voom, edgeR=pvalsIslam$pval$edgeR,  edgeREMLibSize=pvalsIslam$pval$edgeREMLibSizeOldF, DESeq2=pvalsIslam$pval$DESeq2, scde=pvalsIslam$pval$scde, MAST=pvalsIslam$pval$MAST, row.names = rownames(dataIslam)), 
		   padj = data.frame(edgeRRobust=pvalsIslam$padj$edgeR_robust, limma_voom=pvalsIslam$padj$limma_voom, edgeR=pvalsIslam$padj$edgeR, edgeREMLibSize=pvalsIslam$padj$edgeREMLibSizeOldF, DESeq2=pvalsIslam$padj$DESeq2, scde=pvalsIslam$padj$scde, MAST=pvalsIslam$padj$MAST,  row.names = rownames(dataIslam)),
                   truth = truthIslam)
cobraperf <- calculate_performance(cobraIslam, binary_truth = "status")
cobraplot <- prepare_data_for_plot(cobraperf)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))

##################################################################
########################### TRAPNELL #############################
##################################################################
trapnellAssay <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE52529-GPL11154.rds")
#remove wells containing debris
trapnellAssay <- trapnellAssay[,!pData(trapnellAssay)[,"characteristics_ch1.2"]=="debris: TRUE"]
#remove wells that did not contain one cell
trapnellAssay <- trapnellAssay[,!pData(trapnellAssay)[,"characteristics_ch1.4"]!="cells in well: 1"]
countsTrapnell <- round(assay(experiments(trapnellAssay)$gene,"count"))
countsTrapnell <- countsTrapnell[!rowSums(countsTrapnell)==0,]

#paramsTrapnellDropDisp10=getDatasetNB(countsTrapnell,drop.low.lambda=FALSE,drop.extreme.dispersion=.10)
#save(paramsTrapnellDropDisp10,file="/Users/koenvandenberge/PhD_Data/singleCell/conquer/paramsTrapnellDropDisp10Percent.RData")
load("/Users/koenvandenberge/PhD_Data/singleCell/conquer/paramsTrapnellDropDisp10Percent.RData")
#paramsTrapnellDropDisp15=getDatasetNB(countsTrapnell,drop.low.lambda=FALSE,drop.extreme.dispersion=.15)
#save(paramsTrapnellDropDisp15,file="/Users/koenvandenberge/PhD_Data/singleCell/conquer/paramsTrapnellDropDisp15Percent.RData")
#load("/Users/koenvandenberge/PhD_Data/singleCell/conquer/paramsTrapnellDropDisp15Percent.RData")
#load("/Users/koenvandenberge/PhD_Data/singleCell/conquer/paramsTrapnellDropDisp10Percent.RData")


#paramsTrapnellDropDisp20DropLowLambda=getDatasetNB(countsTrapnell,drop.low.lambda=TRUE,drop.extreme.dispersion=.2)
#save(paramsTrapnellDropDisp20DropLowLambda,file="/Users/koenvandenberge/PhD_Data/singleCell/conquer/paramsTrapnellDropDisp20DropLowLambda.Rda")
#load("/Users/koenvandenberge/PhD_Data/singleCell/conquer/paramsTrapnellDropDisp20DropLowLambda.Rda")


set.seed(9)
nSamp <- 80
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 20e3
DEind = sample(1:nTags,floor(nTags*.05),replace=FALSE) #5% differentially expressed
fcSim=(3 + rexp(length(DEind), rate = 1/2)) #adapted from Soneson et al. 2016, Genome Biology
libSizes=round(runif(n=nSamp,min=1e5,max=3e6))
dataTrapnell <- NBsimSingleCell(foldDiff = fcSim, ind=DEind, dataset = countsTrapnell, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, params=paramsTrapnellDropDisp10, noise=TRUE, lib.size=libSizes, randomZero=0.05)

selectedMethods <- c("edgeREstDisp", "DESeq2", "edgeREMLibSizeFastOldF")
pvalsTrapnellNoise <- pval(dataTrapnell, method=selectedMethods, count.type="counts", mc.cores=2, niter=30)

#empirical BCV
design=matrix(rep(1,ncol(countsTrapnell),ncol=1)) #only estimate intercept
dEmp=DGEList(countsTrapnell)
dEmp=calcNormFactors(dEmp)
dEmp=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dEmp,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dEmp,main="Empirical")
plotBCV(dEmp,main="Empirical", col.tagwise=rowSums(countsTrapnell>0))
#simulated BCV
design=model.matrix(~grp)
dSim=DGEList(dataTrapnell$counts)
dSim=calcNormFactors(dSim)
dSim=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dSim,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dSim, main="Simulated")
## density
plot(density(log(countsTrapnell+1)))
lines(density(log(dataTrapnell$counts+1)),col=2)
legend("topright",c("empirical","simulated"),lty=1,col=c("black","red"))
## zero ~ libSize
plot(x=log(colSums(countsTrapnell)),y=colMeans(countsTrapnell==0))
plot(x=log(colSums(dataTrapnell$counts)),y=colMeans(dataTrapnell$counts==0))
## zero ~ cpm
plot(x=aveLogCPM(countsTrapnell),y=rowMeans(countsTrapnell==0))
plot(x=aveLogCPM(dataTrapnell$counts),y=rowMeans(dataTrapnell$counts==0))
library(iCOBRA)

## FDR-TPR
library(iCOBRA)
truth=rep(0,length(pvalsTrapnellNoise$padj$edgeREstDisp))
truth[pvalsTrapnellNoise$indDE]=1
truth=as.data.frame(truth)
rownames(truth)=1:nrow(truth)
cobra <- COBRAData(padj = data.frame(DESeq2 = pvalsTrapnellNoise$padj$DESeq2,
				     edgeREMLibSizeFast = pvalsTrapnellNoise$padj$edgeREMLibSizeFastOldF,
				     edgeR = pvalsTrapnellNoise$padj$edgeREstDisp,
                                     row.names = 1:length(pvalsTrapnellNoise$padj$edgeREstDisp),
                                     stringsAsFactors = FALSE))
cobra <- COBRAData(truth = truth, object_to_extend = cobra)
cobraperf <- calculate_performance(cobra, binary_truth = "truth")
cobraplot <- prepare_data_for_plot(cobraperf, incltruth = TRUE, 
                                   facetted = FALSE)
plot_fdrtprcurve(cobraplot)
plot_roc(cobraplot)














