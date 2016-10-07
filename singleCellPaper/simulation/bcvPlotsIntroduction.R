source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")
setwd("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/")

### make BCV plot islam / trapnell
## BCV plot islam
library(GEOquery)
data = read.delim("expressionTabAdapted_kvdb.txt")
seriesMatrix=getGEO(filename="GSE29087_series_matrix.txt")
countData = data[,8:ncol(data)]
rownames(countData)=data[,1]
#(pData(seriesMatrix)$title)[93:96] #negative controls added by the authors: can be discarded.
countData=countData[,1:92]
fibroID <- grep(x=pData(seriesMatrix)$title,pattern="fibroblast")
stemCellID <- grep(x=pData(seriesMatrix)$title,pattern="stem cell")
colnames(countData)[fibroID] = paste("fibro",1:length(fibroID),sep="_")
colnames(countData)[stemCellID] = paste("stemCell",1:length(stemCellID),sep="_")
cellType=vector(length=ncol(countData))
cellType[fibroID] = "fibro"
cellType[stemCellID] <- "stemCell"
islam = as.matrix(countData)
islam = islam[!rowSums(islam)==0,]

design=model.matrix(~cellType)
d=DGEList(islam)
d=calcNormFactors(d)
d=estimateGLMCommonDisp(d,design, interval=c(0,10))
d=estimateGLMTagwiseDisp(d,design, prior.df=0)
dIslam=d
plotBCV(dIslam)

## BCV plot trapnell
library(MultiAssayExperiment)
dataTrapnell <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE52529-GPL11154.rds") #downloaded from conquer tool from Charlotte Soneson
#remove wells containing debris
dataTrapnell <- dataTrapnell[,!pData(dataTrapnell)[,"characteristics_ch1.2"]=="debris: TRUE"]
#remove wells that did not contain one cell
dataTrapnell <- dataTrapnell[,!pData(dataTrapnell)[,"characteristics_ch1.4"]!="cells in well: 1"]
countsTrapnell <- round(assay(experiments(dataTrapnell)$gene,"count"))
countsTrapnell <- countsTrapnell[!rowSums(countsTrapnell)==0,]
design=matrix(rep(1,ncol(countsTrapnell),ncol=1)) #only estimate intercept
dEmp=DGEList(countsTrapnell)
dEmp=calcNormFactors(dEmp)
dEmp=estimateGLMCommonDisp(dEmp,design, interval=c(0,10))
dEmp=estimateGLMTagwiseDisp(dEmp,design, prior.df=0)
plotBCV(dEmp)


### simulate RNA-Seq data with high number of samples and plot BCV
library(tweeDEseqCountData)
data(pickrell)
pickrell <- as.matrix(exprs(pickrell.eset))
nSamp <- 10
grp <- as.factor(rep(0:1, each = nSamp/2))
nTags = 20e3
design=model.matrix(~grp)
#libSize = sample(round(seq(25e5,5e6,length.out=nSamp)))
libSize <- runif(n=nSamp, 5e6,8e6)
dataNoZI <- NBsim(foldDiff = 1, ind=1, dataset = pickrell, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = 0.1, lib.size=libSize, drop.low.lambda=TRUE)
dSim=DGEList(dataNoZI$counts)
dSim=calcNormFactors(dSim)
dSim=estimateGLMCommonDisp(dSim,design,interval=c(0,10))
dSim=estimateGLMTagwiseDisp(dSim,design,prior.df=0)
plotBCV(dSim)

## add zeroes and plot BCV
dataZeroes = dataNoZI
propZeroes=0.1
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
zeroId[sample(1:length(zeroId),floor(length(zeroId)*propZeroes))]=0
dataZeroes$counts = dataZeroes$counts*zeroId
genesWithAddedZero <- which(rowSums(zeroId)<nSamp)
dZero=DGEList(dataZeroes$counts)
dZero=calcNormFactors(dZero)
dZero=estimateGLMCommonDisp(dZero,design,interval=c(0,10))
dZero=estimateGLMTagwiseDisp(dZero,design,prior.df=0)
plotBCV(dZero)

## plot BCV after zeroWeights
dSimZero <- DGEList(dataZeroes$counts)
dSimZero <- calcNormFactors(dSimZero)
zeroWeights <- zeroWeightsLibSize(counts=dSimZero$counts, design=design, niter=30, initialWeightAt0=TRUE, plotW=TRUE)
dSimZero$weights <- zeroWeights
dSimZero <- estimateGLMCommonDisp(dSimZero,design,interval=c(0,10))
dSimZero <- estimateGLMTagwiseDisp(dSimZero,design,prior.df=0)
plotBCV(dSimZero)


## plot for in paper
#bcvPlotsIntro.pdf
par(mfrow=c(2,2), mar=c(5,4,1,1))
plotBCV(dEmp)
plotBCV(dSim, ylim=c(0,2.5))

plotBCV(dZero)
plotBCV(dSimZero,ylim=c(0,2.5))


##### sample comparison smooth scatter plots
par(mfrow=c(1,2),mar=c(5,5,3,1))
smoothScatter(x=log10(islam[,1]),y=log10(islam[,2]),xlab="log10 counts", ylab="log10 counts", main="", cex.main=2, cex.lab=2, cex.axis=1.5, xlim=c(0,5), ylim=c(0,5))
smoothScatter(x=log10(pickrell[,1]),y=log10(pickrell[,2]),xlab="log10 counts", ylab="log10 counts", main="", cex.main=2, cex.lab=2, cex.axis=1.5, xlim=c(0,5), ylim=c(0,5))



