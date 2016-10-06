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
dataTrapnell <- readRDS("/Users/koenvandenberge/PhD_Data/singleCell/conquer/GSE52529-GPL11154.rds")
countsTrapnell <- round(assay(experiments(dataTrapnell)$gene,"count"))
countsTrapnell <- countsTrapnell[!rowSums(countsTrapnell)==0,]
condition = factor(as.numeric(pData(dataTrapnell)[,"characteristics_ch1.2"]))
design = model.matrix(~condition)
dEmp=DGEList(countsTrapnell)
dEmp=calcNormFactors(dEmp)
dEmp=estimateGLMCommonDisp(dEmp,design, interval=c(0,10))
dEmp=estimateGLMTagwiseDisp(dEmp,design, prior.df=0)
plotBCV(dEmp)


### simulate RNA-Seq data with high number of samples



## add zeroes and plot BCV

## plot BCV after zeroWeights
