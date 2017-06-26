source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v6.R")
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
countsTrapnell <- countsTrapnell[!rowSums(countsTrapnell>0)<5,]
design=matrix(rep(1,ncol(countsTrapnell),ncol=1)) #only estimate intercept
dEmp=DGEList(countsTrapnell)
dEmp=calcNormFactors(dEmp)
dEmp=estimateGLMCommonDisp(dEmp,design, interval=c(0,10))
dEmp=estimateGLMTagwiseDisp(dEmp,design, prior.df=0)
plotBCV(dEmp)
plotBCV(dEmp, col.tagwise=rowSums(dEmp$counts>0))


##### sample comparison smooth scatter plots
library(tweeDEseqCountData)
data(pickrell)
pickrell <- as.matrix(exprs(pickrell.eset))
par(mfrow=c(1,2),mar=c(5,5,3,1))
smoothScatter(x=log10(islam[,1]),y=log10(islam[,2]),xlab="log10(counts+1)", ylab="log10(counts+1)", main="", cex.main=2, cex.lab=2, cex.axis=1.5, xlim=c(0,5), ylim=c(0,5))
smoothScatter(x=log10(pickrell[,1]),y=log10(pickrell[,2]),xlab="log10(counts+1)", ylab="log10(counts+1)", main="", cex.main=2, cex.lab=2, cex.axis=1.5, xlim=c(0,5), ylim=c(0,5))


### conquer EDA over several datasets
library(MultiAssayExperiment) ; library(edgeR)
files=list.files("/Users/koenvandenberge/PhD_Data/singleCell/conquer/")
rdsFiles=files[grep(x=files,pattern=".rds$")]
# discard trimmed data
rdsFiles=rdsFiles[-grep(rdsFiles,pattern="trimmed.rds")]
# discard small datasets
rdsFiles=rdsFiles[-c(2:4)]
rdsNames=c("Buettner, 2015", "Deng, 2014", "Shalek, 2014", "Shalek, 2014", "Trapnell, 2014", "Trapnell, 2014", "Patel, 2014", "Kumar, 2014", "Kumar, 2014", "Guo, 2015", "Engel, 2016", "Meyer, 2016")
rdsNamesSub=rdsNames[-c(4,6,9)]

### P(zero) ~ libsize over datasets in conquer tool.
zeroLibsizeList <- list()
for(i in 1:length(rdsFiles)){
    if(i %in% c(4,6,9)) next
    data=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i]))
    countData <- round(assay(experiments(data)$gene,"count"))    
    if(i==3){ #Shalek was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==5){ #Trapnell was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()	
    }
    if(i==8){ #Kumar was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()	
    }
    dat <- data.frame(logLibSize=log(colSums(countData)), zeroFraction=colMeans(countData==0))
    zeroLibsizeList[[i]]=dat
}
zeroLibsizeList <- zeroLibsizeList[!unlist(lapply(zeroLibsizeList,function(x) is.null(x)))]

#png(filename="~/Desktop/zeroLibSizeConquer.png",width=1000,height=800, res=80)
par(mfrow=c(3,3), mar=c(3,3,1,1), oma=c(3,3,3,3))
for(i in 1:9){ 
    if(i%in%c(1,4,7)) par(mar=c(5,5,1,1)) else par(mar=c(5,4,1,1))
    plot(x=zeroLibsizeList[[i]][,1], y=zeroLibsizeList[[i]][,2], xlab="", ylab="", bty="l", main=rdsNamesSub[i], pch=16, cex=1/3)
    mtext("Fraction of zeros",side=2, outer=TRUE, cex=2)
    mtext("Log library size",side=1, outer=TRUE, cex=2)    
}
#dev.off()

### P(zero) ~ aveLogCPM over datasets in conquer tool.
zeroCpmList <- list()
for(i in 1:length(rdsFiles)){
    if(i %in% c(4,6,9)) next
    data=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i]))
    countData <- round(assay(experiments(data)$gene,"count"))    
    if(i==3){ #Shalek was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==5){ #Trapnell was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()	
    }
    if(i==8){ #Kumar was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()	
    }
    dat <- data.frame(avCpm=aveLogCPM(countData), zeroFraction=rowMeans(countData==0))
    zeroCpmList[[i]]=dat
}
zeroCpmList <- zeroCpmList[!unlist(lapply(zeroCpmList,function(x) is.null(x)))]

png(filename="~/Desktop/zeroCpmConquer.png",width=1000,height=800, res=80)
par(mfrow=c(3,3), mar=c(3,3,1,1), oma=c(3,3,3,3))
for(i in 1:9){ 
    if(i%in%c(1,4,7)) par(mar=c(5,5,1,1)) else par(mar=c(5,4,1,1))
    plot(x=zeroCpmList[[i]][,1], y=zeroCpmList[[i]][,2], xlab="", ylab="", bty="l", main=rdsNamesSub[i], pch=16, cex=1/3)
    mtext("Fraction of zeros",side=2, outer=TRUE, cex=2)
    mtext("Average log CPM",side=1, outer=TRUE, cex=2)    
}
dev.off()

### BCV plot based on 10 samples

### BCV over datasets in conquer tool.
bcvList <- list()
for(i in 1:length(rdsFiles)){
    cat(i)
    if(i %in% c(4,6,9)) next
    data=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i]))
    countData <- round(assay(experiments(data)$gene,"count"))    
    if(i==3){ #Shalek was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()
    }
    if(i==5){ #Trapnell was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()	
    }
    if(i==8){ #Kumar was split up
	data2=readRDS(paste0("/Users/koenvandenberge/PhD_Data/singleCell/conquer/",rdsFiles[i+1]))
	countData1 <- round(assay(experiments(data)$gene,"count"))
	countData2 <- round(assay(experiments(data)$gene,"count"))
	countData = cbind(countData1,countData2)
	rm(data,data2,countData1,countData2); gc()	
    }
    d=DGEList(countData[,1:10])
    d=edgeR::calcNormFactors(d)
    d=estimateDisp(d, prior.df=0)	     
    bcvList[[i]]=d
}
bcvList <- bcvList[!unlist(lapply(bcvList,function(x) is.null(x)))]

png(filename="~/Dropbox/phdKoen/singleCell/figures/supplementary/bcvConquer.png",width=1000,height=800, res=80)
par(mfrow=c(3,3), mar=c(3,3,1,1), oma=c(3,3,3,3))
for(i in 1:9){ 
    if(i%in%c(1,4,7)) par(mar=c(5,5,1,1)) else par(mar=c(5,4,1,1))
    #plot(x=zeroLibsizeList[[i]][,1], y=zeroLibsizeList[[i]][,2], xlab="", ylab="", bty="l", main=rdsNamesSub[i], pch=16, cex=1/3)
    plotBCV(bcvList[[i]], xlab="", ylab="", main=rdsNamesSub[i])
    mtext("Biological coefficient of variation",side=2, outer=TRUE, cex=2)
    mtext("Average Log CPM",side=1, outer=TRUE, cex=2)    
}
dev.off()







