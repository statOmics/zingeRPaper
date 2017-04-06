source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v5.R")
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



### bottomly
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/bottomly_eset.RData")
strain=pData(bottomly.eset)$strain
design=model.matrix(~strain)
bottomly=exprs(bottomly.eset)
dB=DGEList(bottomly)
dB=calcNormFactors(dB)
dBOrig=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dB,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dBOrig)

addZero <- function(idVector,p,nGenes){
    samp <- sample(idVector,floor(p*nGenes))
    idVector2 <- idVector[!idVector%in%samp]
    return(list(samp,idVector2))
}
pZeroTrapnell <- sapply(1:60, function(x) mean(rowSums(countsTrapnell>0)==x))
a=aveLogCPM(dBOrig$counts,offset=getOffset(dBOrig))
a=a-min(a)+0.1
prob=a/sum(a)
prob[a<4.5]=prob[a<4.5]*20
id=sample(1:nrow(bottomly),7500,prob=prob)
start=5
end=10
bottomly2=bottomly[id,]
nGenes=nrow(bottomly2)
nCol=ncol(bottomly2)
idVector=1:nrow(bottomly2)
for(i in start:end){
    hlp <- addZero(idVector=idVector, p=pZeroTrapnell[i-start+1], nGenes=nGenes)
    idVector=hlp[[2]]
    bottomly2[hlp[[1]],sample(1:nCol,nCol-i)]=0
}
addedZeroId=!(bottomly2==bottomly[id,])
bottomly[id,]=bottomly2
dB=DGEList(bottomly)
dB=calcNormFactors(dB)
dB=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dB,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dB, ylim=c(0,6))

##downweighted
dBW=dB
w=matrix(1,nrow=nrow(dB),ncol=ncol(dB))
wSub=w[id,]
wSub[which(addedZeroId)]=0
w[id,]=wSub
dBW$weights=w
keep=rowSums(dBW$counts)>0
dBW=dBW[keep,]
dBW=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dBW,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dBW)

#bcvPlotsIntro.pdf
png("~/Dropbox/phdKoen/singleCell/figures/bcvPlotsIntro_v2.png",width=12,height=8,units="in",res=300)
par(mfrow=c(2,2), mar=c(3,4,0.5,0.5), oma=c(.25,.25,0,0))
plotBCVNoLegend(dEmp,xlab="",ylab="", xlim=c(0,15))
plotBCVNoLegend(dBOrig,xlab="",ylab="", xlim=c(-1,15))
par(mar=c(4,4,0.5,0.5))
plotBCVNoLegend(dB, ylim=c(0,6),xlab="",ylab="", xlim=c(-1,15))
plotBCVNoLegend(dBW,xlab="",ylab="", xlim=c(0,15))
mtext(text="Biological coefficient of variation",side=2,adj=.5,outer=TRUE,line=-2)
mtext(text="Average Log CPM",side=1,adj=.5,outer=TRUE,line=-2)
dev.off()

##### sample comparison smooth scatter plots
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



















#### OLD CODE
### Pickrell BCV plot
library(tweeDEseqCountData)
data(pickrell)
pickrell <- as.matrix(exprs(pickrell.eset))
gender=pData(pickrell.eset)$gender
design=model.matrix(~gender)
dPickrell <- DGEList(pickrell)
dPickrell <- calcNormFactors(dPickrell)
dPickrellOrig <- estimateGLMTagwiseDisp(estimateGLMCommonDisp(dPickrell,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dPickrellOrig)


A=aveLogCPM(dPickrellOrig$counts, offset = getOffset(dPickrellOrig))
d=sqrt(getDispersion(dPickrellOrig))
plot(x=A,y=d, ylim=c(0,4), type="n" , xlab = "Average log CPM", ylab = "Biological coefficient of variation")
points(A, d, pch = 16, cex = 0.2)

### Pickrell BCV plot with augmented zero counts
data(pickrell)
pickrell <- as.matrix(exprs(pickrell.eset))
addZero <- function(idVector,p,nGenes){
    samp <- sample(idVector,floor(p*nGenes))
    idVector2 <- idVector[!idVector%in%samp]
    return(list(samp,idVector2))
}
pZeroTrapnell <- sapply(1:60, function(x) mean(rowSums(countsTrapnell>0)==x))
idVector=1:nrow(pickrell)
nGenes=nrow(pickrell)
nCol=ncol(pickrell)
start=1
end=60
for(i in start:end){
    hlp <- addZero(idVector=idVector, p=pZeroTrapnell[i-start+1], nGenes=nGenes)
    idVector=hlp[[2]]
    pickrell[hlp[[1]],sample(1:nCol,nCol-i)]=0
}
dPickrell <- DGEList(pickrell)
dPickrell <- calcNormFactors(dPickrell)
dPickrell <- estimateGLMTagwiseDisp(estimateGLMCommonDisp(dPickrell,design,interval=c(0,10)),design,prior.df=0)
plotBCV(dPickrell)


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
par(mfrow=c(2,3), mar=c(5,4,1,1))
plotBCV(dEmp)
plotBCV(dSim, ylim=c(0,2.5))

plotBCV(dZero)
plotBCV(dSimZero,ylim=c(0,2.5))



