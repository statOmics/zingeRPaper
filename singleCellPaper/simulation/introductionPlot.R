source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v5.R")
library(scales)
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


design=model.matrix(~cellType)
dIslam=DGEList(islam)
dIslam=calcNormFactors(dIslam)
dIslam=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dIslam,design, interval=c(0,10)),design,prior.df=0)
dIslamTrend=estimateDisp(dIslam,design,prior.df=0)

#zeroWeights=zeroWeightsLibSize(counts=islam,niter=30,design=model.matrix(~cellType))
zeroWeights=zeroWeightsLibSizeDispFast(counts=islam,maxit=100,design=model.matrix(~cellType))
dW=DGEList(islam)
dW=calcNormFactors(dW)
dW$weights=zeroWeights
dW=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dW,design, interval=c(0,10)),design,prior.df=0)
dWTrend=estimateDispWeighted(dW,design,prior.df=0,weights=dW$weights)

###############################
### plot for paper:scRNA-seq ##
###############################
par(mar=c(4.1,4.25,3,1),bty="l", mfrow=c(1,4), cex.lab=1.5, cex.axis=1.5)
plot(dIslam$AveLogCPM,sqrt(dIslam$tagwise.dispersion),pch=16,cex=.2,col=alpha("black",1/3),xlim=c(1.8,12), xlab="Average Log CPM", ylab="BCV")
o <- order(dIslam$AveLogCPM)
lines(dIslam$AveLogCPM[o], sqrt(dIslamTrend$trended.dispersion)[o], col = "red",lwd = 2)
#lines(lowess(x=dIslam$AveLogCPM,y=sqrt(dIslam$tagwise.dispersion),f=1/6),col=2,lwd=2)

plot(x=log(colSums(islam)),y=colMeans(islam==0),pch=19,cex=2/3, xlab="Log library size", ylab="Fraction of zeros", ylim=c(0.05,0.95))
## use weights to extract model
w=zeroWeights
successes <- colSums(1-w) #P(zero)
failures <- colSums(w) #1-P(zero)
counts=islam
counts <- DGEList(counts)
counts <- edgeR::calcNormFactors(counts)
effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
logEffLibSize <- log(effLibSize)
zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")
grid=seq(min(logEffLibSize),max(logEffLibSize),length.out=100)
yHatZero=predict(zeroFit,newdata=data.frame(logEffLibSize=grid),type="response")
lines(x=grid,y=yHatZero,col="salmon",lwd=2)

hist(zeroWeights[islam==0],main="",xlab="Posterior probability")

plot(dW$AveLogCPM,sqrt(dW$tagwise.dispersion),pch=16,cex=.2,xlab="Average Log CPM", ylab="BCV", xlim=c(1.8,12), ylim=c(0,12),col=alpha("black",1/3))
o <- order(dW$AveLogCPM)
lines(dIslam$AveLogCPM[o], sqrt(dIslamTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dW$AveLogCPM[o], sqrt(dWTrend$trended.dispersion)[o], col = "steelblue1",lwd = 2)
#lines(lowess(x=dIslam$AveLogCPM,y=sqrt(getDispersion(dIslam)),f=1/6),col=2,lwd=2)
#lines(lowess(x=dW$AveLogCPM,y=sqrt(getDispersion(dW)),f=1/6),col="steelblue",lwd=2)
legend("topright",c("NB model, scRNA-seq","ZINB model, scRNA-seq"),lty=1,col=c("red","steelblue1"),lwd=2,bty="n")



#########################################################
############# RNA-Seq data simulation ###################
#########################################################
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v5.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

#### no zero inflation simulation
library(Biobase)
load("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/data/bottomly_eset.RData")
bottomly=exprs(bottomly.eset)
bottomly=bottomly[!rowSums(bottomly)==0,]
nSamp <- 10
nTags <- 20e3
set.seed(2)
grp <- as.factor(rep(0:1, each = nSamp/2))
libSize = sample(round(seq(8e6,10e6,length.out=nSamp)))
DEind = sample(1:nTags,floor(nTags/20),replace=FALSE) #5% differentially expressed
fcSim <- (2 + rexp(length(DEind), rate = 1)) #adapted from Soneson et al. 2016, Genome Biology
set.seed(1)
dataNoZI <- NBsim(foldDiff = fcSim, ind=DEind, dataset = bottomly, nTags = nTags, group = grp, verbose = TRUE, add.outlier = FALSE, drop.extreme.dispersion = FALSE, lib.size=libSize, drop.low.lambda=TRUE)

data=dataNoZI$counts
d=DGEList(data)
d=calcNormFactors(d)
design=model.matrix(~grp)
d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
dOriginal=d
dOriginalTrend=estimateGLMTrendedDisp(d,design)
plotBCV(d, ylim=c(0,3))


# add zeroes: here we see an obvious gain
dataZeroes = dataNoZI
propZeroes=0.05
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
set.seed(3)
samp=sample(1:length(zeroId),floor(length(zeroId)*propZeroes))
zeroId[samp]=0
zeroId[dataNoZI$counts==0]=1 #if it already was a zero it is not zero-inflated.
samp=samp[!samp%in%which(dataNoZI$counts==0)] #same
dataZeroes$counts = dataZeroes$counts*zeroId
data=dataZeroes$counts
d=DGEList(data)
d=calcNormFactors(d)
design=model.matrix(~grp)
d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
#d=estimateGLMTagwiseDisp(estimateGLMTrendedDisp(estimateGLMCommonDisp(d,design),design),design,prior.df=0)
dZeroes=d
dZeroesTrend=estimateGLMTrendedDisp(dZeroes,design)
plotBCV(d, ylim=c(0,3))

#zeroWeightsSim=zeroWeightsLibSize(counts=dataZeroes$counts,design=model.matrix(~grp),niter=30)
zeroWeightsSim=zeroWeightsLibSizeDispFast(counts=dataZeroes$counts,design=model.matrix(~grp),maxit=100)
hist(zeroWeightsSim[samp],xlab="post. prob. on count component")

## downweighting
d=DGEList(data)
d=calcNormFactors(d)
d$weights=zeroWeightsSim
d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
dWeighted=d
dWeightedTrend=estimateGLMTrendedDisp(dWeighted,design)
plotBCV(d, ylim=c(0,3))


##### plot for paper: RNA-seq
par(mar=c(4.1,4.25,3,1),bty="l", mfrow=c(1,4), cex.lab=1.5, cex.axis=1.5)
plot(x=dOriginal$AveLogCPM,y=sqrt(dOriginal$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dOriginal$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)

plot(x=dZeroes$AveLogCPM,y=sqrt(dZeroes$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dZeroes$AveLogCPM)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)

hist(zeroId[zeroId==0],xlim=c(0,.95),breaks=seq(0,1,.05),main="",xlab="Posterior probability")
hist(zeroWeightsSim[samp],add=TRUE,breaks=seq(0,1,.05),col=rgb(0.1,0.8,0.1,.2))
legend("topright",c("truth","zingeR"),fill=c(0,rgb(0.1,0.8,0.1,.2)), bty="n",cex=1.75)

#hist(zeroWeightsSim[dataZeroes$counts==0],xlab="Posterior probability",main="")
#hist((zeroId==0)+0, add=TRUE,col=rgb(0.1,0.8,0.1,.2))
#hist(zeroWeightsSim[samp],xlab="Posterior probability",main="")

plot(x=dWeighted$AveLogCPM,y=sqrt(dWeighted$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dWeighted$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
lines(dWeighted$AveLogCPM[o], sqrt(dWeightedTrend$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, RNA-seq","NB model, ZI RNA-seq","ZINB model, ZI RNA-seq"),lty=1,lwd=2,col=c("red","blue","steelblue1"),bty="n", cex=1.5)




#######################################
### one composite plot ################
#######################################
wMean=sapply(1:nrow(zeroWeightsSim), function(i){
	  mean(zeroWeightsSim[i,dataZeroes$counts[i,]==0])
})
wMean[is.na(wMean)]=1
library(Hmisc)
cuts=cut2(wMean,cuts=seq(0,1,by=0.25))


dev.new(width=10,height=5)
##### plot for paper: RNA-seq
library(scales)
par(mar=c(4.1,4.25,3,1),bty="l", mfrow=c(2,4), cex.lab=1.5, cex.axis=1.5)
plot(x=dOriginal$AveLogCPM,y=sqrt(dOriginal$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dOriginal$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)
mtext("a" ,at=-5, font=2, cex=4/3)

#plot(x=dZeroes$AveLogCPM,y=sqrt(dZeroes$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
#cols=colorRampPalette(c("red","yellow","springgreen","royalblue"))(12)


cols=c("red","orange","salmon","black")
plot(x=dZeroes$AveLogCPM,y=sqrt(dZeroes$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha(cols[as.numeric(cuts)],1/3), type="n")
sapply(4:1,function(i){
	points(x=dZeroes$AveLogCPM[cuts==levels(cuts)[i]],y=sqrt(dZeroes$tagwise.dispersion)[cuts==levels(cuts)[i]],pch=16,cex=.2, col=alpha(cols[i]))
	   
})
o<- order(dZeroes$AveLogCPM)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
legend("topleft",legend=c("[0,0.25)", "[0.25,0.5)", "[0.5,0.75)", "[0.75,1]" ),col=cols[1:5],lty=1, bty="n",lwd=2, cex=.8,inset=c(0,-0.045))
mtext("b" ,at=-5, font=2, cex=4/3)

hist(zeroId[zeroId==0],xlim=c(0,.95),breaks=seq(0,1,.05),main="",xlab="Posterior probability")
hist(zeroWeightsSim[samp],add=TRUE,breaks=seq(0,1,.05),col=rgb(0.1,0.8,0.1,.2))
legend("topright",c("truth","zingeR"),fill=c(0,rgb(0.1,0.8,0.1,.2)), bty="n", cex=1.5)
mtext("c" ,at=-0.22, font=2, cex=4/3)

plot(x=dWeighted$AveLogCPM,y=sqrt(dWeighted$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dWeighted$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
lines(dWeighted$AveLogCPM[o], sqrt(dWeightedTrend$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, RNA-seq","NB model, ZI RNA-seq","ZINB model, ZI RNA-seq"),lty=1,lwd=2,col=c("red","blue","steelblue1"),bty="n")
mtext("d" ,at=-5, font=2, cex=4/3)


### scRNA-seq
plot(dIslam$AveLogCPM,sqrt(dIslam$tagwise.dispersion),pch=16,cex=.2,col=alpha("black",1/3),xlim=c(1.8,12), xlab="Average Log CPM", ylab="BCV")
o <- order(dIslam$AveLogCPM)
lines(dIslam$AveLogCPM[o], sqrt(dIslamTrend$trended.dispersion)[o], col = "red",lwd = 2)
mtext("e" ,at=-0.5, font=2, cex=4/3)

plot(x=log(colSums(islam)),y=colMeans(islam==0),pch=19,cex=2/3, xlab="Log library size", ylab="Fraction of zeros", ylim=c(0.05,0.95))
## use weights to extract model
w=zeroWeights
successes <- colSums(1-w) #P(zero)
failures <- colSums(w) #1-P(zero)
counts=islam
counts <- DGEList(counts)
counts <- edgeR::calcNormFactors(counts)
effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
logEffLibSize <- log(effLibSize)
zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")
grid=seq(min(logEffLibSize),max(logEffLibSize),length.out=100)
yHatZero=predict(zeroFit,newdata=data.frame(logEffLibSize=grid),type="response")
lines(x=grid,y=yHatZero,col="salmon",lwd=2)
mtext("f" ,at=9, font=2, cex=4/3)

hist(zeroWeights[islam==0],main="",xlab="Posterior probability")
mtext("g" ,at=-.22, font=2, cex=4/3)

plot(dW$AveLogCPM,sqrt(dW$tagwise.dispersion),pch=16,cex=.2,xlab="Average Log CPM", ylab="BCV", xlim=c(1.8,12), ylim=c(0,12),col=alpha("black",1/3))
o <- order(dW$AveLogCPM)
lines(dIslam$AveLogCPM[o], sqrt(dIslamTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dW$AveLogCPM[o], sqrt(dWTrend$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, scRNA-seq","ZINB model, scRNA-seq"),lty=1,col=c("red","steelblue1"),lwd=2,bty="n")
mtext("h" ,at=-0.5, font=2, cex=4/3)




