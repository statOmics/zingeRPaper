source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v6.R")
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

#zeroWeights=zeroWeightsLibSize(counts=islam, niter=30, design=model.matrix(~cellType))
zeroWeights=zeroWeightsLibSizeDispFast(counts=islam,maxit=100,design=model.matrix(~cellType))
dW=DGEList(islam)
dW=calcNormFactors(dW)
dW$weights=zeroWeights
#dW=estimateGLMTagwiseDisp(estimateGLMCommonDisp(dW,design, interval=c(0,10)),design,prior.df=0)
#dWTrend=estimateDispWeighted(dW,design,prior.df=0,weights=dW$weights)
dW=estimateDisp(dW,design,prior.df=0)

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
lines(dW$AveLogCPM[o], sqrt(dW$trended.dispersion)[o], col = "steelblue1",lwd = 2)
#lines(lowess(x=dIslam$AveLogCPM,y=sqrt(getDispersion(dIslam)),f=1/6),col=2,lwd=2)
#lines(lowess(x=dW$AveLogCPM,y=sqrt(getDispersion(dW)),f=1/6),col="steelblue",lwd=2)
legend("topright",c("NB model, scRNA-seq","ZINB model, scRNA-seq"),lty=1,col=c("red","steelblue1"),lwd=2,bty="n")



#########################################################
############# RNA-Seq data simulation ###################
#########################################################
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v6.R")
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

### ROC curve for identifying excess zeros
#pvalSeq = c(1e-15,1e-14,1e-13,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
#falses=which(dataNoZI$counts==0)
#tpr=fpr=vector(length=length(pvalSeq))
#for(i in 1:length(pvalSeq)){
#    excessID <- which(zeroWeightsSim<=pvalSeq[i])
#    tpr[i] <- mean(samp%in%excessID)
#    fpr[i] <- mean(falses%in%excessID)
#}
#plot(x=fpr,y=tpr,type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")
#points(x=fpr[pvalSeq==1e-5],y=tpr[pvalSeq==1e-5],col=2,pch=19)
#points(x=fpr[pvalSeq==1e-2],y=tpr[pvalSeq==1e-2],col=2,pch=19) #w=0.01
#points(x=fpr[519],y=tpr[519],col=2,pch=19) #w=0.05

weights=zeroWeightsLibSizeDispFast(counts=dataZeroes$counts, design=model.matrix(~grp), maxit=200, plotW=TRUE)
pvalSeq = c(1e-15,1e-14,1e-13,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
## ROC curve stratified by average expression
d=DGEList(dataZeroes$counts)
d=calcNormFactors(d)
acpm=aveLogCPM(d)
library(Hmisc)
cuts=cut2(acpm,g=3)
tpr2=fpr2=list()
for(k in 1:length(levels(cuts))){
    tpr2[[k]]=vector(length=length(pvalSeq))
    fpr2[[k]]=vector(length=length(pvalSeq))
}
for(j in 1:length(levels(cuts))){
    groupID <- cuts==levels(cuts)[j]
    wSub=weights[groupID,]
    falses=which(dataNoZI$counts[cuts==levels(cuts)[j]]==0)
    trues=which(zeroId[cuts==levels(cuts)[j],]==0)
    for(i in 1:length(pvalSeq)){
        excessID <- which(wSub<=pvalSeq[i])
	tpr2[[j]][i] <- mean(trues%in%excessID)
    	fpr2[[j]][i] <- mean(falses%in%excessID)
    }
}
cols=colorRampPalette(c("skyblue","darkblue"))(length(levels(cuts)))
plot(x=fpr2[[1]],y=tpr2[[1]],col="steelblue",type="n", xlab="False positive rate", ylab="True positive rate")
for(k in 1:length(levels(cuts))) lines(x=fpr2[[k]],y=tpr2[[k]],col=cols[k], lty=k, lwd=2)
legend("bottomright",levels(cuts),col=cols,lty=1:length(levels(cuts)), cex=4/3, lwd=2)






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
#cuts=cut2(wMean,cuts=seq(0,1,by=0.25))
cuts=cut2(wMean,cuts=c(0,0.25,0.75,1))


#dev.new(width=10,height=5)
png("introPlotZingeR.png", width=10,height=5, units="in", res=330)
##### plot for paper: RNA-seq
library(scales)
par(mar=c(4.1,4.25,3,1),bty="l", mfrow=c(2,4), cex.lab=1.5, cex.axis=1.5)
plot(x=dOriginal$AveLogCPM,y=sqrt(dOriginal$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dOriginal$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)
mtext("a" ,at=-5, font=2, cex=4/3)

#plot(x=dZeroes$AveLogCPM,y=sqrt(dZeroes$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
#cols=colorRampPalette(c("red","yellow","springgreen","royalblue"))(12)


#cols=c("red","orange","salmon","black")
cols=c("red","gold","black")
plot(x=dZeroes$AveLogCPM,y=sqrt(dZeroes$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha(cols[as.numeric(cuts)],1/3), type="n")
#sapply(c(4,1,2,3),function(i){
sapply(c(3,1,2),function(i){
	points(x=dZeroes$AveLogCPM[cuts==levels(cuts)[i]],y=sqrt(dZeroes$tagwise.dispersion)[cuts==levels(cuts)[i]],pch=16,cex=.2, col=alpha(cols[i]))
})
o<- order(dZeroes$AveLogCPM)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
#legend("topleft",legend=c("[0,0.25)", "[0.25,0.5)", "[0.5,0.75)", "[0.75,1]" ),col=cols[1:5],lty=1, bty="n",lwd=2, cex=.8,inset=c(0,-0.045))
legend("topleft",legend=c("[0,0.25)", "[0.25,0.75)", "[0.75,1]" ),col=cols,lty=1, bty="n",lwd=2, cex=.9,inset=c(0,-0.045))
mtext("b" ,at=-5, font=2, cex=4/3)

#hlpHist = hist(zeroId[zeroId==0], breaks=seq(0,1,.05), plot=FALSE)
#hlpHist$counts[hlpHist$counts>0] = log(hlpHist$counts[hlpHist$counts>0])
#plot(hlpHist)
#hlpHistSim = hist(zeroWeightsSim[samp],breaks=seq(0,1,.05), plot=FALSE)
#hlpHistSim$counts[hlpHistSim$counts>0] = log(hlpHistSim$counts[hlpHistSim$counts>0])
#plot(hlpHistSim, add=TRUE, col=rgb(0.1,0.8,0.1,.2))
hist(zeroId[zeroId==0],xlim=c(0,.95),breaks=seq(0,1,.05),main="",xlab="Posterior probability")
hist(zeroWeightsSim[samp],add=TRUE,breaks=seq(0,1,.05),col=rgb(0.1,0.8,0.1,.2))
legend("topright",c("nr. true excess zeros","zingeR probabilities"),fill=c(0,rgb(0.1,0.8,0.1,.2)), bty="n", cex=1.25)
mtext("c" ,at=-0.22, font=2, cex=4/3)

plot(x=dWeighted$AveLogCPM,y=sqrt(dWeighted$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dWeighted$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
lines(dWeighted$AveLogCPM[o], sqrt(dWeightedTrend$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, NB simul.","NB model, ZINB simul.","ZINB model, ZINB simul."),lty=1,lwd=2,col=c("red","blue","steelblue1"),bty="n")
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
lines(dW$AveLogCPM[o], sqrt(dW$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, scRNA-seq","ZINB model, scRNA-seq"),lty=1,col=c("red","steelblue1"),lwd=2,bty="n")
mtext("h" ,at=-0.5, font=2, cex=4/3)
dev.off()




#######################################
### one composite plot, version 2 ################
#######################################
wMean=sapply(1:nrow(zeroWeightsSim), function(i){
	  mean(zeroWeightsSim[i,dataZeroes$counts[i,]==0])
})
wMean[is.na(wMean)]=1
library(Hmisc)
cuts=cut2(wMean,cuts=c(0,0.25,0.75,1))


#png("~/Dropbox/phdKoen/singleCell/figures/introBCV_v2.png", width=10,height=5, units="in", res=330)
pdf("~/Dropbox/phdKoen/singleCell/figures/introBCV_v2.pdf", width=10,height=5)
#dev.new(width=10,height=5)
##### plot for paper: RNA-seq
library(scales)
par(mar=c(4.1,4.25,3,1),bty="l", mfrow=c(2,4), cex.lab=1.5, cex.axis=1.5)

cols=c("red","gold","black")
plot(x=dZeroes$AveLogCPM,y=sqrt(dZeroes$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha(cols[as.numeric(cuts)],1/3), type="n")
#sapply(c(4,1,2,3),function(i){
sapply(c(3,1,2),function(i){
	points(x=dZeroes$AveLogCPM[cuts==levels(cuts)[i]],y=sqrt(dZeroes$tagwise.dispersion)[cuts==levels(cuts)[i]],pch=16,cex=.2, col=alpha(cols[i]))
})
o<- order(dZeroes$AveLogCPM)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
#legend("topleft",legend=c("[0,0.25)", "[0.25,0.5)", "[0.5,0.75)", "[0.75,1]" ),col=cols[1:5],lty=1, bty="n",lwd=2, cex=.8,inset=c(0,-0.045))
legend("topleft",legend=c("[0,0.25)", "[0.25,0.75)", "[0.75,1]" ),col=cols,lty=1, bty="n",lwd=2, cex=.9,inset=c(0,-0.045))
mtext("a" ,at=-5, font=2, cex=4/3)

cols=colorRampPalette(c("skyblue","darkblue"))(length(levels(cuts)))
plot(x=fpr2[[1]],y=tpr2[[1]],col="steelblue",type="n", xlab="False positive rate", ylab="True positive rate")
for(k in 1:length(levels(cuts))) lines(x=fpr2[[k]],y=tpr2[[k]],col=cols[k], lty=k, lwd=2)
legend("bottomright",levels(cuts),col=cols,lty=1:length(levels(cuts)), cex=1.1, lwd=2, bty="n")
mtext("b" ,at=-0.22, font=2, cex=4/3)


hist(zeroId[zeroId==0],xlim=c(0,.95),breaks=seq(0,1,.05),main="",xlab="Posterior probability")
hist(zeroWeightsSim[samp],add=TRUE,breaks=seq(0,1,.05),col=rgb(0.1,0.8,0.1,.2))
legend("topright",c("nr. true excess zeros","zingeR probabilities"),fill=c(0,rgb(0.1,0.8,0.1,.2)), bty="n", cex=1.25)
mtext("c" ,at=-0.22, font=2, cex=4/3)


plot(x=dWeighted$AveLogCPM,y=sqrt(dWeighted$tagwise.dispersion),pch=16,cex=.2,ylim=c(0,2.5),xlab="Average Log CPM", ylab="BCV",col=alpha("black",1/3))
o <- order(dWeighted$AveLogCPM)
lines(dOriginal$AveLogCPM[o], sqrt(dOriginalTrend$trended.dispersion)[o], col = "red",lwd = 2)
lines(dZeroes$AveLogCPM[o], sqrt(dZeroesTrend$trended.dispersion)[o], col = "blue",lwd = 2)
lines(dWeighted$AveLogCPM[o], sqrt(dWeightedTrend$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, NB simul.","NB model, ZINB simul.","ZINB model, ZINB simul."),lty=1,lwd=2,col=c("red","blue","steelblue1"),bty="n")
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
lines(dW$AveLogCPM[o], sqrt(dW$trended.dispersion)[o], col = "steelblue1",lwd = 2)
legend("topright",c("NB model, scRNA-seq","ZINB model, scRNA-seq"),lty=1,col=c("red","steelblue1"),lwd=2,bty="n")
mtext("h" ,at=-0.5, font=2, cex=4/3)
dev.off()


################### Introduction BCV plot
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v6.R")
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
#d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
#dOriginal=d
#dOriginalTrend=estimateGLMTrendedDisp(d,design)
d=estimateDisp(d,design,prior.df=0)
dNoZero=d
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
#d=estimateGLMTagwiseDisp(estimateGLMCommonDisp(d,design),design,prior.df=0)
#d=estimateGLMTagwiseDisp(estimateGLMTrendedDisp(estimateGLMCommonDisp(d,design),design),design,prior.df=0)
#dZeroes=d
#dZeroesTrend=estimateGLMTrendedDisp(dZeroes,design)
d=estimateDisp(d,design,prior.df=0)
dZero=d
plotBCV(d, ylim=c(0,3))

## downweight introduced zeros
dW=DGEList(data)
dW=calcNormFactors(dW)
zeroId[zeroId==0]=1e-6
zeroId[zeroId==1]=.9999
dW$weights=zeroId
dW=estimateDisp(dW,design,prior.df=0)
plotBCV(dW)

## performance of edgeR or downweighted edgeR on ZI RNA-seq data
edgeRPerf = edgeR.pfun(dataZeroes$counts, group=grp, design)
edgeRWeightedPerf = edgeRWeightedOldF.pfun(dataZeroes$counts, group=grp, design, weights=zeroId)
library(iCOBRA)
truth=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truth[dataZeroes$indDE,"status"]=1
cobra = COBRAData(pval=data.frame(edgeR=edgeRPerf[,"pval"],
				  edgeRWeighted=edgeRWeightedPerf[,"pval"],
				  row.names=rownames(dataZeroes)),
		  padj=data.frame(edgeR=p.adjust(edgeRPerf[,"pval"],"fdr"),
				  edgeRWeighted=p.adjust(edgeRWeightedPerf[,"pval"],"fdr"),
				  row.names=rownames(dataZeroes)),
		  truth=truth)
cobraperf = calculate_performance(cobra, binary_truth="status")
colors = c(edgeR="red", edgeRWeighted="chocolate1")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
hlp = plot_fdrtprcurve(cobraplot, pointsize=2, xaxisrange = c(0, 0.4))$data


plot(x=hlp$FDR[hlp$method=="edgeR"], y=hlp$TPR[hlp$method=="edgeR"], type="l", xlim=c(0,0.4), col="red", lwd=2, xlab="False discovery rate", ylab="True positive rate")
lines(x=hlp$FDR[hlp$method=="edgeRWeighted"], y=hlp$TPR[hlp$method=="edgeRWeighted"], type="l", xlim=c(0,0.4), col="chocolate1", lwd=2, xlab="False discovery rate", ylab="True positive rate")

#cutOffPointsEdgeR = c(which.min(abs(hlp$CUTOFF[hlp$method=="edgeR"] - 0.01)),
#		     which.min(abs(hlp$CUTOFF[hlp$method=="edgeR"] - 0.05)),
#		     which.min(abs(hlp$CUTOFF[hlp$method=="edgeR"] - 0.1)))
#points(x=hlp$FDR[cutOffPointsEdgeR], y=hlp$TPR[cutOffPointsEdgeR], pch=16,col=c("white","red")[(hlp$FDR[cutOffPointsEdgeR]<c(0.01,0.05,0.1))+1])

plotBCVIk = function (y, xlab = "Average log CPM", ylab = "Biological coefficient of variation",
    pch = 16, cex = 0.2, col.common = "red", col.trend = "blue",
    col.tagwise = "black", ...)
{
    ## copied from plotBCV function from edgeR package.
    if (!is(y, "DGEList"))
        stop("y must be a DGEList.")
    A <- y$AveLogCPM
    if (is.null(A))
        A <- aveLogCPM(y$counts, offset = getOffset(y))
    disp <- getDispersion(y)
    if (is.null(disp))
        stop("No dispersions to plot")
    if (attr(disp, "type") == "common")
        disp <- rep(disp, length = length(A))
    plot(A, sqrt(disp), xlab = xlab, ylab = ylab, type = "n",
        ...)
    labels <- cols <- lty <- pt <- NULL
    if (!is.null(y$tagwise.dispersion)) {
        points(A, sqrt(y$tagwise.dispersion), pch = pch, cex = cex,
            col = col.tagwise)
        labels <- c(labels, "Tagwise")
        cols <- c(cols, col.tagwise)
        lty <- c(lty, -1)
        pt <- c(pt, pch)
    }
    if (!is.null(y$common.dispersion)) {
        abline(h = sqrt(y$common.dispersion), col = col.common,
            lwd = 2)
        labels <- c(labels, "Common")
        cols <- c(cols, col.common)
        lty <- c(lty, 1)
        pt <- c(pt, -1)
    }
    if (!is.null(y$trended.dispersion)) {
        o <- order(A)
        lines(A[o], sqrt(y$trended.dispersion)[o], col = col.trend,
            lwd = 2)
        labels <- c(labels, "Trend")
        cols <- c(cols, col.trend)
        lty <- c(lty, 1)
        pt <- c(pt, -1)
    }
    #legend("topright", legend = labels, lty = lty, pch = pt,
     #   pt.cex = cex, lwd = 2, col = cols)
    invisible()
}



### composite plot
#png("~/Dropbox/phdKoen/singleCell/figures/introBCVRNAseq.png",width=10,height=8,units="in",res=100)
pdf("~/Dropbox/phdKoen/singleCell/figures/introBCVRNAseq.pdf",width=10,height=8)
par(mfrow=c(2,2), mar=c(5,5,3,1))
plotBCVIk(dNoZero, col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l")
mtext("a",side=3, at=-4.9, cex=4/3,font=2)
plotBCVIk(dZero, col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l")
mtext("b",side=3, at=-5.1, cex=4/3,font=2)
plotBCVIk(dW, col.common=NULL, col.trend="red", cex.axis=1.5, cex.lab=1.25, bty="l")
mtext("c",side=3, at=-5.1, cex=4/3,font=2)
plot(x=hlp$FDR[hlp$method=="edgeR"], y=hlp$TPR[hlp$method=="edgeR"], type="l", xlim=c(0,0.4), col="red", lwd=2, xlab="False discovery proportion", ylab="True positive rate", cex.axis=1.5, cex.lab=1.5, bty="l")
lines(x=hlp$FDR[hlp$method=="edgeRWeighted"], y=hlp$TPR[hlp$method=="edgeRWeighted"], type="l", xlim=c(0,0.4), col="chocolate1", lwd=2, xlab="False discovery rate", ylab="True positive rate")
legend("bottomright",c("edgeR","weighted edgeR"),bty="n",lty=1,lwd=2,col=c("red","chocolate1"), cex=1.33)
mtext("d",side=3,at=-0.1,cex=4/3,font=2)
dev.off()
