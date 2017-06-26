source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/simulation/simulationHelpFunctions_v6.R")
source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

#########################################################
############# RNA-Seq data simulation ###################
#########################################################
#### no zero inflation simulation
library(Biobase) ; library(genefilter)
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
selectedMethods <- c("DESeq2Zero_adjustedDf_posCountsNormZeroWeights", "edgeREMLibSizeDispFastOldFFilteredEdgeR", "MAST", "limma_voom","edgeR","DESeq2", "DESeq2_poscounts" , "NODES")
group=grp
pvalsNoZI <- pval(dataNoZI, method=selectedMethods, count.type="counts", mc.cores=1, niter=200)
#save(pvalsNoZI,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsNoZI.rda")
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsNoZI.rda")
scdeP=scde.pfun(dataNoZI$counts,grp)


library(iCOBRA)
truthNoZI=data.frame(status=rep(0,nTags), row.names=rownames(dataNoZI))
truthNoZI[dataNoZI$indDE,"status"]=1
cobraNoZI <- COBRAData(pval =data.frame(
					limma_voom=pvalsNoZI$pval$limma_voom,
					edgeR=pvalsNoZI$pval$edgeR,
					zingeR_edgeR=pvalsNoZI$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					DESeq2=pvalsNoZI$pval$DESeq2,
					DESeq2_poscounts=pvalsNoZI$pval$DESeq2_poscounts,
					zingeR_DESeq2=pvalsNoZI$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					MAST=pvalsNoZI$pval$MAST,
					NODES=pvalsNoZI$pval$NODES,
					scde=scdeP[,"pval"],
					row.names = rownames(dataNoZI)),
		   padj = data.frame(
				     limma_voom=pvalsNoZI$padj$limma_voom,
				     edgeR=pvalsNoZI$padj$edgeR,
				     zingeR_edgeR=pvalsNoZI$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
				     DESeq2=pvalsNoZI$padj$DESeq2,
				     DESeq2_poscounts=pvalsNoZI$padj$DESeq2_poscounts,
				     zingeR_DESeq2=pvalsNoZI$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
				     MAST=pvalsNoZI$padj$MAST,
				     NODES=pvalsNoZI$padj$NODES,
				     scde=scdeP[,"padj"],
				     row.names = rownames(dataNoZI)),
                   truth = truthNoZI)
cobraperf <- calculate_performance(cobraNoZI, binary_truth = "status")
colors=c(limma_voom="blue", limma_voomZero="steelblue", edgeRTruth="chocolate1", edgeR="red", zingeR_edgeR="salmon", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", MAST="darkturquoise", metagenomeSeq="green", scde="grey", NODES="black")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
plot_fdrtprcurve(cobraplot, pointsize=2)
#plot_roc(cobraplot,xaxisrange=c(0,0.1), yaxisrange=c(0,.9))

# add zeroes: here we see an obvious gain
dataZeroes = dataNoZI
propZeroes=0.05
zeroId = matrix(1,nrow=nrow(dataNoZI),ncol=ncol(dataNoZI))
set.seed(46)
samp=sample(1:length(zeroId),floor(length(zeroId)*propZeroes))
zeroId[samp]=0
zeroId[dataNoZI$counts==0]=1 #if it already was a zero it is not zero-inflated.
samp=samp[!samp%in%which(dataNoZI$counts==0)] #same
dataZeroes$counts = dataZeroes$counts*zeroId
selectedMethods <- c("DESeq2Zero_adjustedDf_posCountsNormZeroWeights", "edgeREMLibSizeDispFastOldFFilteredEdgeR", "MAST", "limma_voom","edgeR","DESeq2", "DESeq2_poscounts" , "NODES")
pvalsZeroes = pval(dataZeroes, method=selectedMethods, count.type="counts", mc.cores=1, niter=200)
#save(pvalsZeroes,file="~/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsZeroes.rda")
load("~/Dropbox/PhD/Research/zeroInflation/singleCell/pvalsZeroes.rda")
## edgeR with ground truth
weightHlp = zeroId
weightHlp[weightHlp==0]=1e-6 #does not really allow zero weights.
pvalEdgeRGroundTruth <- edgeRWeightedOldF.pfun(counts=dataZeroes$counts, group=grp, weights=weightHlp)
## DESeq2 with ground truth
pvalDESeq2GroundTruth <- DESeq2_weightedT.pfun(counts=dataZeroes$counts, group=grp, weights=zeroId)
scdePZero=scde.pfun(dataZeroes$counts,grp)

##performance curves
#rnaSeqPerformanceWithZeroes.pdf
truthZero=data.frame(status=rep(0,nTags), row.names=rownames(dataZeroes))
truthZero[dataZeroes$indDE,"status"]=1
cobraZero <- COBRAData(pval =data.frame(
					limma_voom=pvalsZeroes$pval$limma_voom,
					DESeq2=pvalsZeroes$pval$DESeq2,
					DESeq2_poscounts=pvalsZeroes$pval$DESeq2_poscounts,
					DESeq2Truth=pvalDESeq2GroundTruth[,"pval"],
					zingeR_DESeq2=pvalsZeroes$pval$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
					edgeR=pvalsZeroes$pval$edgeR,
					zingeR_edgeR=pvalsZeroes$pval$edgeREMLibSizeDispFastOldFFilteredEdgeR,
					edgeRTruth = pvalEdgeRGroundTruth[,1],
					MAST=pvalsZeroes$pval$MAST,
					NODES=pvalsZeroes$pval$NODES,
					scde=scdePZero[,"pval"],
					row.names = rownames(dataZeroes)),
		   padj = data.frame(
				     limma_voom=pvalsZeroes$padj$limma_voom,
				     DESeq2=pvalsZeroes$padj$DESeq2,
				     DESeq2_poscounts=pvalsZeroes$padj$DESeq2_poscounts,
				     DESeq2Truth=pvalDESeq2GroundTruth[,"padj"],
				     zingeR_DESeq2=pvalsZeroes$padj$DESeq2Zero_adjustedDf_posCountsNormZeroWeights,
				     edgeR=pvalsZeroes$padj$edgeR,
				     zingeR_edgeR=pvalsZeroes$padj$edgeREMLibSizeDispFastOldFFilteredEdgeR,
				     edgeRTruth = pvalEdgeRGroundTruth[,2],
				     MAST=pvalsZeroes$padj$MAST,
				     NODES=pvalsZeroes$padj$NODES,
				     scde=scdePZero[,"padj"],
				     row.names = rownames(dataZeroes)),
                   truth = truthZero)
cobraperf <- calculate_performance(cobraZero, binary_truth = "status")
colors=c(limma_voom="blue", edgeRTruth="chocolate1", edgeR="red", zingeR_edgeR="salmon", DESeq2="brown", DESeq2_poscounts="navajowhite2", zingeR_DESeq2="darkseagreen", DESeq2Truth="forestgreen", MAST="darkturquoise", metagenomeSeq="green", scde="grey", NODES="black")
colsCobra=colors[match(sort(names(cobraperf@overlap)[1:(ncol(cobraperf@overlap)-1)]),names(colors))]
cobraplotZero <- prepare_data_for_plot(cobraperf, colorscheme=colsCobra)
plot_fdrtprcurve(cobraplotZero, pointsize=2)


## two-panel plot
library(cowplot)
prow <- plot_grid( plot_fdrtprcurve(cobraplotZero, pointsize=2) + theme(legend.position="none") + xlab("FDP"),
           plot_fdrtprcurve(cobraplot, pointsize=2) + xlab("FDP") + theme(legend.position="none"),
           align = 'vh',
           labels = c("a", "b"),
           hjust = -1,
           nrow = 1
           )
legend_b <- get_legend(plot_fdrtprcurve(cobraplotZero, pointsize=2) + theme(legend.position="bottom"))
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
png("~/Dropbox/phdKoen/singleCell/figures/supplementary/RNASeq_composite.png", width=7,height=8, units="in", res=300)
p
dev.off()


## plot a histogram of the weights for the introduced zeroes
weights=zeroWeightsLibSizeDispFast(counts=dataZeroes$counts, design=model.matrix(~grp), maxit=100, plotW=TRUE)
hist(weights[samp[dataNoZI$counts[samp]!=0]], xlab="weight", xlim=c(0,1), main="")


## ROC curve for identifying excess zeros
weights=zeroWeightsLibSizeDispFast(counts=dataZeroes$counts, design=model.matrix(~grp), maxit=200, plotW=TRUE)
hist(weights[samp], xlab="weight", xlim=c(0,1), main="")

png("~/Dropbox/phdKoen/singleCell/figures/supplementary/rocExcessZerosRNASeqSimulation_composite.png", width=8,height=6, units="in", res=300)
pvalSeq = c(1e-15,1e-14,1e-13,1e-12,1e-10,1e-9,1e-8,1e-7,1e-6,seq(.00001,.005,by=.00001),seq(.005,1,by=.005))
falses=which(dataNoZI$counts==0)
tpr=fpr=vector(length=length(pvalSeq))
for(i in 1:length(pvalSeq)){
    excessID <- which(weights<=pvalSeq[i])
    tpr[i] <- mean(samp%in%excessID)
    fpr[i] <- mean(falses%in%excessID)
}
par(mfrow=c(1,2))
plot(x=fpr,y=tpr,type="l", xlab="False positive rate", ylab="True positive rate", lwd=2, col="steelblue")
points(x=fpr[pvalSeq==1e-5],y=tpr[pvalSeq==1e-5],col=2,pch=19)
points(x=fpr[pvalSeq==1e-2],y=tpr[pvalSeq==1e-2],col=2,pch=19)
points(x=fpr[519],y=tpr[519],col=2,pch=19) #w=0.05
text(x=fpr[pvalSeq==1e-5]+0.13,y=tpr[pvalSeq==1e-5],labels="w=1e-5")
text(x=fpr[pvalSeq==1e-2]+0.13,y=tpr[pvalSeq==1e-2],labels="w=1e-2")
text(x=fpr[519]+0.13,y=tpr[519],labels="w=5e-2")

## ROC curve stratified by average expression
d=DGEList(dataZeroes$counts)
d=calcNormFactors(d)
acpm=aveLogCPM(d)
library(Hmisc)
cuts=cut2(acpm,g=5)
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
for(k in 1:length(levels(cuts))) lines(x=fpr2[[k]],y=tpr2[[k]],col=cols[k], lty=k)
legend("bottomright",levels(cuts),col=cols,lty=1:length(levels(cuts)))
dev.off()















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
