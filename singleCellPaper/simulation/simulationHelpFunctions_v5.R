library(edgeR)
library(gamlss.dist)
library(mgcv)
library(gamlss)
library(gamlss.tr)

gen.trun(par=0, family="NBI", name="ZeroTruncated", type="left", varying=FALSE)
getParamsZTNB <- function(counts, offset, design=NULL) {  
    require(MASS)
    libSize=offset
    #fit a ZTNB model only on positive counts part
    countsModel = counts[counts>0]
    if(length(countsModel)<2) stop("Need at least two positive counts")
    libSizeModel = libSize[counts>0]
    if(is.null(design)){designFit=matrix(1,nrow=length(countsModel),ncol=1)} else {designFit=design[counts>0,]}
    fit=try(gamlss(formula=countsModel~-1+designFit+offset(log(libSizeModel)), family="NBIZeroTruncated", control=gamlss.control(trace=FALSE, n.cyc=300)),silent=TRUE)

    if(class(fit)[1]=="try-error") return(c(dispersion=NA, lambda=NA))
    #lambda=exp(mean(fit$mu.coefficients)) #geometric mean
    lambda=mean(countsModel/libSizeModel)
    dispersion=exp(fit$sigma.coefficients)
    return(c(dispersion=dispersion,lambda=lambda))
}

getDatasetZTNB = function(counts, design, drop.extreme.dispersion = FALSE, cpm= "AveLogCPM"){
    
        #### estimate lambda and overdispersion based on ZTNB.
	d <- DGEList(counts)
	cp <- cpm(d,normalized.lib.sizes=TRUE)
	dFiltered=d
	dFiltered <- edgeR::calcNormFactors(dFiltered)	
        dFiltered$AveLogCPM <- aveLogCPM(dFiltered)
	params=t(apply(dFiltered$counts,1,function(x) getParamsZTNB(counts=x,offset=dFiltered$samples$lib.size, design=design)))
	rmRows = which(params[,2]>1) #impossibly high lambda
	rmRows2 = which(params[,2]==0) #zero lambda
	naRows = which(apply(params,1, function(row) any(is.na(row)))) #not fitted
	nonZeroDispRows = which(params[,1]<0 | params[,1]==0) #negative dispersion
	throwRows = c(rmRows,rmRows2,naRows,nonZeroDispRows)
	params = params[-throwRows,]
	
	### estimate logistic GAM P(zero) ~ s(aveLogCPM) + logLibSize
	### use unfiltered data for this model.
	propZero = colMeans(counts==0)
	propZeroGene = rowMeans(counts==0)	
	d <- DGEList(counts)
	d <- edgeR::calcNormFactors(d)
	if(cpm=="AveLogCPM"){ avCpm <- aveLogCPM(d)} else if(cpm=="aCpm"){ avCpm <- aCPM(d$counts)} else {stop("cpm must be either AveLogCPM or aCPM")}
	cpmHist = hist(avCpm, breaks=150, plot=FALSE)
    	breaks = cpmHist$breaks
    	mids = cpmHist$mids
    	midsHlp=rep(mids,ncol(d$counts))
	logLibSize = log(colSums(counts))	
    	logLibHlp=rep(logLibSize,each=length(mids))
	binHlp=sapply(breaks[-length(breaks)],function(x) avCpm>x)
  	binId=apply(binHlp,1,function(x) max(which(x)))
	nonNullCounts = t(sapply(1:length(mids), function(bin){
			    binRows <- binId==bin
			    if(sum(binRows)==0) rep(0,ncol(counts)) else
			    if(sum(binRows)==1) (counts[which(binRows),]!=0)*1 else
				colSums(counts[which(binRows),]!=0) 
	    }))
	nullCounts = t(sapply(1:length(mids), function(bin){
		    	binRows <- binId==bin
		    	if(sum(binRows)==0) rep(0,ncol(counts)) else
		    	if(sum(binRows)==1) (counts[which(binRows),]==0)*1 else 
			    colSums(counts[which(binRows),]==0)  	   
	    }))
	expectCounts=cbind(c(nullCounts),c(nonNullCounts))
	#zeroFit=mgcv::gam(expectCounts~s(midsHlp)+logLibHlp,family=binomial)
	zeroFit=mgcv::gam(expectCounts~s(midsHlp,by=logLibHlp),family=binomial)

	### drop extreme dispersions
        dFiltered$AveLogCPM <- aveLogCPM(dFiltered)	
	dFiltered$AveLogCPM <- dFiltered$AveLogCPM[-throwRows]
	propZeroGene = propZeroGene[-throwRows]
	params=data.frame(dispersion=params[,1], lambda=params[,2], aveLogCpm=dFiltered$AveLogCPM, propZeroGene=propZeroGene)
	dispersion <- params$dispersion
	AveLogCPM <- params$aveLogCpm
	lambda <- params$lambda
	propZeroGene <- params$propZeroGene
	
	if(is.numeric(drop.extreme.dispersion))
	{   
		bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE, na.rm=TRUE)
		ids <- dispersion <= bad
		AveLogCPM <- AveLogCPM[ids]
		dispersion <- dispersion[ids]
		lambda <- lambda[ids]
		propZeroGene <- propZeroGene[ids]
		params <- params[ids,]
		dFiltered <- dFiltered[ids,]
	}
	#lambda=lambda/sum(lambda) #make sure they sum to 1
	dataset.AveLogCPM <- AveLogCPM
	dataset.dispersion <- dispersion
	dataset.lambda <- lambda
	dataset.propZeroGene <- propZeroGene
	dataset.lib.size <- d$samples$lib.size
	dataset.nTags <- nrow(d) 
	list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags, dataset.propZeroFit=zeroFit, dataset.lambda=lambda, dataset.propZeroGene=propZeroGene, dataset.breaks = breaks, dataset.cpm=cpm)
}

NBsimSingleCell <- function(dataset, group, nTags = 10000, nlibs = length(group), lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1, pUp=.5, foldDiff=3, verbose=TRUE, ind=NULL, params=NULL, noiseCell=0, noiseGene=0, randomZero=0.025, cpm=c("aCpm","AveLogCPM"), max.dispersion=400, min.dispersion=0.1)
{
	require(edgeR)
	group = as.factor(group)
	expit=function(x) exp(x)/(1+exp(x))
	logit=function(x) log(x/(1-x))

	sample.fun <- function(object)
	{
		nlibs <- object$nlibs
		nTags <- object$nTags
		AveLogCPM <-object$dataset$dataset.AveLogCPM
		dispersion <- object$dataset$dataset.dispersion
		lambda <- object$dataset$dataset.lambda
		#lambda <- (2^AveLogCPM)/1e6
		propZeroGene <- dat$dataset$dataset.propZeroGene
                id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
		object$AveLogCPM <- AveLogCPM[id_r]
		Lambda <- lambda[id_r]
		#Lambda <- Lambda/sum(Lambda) #normalize so they all sum to 1
		Dispersion <- dispersion[id_r]
		Dispersion[Dispersion>max.dispersion] = max.dispersion
		Dispersion[Dispersion<min.dispersion] = min.dispersion
		propZeroGene <- propZeroGene[id_r]	
		Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
		object$Lambda <- Lambda
		Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
		object$Dispersion <- Dispersion
		object$propZeroGene <- propZeroGene
		object
	}
	diff.fun <- function(object)
	{
		group <- object$group
		pUp <-  object$pUp 
		foldDiff <- object$foldDiff
		Lambda <- object$Lambda
		nTags <- object$nTags
		g <- group == levels(group)[1]
		if(length(ind)>0 & !all(foldDiff==1)) {
			fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
			Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
			Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir)) 
			object$Lambda <- Lambda
			object$indDE <- ind
			object$indNonDE <- (1:nTags)[-ind]
			foldDiff[fcDir==1] <- 1/foldDiff[fcDir==1]
			object$foldDiff <- foldDiff #group2 / group1
			}
		if(all(foldDiff==1)) object$indDE <- NA
		object
	}
	sim.fun <- function(object)
	{   
		Lambda <- object$Lambda
		Dispersion <- object$Dispersion
		nTags <- object$nTags
		nlibs <- object$nlibs
		lib.size <- object$lib.size
		zeroFit <- dat$dataset$dataset.propZeroFit
		propZeroGene <- dat$propZeroGene
		propZeroGene[propZeroGene==1] <- 1-1e-4
		propZeroGene[propZeroGene==0] <- 1e-4
		design <- object$design
		avLogCpm <- object$AveLogCPM
		mids <- object$dataset$dataset.mids
		breaks <- object$dataset$dataset.breaks
		## get matrix of zero probabilities
		libPredict=rep(log(lib.size),each=length(avLogCpm))
		cpmPredict=rep(avLogCpm,length(lib.size))
		## no noise
		if((noiseCell+noiseGene)==0){
		    #zeroProbMat = matrix(predict(zeroFit, newdata=data.frame(logLibHlp=libPredict, midsHlp=cpmPredict), type="response"), byrow=FALSE, ncol=nlibs)
		    zeroProbMatLink = matrix(predict(zeroFit, newdata=data.frame(logLibHlp=libPredict, midsHlp=cpmPredict), type="link"), byrow=FALSE, ncol=nlibs)
		    meanDiff = rowMeans(zeroProbMatLink)-logit(propZeroGene)
		    zeroProbMat = expit(sweep(zeroProbMatLink,1,meanDiff,"-"))
		    #zeroProbMat=expit(zeroProbMatLink) #no empirical adjustment
		} else if((noiseCell+noiseGene)>0){
		    ## TO DO: code for adding noise is not up to date
		    zeroProbMat = matrix(predict(zeroFit, newdata=data.frame(logLibHlp=libPredict, midsHlp=cpmPredict), type="link"), byrow=FALSE, ncol=nlibs)
		    noiseCol=rnorm(n=ncol(zeroProbMat),mean=0,sd=noiseCell)
		    noise=rnorm(n=nrow(zeroProbMat),sd=noiseGene)
		    zeroProbMat = expit(sweep(zeroProbMat+noise,2,FUN="+",STATS=noiseCol))
		    #zeroProbMat = expit(sweep(zeroProbMat,2,FUN="+",STATS=noiseCol))
		}
		## introduce random zeroes for lower count genes
		zeroProbMat[sample(1:length(zeroProbMat), floor(randomZero*length(zeroProbMat)))]=1-1e-5		
		#lowCountGenesID <- which(rowMeans(zeroProbMat)>0.15)
		#hlp=zeroProbMat[lowCountGenesID,]
		#hlp[sample(1:length(hlp), floor(randomZero*length(hlp)))]=1-1e-5
		#zeroProbMat[lowCountGenesID,]=hlp

		## adjust lib size for adding zeroes by calculating expected loss
		#avCount=sweep(Lambda,2,lib.size,"*")
		#expectedLoss=avCount*zeroProbMat
		#libSizeCountSim = lib.size + colSums(expectedLoss)
		
		## adjust mu for adding zeroes
		mu=sweep(Lambda,2,lib.size,"*")
		mu[mu<1e-16] = 1
		adjustment = zeroProbMat*mu
		mu=mu+adjustment
		
		## simulate counts acc to a zero-adjusted NB model
		#mu=sweep(Lambda,2,libSizeCountSim,"*")
		#mu[mu<1e-16] = 0.5
		counts = rZANBI(n=nTags*nlibs, mu=mu, sigma=Dispersion, nu=zeroProbMat)

		## the rZANBI function rarely simulates Inf values for very low mu estimates. Resimulate for these genes using same params, if present	
		## also, resample features with all zero counts
		zeroCountsId <- which(rowSums(counts)==0)
		infId <- which(apply(counts,1,function(row) any(is.infinite(row))))
		while(length(zeroCountsId)>0 | length(infId)>0){
		    if(length(zeroCountsId)>0){ #resimulate all zero gene
		    	counts[zeroCountsId,] = rZANBI(n=length(zeroCountsId)*nlibs, mu=mu[zeroCountsId,], sigma=Dispersion[zeroCountsId,], nu=zeroProbMat[zeroCountsId,])
		    }
		    if(length(infId)>0){ #resimulate Inf
			counts[infId,] <- rZANBI(n=length(infId)*nlibs, mu=mu[infId,], sigma=Dispersion[infId,], nu=zeroProbMat[infId,])
		    }
		    zeroCountsId <- which(rowSums(counts)==0)
		    infId <- which(apply(counts,1,function(row) any(is.infinite(row))))
		}
		
		rownames(counts) <- paste("ids", 1:nTags, sep = "")
		object$counts <- counts
		object
	}
			
        if(verbose) message("Preparing dataset.\n")	
	  if(is.null(params)){ 
	      dataset <- getDatasetZTNB(counts = dataset, drop.extreme.dispersion = drop.extreme.dispersion, drop.low.lambda = drop.low.lambda)
	  } else {
	  dataset <- params
	  }
	  dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pUp = pUp, foldDiff = foldDiff))
	if(cpm=="aCpm") dat$dataset$dataset.AveLogCPM = dat$dataset$dataset.aCpm
	

	if(is.null(dat$lib.size)){
	  dat$lib.size <- sample(dataset$dataset.lib.size, nlibs, replace=TRUE)}
	if(is.null(nTags)) dat$nTags <- dat$dataset$dataset.nTags 
        if(verbose) message("Sampling.\n")	
	      dat <- sample.fun(dat)
        if(verbose) message("Calculating differential expression.\n")	
	      dat <- diff.fun(dat)
        if(verbose) message("Simulating data.\n")	
	      dat <- sim.fun(dat)
	dat 	
}


NBsimFold <-
function(fold_seq = seq(1, 3, by = 0.25), dataset, group, nTags = 10000, nlibs = length(group), fix.dispersion = NA, lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1,  add.outlier = FALSE, outlierMech = c("S", "R", "M"), pOutlier = 0.1, min.factor = 1.5, max.factor = 8, pDiff=.1, pUp=.5, name = NULL)
{  
    ## NBsimFold generates a series of simulated counts followed by the NB model ##
	o <- order(fold_seq)
	fold_seq <- fold_seq[o]
	fold_seq <- unique(round(fold_seq, 4))
	if(is.null(name)) name <- deparse(substitute(dataset))
	if(is.matrix(dataset)) dataset <- getDataset(counts =dataset, drop.extreme.dispersion = drop.extreme.dispersion)
	if(all(fold_seq[1] == fold_seq) ) names(fold_seq) <- paste("sample_", seq(fold_seq), sep = "" )
	else names(fold_seq) <- paste("fold_", fold_seq, sep = "" )
	re_split <- function(object)
    {out <- list()
		for (what in names(object[[1L]])) {
			out[[what]] <- lapply(object, .subset2, what)
			if(identical(out[[what]][[1L]], out[[what]][[2L]]))
			out[[what]] <- out[[what]][[1L]]}
		out}
	
	rere_split <- function(object, x)
	{ object[[x]] <- re_split(object[[x]])
	  object}
		
	
	simFold <- lapply(fold_seq, function(x) NBsim(foldDiff = x, dataset = dataset, group = group, nTags = nTags, nlibs = nlibs, fix.dispersion = fix.dispersion, lib.size = lib.size, drop.low.lambda = drop.low.lambda, drop.extreme.dispersion = drop.extreme.dispersion,  add.outlier = add.outlier, outlierMech = outlierMech, pOutlier = pOutlier, min.factor = min.factor, max.factor = max.factor, pDiff = pDiff, pUp=pUp, name = name))
	
    simFold <- re_split(simFold)
    for(i in seq_along(outlierMech))
       simFold <- rere_split(simFold, outlierMech[i])
	simFold$fold_seq <- names(fold_seq)
    simFold <- new("FoldList", simFold)
    simFold
}


getVersion <-
function(x)
{
  ## it is low-level function of pval ##
  x1 <- gsub("(\\_)(\\w+)", "", x)	
  v <- unlist(lapply(x1, function(z) {options(warn = -1)
                                      desp<- packageDescription(z)
                                      if(length(desp) == 1)
                                        return("unknown")
                                      else desp$Version
                                      }))
  paste0(x,"_", v )	
}	

rmVersion <-
function(x)
{
 ## it is low-level function of pval ##
 x1 <- strsplit(x, "\\_")
 x1 <- lapply(x1, function(x) x[-length(x)])	
 sapply(x1, paste0, collapse = "_")	
	
}

odd <- function(x) 
{   
	## it is low-level function of pval ##
	y <- seq(x)
	idx <- y %% 2 != 0
	x[idx]
}

mainShow <-
function(count.type, count.name, group, pOutlier)
{ 
	## it is low-level function of pval ##
	pOutlier <- paste(100*pOutlier, "% ", "outliers", sep = "")
	group <- as.factor(group)
	group <- paste0(sum(group == levels(group)[1]), "vs", sum(group == levels(group)[2]))
	if(count.type == "counts")
	paste0("No outliers", "/", count.name, "/", group)
	else
	paste0(pOutlier, "/", count.type, "/",  count.name, "/", group)
	
}


resetPar <- function() {
    ## this re-set args of par for plot ##
	dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}


DESeq2.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  dse = estimateSizeFactors(dse)
	  dse = estimateDispersions(dse)
	  dse = nbinomWaldTest(dse, betaPrior=TRUE)
	  res <- results(dse)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
          out
}

DESeq2_poscounts.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  dse <- estimateSizeFactors(dse,type="poscounts")
	  dse <- estimateDispersions(dse)
	  dse <- nbinomWaldTest(dse, betaPrior=TRUE)
	  res <- results(dse)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
      out
}

DESeq2_poscounts_DESeq.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  dse <- estimateSizeFactors(dse,type="poscounts")
	  dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
      out
}

DESeq2_poscounts_noShrink.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  dse <- estimateSizeFactors(dse,type="poscounts")
	  dse <- estimateDispersions(dse)
	  dse <- nbinomWaldTest(dse, betaPrior=FALSE)
	  res <- results(dse)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
      out
}

DESeq2_poscounts_noFiltering.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  dse <- estimateSizeFactors(dse,type="poscounts")
	  dse <- estimateDispersions(dse)
	  dse <- nbinomWaldTest(dse, betaPrior=TRUE)
	  res <- results(dse, independentFiltering=FALSE)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
      out
}

DESeq2_poscounts_noImputation.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  dse <- estimateSizeFactors(dse,type="poscounts")
	  dse <- estimateDispersions(dse)
	  dse <- nbinomWaldTest(dse, betaPrior=TRUE)
	  res <- results(dse, minReplicatesForReplace=Inf)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
      out
}


DESeq2Zero_adjustedDf_DEseq2normZeroWeights.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2) ; library(genefilter)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~group), plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2")
	  dimnames(zeroWeights) = NULL
	  assays(dse)[["weights"]] = zeroWeights
	  #dse = estimateSizeFactors(dse)
	  #dse = estimateDispersions(dse)
	  #dse=nbinomWaldTest(dse)
	  dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse)
	  baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
	  pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
	  padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
	  out <- cbind(pval = pvalDesZero, padj = padjusted$padj, lfc = res$log2FoldChange)
	  #resDESeq2Zero=out
      out
}


DESeq2Zero_adjustedDf_posCountsNormZeroWeights.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2) ; library(genefilter)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~group), plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos")
	  dimnames(zeroWeights) = NULL
	  assays(dse)[["weights"]] = zeroWeights
	  dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
	  dse = estimateDispersions(dse)
	  dse=nbinomWaldTest(dse, betaPrior=TRUE)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse)
	  baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
	  pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
	  padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
	  out <- cbind(pval = pvalDesZero, padj = padjusted$padj, lfc = res$log2FoldChange)
      out
}

DESeq2Zero_wald_posCountsNormZeroWeights.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2) ; library(genefilter)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~group), plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos")
	  dimnames(zeroWeights) = NULL
	  assays(dse)[["weights"]] = zeroWeights
	  dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
	  dse = estimateDispersions(dse)
	  dse=nbinomWaldTest(dse, betaPrior=TRUE)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse)
	  #baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
	  #pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
	  #padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
      out
}

DESeq2_weightedT.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL, weights)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  dimnames(weights) = NULL
	  assays(dse)[["weights"]] = weights
	  dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse)
	  baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
	  pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(weights)-2))
	  padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
	  out <- cbind(pval = pvalDesZero, padj = padjusted$padj, lfc = res$log2FoldChange)
          out
}

DESeq2_noShrink.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  dse <- DESeq(dse, betaPrior=FALSE)
	  res <- results(dse)
	  out <- cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
          out
}

DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noShrink.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2) ; library(genefilter)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~group), plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos")
	  dimnames(zeroWeights) = NULL
	  assays(dse)[["weights"]] = zeroWeights
	  dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
	  dse = estimateDispersions(dse)
	  dse=nbinomWaldTest(dse, betaPrior=FALSE)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse)
	  baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
	  pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
	  padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
	  out <- cbind(pval = pvalDesZero, padj = padjusted$padj, lfc = res$log2FoldChange)
	  #resDESeq2Zero=out
      out
}

DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noFiltering.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
	  library(DESeq2) ; library(genefilter)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~group), plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos")
	  dimnames(zeroWeights) = NULL
	  assays(dse)[["weights"]] = zeroWeights
	  dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
	  dse = estimateDispersions(dse)
	  dse=nbinomWaldTest(dse, betaPrior=TRUE)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse, independentFiltering=FALSE)
	  baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
	  pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
	  padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
	  out <- cbind(pval = pvalDesZero, padj = padjusted$padj, lfc = res$log2FoldChange)
	  #resDESeq2Zero=out
      out
}

DESeq2Zero_adjustedDf_posCountsNormZeroWeights_noImputation.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
  {   
      ## implement DESeq2 ##
      ###### UNFINISHED ########
	  library(DESeq2) ; library(genefilter)
	  colData <- data.frame(group)
	  dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	  colData(dse)$group <- as.factor(colData(dse)$group)
	  zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~group), plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos")
	  dimnames(zeroWeights) = NULL
	  assays(dse)[["weights"]] = zeroWeights
	  dse = DESeq2::estimateSizeFactors(dse, type = "poscounts")
	  dse = estimateDispersions(dse)
	  dse=nbinomWaldTest(dse, betaPrior=TRUE)
	  #dse <- DESeq(dse, betaPrior=TRUE)
	  res <- results(dse, minReplicatesForReplace=Inf)
	  baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
	  pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
	  padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
	  out <- cbind(pval = pvalDesZero, padj = padjusted$padj, lfc = res$log2FoldChange)
	  #resDESeq2Zero=out
      out
}


DESeq2_noCook.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
{   
	## implement DESeq2 turning off cooksCutoff ##
	library(DESeq2)
	colData <- data.frame(group)
	dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	colData(dse)$group <- as.factor(colData(dse)$group)
	dse <- DESeq(dse, betaPrior=TRUE)
	res <- results(dse, cooksCutoff = Inf)
	cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
}

DESeq2_noImputation.pfun <-
function(counts, group, design = NULL, mc.cores = 4, niter=NULL)
{   
	## implement DESeq2 turning off imputation for outliers ##
	library(DESeq2)
	colData <- data.frame(group)
	dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	colData(dse)$group <- as.factor(colData(dse)$group)
	dse <- DESeq(dse,minReplicatesForReplace = Inf, betaPrior=TRUE)
	res <- results(dse)
	cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
}

DESeq2_noFiltering.pfun <- function(counts, group, design = NULL, mc.cores = 4, niter=NULL){   
	## implement DESeq2 turning off cooksCutoff ##
	library(DESeq2)
	colData <- data.frame(group)
	dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	colData(dse)$group <- as.factor(colData(dse)$group)
	dse <- DESeq(dse, betaPrior=TRUE)
	res <- results(dse, independentFiltering=FALSE)
	cbind(pval = res$pvalue, padj = res$padj, lfc = res$log2FoldChange)
}

edgeR.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
{
    ## edgeR standard pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)
	d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
	d <- estimateGLMTrendedDisp(d,design=design)
	d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
	f <- glmFit(d, design = design)
	lr <- glmLRT(f, coef=2)
	lfc <- lr$table$logFC
	pval = lr$table$PValue
	padj = p.adjust(pval, "BH")
	out = cbind(pval = pval, padj = padj, lfc = lfc)
	return(out)
}

edgeRHurdle.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
{
    ## edgeR standard pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	zeroId=counts==0
	weights=1-zeroId
	d$weights=weights
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)
	d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
	d <- estimateGLMTrendedDisp(d,design=design)
	d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
	f <- glmFit(d, design = design)
	lr <- glmLRT(f, coef=2)
	lfc <- lr$table$logFC
	pval = lr$table$PValue
	padj = p.adjust(pval, "BH")
	out = cbind(pval = pval, padj = padj, lfc = lfc)
	return(out)
}



edgeRFiltered.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
{
    ## edgeR standard pipeline ##
	library(edgeR) ; library(genefilter)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)
	d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
	d <- estimateGLMTrendedDisp(d,design=design)
	d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
	f <- glmFit(d, design = design)
	lr <- glmLRT(f, coef=2)
	lfc <- lr$table$logFC
	pval = lr$table$PValue
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	out = cbind(pval = pval, padj = padj, lfc = lfc)
	return(out)
}

edgeRWeightedOldF.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL, weights=matrix(1,nrow=nrow(counts),ncol=ncol(counts)))
{
    ## edgeR standard pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)
	d$weights <- weights
	d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
	d <- estimateGLMTrendedDisp(d,design=design)
	d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	lr <- glmLRTOld(edger.fit,coef=2,test="F")
	pval = lr$table$PValue
	padj = p.adjust(pval, "BH")
	out = cbind(pval = pval, padj = padj)
	out[is.na(out)]=1
	return(out)
}

edgeROldF.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
{
    ## edgeR standard pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)
	d <- estimateGLMCommonDisp(d,design=design, interval=c(0,10))
	d <- estimateGLMTrendedDisp(d,design=design)
	d <- estimateGLMTagwiseDisp(d, design = design, prior.df = prior.df)
	edger.fit <- glmFit(d, design) #uses weights
  	lr <- glmLRTOld(edger.fit,coef=2,test="F", ZI=FALSE)
	pval = lr$table$PValue
	padj = p.adjust(pval, "BH")
	lfc <- lr$table$logFC
	out = cbind(pval = pval, padj = padj, lfc=lfc)
	out[is.na(out)]=1
	return(out)
}



edgeREstDisp.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
{
    ## edgeR standard pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)
	d <- estimateDisp(d,design)
	f <- glmFit(d, design = design)
	lr <- glmLRT(f, coef=2)
	pval = lr$table$PValue
	padj = p.adjust(pval, "BH")
	out = cbind(pval = pval, padj = padj)
	out[is.na(out)]=1
	return(out)
}


edgeR_robust.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10, niter=NULL)
{   
    ## edgeR-robsut pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)	
	dw <- estimateGLMRobustDisp(d,design=design, prior.df=prior.df, maxit = 6)
	fw <- glmFit(dw, design=design)
	lrw <- glmLRT(fw,coef=2)
	lfc <- lrw$table$logFC
   	pval = lrw$table$PValue
	padj = p.adjust(pval, "BH")
	cbind(pval = pval, padj = padj, lfc = lfc)
}

edgeR_rdev.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10)
{   
    ## edgeR-robust by deviance residual ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
 	design = model.matrix(~group)	
	dw <- estimateGLMRobustDisp(d,design=design, prior.df=prior.df, maxit = 6, residual.type = "deviance")
	fw <- glmFit(dw, design=design)
	lrw <- glmLRT(fw,coef=2)
   	pval = lrw$table$PValue
	padj = p.adjust(pval, "BH")
	cbind(pval = pval, padj = padj)
}

edgeR_rans.pfun <-
function(counts, group, design = NULL, mc.cores = 4, prior.df=10)
{   
	## edgeR-robust by anscombe residual ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
	dw <- estimateGLMRobustDisp(d,design=design, prior.df=prior.df, maxit = 6, residual.type = "anscombe")
	fw <- glmFit(dw, design=design)
	lrw <- glmLRT(fw,coef=2)
   	pval = lrw$table$PValue
	padj = p.adjust(pval, "BH")
	cbind(pval = pval, padj = padj)
}

limma_voom.pfun <-
function(counts, group, design = NULL, mc.cores = 2, niter=NULL) 
{   
	## limma voom pipeline ##
	library(limma)
 	design = model.matrix(~group)
	nf <- edgeR::calcNormFactors(counts)
	y <- voom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf)
	fit <- lmFit(y, design)
	fit <- eBayes(fit)
	tt <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")
	pval <- tt$P.Value
	padj <- tt$adj.P.Val
	lfc <- tt$logFC
	cbind(pval = pval, padj = padj, lfc=lfc)
}

limma_voomHurdle.pfun <-
function(counts, group, design = NULL, mc.cores = 2, niter=NULL) 
{   
	## limma voom pipeline ##
	library(limma)
	zeroId=counts==0
 	design = model.matrix(~group)
	nf <- edgeR::calcNormFactors(counts)
	y <- voom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf, weights=1-zeroId)
	y$weights=1-zeroId
	fit <- lmFit(y, design)
	fit <- eBayes(fit)
	tt <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")
	pval <- tt$P.Value
	padj <- tt$adj.P.Val
	lfc <- tt$logFC
	cbind(pval = pval, padj = padj, lfc=lfc)
}

limma_voomHurdleHeteroscedastic.pfun <-
function(counts, group, design = NULL, mc.cores = 2, niter=NULL) 
{   
	## limma voom pipeline ##
	library(limma)
	zeroId=counts==0
 	design = model.matrix(~group)
	nf <- edgeR::calcNormFactors(counts)
	y <- voom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf, weights=1-zeroId)
	y$weights=(1-zeroId)*y$weights
	fit <- lmFit(y, design)
	fit <- eBayes(fit)
	tt <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")
	pval <- tt$P.Value
	padj <- tt$adj.P.Val
	lfc <- tt$logFC
	cbind(pval = pval, padj = padj, lfc=lfc)
}


limma_voomFiltered.pfun <-
function(counts, group, design = NULL, mc.cores = 2, niter=NULL) 
{   
	## limma voom pipeline ##
	library(limma) ; library(genefilter)
 	design = model.matrix(~group)
	nf <- edgeR::calcNormFactors(counts)
	y <- voom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf)
	fit <- lmFit(y, design)
	fit <- eBayes(fit)
	tt <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")
	pval <- tt$P.Value
	baseMean = unname(rowMeans(sweep(counts,2,nf,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	lfc <- tt$logFC
	cbind(pval = pval, padj = padj, lfc=lfc)
}



limma_voomZeroFiltered.pfun <-
function(counts, group, design = NULL, mc.cores = 2, niter=100) 
{   
	## limma voom pipeline ##
	library(limma) ; library(genefilter)
 	design = model.matrix(~group)
	nf <- edgeR::calcNormFactors(counts)
	zeroWeights <- zeroWeightsLibSizeDispFast(counts=counts, design=design, maxit=niter)
	y <- voom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf, weights=zeroWeights)
	y$weights=y$weights*zeroWeights
	fit <- lmFit(y, design, weights=y$weights)
	fit$df.residual=rowSums(zeroWeights)-ncol(design)
	fit <- eBayes(fit)
	tt <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")
	pval <- tt$P.Value
	baseMean = unname(rowMeans(sweep(counts,2,nf,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	lfc <- tt$logFC
	cbind(pval = pval, padj = padj, lfc=lfc)
}



limma_voomZero2Filtered.pfun <-
function(counts, group, design = NULL, mc.cores = 2, niter=100) 
{   
	## limma voom pipeline ##
	library(limma) ; library(genefilter)
 	design = model.matrix(~group)
	nf <- edgeR::calcNormFactors(counts)
	zeroWeights <- zeroWeightsLibSizeDispFast(counts=counts, design=design, maxit=niter)
	y <- weightedVoom(counts, design, plot=FALSE, lib.size = colSums(counts)*nf, weights=zeroWeights)
	y$weights=y$weights*zeroWeights
	fit <- weightedLmFit(y, design, weights=y$weights)
	fit$df.residual=rowSums(zeroWeights)-ncol(design)
	fit <- eBayes(fit)
	tt <- topTable(fit,coef=2,n=nrow(counts), sort.by = "none")
	pval <- tt$P.Value
	baseMean = unname(rowMeans(sweep(counts,2,nf,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	lfc <- tt$logFC
	cbind(pval = pval, padj = padj, lfc=lfc)
}


baySeq.pfun <-
function(counts, group, design = NULL, p.out = "pvalue", mc.cores = 4)
{   
    ## baySeq pipeline ## 
	library(baySeq)
	library(snow)
	cl <- snow:::makeCluster(mc.cores, "SOCK")
	group <- as.numeric(as.character(group))
	cd <- new("countData", data = counts, replicates = group, groups = list(NDE = rep(1, length(group)), DE = group))
	cd@libsizes <- getLibsizes(cd)
	cd <- getPriors.NB(cd, equalDispersions = TRUE, estimation = "QL", cl = cl)
	cd <- getLikelihoods.NB(cd,pET = "BIC", cl = cl)
	try(snow:::stopCluster(cl), silent = TRUE)
	cd.table <- topCounts(cd, group = "DE", number = nrow(counts))
	id <- match( row.names(counts), row.names(cd.table))
	padj <- cd.table[id,]$FDR
	cbind(pval = pval, padj = padj)	#pval, padj are identical
}


samr_SAMseq.pfun <- 
function(counts, group, design = NULL, mc.cores = 4)
{
  ## SAMseq pipeline ##
  library(samr)
  f <- SAMseq(counts, group, resp.type = "Two class unpaired", fdr.output = 1)
  f.table = rbind(f$siggenes.table$genes.up, f$siggenes.table$genes.lo)
	fdr = rep(NA, nrow(counts)) #contains NA value
  fdr[as.numeric(f.table[, "Gene Name"])] <- as.numeric(f.table[, "q-value(%)"])/100  
  padj <- fdr
  cbind(pval = pval, padj = padj) #pval, padj are identical
}


#ShrinkBayes.pfun <- 
#function(counts, group, design, p.out = "pvalue", mc.cores)
#{   library(ShrinkBayes)
#	try(sfStop(), silent = TRUE)
#	try(sfRemoveAll(), silent = TRUE)
#	try(rm(g, envir = .GlobalEnv), silent = TRUE)
#	g <- as.factor(group)
#    assign("g", g, .GlobalEnv)
#	form <- y ~ 1 + g
#	form0 <- y ~ 1
#	shrinksimul <- ShrinkSeq(form = form, dat = counts, shrinkfixed = "g", fams = "zinb", ncpus = mc.cores)
#	fitall <- FitAllShrink(form, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)
#	fitall0 <- FitAllShrink(form0, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)
#	ptm <- proc.time()
#	npprior <-MixtureUpdatePrior(fitall = fitall, fitall0 = fitall0, shrinkpara="g", ncpus = mc.cores)
#	print(proc.time() - ptm)
#	nppostshr <- MixtureUpdatePosterior(fitall, npprior, fitall0, ncpus = mc.cores)
#	lfdr <- SummaryWrap(nppostshr)
#	fdr <- BFDR(lfdr)
#    rm(g, envir = .GlobalEnv)
#    pval <- as.vector(fdr)
#}


ShrinkBayes.pscript <-  expression(
   {
    ## ShrinkBayes pipeline implemented on special envir ##  
	ptm <- proc.time()
    library(ShrinkBayes)
	g <- as.factor(group)
	form <- y ~ 1 + g
	form0 <- y ~ 1
	shrinksimul <- ShrinkSeq(form = form, dat = counts, shrinkfixed = "g", fams = "zinb", ncpus = mc.cores)
	fitall <- FitAllShrink(form, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)
	fitall0 <- FitAllShrink(form0, dat = counts, fams = "zinb",shrinksimul = shrinksimul, ncpus = mc.cores)  
	npprior <-MixtureUpdatePrior(fitall = fitall, fitall0 = fitall0, shrinkpara="g", ncpus = mc.cores)
	nppostshr <- MixtureUpdatePosterior(fitall, npprior, fitall0, ncpus = mc.cores)
	lfdr <- SummaryWrap(nppostshr)
	print(proc.time() - ptm)
	try(sfStop(), silent = TRUE)
	try(sfRemoveAll(), silent = TRUE)
	fdr <- BFDR(lfdr)
	pGlobal <- as.vector(fdr)
	rm(list = c("ptm", "g", "form", "form0", "shrinksimul", "fitall", "fitall0", "npprior" ,"nppostshr", "lfdr", "fdr"))
	padj <- pGlobal
	cbind(pval = pval, padj = padj) #pval, padj are identical
	})

EBSeq.pfun <- 
function(counts, group, design = NULL, mc.cores = 4)
{
  ## EBSeq pipeline ##
  library(EBSeq)
  group <- as.factor(group)
  sizes = MedianNorm(counts)
  f <- EBTest(Data = counts, Conditions = group, sizeFactors = sizes, maxround = 5)
  pp = GetPPMat(f)
  padj <- fdr <- 1 - pp[, "PPDE"]
	if(!length(pval) == nrow(counts)) #check rm 0 counts
	{
	 i <- match(names(padj), rownames(counts))
	 Padj <- rep(NA, nrow(counts))
	 Padj[i] <- padj
	 padj <- Padj	
	}	
  cbind(pval = pval, padj = padj) #pval, padj are identical
}	


###### extra functions added by Koen

.comboGroups <- function(truths) 
# Function that returns a list of vectors of indices,
# where each vector refers to the rows with the same
# combination of TRUE/FALSE values in 'truths'.
# 
# written by Aaron Lun
# Created 24 October 2014
{
#	Integer packing will only work for 31 libraries at a time.	
	assembly <- list()
	collected <- 0L
	step <- 31L
	bits <- as.integer(2^(1:step-1L))

	while (collected < ncol(truths)) {
		upper <- pmin(ncol(truths) - collected, step)
		keys <- t(truths[,collected+1:upper,drop=FALSE]) * bits[1:upper]
		assembly[[length(assembly)+1L]] <- as.integer(colSums(keys))
		collected <- collected + step
	}

#	Figuring out the unique components.
	o <- do.call(order, assembly)
	nr <- nrow(truths)
	is.different <- logical(nr)
	for (i in 1:length(assembly)) { 
		is.different <- is.different | c(TRUE, diff(assembly[[i]][o])!=0L)
	}
	first.of.each <- which(is.different)
	last.of.each <- c(first.of.each[-1]-1L, nr)

#	Returning the groups.
	output <- list()
	for (u in 1:length(first.of.each)) {
		output[[u]] <- o[first.of.each[u]:last.of.each[u]]
	}
	return(output)
}



.residDF <- function(zero, design)
#	Effective residual degrees of freedom after adjusting for exact zeros
#	Gordon Smyth and Aaron Lun
#	Created 6 Jan 2014.  Last modified 2 Sep 2014
{
	nlib <- ncol(zero)
	ncoef <- ncol(design)
	nzero <- as.integer(rowSums(zero))

#	Default is no zero
	DF <- rep(nlib-ncoef,length(nzero))

#	All zero case
	DF[nzero==nlib] <- 0L

#	Anything in between?
	somezero <- nzero>0L & nzero<nlib
	if(any(somezero)) {
		zero2 <- zero[somezero,,drop=FALSE]
		groupings <- .comboGroups(zero2)

#		Identifying the true residual d.f. for each of these rows.			
		DF2 <- nlib-nzero[somezero]
		for (u in 1:length(groupings)) {
			i <- groupings[[u]]
			zeroi <- zero2[i[1],]
			DF2[i] <- DF2[i]-qr(design[!zeroi,,drop=FALSE])$rank
		}
		DF2 <- pmax(DF2, 0L)
		DF[somezero] <- DF2
	}

	DF
}



estimateDispWeighted = function (y, design = NULL, prior.df = NULL, trend.method = "locfit", tagwise = TRUE, span = NULL, min.row.sum = 5, grid.length = 21, 
    grid.range = c(-10, 10), robust = FALSE, winsor.tail.p = c(0.05, 
        0.1), tol = 1e-06, weights=NULL) 
{
    #adjusted by Koen VdB on 04 March 2016
    if (!is(y, "DGEList")) 
        stop("y must be a DGEList")
    trend.method <- match.arg(trend.method, c("none", "loess", 
        "locfit", "movingave"))
    ntags <- nrow(y$counts)
    nlibs <- ncol(y$counts)
    offset <- getOffset(y)
    AveLogCPM <- aveLogCPM(y)
    offset <- expandAsMatrix(offset, dim(y))
    sel <- rowSums(y$counts) >= min.row.sum
    spline.pts <- seq(from = grid.range[1], to = grid.range[2], 
        length = grid.length)
    spline.disp <- 0.1 * 2^spline.pts
    grid.vals <- spline.disp/(1 + spline.disp)
    l0 <- matrix(0, sum(sel), grid.length)
    if (is.null(design)) {
        cat("Design matrix not provided. Switch to the classic mode.\n")
        group <- y$samples$group <- as.factor(y$samples$group)
        if (length(levels(group)) == 1) 
            design <- matrix(1, nlibs, 1)
        else design <- model.matrix(~group)
        if (all(tabulate(group) <= 1)) {
            warning("There is no replication, setting dispersion to NA.")
            y$common.dispersion <- NA
            return(y)
        }
        pseudo.obj <- y[sel, ]
        q2q.out <- equalizeLibSizes(y[sel, ], dispersion = 0.01)
        pseudo.obj$counts <- q2q.out$pseudo
        ysplit <- splitIntoGroups(pseudo.obj)
        delta <- optimize(commonCondLogLikDerDelta, interval = c(1e-04, 
            100/(100 + 1)), tol = tol, maximum = TRUE, y = ysplit, 
            der = 0)
        delta <- delta$maximum
        disp <- delta/(1 - delta)
        q2q.out <- equalizeLibSizes(y[sel, ], dispersion = disp)
        pseudo.obj$counts <- q2q.out$pseudo
        ysplit <- splitIntoGroups(pseudo.obj)
        for (j in 1:grid.length) for (i in 1:length(ysplit)) l0[, 
            j] <- condLogLikDerDelta(ysplit[[i]], grid.vals[j], 
            der = 0) + l0[, j]
    }
    else {
        design <- as.matrix(design)
        if (ncol(design) >= ncol(y$counts)) {
            warning("No residual df: setting dispersion to NA")
            y$common.dispersion <- NA
            return(y)
        }
        glmfit <- glmFit(y$counts[sel, ], design, offset = offset[sel, 
            ], dispersion = 0.05, prior.count = 0, weights=weights[sel,]) ###www 
        zerofit <- (glmfit$fitted.values < 1e-04) & (glmfit$counts < 
            1e-04)
        by.group <- .comboGroups(zerofit)
        for (subg in by.group) {
            cur.nzero <- !zerofit[subg[1], ]
            if (!any(cur.nzero)) {
                next
            }
            if (all(cur.nzero)) {
                redesign <- design
            }
            else {
                redesign <- design[cur.nzero, , drop = FALSE]
                QR <- qr(redesign)
                redesign <- redesign[, QR$pivot[1:QR$rank], drop = FALSE]
                if (nrow(redesign) == ncol(redesign)) {
                  next
                }
            }
            last.beta <- NULL
            for (i in 1:grid.length) {
                out <- adjustedProfileLik(spline.disp[i], y = y$counts[sel, 
                  ][subg, cur.nzero, drop = FALSE], design = redesign, 
                  offset = offset[sel, ][subg, cur.nzero, drop = FALSE], 
                  start = last.beta, get.coef = TRUE, weights=weights[sel,][subg, cur.nzero, drop = FALSE]) ###www
                l0[subg, i] <- out$apl
                last.beta <- out$beta
            }
        }
    }
    out.1 <- WLEB(theta = spline.pts, loglik = l0, covariate = AveLogCPM[sel], 
        trend.method = trend.method, span = span, individual = FALSE, 
        m0.out = TRUE)
    y$common.dispersion <- 0.1 * 2^out.1$overall
    disp.trend <- 0.1 * 2^out.1$trend
    y$trended.dispersion <- rep(disp.trend[which.min(AveLogCPM[sel])], 
        ntags)
    y$trended.dispersion[sel] <- disp.trend
    y$trend.method <- trend.method
    y$AveLogCPM <- AveLogCPM
    y$span <- out.1$span
    if (!tagwise) 
        return(y)
    if (is.null(prior.df)) {
        glmfit <- glmFit(y$counts[sel, ], design, offset = offset[sel, 
            ], dispersion = disp.trend, prior.count = 0, weights=weights[sel,]) ###www
        df.residual <- glmfit$df.residual
        zerofit <- (glmfit$fitted.values < 1e-04) & (glmfit$counts < 
            1e-04)
        df.residual <- .residDF(zerofit, design)
        s2 <- glmfit$deviance/df.residual
        s2[df.residual == 0] <- 0
        s2 <- pmax(s2, 0)
        s2.fit <- squeezeVar(s2, df = df.residual, covariate = AveLogCPM[sel], 
            robust = robust, winsor.tail.p = winsor.tail.p)
        prior.df <- s2.fit$df.prior
    }
    ncoefs <- ncol(design)
    prior.n <- prior.df/(nlibs - ncoefs)
    if (trend.method != "none") {
        y$tagwise.dispersion <- y$trended.dispersion
    }
    else {
        y$tagwise.dispersion <- rep(y$common.dispersion, ntags)
    }
    too.large <- prior.n > 1e+06
    if (!all(too.large)) {
        temp.n <- prior.n
        if (any(too.large)) {
            temp.n[too.large] <- 1e+06
        }
        out.2 <- WLEB(theta = spline.pts, loglik = l0, prior.n = temp.n, 
            covariate = AveLogCPM[sel], trend.method = trend.method, 
            span = span, overall = FALSE, trend = FALSE, m0 = out.1$shared.loglik)
        if (!robust) {
            y$tagwise.dispersion[sel] <- 0.1 * 2^out.2$individual
        }
        else {
            y$tagwise.dispersion[sel][!too.large] <- 0.1 * 2^out.2$individual[!too.large]
        }
    }
    if (!robust) {
        y$prior.df <- prior.df
        y$prior.n <- prior.n
    }
    else {
        y$prior.df <- y$prior.n <- rep(Inf, ntags)
        y$prior.df[sel] <- prior.df
        y$prior.n[sel] <- prior.n
    }
    y
}


zeroWeightsLibSizeCpm <- function(counts, design, initialWeightAt0=TRUE, niter=30, plot=FALSE, plotW=FALSE){
    require(edgeR)
    require(mgcv)
    counts <- DGEList(counts)
    counts <- edgeR::calcNormFactors(counts)
    effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
    logEffLibSize <- log(effLibSize)
    avCpm <- aveLogCPM(counts)
    histCpm <- hist(avCpm,breaks=200, plot=FALSE)
    cpmMids <- histCpm$mids
    cpmBreaks <- histCpm$breaks
    cpmMidsFit <- rep(cpmMids,each=ncol(counts))
    logEffLibSizeFit <- rep(logEffLibSize,length(cpmMids))
    binId <- sapply(avCpm, function(x) sum(cpmBreaks<x))

    ## E-step initialization
    zeroId <- counts$counts==0
    w <- matrix(1,nrow=nrow(counts),ncol=ncol(counts))
    if(initialWeightAt0) w[zeroId] <- 0.1 else w[zeroId] <- 0.9

    for(i in 1:niter){
	counts$weights <- w
	
	### M-step counts
	counts <- estimateDispWeighted(counts, design, prior.df=0, trend.method="none", weights=counts$weights)
	if(plot) plotBCV(counts)
	fit <- glmFit(counts, design)
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)
	
	### M-step mixture parameter: model zero probability
	expZeroCpmBin <- sapply(1:length(cpmMids), function(i){
		   binGenes <- which(binId==i)
		   if(length(binGenes)==0){return(cbind(rep(0,ncol(counts)),rep(0,ncol(counts))))} else
		       if(length(binGenes)==1){ return(cbind(1-w[binGenes,],w[binGenes,]))} else
			   if(length(binGenes)>1) return(cbind(colSums(1-w[binGenes,]),colSums(w[binGenes,])))
	}, simplify=FALSE) #first column is P(zero) hence successes
	successes <- unlist(lapply(expZeroCpmBin, function(x) x[,1]))
	failures <- unlist(lapply(expZeroCpmBin, function(x) x[,2]))
	zeroFit <- gam(cbind(successes,failures) ~ logEffLibSizeFit + s(cpmMidsFit), family="binomial")
	pi0Hat <- t(sapply(avCpm,function(x) predict(zeroFit,newdata=data.frame(cpmMidsFit=x,logEffLibSizeFit=logEffLibSize), type="response")))
		
	## E-step: Given estimated parameters, calculate expected value of weights
	w <- 1-pi0Hat*zeroId/(pi0Hat*zeroId+(1-pi0Hat)*likC+1e-15) 
	if(plotW) hist(w[zeroId])
    }
    return(w)
}



zeroWeightsLibSize <- function(counts, design, initialWeightAt0=TRUE, niter=30, plot=FALSE, plotW=FALSE, designZI=NULL){
    require(edgeR)
    counts <- DGEList(counts)
    counts <- edgeR::calcNormFactors(counts)
    effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
    logEffLibSize <- log(effLibSize)
    zeroId <- counts$counts==0
    w <- matrix(1,nrow=nrow(counts),ncol=ncol(counts))
    for(k in 1:ncol(w)) w[counts$counts[,k]==0,k] <- 1-mean(counts$counts[,k]==0)

    for(i in 1:niter){
	counts$weights <- w
	
	### M-step counts
	counts <- estimateDispWeighted(counts, design, prior.df=0, trend.method="none", weights=counts$weights)
	if(plot) plotBCV(counts)
	fit <- glmFit(counts, design)
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)
	
	### M-step mixture parameter: model zero probability
	successes <- colSums(1-w) #P(zero)
	failures <- colSums(w) #1-P(zero)
	if(is.null(designZI)){
	zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")} else{
	zeroFit <- glm(cbind(successes,failures) ~-1+designZI, family="binomial")}
	pi0Hat <- predict(zeroFit,type="response") 
	
	## E-step: Given estimated parameters, calculate expected value of weights
	pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
	w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)
	if(plotW) hist(w[zeroId])
    }
    return(w)
}

zeroWeightsLibSizeFast <- function(counts, design, initialWeightAt0=TRUE, maxit=100, plot=FALSE, plotW=FALSE, designZI=NULL, wTol=1e-4){
    require(edgeR)
    if(plot | plotW) par(mfrow=c(1,plot+plotW))    
    counts <- DGEList(counts)
    counts <- edgeR::calcNormFactors(counts)
    effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
    logEffLibSize <- log(effLibSize)
    zeroId <- counts$counts==0
    w <- matrix(1,nrow=nrow(counts),ncol=ncol(counts), dimnames=list(c(1:nrow(counts)), NULL))
    #if(initialWeightAt0) w[zeroId] <- 0.1 else w[zeroId] <- 0.9
      ## starting values based on P(zero) in the library
    for(k in 1:ncol(w)) w[counts$counts[,k]==0,k] <- 1-mean(counts$counts[,k]==0)
    
    wFinal <- matrix(NA,nrow=nrow(counts),ncol=ncol(counts))
    active <- rowSums(counts$counts==0)>0 #work with genes with at least 1 zero
    wFinal[!active,]=w[!active,]
    counts <- counts[active,]
    w <- w[active,]
    llOld <- matrix(-1e4,nrow=nrow(counts),ncol=ncol(counts))
    likCOld <- matrix(0,nrow=nrow(counts),ncol=ncol(counts))

    for(i in 1:maxit){
        zeroId <- counts$counts==0	
	counts$weights <- w
	
	### M-step counts
	counts <- estimateGLMCommonDisp(counts, design, interval=c(0,10))
	counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0, min.row.sum=1)
	if(plot) plotBCV(counts)
	fit <- glmFit(counts, design)
	if(i>1) likCOld <- likC[!converged,]	
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)
	
	### M-step mixture parameter: model zero probability
	successes <- colSums(1-w) #P(zero)
	failures <- colSums(w) #1-P(zero)
	if(is.null(designZI)){
	zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")} else{
	zeroFit <- glm(cbind(successes,failures) ~-1+designZI, family="binomial")}
	pi0Hat <- predict(zeroFit,type="response") 
	
	## E-step: Given estimated parameters, calculate expected value of weights
	pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
	wOld <- w
	w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)
	rownames(w) <- rownames(wOld)
	if(plotW) hist(w[zeroId])

	## expected complete data log-likelihood
	if(i>1) llOld <- ll[!converged,]
	ll <- w*log(pi0HatMat) + (1-w)*log(1-pi0HatMat) + (1-w)*log(likC)

	## filter genes with converged weights
	#converged <- rowSums( (likC-likCOld)/likC )<1e-3

	   ## weights
	#converged <- rowSums(abs(w-wOld)<wTol)==ncol(counts$counts)
	  ## metagenomeSeq's stopping rule
	#nll <- -rowSums(ll)
	#nllOld <- -rowSums(llOld)
	#epsilon <- (nllOld-nll)/nllOld
	#converged <- epsilon<1e-4 #claims convergence way too fast.
	converged=FALSE
	if(any(converged)){
	    wFinal[as.numeric(rownames(w)[converged]),] = w[converged,]
	    w <- matrix(w[!converged,],ncol=ncol(counts), dimnames=list(c(rownames(w)[!converged]),NULL))
	    counts <- counts[!converged,]
	}
	#cat(paste0("mean diff in L: ",round(mean(rowSums(exp(ll)-exp(llOld))),2),". active features: ",nrow(counts),"\n"))
	cat(paste0("iteration ",i))
	#metagenomeSeq for example seems to report mean(ll) instead of difference.
	if(all(converged)) break
	if(i==maxit){
	    wFinal[apply(wFinal,1,function(row) any(is.na(row))),] = w
	    break
	}
    }
    return(wFinal)
}



zeroWeightsLibSizeDispFast <- function(counts, design, initialWeightAt0=TRUE, maxit=100, plot=FALSE, plotW=FALSE, designZI=NULL, llTol=1e-4, normalization="TMM"){
    require(edgeR) ; require(DESeq2)
    if(plot | plotW) par(mfrow=c(1,plot+plotW))    
    counts <- DGEList(counts)
    #counts <- edgeR::calcNormFactors(counts)
    if(normalization=="TMM"){
    	counts = edgeR::calcNormFactors(counts)
    } else if(normalization=="DESeq2"){
    	 dse = DESeqDataSetFromMatrix(counts$counts, colData=data.frame(group=group), design=~group)
    	 dse = DESeq2::estimateSizeFactors(dse)
    	 counts$samples$norm.factors = 1/dse$sizeFactor
    } else if(normalization=="DESeq2_pos"){
    	dse = DESeqDataSetFromMatrix(counts$counts, colData=data.frame(group=group), design=~group)
    	 dse = DESeq2::estimateSizeFactors(dse, type="poscounts")
    	 counts$samples$norm.factors = 1/dse$sizeFactor
    } else stop("normalization must be either TMM, DESeq2 or DESeq2_pos")
    effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
    logEffLibSize <- log(effLibSize)
    zeroId <- counts$counts==0
    w <- matrix(1,nrow=nrow(counts),ncol=ncol(counts), dimnames=list(c(1:nrow(counts)), NULL))
      ## starting values based on P(zero) in the library
    for(k in 1:ncol(w)) w[counts$counts[,k]==0,k] <- 1-mean(counts$counts[,k]==0)
    
    llOld <- matrix(-1e4,nrow=nrow(counts),ncol=ncol(counts))
    likCOld <- matrix(0,nrow=nrow(counts),ncol=ncol(counts))
    converged=FALSE
    j=0

    for(i in 1:maxit){
	j=j+1
        zeroId <- counts$counts==0	
	counts$weights <- w
	
	### M-step counts
	#only estimate dispersions every 5 iterations
	#if(i==1 | is.wholenumber(i/10)){
	if(i==1 | converged){
	counts <- estimateGLMCommonDisp(counts, design, interval=c(0,10))
	counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0, min.row.sum=1)
	}
	if(plot) plotBCV(counts)	
	fit <- glmFit(counts, design)
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)
	
	### M-step mixture parameter: model zero probability
	successes <- colSums(1-w) #P(zero)
	failures <- colSums(w) #1-P(zero)
	if(is.null(designZI)){
	zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")} else{
	zeroFit <- glm(cbind(successes,failures) ~-1+designZI, family="binomial")}
	pi0Hat <- predict(zeroFit,type="response") 
	
	## E-step: Given estimated parameters, calculate expected value of weights
	pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
	w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)

	## data log-likelihood
	if(i>1) llOld=ll
	ll <- log(pi0HatMat*zeroId + (1-pi0HatMat)*likC)

	delta <- (rowSums(ll)-rowSums(llOld))/(rowSums(llOld)+llTol)
	if(mean(abs(delta) < llTol)>.999){ #if 99.9% has converged
	    if(j==1 & mean(abs(delta) < llTol)>.999){ #final convergence?
	    	cat(paste0("converged. \n")) ; return(w)}
	    j=0 
	    converged=TRUE} else {converged=FALSE}
	cat(paste0("iteration: ",i,". mean conv.: ",mean(abs(delta) < llTol),"\n"))
	if(plotW) hist(w[zeroId],main=paste0("iteration: ",i,". mean conv.: ",mean(abs(delta) < llTol)))	
    }
    return(w)
}



zeroWeightsLibSizeScran <- function(counts, design, initialWeightAt0=TRUE, niter=30, plot=FALSE, plotW=FALSE, designZI=NULL){
    require(edgeR) ; require(scran)
    counts <- DGEList(counts)
    sf <- computeSumFactors(counts$counts, positive=TRUE)
    if(any(sf==0)){
	keep <- sf>0
	sf <- sf[keep]
	counts <- counts[,keep]
	design <- design[keep,]
    }
    counts$samples$norm.factors <- 1/sf
    effLibSize <- counts$samples$lib.size*counts$samples$norm.factors
    logEffLibSize <- log(effLibSize)
    zeroId <- counts$counts==0
    w <- matrix(1,nrow=nrow(counts),ncol=ncol(counts))
    if(initialWeightAt0) w[zeroId] <- 0.1 else w[zeroId] <- 0.9

    for(i in 1:niter){
	print(i)
	counts$weights <- w
	
	### M-step counts
	#counts <- estimateDisp(counts, design, prior.df=0, trend.method="none")
	counts <- estimateGLMCommonDisp(counts, design)
	counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0)
	if(plot) plotBCV(counts)
	fit <- glmFit(counts, design)
	likC <- dnbinom(counts$counts, mu=fit$fitted.values, size=1/counts$tagwise.dispersion)
	
	### M-step mixture parameter: model zero probability
	successes <- colSums(1-w) #P(zero)
	failures <- colSums(w) #1-P(zero)
	if(is.null(designZI)){
	zeroFit <- glm(cbind(successes,failures) ~ logEffLibSize, family="binomial")} else{
	    zeroFit <- glm(cbind(successes,failures) ~-1+designZI, family="binomial")}
	pi0Hat <- predict(zeroFit,type="response") 
	
	## E-step: Given estimated parameters, calculate expected value of weights
	pi0HatMat <- expandAsMatrix(pi0Hat,dim=dim(counts),byrow=TRUE)
	w <- 1-pi0HatMat*zeroId/(pi0HatMat*zeroId+(1-pi0HatMat)*likC*zeroId+1e-15)
	if(plotW) hist(w[zeroId])
    }
    return(w)
}



pvalueAdjustment_kvdb <- function(baseMean, filter, pValue,
                             theta, alpha=0.05, pAdjustMethod="BH") {
  # perform independent filtering
    if (missing(filter)) {
      filter <- baseMean
    }
    if (missing(theta)) {
      lowerQuantile <- mean(filter == 0)
      if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
      theta <- seq(lowerQuantile, upperQuantile, length=50)
    }

    # do filtering using genefilter
    stopifnot(length(theta) > 1)
    filtPadj <- filtered_p(filter=filter, test=pValue,
                           theta=theta, method=pAdjustMethod) 
    numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
    # prevent over-aggressive filtering when all genes are null,
    # by requiring the max number of rejections is above a fitted curve.
    # If the max number of rejection is not greater than 10, then don't
    # perform independent filtering at all.
    lo.fit <- lowess(numRej ~ theta, f=1/5)
    if (max(numRej) <= 10) {
      j <- 1
    } else { 
      residual <- if (all(numRej==0)) {
        0
      } else {
        numRej[numRej > 0] - lo.fit$y[numRej > 0]
      }
      thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
      j <- if (any(numRej > thresh)) {
        which(numRej > thresh)[1]
      } else {
        1  
      }
    }
    padj <- filtPadj[, j, drop=TRUE]
    cutoffs <- quantile(filter, theta)
    filterThreshold <- cutoffs[j]
    filterNumRej <- data.frame(theta=theta, numRej=numRej)
    filterTheta <- theta[j]

    return(list(padj=padj, filterThreshold=filterThreshold, filterTheta=filterTheta, filterNumRej = filterNumRej, lo.fit=lo.fit, alpha=alpha))

}


edgeREM.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	zeroWeights = zeroWeightsLibSizeCpm(d, design, plot=FALSE, niter=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights)
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
  	edger.lrt <- glmLRT(edger.fit,coef=2)
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}


edgeREMOldF.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	zeroWeights = zeroWeightsLibSizeCpm(d, design, plot=FALSE, niter=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights)
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}

edgeREMLibSize.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	#not adding a design matrix models the zeroes with the library size automatically
	zeroWeights = zeroWeightsLibSize(d, design, plot=FALSE, niter=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
  	edger.lrt <- glmLRT(edger.fit,coef=2)
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}

edgeREMLibSizeOldF.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	#not adding a design matrix models the zeroes with the library size automatically
	zeroWeights = zeroWeightsLibSize(d, design, plot=FALSE, niter=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	lfc <- edger.lrt$table$logFC
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj,lfc)
	out[is.na(out)] <- 1
	return(out)
}

edgeREMLibSizeOldFFiltered.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	#not adding a design matrix models the zeroes with the library size automatically
	zeroWeights = zeroWeightsLibSize(d, design, plot=FALSE, niter=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	lfc <- edger.lrt$table$logFC
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))	
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj,lfc)
	out[is.na(out)] <- 1
	return(out)
}

edgeREMLibSizeFastOldF.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	#not adding a design matrix models the zeroes with the library size automatically
	zeroWeights = zeroWeightsLibSizeFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}


edgeREMLibSizeFastOldFFilteredEdgeR.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR) ; library(genefilter)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	#not adding a design matrix models the zeroes with the library size automatically
	zeroWeights = zeroWeightsLibSizeFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	pval <- edger.lrt$table$PValue
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))	
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	#padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}

edgeREMLibSizeDispFastOldFFilteredEdgeR.pfun=function(counts, group, design=NULL, mc.cores=2, niter=200){
	library(edgeR) ; library(genefilter)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
    design = model.matrix(~ group)
	zeroWeights = zeroWeightsLibSizeDispFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
	lfc <- edger.lrt$table$logFC
  	pval <- edger.lrt$table$PValue
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	out=cbind(pval,padj,lfc)
	return(out)
}

WaldTest <- function(y)
{
    X <- y$design
    coef <- y$coef
    phi <- y$dispersion
    mu <- y$fitted
    mu <- mu+1e-6
    nr <- nrow(mu)
    vb <- coef
    for (i in seq(nr))
    { 
        W <- diag((mu[i,]^-1+phi[i])^-1)
        xtwx <- t(X) %*% W %*% X
        xtwxInv <- solve(xtwx)
    vb[i,] <- diag(xtwxInv %*% xtwx %*% xtwxInv)
    }
    coef/sqrt(vb)
}


edgeREMLibSizeDispFastWaldFilteredEdgeR.pfun=function(counts, group, design=NULL, mc.cores=2, niter=200){
	library(edgeR) ; library(genefilter)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	zeroWeights = zeroWeightsLibSizeDispFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	edger.fit <- glmFit(d, design) #uses weights
	waldStats=WaldTest(edger.fit)
	#pval=(1-pnorm(abs(waldStats[,"group1"])))*2
	pval=(1-pt(abs(waldStats[,"group1"]),df=rowSums(zeroWeights)-ncol(design)))*2	
	lfc=NA
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	out=cbind(pval,padj,lfc)
	return(out)
}


edgeREMDESeq2Norm.pfun=function(counts, group, design=NULL, mc.cores=2, niter=200){
	library(edgeR) ; library(genefilter)

  	colData <- data.frame(group)
	dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	hlp=estimateSizeFactors(dse)
    ## edgeR standard pipeline ##
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
 	d$samples$norm.factors=1/hlp$sizeFactor
        design = model.matrix(~ group)
	zeroWeights = zeroWeightsLibSizeDispFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
	lfc <- edger.lrt$table$logFC
  	pval <- edger.lrt$table$PValue
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	out=cbind(pval,padj,lfc)
	return(out)
}



edgeREMLibSizeDispFastLRTFilteredEdgeR.pfun=function(counts, group, design=NULL, mc.cores=2, niter=200){
	library(edgeR) ; library(genefilter)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	zeroWeights = zeroWeightsLibSizeDispFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	edger.fit <- glmFit(d, design) #uses weights
  	edger.lrt <- glmLRT(edger.fit,coef=2)
	lfc <- edger.lrt$table$logFC
  	pval <- edger.lrt$table$PValue
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pval)
  	padj <- hlp$padj
	out=cbind(pval,padj,lfc)
	return(out)
}




edgeREMLibSizeFastOldFFilteredDESeq2.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
    	## DESeq2 to check filtering
	library(DESeq2)
	colData <- data.frame(group)
	dse <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ group)
	colData(dse)$group <- as.factor(colData(dse)$group)
	dse <- DESeq(dse)
	res <- results(dse)
	filtered <- is.na(res$padj)

	## edgeR
	library(edgeR)
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
	keep <- !filtered
	#d <- d[keep,]
	#d$samples$lib.size <- colSums(d$counts)
        design = model.matrix(~ group)
	#not adding a design matrix models the zeroes with the library size automatically
	zeroWeights = zeroWeightsLibSizeFast(d, design, plot=FALSE, maxit=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
	pval <- vector(length=nrow(d))
	pval[]=NA
	pval[keep]=edger.lrt$table$PValue[keep]
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}

edgeREMLibSizeOldFScran.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR) ; require(scran)
	d <- DGEList(counts = counts, group = group )
        if(is.null(design)) design = model.matrix(~ group)
	sf <- computeSumFactors(d$counts, positive=TRUE)
	if(any(sf==0)){
		keep <- sf>0
		sf <- sf[keep]
		d <- d[,keep]
		design <- design[keep,]
	}
	#not adding a design matrix models the zeroes with the library size automatically
	zeroWeights = zeroWeightsLibSize(d, design, plot=FALSE, niter=niter, initialWeightAt0=TRUE, plotW=FALSE)
	d$samples$norm.factors <- 1/sf
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}


edgeREMLibSizeRichnessOldF.pfun=function(counts, group, design=NULL, mc.cores=2, niter=50){
	library(edgeR) ; library(vegan)
	expRichness <- log(rarefy(t(counts),sample=min(colSums(counts))))
	d <- DGEList(counts = counts, group = group )
	d <- edgeR::calcNormFactors(d)
        design = model.matrix(~ group)
	effLibSize <- d$samples$lib.size*d$samples$norm.factors
	logEffLibSize <- log(effLibSize)
	zeroWeights = zeroWeightsLibSize(d, design, plot=FALSE, niter=niter, initialWeightAt0=TRUE, plotW=FALSE, designZI=model.matrix(~logEffLibSize+expRichness))
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	pval <- edger.lrt$table$PValue
  	pval[rowSums(counts) == 0] <- NA
  	padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=cbind(pval,padj)
	out[is.na(out)] <- 1
	return(out)
}

scde.pfun <- function(counts, group, design=NULL, mc.cores=1, niter=NULL){
    counts = matrix(as.integer(counts),nrow=nrow(counts),ncol=ncol(counts))
    if(is.null(colnames(counts))) colnames(counts)=paste0("sample",1:ncol(counts))
    require(scde)
    
    # calculate error models
    o.ifm <- scde.error.models(counts = counts, groups = group, n.cores = mc.cores, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 0)
    # estimate gene expression prior
    o.prior <- scde.expression.prior(models = o.ifm, counts = counts, length.out = 400, show.plot = FALSE)
    # calculate differential expression
    ediff <- scde.expression.difference(o.ifm, counts, o.prior, groups  =  group, n.randomizations  =  150, n.cores  =  mc.cores, verbose  =  0)
    lfc <- ediff$mle
    pval=(1-pnorm(abs(ediff$Z)))*2
    padj=(1-pnorm(abs(ediff$cZ)))*2 
    out = cbind(pval,padj,lfc)
    return(out)
}

MAST_oldWrong.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    require(MAST)
    #convert to tpm
    tpm <- counts*1e6/colSums(counts) #consider equal length across all features since we did not take length into account while simulating.
    sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(group=group))
    ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
    CD <- cData(sca)
    CD$ngeneson <- ngeneson
    CD$cngeneson <- CD$ngeneson-mean(ngeneson)
    cData(sca) <- CD  
    ## differential expression 
    fit <- zlm.SingleCellAssay(~cngeneson+group,sca=sca,method="bayesglm",ebayes=TRUE)
    L=matrix(0,nrow=ncol(coef(fit,"D")))
    rownames(L)=colnames(coef(fit,"D"))
    L["group1",]=1
    lrFit <- lrTest(fit, hypothesis=L)
    pval=lrFit[,'hurdle','Pr(>Chisq)']
    padj=p.adjust(pval,method="BH")
    lfc=NA
    out=cbind(pval,padj,lfc)
    return(out)
}


MAST.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    require(MAST)
    #convert to tpm
    tpm <- log2((counts+1)*1e6/colSums(counts)) #consider equal length across all features since we did not take length into account while simulating.
    sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(group=group))
    ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
    CD <- cData(sca)
    CD$ngeneson <- ngeneson
    CD$cngeneson <- CD$ngeneson-mean(ngeneson)
    cData(sca) <- CD  
    ## differential expression 
    fit <- zlm.SingleCellAssay(~cngeneson+group,sca=sca,method="bayesglm",ebayes=TRUE)
    L=matrix(0,nrow=ncol(coef(fit,"D")))
    rownames(L)=colnames(coef(fit,"D"))
    L["group1",]=1
    lrFit <- lrTest(fit, hypothesis=L)
    pval=lrFit[,'hurdle','Pr(>Chisq)']
    padj=p.adjust(pval,method="BH")
    lfc=NA
    out=cbind(pval,padj,lfc)
    return(out)
}



MAST_count_oldWrong.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    require(MAST)
    #convert to tpm
    tpm <- counts*1e6/colSums(counts) #consider equal length across all features since we did not take length into account while simulating.
    sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(group=group))
    ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
    CD <- cData(sca)
    CD$ngeneson <- ngeneson
    CD$cngeneson <- CD$ngeneson-mean(ngeneson)
    cData(sca) <- CD  
    ## differential expression 
    fit <- zlm.SingleCellAssay(~cngeneson+group,sca=sca,method="bayesglm",ebayes=TRUE)
    L=matrix(0,nrow=ncol(coef(fit,"D")))
    rownames(L)=colnames(coef(fit,"D"))
    L["group1",]=1
    lrFit <- lrTest(fit, hypothesis=L)
    pval=lrFit[,'cont','Pr(>Chisq)'] 
    padj=p.adjust(pval,method="BH")
    lfc=NA
    out=cbind(pval,padj,lfc)
    return(out)
}


MAST_count.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    require(MAST)
    #convert to tpm
    tpm <- log2((counts+1)*1e6/colSums(counts)) #consider equal length across all features since we did not take length into account while simulating.
    sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(group=group))
    ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
    CD <- cData(sca)
    CD$ngeneson <- ngeneson
    CD$cngeneson <- CD$ngeneson-mean(ngeneson)
    cData(sca) <- CD  
    ## differential expression 
    fit <- zlm.SingleCellAssay(~cngeneson+group,sca=sca,method="bayesglm",ebayes=TRUE)
    L=matrix(0,nrow=ncol(coef(fit,"D")))
    rownames(L)=colnames(coef(fit,"D"))
    L["group1",]=1
    lrFit <- lrTest(fit, hypothesis=L)
    pval=lrFit[,'cont','Pr(>Chisq)'] 
    padj=p.adjust(pval,method="BH")
    lfc=NA
    out=cbind(pval,padj,lfc)
    return(out)
}



positiveHurdleTPM.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    normFactors <- edgeR::calcNormFactors(counts)
    #counts <- sweep(counts[1:3,],2,FUN="*",STATS=normFactors)
    CDR <- apply(counts,2,function(x) mean(x==0))
    tpm <- counts*1e6/colSums(counts) #consider equal length across all features since we did not take length into account while simulating.
    response <- log2(tpm+1)
    zeroID <- counts==0
    weights <- 1-zeroID
    pval <- sapply(1:nrow(response), function(i) try(summary(lm(response[i,]~group+CDR,weights=weights[i,]))[[5]][2,4]))
    pval <- as.numeric(pval) #gives NA for genes only expressed in one condition
    padj <- p.adjust(pval,"BH")
    lfc <- NA
    out=cbind(pval,padj,lfc)
    return(out)
}

positiveHurdleTPMCDR.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    normFactors <- edgeR::calcNormFactors(counts)
    #counts <- sweep(counts[1:3,],2,FUN="*",STATS=normFactors)
    CDR <- apply(counts,2,function(x) mean(x==0))
    tpm <- counts*1e6/colSums(counts) #consider equal length across all features since we did not take length into account while simulating.
    response <- log2(tpm+1)
    zeroID <- counts==0
    weights <- 1-zeroID
    pval <- sapply(1:nrow(response), function(i) try(summary(lm(response[i,]~group+CDR,weights=weights[i,]))[[5]][2,4]))
    pval <- as.numeric(pval) #gives NA for genes only expressed in one condition
    padj <- p.adjust(pval,"BH")
    lfc <- NA
    out=cbind(pval,padj,lfc)
    return(out)
}

positiveHurdleLimmaTPM.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
	#this function gives identical point estimates compared to a positive hurdle limma analysis of the like
	#d=DGEList(counts)
	#d=calcNormFactors(d)
	#weights=1-(counts==0)
	#v=voom(d,design,weights=weights)
	#v$weights=weights
	#fit=lmFit(v,design)
    normFactors <- edgeR::calcNormFactors(counts)
    lib.size=normFactors*colSums(counts)
    response <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06)) #consider equal length across all features since we did not take length into account while simulating.
    zeroID <- counts==0
    weights <- 1-zeroID
    pval <- sapply(1:nrow(response), function(i) try(summary(lm(response[i,]~group,weights=weights[i,]))[[5]][2,4]))
    pval <- as.numeric(pval) #gives NA for genes only expressed in one condition
    padj <- p.adjust(pval,"BH")
    lfc <- NA
    out=cbind(pval,padj,lfc)
    return(out)
}


regularHurdleTPM.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    CDR <- apply(counts,2,function(x) mean(x==0))
    zeroID <- (counts==0)+0
    weights <- 1-zeroID
    lib.size=colSums(counts)*edgeR::calcNormFactors(counts)
    tpm <- counts*1e6/lib.size #consider equal length across all features since we did not take length into account while simulating.
    response <- log2(tpm+1)
    pval=vector(length=nrow(counts))

    for(g in 1:nrow(counts)){
    ### binomial model
    mFull=glm(zeroID[g,]~group+CDR,family="binomial")
    mNull=glm(zeroID[g,]~CDR,family="binomial")
    lrtStatBinom=as.numeric(logLik(mFull)-logLik(mNull))*2

    ### count model
    mContFull=lm(response[g,]~group+CDR,weights=weights[g,])
    mContNull=lm(response[g,]~CDR,weights=weights[g,])
    lrStatCont=as.numeric(logLik(mContFull)-logLik(mContNull))*2

    lrStatHurdle=lrtStatBinom+lrStatCont
    pval[g]=1-pchisq(lrStatHurdle,df=2)
    }
    pval <- as.numeric(pval) #gives NA for genes only expressed in one condition
    padj <- p.adjust(pval,"BH")
    lfc <- NA
    out=cbind(pval,padj,lfc)
    return(out)
}



monocle.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){	
    require(monocle)
    pheno=AnnotatedDataFrame(data.frame(group=group))
    rownames(pheno)=colnames(counts)
    fData=AnnotatedDataFrame(data.frame(geneName=rownames(counts)))
    rownames(fData)=rownames(counts)
    cd = newCellDataSet(counts, phenoData=pheno, featureData=fData, expressionFamily=negbinomial(), lowerDetectionLimit=1)
    cd = estimateSizeFactors(cd)
    cd = estimateDispersions(cd)
    res = differentialGeneTest(cd, fullModelFormulaStr = "~group", reducedModelFormulaStr="~1", relative_expr=FALSE, verbose=TRUE)
    pval = res$pval
    padj = p.adjust(pval,"BH")
    out=cbind(pval,padj)
    return(out)
}

monocle_tpm.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){	
    require(monocle)
    pheno=AnnotatedDataFrame(data.frame(group=group))
    rownames(pheno)=colnames(counts)
    fData=AnnotatedDataFrame(data.frame(geneName=rownames(counts)))
    rownames(fData)=rownames(counts)
    tpm <- counts*1e6/colSums(counts) #consider equal length across all features since we did not take length into account while simulating.    
    cd = newCellDataSet(tpm, phenoData=pheno, featureData=fData, expressionFamily=negbinomial.size(), lowerDetectionLimit=0.1)
    cd = DESeq2::estimateSizeFactors(cd)
    cd = DESeq2::estimateDispersions(cd)
    res = differentialGeneTest(cd, fullModelFormulaStr = "~group", reducedModelFormulaStr="~1", relative_expr=FALSE)
    pval = res$pval
    padj = p.adjust(pval,"BH")
    out=cbind(pval,padj)
    return(out)
}

metagenomeSeq.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
    require(metagenomeSeq)
    design <- model.matrix(~group)
    pheno <- AnnotatedDataFrame(data.frame(group=group))
    rownames(pheno) <- colnames(counts)
    p <- cumNormStatFast(counts)
    dat <- newMRexperiment(counts=counts, phenoData=pheno, featureData = NULL, libSize = colSums(counts), normFactors = metagenomeSeq::calcNormFactors(counts, p=p))
    fit <- fitZig(dat,design)
    lfc <- fit$eb$coefficients[,"group1"]
    pval <- fit$eb$p.value[,"group1"]
    padj <- p.adjust(pval)
    out <- cbind(pval,padj,lfc)
    return(out)
}

NODES.pfun <- function(counts, group, design=NULL, mc.cores=2, niter=NULL){
	require(NODES)
	g=ifelse(group==0,"A","B")
	colnames(counts)=g
	normCounts=pQ(counts)
	res=NODES(data=normCounts,group=colnames(normCounts))
	pval=vector(length=nrow(counts))
	names(pval)=rownames(counts)
	pval[rownames(normCounts)]=res$Fisher
	pval[is.na(pval)]=1
	padj=p.adjust(pval,"BH")
	lfc=NA
	out=cbind(pval,padj,lfc)
	return(out)
}



### end of added functions by Koen



getDataset <- function(counts, drop.extreme.dispersion = 0.1, drop.low.lambda = TRUE) {  
## this function generates NB parameters from real dataset ##
## it is low-level function of NBsim ##
	d <- DGEList(counts)
	d <- edgeR::calcNormFactors(d)
	cp <- round(cpm(d,normalized.lib.sizes=TRUE),1)
	if(drop.low.lambda){
	    d <- d[rowSums(cp>1) >= 2, ]
	}
	d$AveLogCPM <- log2(rowMeans(cpm(d, prior.count = 1e-5)))
	d <- estimateGLMCommonDisp(d)
	d <- estimateGLMTrendedDisp(d)
	d <- estimateGLMTagwiseDisp(d)
	dispersion <- d$tagwise.dispersion
	AveLogCPM <- d$AveLogCPM
	if(is.numeric(drop.extreme.dispersion))
	{   
		bad <- quantile(dispersion, 1-drop.extreme.dispersion, names = FALSE)
		ids <- dispersion <= bad
		AveLogCPM <- AveLogCPM[ids]
		dispersion <- dispersion[ids]
	}
	dataset.AveLogCPM <- AveLogCPM
	dataset.dispersion <- dispersion
	dataset.lib.size <- d$samples$lib.size
	dataset.nTags <- nrow(d)
	list(dataset.AveLogCPM = dataset.AveLogCPM, dataset.dispersion = dataset.dispersion, dataset.lib.size = dataset.lib.size, dataset.nTags = dataset.nTags)
}



NBsim <-
function(dataset, group, nTags = 10000, nlibs = length(group), fix.dispersion = NA, lib.size = NULL, drop.low.lambda = TRUE, drop.extreme.dispersion = 0.1,  add.outlier = FALSE, outlierMech = c("S", "R", "M"), pOutlier = 0.1, min.factor = 1.5, max.factor = 10, pDiff=.1, pUp=.5, foldDiff=3, name = NULL, save.file = FALSE, file = NULL, only.add.outlier = FALSE, verbose=TRUE, ind=NULL)

{   
## NBsim generate simulated count from the real dataset followed by the NB model ##		
	require(edgeR)
	group = as.factor(group)

	sample.fun <- function(object)
	{
		## it is low-level function of NBsim ##
        ## it samples from the real dataset ## 
		nlibs <- object$nlibs
		nTags <- object$nTags
		AveLogCPM <-object$dataset$dataset.AveLogCPM
		dispersion <- object$dataset$dataset.dispersion
		
                id_r <- sample(length(AveLogCPM), nTags, replace = TRUE)
		object$AveLogCPM <- AveLogCPM[id_r] #added by Koen to use for adding zeroes
		Lambda <- 2^(AveLogCPM[id_r])
		Lambda <- Lambda/sum(Lambda)
		Dispersion <- dispersion[id_r]
		id_0<- Lambda == 0
		Lambda <- Lambda[!id_0]
		Dispersion <- Dispersion[!id_0]
		Lambda <- expandAsMatrix(Lambda, dim = c(nTags, nlibs))
		object$Lambda <- Lambda
		if(!is.na(fix.dispersion))
		Dispersion <- expandAsMatrix(fix.dispersion, dim = c(nTags, nlibs))
		else Dispersion <- expandAsMatrix(Dispersion, dim = c(nTags, nlibs))
		object$Dispersion <- Dispersion
		object
		
	}
	diff.fun <- function(object)
	{ 
		
        ## it is low-level function of NBsim ##
        ## it creates diff genes according to foldDiff ## 
		group <- object$group
		pDiff <- object$pDiff
		pUp <-  object$pUp 
		foldDiff <- object$foldDiff
		Lambda <- object$Lambda
		nTags <- object$nTags
		g <- group == levels(group)[1]
		## added by Koen to specify DE index yourself and allows to specify foldDiff as matrix
		if(is.null(ind)) ind <- sample(nTags, floor(pDiff*nTags))
		##
		if(length(ind)>0 & !mean(foldDiff==1)==1 ) {
			fcDir <- sample(c(-1,1), length(ind), prob=c(1-pUp,pUp), replace=TRUE)
			Lambda[ind,g] <- Lambda[ind,g]*exp(log(foldDiff)/2*fcDir)
			Lambda[ind,!g] <- Lambda[ind,!g]*exp(log(foldDiff)/2*(-fcDir)) 
            #Lambda <- t(t(Lambda)/colSums(Lambda))
			object$Lambda <- Lambda
			object$indDE <- ind
			object$indNonDE <- (1:nTags)[-ind]
			object$mask_DEup <- object$mask_DEdown <- object$mask_nonDE <- expandAsMatrix(FALSE, dim = dim(Lambda))
			object$mask_DEup[ind[fcDir == 1], g] <- TRUE
			object$mask_DEup[ind[fcDir == -1], !g] <- TRUE
			object$mask_DEdown[ind[fcDir == -1], g] <- TRUE
			object$mask_DEdown[ind[fcDir == 1], !g] <- TRUE
			object$mask_nonDE[-ind,] <- TRUE
			object$mask_DE <- object$mask_DEup | object$mask_DEdown}
		if(mean(foldDiff==1)==1 | pDiff == 0)
		object$indDE <- NA
		object
	}
	sim.fun <- function(object)
	{   
        ## it is low-level function of NBsim ##
        ## it simulate counts using rnbinom ## 
		Lambda <- object$Lambda
		Dispersion <- object$Dispersion
		nTags <- object$nTags
		nlibs <- object$nlibs
		lib.size <- object$lib.size
		counts <- matrix(rnbinom(nTags*nlibs, mu = t(t(Lambda)*lib.size), size = 1/Dispersion), nrow = nTags, ncol = nlibs) 
		rownames(counts) <- paste("ids", 1:nTags, sep = "")
		object$counts <- counts
		object
	}
			
	outlier.fun <- function(object, outlierMech, pOutlier, min.factor = 2, max.factor = 5)
        {   
        ## it is low-level function of NBsim ##
        ## it makes outlier ## 
	        outlierMech <- match.arg(outlierMech, c("S", "M", "R"))
	        dim <- dim(object$counts)
                outlier.factor <- function() runif(1, min.factor, max.factor)
                countAddOut <- object$counts
                LambdaAddOut <- object$Lambda
	        DispersionAddOut <- object$Dispersion	
	        switch(outlierMech, 
	               S = {
	                    mask_outlier <- expandAsMatrix(FALSE, dim = dim)
	                    id_r <- which(runif(dim[1]) < pOutlier)
	                    id_c <- sample(dim[2], length(id_r), replace = TRUE)
		            for(i in seq(id_r)) 
		                   mask_outlier[id_r[i], id_c[i]] <- TRUE
                            countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
                            }, 
	               R = {				   
                            mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
                            countAddOut[mask_outlier] <- sapply(countAddOut[mask_outlier], function(z) round(z*outlier.factor()))
                            },
         
	               M = {
                            mask_outlier <- matrix(runif(dim[1]*dim[2]) < pOutlier, dim[1], dim[2])
                            LambdaAddOut[mask_outlier] <- sapply(LambdaAddOut[mask_outlier], function(z) z*outlier.factor())
	                    countAddOut[mask_outlier] <- rnbinom(sum(mask_outlier), mu = t(t(LambdaAddOut)*object$lib.size)[mask_outlier], size = 1/DispersionAddOut[mask_outlier])
                           }
                       )
                 if(!mean(object$foldDiff == 1)==1 & !pDiff == 0)
	         {
	                indDEupOutlier <- which(apply(object$mask_DEup & mask_outlier, 1, any))
	                indDEdownOutlier <- which(apply(object$mask_DEdown & mask_outlier, 1, any))
	                indDEnoOutlier <- which(apply((object$mask_DE & !mask_outlier) , 1, all))
	                indNonDEOutlier <- which(apply(object$mask_nonDE & mask_outlier, 1, any))
	                indNonDEnoOutlier <- which(apply((object$mask_nonDE & !mask_outlier) , 1, all))
	                indDEbothOutlier <- NA
	                o <- indDEupOutlier %in% indDEdownOutlier
	                q <-  indDEdownOutlier %in% indDEupOutlier
	                if(any(o))
	                {
                              indDEupOutlier <- indDEupOutlier[!o]
		              indDEbothOutlier <- indDEupOutlier[o]	
	                      indDEdownOutlier <- indDEdownOutlier[!q]	
	                }	
	         }
	         else
	         {
                        indDEupOutlier <- indDEdownOutlier <- indDEnoOutlier <- indNonDEOutlier <- indNonDEnoOutlier <- indDEbothOutlier <- NA
                 }
	             out <- list(countAddOut = countAddOut, outlierMech = outlierMech, pOutlier = pOutlier, mask_outlier = mask_outlier, indDEupOutlier = indDEupOutlier, 
                                 indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, 
                                 indNonDEnoOutlier = indNonDEnoOutlier, LambdaAddOut = LambdaAddOut, DispersionAddOut = DispersionAddOut) 

        }
	
	calProb <- function(x, l) round(1 -(1 - x)^(1/l), 2) ## calculate probability to make sure all the outlierMech produce the same amount of outliers ##
	##### hlp = as.matrix(islamFilt)

        if(verbose) message("Preparing dataset.\n")	
	if(class(dataset) == "DGEList")
	{   
		dat <- dataset
		dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
	}
	else if(is.character(dataset)) 
	{
		load(dataset)
		dat <- get(gsub("(\\.)(\\w+)", "", basename(dataset)))
		dat[["R"]] <- dat[["S"]] <- dat[["M"]] <- dat[["pOutlier"]] <- dat[["outlierMech"]]<- NULL
	}
	else if(is.matrix(dataset)) 
	{ 
	  if(is.null(name)) name <- deparse(substitute(dataset))	
	  dataset <- getDataset(counts =dataset, drop.extreme.dispersion = drop.extreme.dispersion, drop.low.lambda = drop.low.lambda)
	  dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))		
	}
	else
	dat <- new("DGEList", list(dataset = dataset, nTags = nTags, lib.size = lib.size, nlibs = nlibs, group = group, design = model.matrix(~group), pDiff= pDiff, pUp = pUp, foldDiff = foldDiff, outlierMech = outlierMech, min.factor = min.factor, max.factor = max.factor, name = name))

	if(!only.add.outlier)
	{
		if(is.null(lib.size)){
		    dat$lib.size <- runif(nlibs, min = 0.7*median(dat$dataset$dataset.lib.size), max = 1.3*median(dat$dataset$dataset.lib.size))
		    #propZeroFit=dat$dataset.propZeroFit
	
		}
	    if(is.null(nTags))
	      dat$nTags <- dat$dataset$dataset.nTags 
        if(verbose) message("Sampling.\n")	
	      dat <- sample.fun(dat)
        if(verbose) message("Calculating differential expression.\n")	
	      dat <- diff.fun(dat)
        if(verbose) message("Simulating data.\n")	
	      dat <- sim.fun(dat)
	}
	if(add.outlier){
		outlierMech <- match.arg(outlierMech,  c("S", "R", "M"), several.ok = TRUE)
		if(length(pOutlier)== 1 & length(outlierMech) > 1 & any(outlierMech == "S"))
		{ 
			prob <- calProb(pOutlier, length(group))	
			pOutlier <- rep(pOutlier, length = length(outlierMech))
			pOutlier[!outlierMech == "S"] <- prob	
		}
		else if(!length(pOutlier) == length(outlierMech))
		stop("pOutlier is not equal to outlierMech")
                if(verbose) message("Adding outliers.\n")	
		dat[outlierMech] <- mapply(function(x, y) outlier.fun(dat, outlierMech = x, pOutlier = y, min.factor = min.factor, max.factor = max.factor), x = outlierMech, y = pOutlier, SIMPLIFY = FALSE)
	    dat$pOutlier <- pOutlier
	}
	if(save.file)
	{ 
		
		## save file for shiny app ##
		if(verbose) message("Saving file.\n")
		if(is.null(file)) 
		{ g <- paste0("g", sum(levels(group)[1] == group), "v", sum(levels(group)[2] == group))
			f <- paste0("f", foldDiff)
			if(add.outlier) o <- paste0("o", sprintf( "%02d",100*pOutlier[1L])) 
			else o <- paste0("o", sprintf( "%02d", 0 ))  
			file <- paste0(dat$name, "/", g, f, o, ".Rdata")
			dir.create(dat$name, showWarnings = FALSE)  
		}
		filenm <- eval(gsub("(\\.)(\\w+)", "", basename(file)))
		assign(filenm, dat)
		save(list = filenm, file = file)		
	}
	dat 	
}





pval <-
function(y, ...) ## evaluate DE methods ##
UseMethod("pval")
pval.default <-
function(y, group, design = NULL, method = "edgeR", mc.cores = 4, globalEnvir = FALSE, niter=NULL, ...)
{   
	## evaluate DE methods ##
    ## return to a list of pvalue and runing time ##
    ## pvalue contains pvalue and p-adjust value ##
	gc(FALSE)
	time <- proc.time()
	group <- as.factor(group)
	if(globalEnvir) method <- paste0(method, ".pscript")
	else method <- paste0(method, ".pfun")
	p <- get(method)
	if(globalEnvir)
	{
		L <- list(counts = y, group = group, design = design, mc.cores = mc.cores, p = p)
		e <- list2env(L, envir = .GlobalEnv)
		pvalue <- with(e, eval(p))
		try(rm(list = names(L), envir = e),  silent = TRUE)
		try(rm(pGlobal, envir = e),  silent = TRUE)
	}	
	else pvalue <- p(y, group, design, mc.cores, niter, ...)
	pvalue
	new.time <- proc.time()
	output <- list(pvalue = pvalue, time = new.time - time)
}



pval.DGEList <-
function(y, method, mc.cores = 4, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), save.file = FALSE, name = deparse(substitute(y)), count.type="counts", niter=NULL)
{   
	## evaluate DE methods ##
    ## return to a DGEList including pvalue and other necessary indexs for re-analysis and plot ## 
	library(parallel)	   
	counts = y$counts
	 pOutlier = mask_outlier = indDEupOutlier = indDEdownOutlier = indDEbothOutlier = indDEnoOutlier = indNonDEOutlier = indNonDEnoOutlier = NA
	names(method) <- method 
	group <- y$group
	design <- y$design
	is.parallel <- method %in% parallel.method
	is.globalEnvir <- method %in% globalEnvir.method
	id.re <- !(is.parallel|is.globalEnvir)
	reduced.method <- method[id.re]
	if(any(id.re)) output <-  parallel:::mclapply(reduced.method, function(x) pval(y = counts, group = group, design = design, method = x, niter=niter), mc.cores = mc.cores, mc.preschedule = FALSE) else output <- list()
	if(any(is.parallel))
	{
	  for( i in names(method[is.parallel]))
         output[[i]] <- pval(y = counts, group = group, design = design, method = i, mc.cores = mc.cores, niter=niter)
	}
	if(any(is.globalEnvir))
	{  		
		for( i in names(method[is.globalEnvir]))
		output[[i]] <- pval(y = counts, group = group, design = design, method = i, mc.cores = mc.cores, globalEnvir = TRUE, niter=niter)
	}
	output <- output[method]
	padj <- lapply(output, function(x) x[["pvalue"]][, "padj"])
	pval <- lapply(output, function(x) x[["pvalue"]][, "pval"])
	lfc <- lapply(output, function(x) x[["pvalue"]][, "lfc"])
        time <- lapply(output, function(x) x[["time"]])
	output <- new("DGEList", list(pval = pval, padj = padj, lfc = lfc,  counts = counts, group = group, design = design, indDE = y$indDE, method = names(method), indDEupOutlier = indDEupOutlier, indDEdownOutlier = indDEdownOutlier, indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, indNonDEnoOutlier = indNonDEnoOutlier, time = time)) 
	output$main <- mainShow(count.name = y$name, group = group, pOutlier = pOutlier, count.type=count.type)
	output$methodVersion <- getVersion(method)
	output
}

pval.character <-
function(y, method, count.type = "counts", mc.cores = 6, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), save.file = FALSE, niter=NULL)
{   
    ## for shiny app ##
	fnm <- y
	load(y)
	name <- gsub("(\\.)(\\w+)", "", basename(y))
	y <- get(name)
	pval.DGEList(y = y, method = method, count.type = count.type, mc.cores = mc.cores, parallel.method = parallel.method, globalEnvir.method = globalEnvir.method, save.file = save.file, name = fnm, niter=niter)
}
pval.FoldList <-
function(y, method, count.type = "counts", mc.cores = 6, parallel.method = c("baySeq"), globalEnvir.method = c("ShrinkBayes"), cut.computing = TRUE)
{   
    ## evaluate DE methods for FoldList ##
	library(parallel)
	count.type <- match.arg(count.type, c("counts", "S", "R", "M"))
	if(count.type == "counts") 
	{counts = y$counts
		pOutlier = mask_outlier = indDEupOutlier = indDEdownOutlier = indDEbothOutlier = indDEnoOutlier = indNonDEOutlier = indNonDEnoOutlier = NA}
	else 
	{counts = y[[count.type]]$countAddOut
		pOutlier = y[[count.type]]$pOutlier
		mask_outlier = y[[count.type]]$mask_outlier
		indDEupOutlier = y[[count.type]]$indDEupOutlier
		indDEdownOutlier = y[[count.type]]$indDEdownOutlier
		indDEbothOutlier = y[[count.type]]$indDEbothOutlier	
		indDEnoOutlier = y[[count.type]]$indDEnoOutlier
		indNonDEnoOutlier = y[[count.type]]$indNonDEnoOutlier
		indNonDEOutlier = y[[count.type]]$indNonDEOutlier}
	
	names(method) <- method 
	group <- y$group
	design <- y$design
	is.parallel <- method %in% parallel.method
	is.globalEnvir <- method %in% globalEnvir.method
	id.re <- !(is.parallel|is.globalEnvir)
	reduced.method <- method[id.re] 
	if(any(id.re)) output <- lapply(reduced.method, function(x) parallel:::mclapply(counts, function(w) pval(y = w, method = x, group = group, design = design), mc.cores = mc.cores))
	else output <- list()
	fold_seq <- fold_seq.keep <- y$fold_seq
	if(cut.computing) fold_seq.keep <- odd(fold_seq)
	if(any(is.parallel))
	{
		for(j in fold_seq)
		{ 
		  is.keep <- j %in% fold_seq.keep
		  for( i in names(method[is.parallel]))
			{
				if(any(is.keep)) output[[i]][[j]] <- pval(y = counts[[j]], group = group, design = design, method = i, mc.cores = mc.cores)
				else 
				{
					output[[i]][[j]][["pavlue"]] <- cbind(pval = NA, padj = NA)
					output[[i]][[j]][["time"]] <- NA
				}
			}	
		}
	}
	if(any(is.globalEnvir))
	{
		for(j in fold_seq)
		{ 
			is.keep <- j %in% fold_seq.keep
			for( i in names(method[is.globalEnvir]))
			{
				if(any(is.keep)) output[[i]][[j]] <- pval(y = counts[[j]], group = group, design = design, method = i, mc.cores = mc.cores, globalEnvir = TRUE)
				else 
				{
					output[[i]][[j]][["pavlue"]] <- cbind(pval = NA, padj = NA)
					output[[i]][[j]][["time"]] <- NA
				}
			}	
		}
	}
	output <- output[method]
	padj <- try(lapply(output, lapply, function(x) x[["pvalue"]][, "padj"]), silent = TRUE)
	pval <- try(lapply(output, lapply, function(x) x[["pvalue"]][, "pval"]), silent = TRUE)
	time <- try(lapply(output, lapply, function(x) x[["time"]]), silent = TRUE)
	output <- new("FoldList", list(fold_seq = y$fold_seq, pval = pval, padj = padj, counts = counts, count.type = count.type, group = group, design = design, indDE = y$indDE, method = names(method), indDEupOutlier = indDEupOutlier, indDEdownOutlier = indDEdownOutlier,indDEbothOutlier = indDEbothOutlier, indDEnoOutlier = indDEnoOutlier, indNonDEOutlier = indNonDEOutlier, indNonDEnoOutlier = indNonDEnoOutlier, time = time))
	output$main <- mainShow(count.type = count.type, count.name = y$name, group = group, pOutlier = pOutlier)
	output$methodVersion <- getVersion(method)
	output
}
getPvalVersion <- function(methodVersion, pout = "pval", count.type = "counts", datanm)
{ 
  ## for shiny app ##
  Type <- switch(count.type, counts = "b", S = "s", M = "m", R = "r")
  filenm <- paste0(pout, "_", Type, "_",  basename(datanm), "_", methodVersion, ".Rdata")
  load(paste0(dirname(datanm),"/", rmVersion(methodVersion),"/", filenm))
  get(methodVersion)	   
}	
getPval <- function(dataset,methodVersion, count.type = c("counts", "S", "R", "M"))
{
    ## for shiny app ##
	load(dataset)
	datanm <- gsub("(\\.)(\\w+)", "", dataset)
	dat <- get(basename(datanm))
	count.type <- match.arg(count.type, c("counts", "S", "R", "M"))
	if(count.type == "counts") Dat <- new("DGEList", dat)
	else 
	{
		Dat <- new("DGEList", dat[[count.type]])
		Dat[["counts"]] <- Dat[["countAddOut"]]
	}
	Dat$method <- Dat$methodVersion <- methodVersion
	Dat$group = dat$group
	Dat$indDE = dat$indDE
	Dat$name = dat$name
	Dat$main <- mainShow(count.type = count.type, count.name = Dat$name, group = Dat$group, pOutlier = Dat$pOutlier)
	index <- c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier")
	names(index) <- index
	indDiff <- lapply(index, function(x) Dat[[x]])
	indDiff <- indDiff[!sapply(indDiff, is.null)]
	indDiff <- indDiff[!is.na(indDiff)]
	Dat$index <- names(indDiff)
	names(methodVersion) <- methodVersion
	Dat[["padj"]] <- lapply(methodVersion, getPvalVersion, pout = "padj", count.type = count.type, datanm = datanm)
	Dat[["pval"]] <- lapply(methodVersion, getPvalVersion, pout = "pval", count.type = count.type, datanm = datanm)
	Dat 
	
}
	

resetPar <- function() {
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    op
}

roPlot <-
function(y, ...) 
UseMethod("roPlot")
## plot ROC curve ##

roPlot.default <-
function(y, indDiff, returnData=FALSE, plot.max.fpr = 0.4, add = FALSE, cex.axis = 2, threshold = 0.05, col = 1, cex.threshold = 3, plot.max.tpr=1, ...)
{   
    ## plot ROC curve ##
	#old.par <- par(c("mar", "mgp", "cex.axis"))
	#par(mar=c(4,5,3,2))
        #par(mgp = c(2.6, 1, 0))
	#par(cex.axis = cex.axis)
	#on.exit(par(old.par))
	library(ROCR)
	if(any(is.na(y)))
	{
      y[is.na(y)] <- 1
	}	
	y = 1 - y
	label <- as.factor(rep("nonDE", length(y)))
	levels(label) <- c("nonDE", "DE")
	label[indDiff] <- "DE"
	pred <- prediction(y, label, label.ordering = c("nonDE", "DE"))
	perf <- performance(pred, "tpr", "fpr")
	if(is.null(plot.max.fpr))
	plot.max.fpr <- 1
	plot(perf, xlim = c(0, plot.max.fpr), ylim=c(0,plot.max.tpr), col = col, add = add, ...)
	if(!is.null(threshold))
	{
		fpr <- approx(y = perf@x.values[[1]], x = perf@alpha.values[[1]], xout = 1- threshold)$y
		tpr <- approx(y = perf@y.values[[1]], x = perf@x.values[[1]], xout = fpr)$y
		points(x = fpr, y = tpr, pch = 4, col = col, cex = cex.threshold, ...)
	}
	#if(returnData) return(perf) #added by Koen Vdb
	#par(resetPar())
}
roPlot.DGEList <-
function(y, plot.max.fpr = 0.4, plot.max.tpr=1, pout = "padj", threshold = 0.05, selected.method = NULL, show.legend = TRUE, short.main = FALSE, col = NULL, lty = 1, lwd = 5, box.lwd = 1.5, cex.main = 2.5, cex.axis=2, cex.lab = 2.1, cex = 1.6, cex.threshold = 8)
{   
	## plot ROC curve for DGEList ##
	library(ROCR)
	main <- y$main
	if(short.main) 
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")

	}
	if(is.null(selected.method)) 
	{
		method <- y$method
		methodVersion <- y$methodVersion
	}
	else 
	{
		method <- match.arg(selected.method, y$method, several.ok = TRUE)
		methodVersion <- y$methodVersion[match(method, y$method)]
	}
	pout <- match.arg(pout, c("pval", "padj"), several.ok = TRUE)
	pvalue <- y[[pout]][method]
	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow", "steelblue", "salmon", "violetred4", "darkred", "skyblue4")
	if(is.null(col)) col <- pre.col[seq(method)]		
	roPlot(y = pvalue[[1L]], indDiff = y$indDE, threshold = threshold, plot.max.fpr = plot.max.fpr, plot.max.tpr = plot.max.tpr, col = col[1L], main = main, lty = lty, lwd = lwd, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, cex.threshold = cex.threshold)
	if(is.list(pvalue)) mapply(function(u, v, w) roPlot(y = u, indDiff = y$indDE, threshold = threshold, plot.max.fpr = plot.max.fpr, plot.max.tpr = plot.max.tpr, col = v, lwd = lwd, lty = lty, cex.main = cex.main, cex.axis=cex.axis, cex.lab = cex.lab, cex.threshold = cex.threshold, add = TRUE), u = pvalue[-1L], v = col[-1L])
	if(show.legend) legend("bottomright", methodVersion, col = col, lty = lty, lwd = lwd, box.lwd = box.lwd, cex = cex)	
	if(!is.null(threshold)&show.legend) legend("topleft", paste0("P_adj_value=", threshold), col = "black", pch = 4, lty = NA, lwd = lwd, box.lwd = box.lwd, cex = cex, pt.cex = 1.5*cex)	
}



fdPlot <-
function(y, ...) 
UseMethod("fdPlot")
## plot False discovery number curve ##

fdPlot.default <-
function(y, indDiff, add=FALSE, xlab="Number of genes selected", 
ylab="Number of false discoveries", lwd=4, type="l", ... ) 
{
	
	## plot False discovery number curve ##	
#	old.par <- par(c("mar", "mgp"))
#	par(mar=c(4,5,3,2))
#    par(mgp = c(2.6, 1, 0))
#	on.exit(par(old.par))	
	x <- 1:length(indDiff)
	o <- order(y)
	w <- !o[x] %in% indDiff
	y1 <- cumsum(w)
	matplot(x, y1, xlab=xlab, ylab=ylab, lwd=lwd, type=type, add=add, ...)
}
fdPlot.DGEList <-
function(y, pout = "padj", selected.method = NULL, short.main = FALSE, show.legend = TRUE, col = NULL, lty = 1, lwd = 5, box.lwd = 1.5, cex.main = 2.5, cex.axis=2, cex.lab = 2.1, cex = 1.3, xlim = NULL)
{   
	## plot False discovery number curve for DGEList ## 
	main <- y$main
	if(short.main) 
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
	}
	if(is.null(selected.method)) 
	{
		method <- y$method
		methodVersion <- y$methodVersion
	}
	else 
	{
		method <- match.arg(selected.method, y$method, several.ok = TRUE)
		methodVersion <- y$methodVersion[match(method, y$method)]
	}
	pout <- match.arg(pout, c("pval", "padj"), several.ok = TRUE)
	pvalue <- y[[pout]][method]

	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]	
	fdPlot(y = pvalue[[1L]], indDiff = y$indDE,  log="y", col = col[1L], main = main, lty = lty, lwd = lwd, cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, xlim = xlim)
	if(is.list(pvalue)) mapply(function(u, v, w) fdPlot(y = u, indDiff = y$indDE, log="y", col = v, lwd = lwd, lty = lty, cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab, add = TRUE), u = pvalue[-1L], v = col[-1L])
	if(show.legend) legend("bottomright", methodVersion, col = col, lty = lty, lwd = lwd, box.lwd = box.lwd, cex = cex)
}

getPower <-
function(y, ...) 
UseMethod("getPower")
## plot Power curve ##
getPower.default <-
function(y, indDiff, threshold)
{ 
	## plot Power curve ##
	if(all(is.na(y)))
    NA
	else if(all(is.na(indDiff)))
    contable(score = y, indDiff = indDiff, threshold = threshold, output = "fpr")
	else
    contable(score = y, indDiff = indDiff, threshold = threshold, output = "power")
}
getPower.DGEList <-
function(y, pout = "padj", threshold = 0.05, index = c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"), byAveLogCPM = FALSE, cutCPM = 4, selected.method = NULL, short.main = FALSE, plot = FALSE, col = NULL, cex.main = 2.5, cex.axis=2, cex.sub = 1.5, cex.lab = 2.1, ylim = NULL, ...)
{
	## plot Power curve for DGEList ##
	index <- match.arg(index, c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"))
	main <- y$main
	if(short.main)
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
		
	}
    
    if(is.null(selected.method))
    {
        method <- y$method
        methodVersion <- y$methodVersion
    }
    else
    {
        method <- match.arg(selected.method, y$method, several.ok = TRUE)
        methodVersion <- y$methodVersion[match(method, y$method)]
    }
	pout <- match.arg(pout, c("pval", "padj"))
	pvalue <- y[[pout]][method]
	indDiff <- y[[index]]
	
	if(byAveLogCPM)
    {
        threshold <- rep(threshold, cutCPM)
        diff <- rep(FALSE, nrow(y$counts))
        diff[indDiff] <- TRUE
        d <- DGEList(counts = y$counts, group = y$group)
        d <- edgeR::calcNormFactors(d)
        AveLogCPM <- aveLogCPM(d)
        o <- order(AveLogCPM)
        l <- length(o)/cutCPM
        oo <- split(o, ceiling(seq_along(o)/l))
        AveLogCPM.list <- lapply(oo, function(x) AveLogCPM[x])
        nm <- lapply(AveLogCPM.list, function(x) round(range(x), 2))
        nm <- unlist(lapply(nm, function(x) paste0("(", paste0(x, collapse = ","), "]")))
        
        power <- list()
        power[["all"]] <- unlist(lapply(pvalue, getPower, indDiff = indDiff, threshold = threshold[1]))
        for( i in 1:length(oo))
        {
            diff_idx <- which(diff[oo[[i]]])
            pvalue_idx <- lapply(pvalue, function(x) x[oo[[i]]])
            power[[nm[i]]] <- unlist(lapply(pvalue_idx, getPower, indDiff = diff_idx, threshold = threshold[i]))
        }
        power <- do.call("rbind", power)
	}
    else
    {
        power <- unlist(lapply(pvalue, getPower, indDiff = indDiff, threshold = threshold))
	}
	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]
    out <- list()
    out$power <- power
	out$index <- gsub("ind", "", index)
    out$main <- main
    out$method <- method
	out$methodVersion <- methodVersion
    out$fold_seq <- y$fold_seq
    if(plot)
    if(byAveLogCPM) matPlot(out, output = "power", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    else powBarPlot(out, col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
	out
}



getFDR <-
function(y, ...)
UseMethod("getFDR")
## plot Power curve ##
getFDR.default <-
function(y, indDiff, threshold)
{
	## plot Power curve ##
	if(all(is.na(y)))
    NA
	else if(all(is.na(indDiff)))
    NA
	else
    contable(score = y, indDiff = indDiff, threshold = threshold, output = "fdr")
}
getFDR.DGEList <-
function(y, pout = "padj", threshold = 0.05, index = c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"), byAveLogCPM = FALSE, cutCPM = 4, selected.method = NULL, short.main = FALSE, plot = FALSE, col = NULL, cex.main = 2.5, cex.axis=2, cex.sub = 1.5, cex.lab = 2.1, ylim = NULL, ...)
{
	## plot Power curve for DGEList ##
	index <- match.arg(index, c("indDE", "indDEupOutlier", "indDEdownOutlier", "indDEbothOutlier", "indDEnoOutlier"))
	main <- y$main
	if(short.main)
	{
		main <- strsplit(main, "/")[[1]]
		main <- main[-c(length(main), length(main)-1)]
		main <- paste0(main, collapse = "/")
		
	}
    
    if(is.null(selected.method))
    {
        method <- y$method
        methodVersion <- y$methodVersion
    }
    else
    {
        method <- match.arg(selected.method, y$method, several.ok = TRUE)
        methodVersion <- y$methodVersion[match(method, y$method)]
    }
	pout <- match.arg(pout, c("pval", "padj"))
	pvalue <- y[[pout]][method]
	indDiff <- y[[index]]
	
	if(byAveLogCPM)
	{
        threshold <- rep(threshold, cutCPM)
        diff <- rep(FALSE, nrow(y$counts))
        diff[indDiff] <- TRUE
        d <- DGEList(counts = y$counts, group = y$group)
        d <- edgeR::calcNormFactors(d)
        AveLogCPM <- aveLogCPM(d)
        o <- order(AveLogCPM)
        l <- length(o)/cutCPM
        oo <- split(o, ceiling(seq_along(o)/l))
        AveLogCPM.list <- lapply(oo, function(x) AveLogCPM[x])
        nm <- lapply(AveLogCPM.list, function(x) round(range(x), 2))
        nm <- unlist(lapply(nm, function(x) paste0("(", paste0(x, collapse = ","), "]")))
        
        fdr <- list()
        fdr[["all"]] <- unlist(lapply(pvalue, getFDR, indDiff = indDiff, threshold = threshold[1]))
        for( i in 1:length(oo))
        {
            diff_idx <- which(diff[oo[[i]]])
            pvalue_idx <- lapply(pvalue, function(x) x[oo[[i]]])
            fdr[[nm[i]]] <- unlist(lapply(pvalue_idx, getFDR, indDiff = diff_idx, threshold = threshold[i]))
        }
        fdr <- do.call("rbind", fdr)
	}
    else
    {
        fdr <- unlist(lapply(pvalue, getFDR, indDiff = indDiff, threshold = threshold))
	}
	pre.col <- c("black", "blue", "purple", "gray", "tan3", "red", "green", "powderblue", "chartreuse4", "yellow")
	if(is.null(col)) col <- pre.col[seq(method)]
    out <- list()
    out$fdr <- fdr
	out$index <- gsub("ind", "", index)
    out$main <- main
    out$method <- method
	out$methodVersion <- methodVersion
    out$fold_seq <- y$fold_seq
    if(plot)
    if(byAveLogCPM) matPlot(out, output = "fdr", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
    else barPlot(out, output = "fdr", col = col, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, cex.sub = cex.sub, ylim = ylim, ...)
	out
}





