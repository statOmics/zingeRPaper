
runEdgeR <- function(e) {
  library(edgeR)
  condition = pData(e)$condition
  pickingSession = pData(e)$pickingSession
  design <- model.matrix(~ condition + pickingSession)
  dgel <- DGEList(exprs(e))
  dgel <- edgeR::calcNormFactors(dgel)
  dgel=estimateDisp(dgel,design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit, coef="conditionB")
  pvals <- edger.lrt$table$PValue
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj)
}

runEdgeRRobust <- function(e) {
    library(edgeR)
  condition = pData(e)$condition
  pickingSession = pData(e)$pickingSession
  design <- model.matrix(~ condition + pickingSession)
  dgel <- DGEList(exprs(e))
  dgel <- edgeR::calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit, coef="conditionB")
  pvals <- edger.lrt$table$PValue
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj)
}


runVoom <- function(e) {
  library(limma)
  condition = pData(e)$condition
  pickingSession = pData(e)$pickingSession
  design <- model.matrix(~ condition + pickingSession)
  dgel <- DGEList(exprs(e))
  dgel <- edgeR::calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef="conditionB",n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value 
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=tt$logFC)
}

runDESeq2 <- function(e, retDDS=FALSE) {
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(exprs(e), colData=DataFrame(pData(e)), design=~ condition + pickingSession)
  #dds <- DESeq(dds, betaPrior=TRUE, quiet=TRUE, modelMatrixType="standard")
  dds = estimateSizeFactors(dds)
  dds = estimateDispersions(dds)
  dds = nbinomWaldTest(dds, betaPrior=TRUE, modelMatrixType="standard")
  res <- results(dds, name="condition_B_vs_A")
  pvals <- res$pvalue
  padj <- res$padj
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

runDESeq2_poscounts <- function(e, retDDS=FALSE) {
    library(DESeq2)
  dds <- DESeqDataSetFromMatrix(exprs(e), colData=DataFrame(pData(e)), design=~ condition + pickingSession)
  dds <- estimateSizeFactors(dds,type="poscounts")
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds, betaPrior=TRUE, modelMatrixType="standard")
  res <- results(dds, name="condition_B_vs_A")
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

zeroWeightsLibSizeDispFast <- function(counts, design, colData=NULL, initialWeightAt0=TRUE, maxit=100, plot=FALSE, plotW=FALSE, designZI=NULL, llTol=1e-4, normalization="TMM", llOffset=1e-6){
    require(edgeR) ; require(DESeq2)
    if(plot | plotW) par(mfrow=c(1,plot+plotW))

    if(normalization=="TMM"){
      counts <- DGEList(counts)
      counts = edgeR::calcNormFactors(counts)
    } else if(normalization=="DESeq2"){
      designFormula=as.formula(paste0("~",paste(names(attr(design,"contrasts")),collapse="+")))
dse = DESeqDataSetFromMatrix(counts, colData=colData, design=designFormula)
       dse = DESeq2::estimateSizeFactors(dse)
       counts <- DGEList(counts)
       counts$samples$norm.factors = 1/dse$sizeFactor
    } else if(normalization=="DESeq2_pos"){
      designFormula=as.formula(paste0("~",paste(names(attr(design,"contrasts")),collapse="+")))
      dse = DESeqDataSetFromMatrix(counts, colData=colData, design=designFormula)
       dse = DESeq2::estimateSizeFactors(dse, type="poscounts")
       counts <- DGEList(counts)
       counts$samples$norm.factors = 1/dse$sizeFactor
    }

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
  #counts <- estimateGLMCommonDisp(counts, design, interval=c(0,10))
  #counts <- estimateGLMTagwiseDisp(counts, design, prior.df=0, min.row.sum=1)
  counts = estimateDisp(counts, design, prior.df=0, min.row.sum=1)
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
  ll <- log(pi0HatMat*zeroId + (1-pi0HatMat)*likC + llOffset)

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

source("~/Dropbox/phdKoen/singleCell/githubPaper/singleCellPaper/method/glmLRTOld.R")

runEdgeREMLibSize=function(e){
    #function(counts, group, design=NULL, mc.cores=2, niter=50){
  library(edgeR) ; library(genefilter)
  condition = pData(e)$condition
  pickingSession = pData(e)$pickingSession
  design <- model.matrix(~ condition + pickingSession)	
	d <- DGEList(exprs(e))
	d <- edgeR::calcNormFactors(d)
	#not adding a design matrix models the zeroes with the library size automatically
  effLogLibSize = log(d$samples$lib.size*d$samples$norm.factors)
  pickingSession = pData(e)[,"Picking sessions"]
  designZI = model.matrix(~effLogLibSize + pickingSession)
	zeroWeights = zeroWeightsLibSizeDispFast(d, design, plot=FALSE, maxit=200, initialWeightAt0=TRUE, plotW=FALSE, designZI=designZI)
	d$weights = zeroWeights
	d=estimateDisp(d,design)
	#plotBCV(d)
	edger.fit <- glmFit(d, design) #uses weights
	edger.fit$df.residual <- rowSums(edger.fit$weights)-ncol(design)
  	edger.lrt <- glmLRTOld(edger.fit,coef=2,test="F")
  	pvals <- edger.lrt$table$PValue
	baseMean = unname(rowMeans(sweep(d$counts,2,d$samples$norm.factors,FUN="*")))	
	hlp <- pvalueAdjustment_kvdb(baseMean=baseMean, pValue=pvals)
  	padj <- hlp$padj
	#padj <- p.adjust(pval,method="BH")
  	padj[is.na(padj)] <- 1
	out=list(pvals=pvals,padj=padj)
	out[is.na(out)] <- 1
	return(out)
}



runScde <- function(e){
    require(scde)
    # calculate models
    o.ifm <- scde.error.models(counts = exprs(e), groups = pData(e)$condition, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 0)

    # filter out cells that don't show positive correlation with
    # the expected expression magnitudes (very poor fits)
    #valid.cells <- o.ifm$corr.a > 0
    #table(valid.cells)
    #o.ifm <- o.ifm[valid.cells, ]

    # estimate gene expression prior
    o.prior <- scde.expression.prior(models = o.ifm, counts = exprs(e), length.out = 400, show.plot = FALSE)
    # run differential expression tests on all genes.
    ediff <- scde.expression.difference(o.ifm, exprs(e), o.prior, groups  =  pData(e)$condition, n.randomizations  =  100, n.cores  =  1, verbose  =  0, batch=pData(e)$pickingSession)
    pvals=(1-pnorm(abs(ediff$batch.adjusted$Z)))*2
    padj=p.adjust(pvals,method="BH")
    list(Z=ediff$Z,pvals=pvals,padj=padj)
}


runMAST <- function(e){
    require(MAST)
    counts=exprs(e)
    tpm <- counts*1e6/colSums(counts)
    sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(group=pData(e)$condition, pickingSession=pData(e)$pickingSession))
    ngeneson <- apply(exprs(e),2,function(x) mean(x>0))
    CD <- cData(sca)
    CD$ngeneson <- ngeneson
    CD$cngeneson <- CD$ngeneson-mean(ngeneson)
    cData(sca) <- CD
     ## differential expression 
    fit <- zlm.SingleCellAssay(~cngeneson+group+pickingSession, sca=sca, method="bayesglm", ebayes=TRUE)
    L=matrix(0,nrow=ncol(coef(fit,"D")))
    rownames(L)=colnames(coef(fit,"D"))
    L["groupB",]=1
    lrFit <- lrTest(fit, hypothesis=L)
    pval=lrFit[,'hurdle','Pr(>Chisq)']
    padj=p.adjust(pval,method="BH")
    list(pvals=pval,padj=padj)
}


runMetagenomeSeq <- function(e){
  require(metagenomeSeq)
  condition = pData(e)$condition
  pickingSession = pData(e)$pickingSession
  design <- model.matrix(~ condition + pickingSession)    
  pheno <- AnnotatedDataFrame(data.frame(group=condition, pickingSession=pickingSession))
  rownames(pheno) <- colnames(exprs(e))
  p <- cumNormStatFast(exprs(e))
  dat <- newMRexperiment(counts=exprs(e), phenoData=pheno, featureData = NULL, libSize = colSums(exprs(e)), normFactors = metagenomeSeq::calcNormFactors(exprs(e), p=p))
  fit <- fitZig(dat,design)
  pvals <- fit$eb$p.value[,"conditionB"]
  padj <- p.adjust(pvals,method="BH")
  list(pvals=pvals, padj=padj)
}


runDESeq2Zero <- function(e){
      ## implement DESeq2 ##
    library(DESeq2) ; library(genefilter)
    condition = pData(e)$condition
  pickingSession = pData(e)$pickingSession
  dse <- DESeqDataSetFromMatrix(exprs(e), colData=DataFrame(pData(e)), design=~ condition + pickingSession)
    dse <- estimateSizeFactors(dse, type="poscounts")
    effLogLibSize <- log(colSums(counts(dse))*(1/sizeFactors(dse)))
    designZI=model.matrix(~effLogLibSize + pickingSession)
    zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~condition + pickingSession), colData=colData(dse), plot=FALSE, maxit=200, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos", designZI=designZI)
    dimnames(zeroWeights) = NULL
    assays(dse)[["weights"]] = zeroWeights
    dse <- estimateDispersions(dse)
    dse <- nbinomWaldTest(dse, betaPrior=TRUE, modelMatrixType="standard")
    #dse <- DESeq(dse, betaPrior=TRUE)
    res <- results(dse, name="condition_B_vs_A")
    baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
    pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
    padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
    list(pvals=pvalDesZero, padj=padjusted$padj)
}

