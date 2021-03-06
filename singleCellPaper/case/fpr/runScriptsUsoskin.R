

runEdgeR <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- edgeR::calcNormFactors(dgel)
  dgel=estimateDisp(dgel,design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  #predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  #predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj)
}

runEdgeRRobust <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- edgeR::calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj)
}


runVoom <- function(e) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- edgeR::calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value 
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=tt$logFC)
}

runDESeq2 <- function(e, retDDS=FALSE) {
    library(DESeq2)
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,betaPrior=TRUE,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

runDESeq2_poscounts <- function(e, retDDS=FALSE) {
    library(DESeq2)
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- estimateSizeFactors(dds,type="poscounts")
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds,betaPrior=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

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

zeroWeightsLibSizeDispFast <- function(counts, design, colData=NULL, initialWeightAt0=TRUE, maxit=100, plot=FALSE, plotW=FALSE, designZI=NULL, llTol=1e-4, normalization="TMM"){
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
	design <- model.matrix(~ pData(e)$condition)
	library(edgeR) ; library(genefilter)
	d <- DGEList(exprs(e))
	d <- edgeR::calcNormFactors(d)
	#not adding a design matrix models the zeroes with the library size automatically
  effLogLibSize = log(d$samples$lib.size*d$samples$norm.factors)
  pickingSession = pData(e)[,"Picking sessions"]
  designZI = model.matrix(~effLogLibSize + pickingSession)
	zeroWeights = zeroWeightsLibSizeDispFast(d, design, plot=FALSE, maxit=200, initialWeightAt0=TRUE, plotW=FALSE, designZI=designZI)
	d$weights = zeroWeights
	d=estimateDispWeighted(d,design,weights=zeroWeights, grid.range=c(-15,15))
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
    ediff <- scde.expression.difference(o.ifm, exprs(e), o.prior, groups  =  pData(e)$condition, n.randomizations  =  100, n.cores  =  1, verbose  =  0)
    pvals=(1-pnorm(abs(ediff$Z)))*2
    padj=p.adjust(pvals,method="BH")
    list(Z=ediff$Z,pvals=pvals,padj=padj)
}


runMAST <- function(e){
    require(MAST)
    counts=exprs(e)
    tpm <- counts*1e6/colSums(counts)
    sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(group=pData(e)$condition))
    ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
    CD <- cData(sca)
    CD$ngeneson <- ngeneson
    CD$cngeneson <- CD$ngeneson-mean(ngeneson)
    cData(sca) <- CD
     ## differential expression 
    fit <- zlm.SingleCellAssay(~cngeneson+group,sca=sca,method="bayesglm",ebayes=TRUE)
    L=matrix(0,nrow=ncol(coef(fit,"D")))
    rownames(L)=colnames(coef(fit,"D"))
    L["groupB",]=1
    lrFit <- lrTest(fit, hypothesis=L)
    pval=lrFit[,'hurdle','Pr(>Chisq)']
    padj=p.adjust(pval,method="BH")
    list(pvals=pval,padj=padj)
}

runMAST_count <- function(e){
    require(MAST)
    counts=exprs(e)
    tpm <- counts*1e6/colSums(counts)
    sca <- FromMatrix('SingleCellAssay', t(tpm), cData=data.frame(group=pData(e)$condition))
    ngeneson <- apply(exprs(sca),1,function(x)mean(x>0))
    CD <- cData(sca)
    CD$ngeneson <- ngeneson
    CD$cngeneson <- CD$ngeneson-mean(ngeneson)
    cData(sca) <- CD
     ## differential expression 
    fit <- zlm.SingleCellAssay(~cngeneson+group,sca=sca,method="bayesglm",ebayes=TRUE)
    L=matrix(0,nrow=ncol(coef(fit,"D")))
    rownames(L)=colnames(coef(fit,"D"))
    L["groupB",]=1
    lrFit <- lrTest(fit, hypothesis=L)
    pval=lrFit[,'cont','Pr(>Chisq)']
    padj=p.adjust(pval,method="BH")
    list(pvals=pval,padj=padj)
}


runLimmaHurdle <- function(e){
    ## limma voom pipeline ##
	library(limma)
	counts=exprs(e)
	group=pData(e)$condition
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
	list(pvals = pval, padj = padj)
}



runEdgeRHurdle <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel$weights <- 1-(exprs(e)==0)
  dgel <- edgeR::calcNormFactors(dgel)
  dgel=estimateDisp(dgel,design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  pvals <- edger.lrt$table$PValue
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj)
}

runMetagenomeSeq <- function(e){
  require(metagenomeSeq)
  design <- model.matrix(~pData(e)$condition)
  pheno <- AnnotatedDataFrame(data.frame(group=pData(e)$condition))
  rownames(pheno) <- colnames(exprs(e))
  p <- cumNormStatFast(exprs(e))
  dat <- newMRexperiment(counts=exprs(e), phenoData=pheno, featureData = NULL, libSize = colSums(exprs(e)), normFactors = metagenomeSeq::calcNormFactors(exprs(e), p=p))
  fit <- fitZig(dat,design)
  pvals <- fit$eb$p.value[,"pData(e)$conditionB"]
  padj <- p.adjust(pvals,method="BH")
  list(pvals=pvals, padj=padj)
}


runDESeq2Zero <- function(e){
      ## implement DESeq2 ##
    library(DESeq2) ; library(genefilter)
    condition=pData(e)$condition
    colData=DataFrame(pData(e))
    dse <- DESeqDataSetFromMatrix(exprs(e), colData, ~ condition)
    dse <- estimateSizeFactors(dse, type="poscounts")
    effLogLibSize <- log(colSums(counts(dse))*(1/sizeFactors(dse)))
    pickingSession = pData(e)[,"Picking sessions"]
    designZI=model.matrix(~effLogLibSize + pickingSession)
    zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~condition), colData=colData, plot=FALSE, maxit=200, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos", designZI=designZI)
    dimnames(zeroWeights) = NULL
    assays(dse)[["weights"]] = zeroWeights
    dse <- estimateDispersions(dse)
    dse <- nbinomWaldTest(dse, betaPrior=TRUE)
    #dse <- DESeq(dse, betaPrior=TRUE)
    res <- results(dse)
    baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
    pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
    padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
    list(pvals=pvalDesZero, padj=padjusted$padj)
}

# try betaPrior=TRUE and betaPrior=FALSE
# for both poscounts in dse and zeroWEights
# poscounts in zeroWEeigts but tmm in dse
# tmm in zeroweights and dse



 # runDESeq2ZeroTest <- function(e){
 #       ## implement DESeq2 ##
 #     library(DESeq2) ; library(genefilter)
 #     condition=pData(e)$condition
 #     colData=DataFrame(pData(e))
 #         pickingSession = pData(e)[,"Picking sessions"]
 #     dse <- DESeqDataSetFromMatrix(exprs(e), colData, ~ condition)
 #     #dse <- estimateSizeFactors(dse, type="poscounts")
 #     #dse = estimateSizeFactors(dse)
 #      sizeFactors(dse) = rep(1,120)
 #      #d=DGEList(exprs(e))
 #      #d=edgeR::calcNormFactors(d)
 #     #sizeFactors(dse)=1/d$samples$norm.factors
 #     effLogLibSize <- log(colSums(counts(dse))*(1/sizeFactors(dse)))
 #     designZI=model.matrix(~effLogLibSize + pickingSession)
 #     zeroWeights = zeroWeightsLibSizeDispFast(counts(dse), design=model.matrix(~condition), colData=colData, plot=FALSE, maxit=200, initialWeightAt0=TRUE, plotW=FALSE, normalization="DESeq2_pos", designZI=designZI)
 #     dimnames(zeroWeights) = NULL
 #     assays(dse)[["weights"]] = zeroWeights
 #     dse <- estimateDispersions(dse)
 #     dse <- nbinomWaldTest(dse, betaPrior=FALSE)
 #     #dse <- DESeq(dse, betaPrior=TRUE)
 #     res <- results(dse, cooksCutoff=Inf)
 #     baseMean=unname(rowMeans(sweep(counts(dse),2,1/sizeFactors(dse),FUN="*")))
 #     pvalDesZero = 2*(1-pt(abs(res$stat),df=rowSums(zeroWeights)-2))
 #     padjusted = pvalueAdjustment_kvdb(pValue=pvalDesZero, filter=baseMean, alpha=0.05)
 #     list(pvals=pvalDesZero, padj=padjusted$padj)
 # }
