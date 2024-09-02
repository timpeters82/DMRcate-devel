sequencing.annotate <- function(obj, methdesign, all.cov=FALSE, contrasts = FALSE, 
                                cont.matrix = NULL, fdr = 0.05, coef, ...){
  
  if(is(obj, "data.frame")){
    if (!(all(c("stat", "chr", "pos", "diff", "fdr") %in% colnames(obj)) | all(c("stat", "chr", "pos", "fdrs") %in% colnames(obj)))) {
      stop("Error: object does not contain all required columns, was it created by DSS::DMLtest() or DSS::DMLtest.multiFactor?")
    }
    if(!"diff" %in% colnames(obj)){
      obj$diff <- 0
    }
    obj$pval <- 2*pnorm(-abs(obj$stat))
    obj$fdr <- p.adjust(obj$pval, method="BH")
    nsig <- sum(obj$fdr < fdr)
    annotated <- GRanges(as.character(obj$chr), IRanges(obj$pos, obj$pos), stat = obj$stat,
                         rawpval = obj$pval, diff = obj$diff, ind.fdr = obj$fdr, is.sig = obj$fdr < fdr)
    names(annotated) <- rownames(obj)
    annotated <- sort(annotated)
    
  } else if(is(obj, "BSseq")){
    phen <- pData(obj)
    obj <- BSseq(M=as.matrix(getCoverage(obj, type = "M")), 
                 Cov=as.matrix(getCoverage(obj, type = "Cov")), 
                 pos=start(obj), chr=seqnames(obj))
    colData(obj) <- phen
    if(any(width(obj) > 1)){
      stop("Error: all ranges in the BSseq object must be single nucleotides with width 1.")
    }
    if(is.null(rownames(colData(obj)))){
      stop("Error: BSseq object must be annotated with colData with sample IDs as rownames of the data.frame.")
    }
    if(all.cov){
      message("Filtering out all CpGs where at least one sample has zero coverage...")
      obj <-obj[apply(getCoverage(obj, type = "Cov"), 1, function (x) all(x > 0)),]
    } else {
      message("Filtering out CpGs where no samples have coverage...")
      obj <- obj[rowSums(getCoverage(obj, type = "Cov")) > 0,]
    }
    message("Processing BSseq object...")
    obj <- orderBSseq(obj)
    meth <- getCoverage(obj, type="M")
    unmeth <- getCoverage(obj, type="Cov") - meth
    countmatrix <- eval(parse(text=paste0("cbind(", paste(sapply(1:ncol(meth),
                                                                 function (x) gsub("idx", x, "meth[,idx], unmeth[,idx]")), collapse=', '), ")")))
    
    colnames(countmatrix) <- paste(rep(rownames(colData(obj)), each=2), 
                                   rep(c("C", "T"), times=ncol(obj)), sep=".")
    message("Transforming counts...")
    stopifnot(is.matrix(methdesign))
    ym <- voom(countmatrix, methdesign, lib.size = rep(colSums(meth+unmeth), each=2))$E
    if (contrasts & is.null(cont.matrix)) {
      stop("Error: a contrast matrix must be specified if contrasts = TRUE")
    } 
    message("Fitting model...")
    fit <- lmFit(ym, methdesign, ...)
    if (contrasts) {
      stopifnot(coef %in% colnames(cont.matrix))
      fit <- contrasts.fit(fit, cont.matrix)
    }
    fit <- eBayes(fit)
    tt <- topTable(fit, coef = coef, number = nrow(ym), sort.by = "none")
    nsig <- sum(tt$adj.P.Val < fdr)
    annotated <- GRanges(as.character(seqnames(obj)), IRanges(start(obj), start(obj)), stat = tt$t,
                         rawpval = tt$P.Value, diff = tt$logFC, ind.fdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < fdr)
    names(annotated) <- paste(seqnames(obj), start(obj), sep=":")
    annotated <- sort(annotated)
    
  } else {
    stop("Error: obj must be a data.frame or BSseq object")
  }
  if (nsig == 0) {
    message("Your contrast returned no individually significant CpGs. Consider increasing the 'fdr' parameter using changeFDR(), but be warned there is an increased risk of Type I errors.")
  }
  if (nsig > 0 & nsig <= 100) {
    message(paste("Your contrast returned", nsig, 
                  "individually significant CpGs; a small but real effect. Consider increasing the 'fdr' parameter using changeFDR(), but be warned there is an increased risk of Type I errors."))
  }
  if (nsig > 100 & nsig <= 100000) {
    message(paste("Your contrast returned", nsig, 
                  "individually significant CpGs. We recommend the default setting of pcutoff in dmrcate()."))
  }
  if (nsig > 100000) {
    message(paste("Your contrast returned", nsig, 
                  "individually significant CpGs; this is plenty. Consider decreasing the 'fdr' parameter using changeFDR(), for more precise DMR definition."))
  }
  return(new("CpGannotated", ranges=annotated, betas=matrix("Betas can be obtained with bsseq::getMeth(obj).", 1, 1)))
}