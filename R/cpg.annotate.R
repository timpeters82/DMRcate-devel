cpg.annotate <- function (datatype = c("array", "sequencing"), object, what = c("Beta", 
                                                                "M"), arraytype = c("EPIC", "450K"), analysis.type = c("differential", 
                                                                                                                       "variability", "ANOVA", "diffVar"), design, contrasts = FALSE, 
          cont.matrix = NULL, fdr = 0.05, coef, 
          varFitcoef=NULL, topVarcoef=NULL, ...) 
{
  analysis.type <- match.arg(analysis.type)
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  if (datatype == "array") {
    stopifnot(class(object)[1] %in% c("matrix", "GenomicRatioSet"))
    if (is(object, "matrix")) {
      if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               what = what)
      }
      if (arraytype == "EPIC") {
        grset <- makeGenomicRatioSetFromMatrix(mat = object, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               what = what)
      }
    }
    else {
      grset <- object
    }
    object <- getM(grset)
    switch(analysis.type, differential = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fit <- lmFit(object, design, ...)
      if (contrasts) {
        stopifnot(coef %in% colnames(cont.matrix))
        fit <- contrasts.fit(fit, cont.matrix)
      }
      fit <- eBayes(fit)
      tt <- topTable(fit, coef = coef, number = nrow(object))
      nsig <- sum(tt$adj.P.Val < fdr)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      betafit <- lmFit(ilogit2(object), design, ...)
      if (contrasts) {
        betafit <- contrasts.fit(betafit, cont.matrix)
      }
      betafit <- eBayes(betafit)
      betatt <- topTable(betafit, coef = coef, number = nrow(object))
      m <- match(rownames(tt), rownames(betatt))
      tt$diff <- betatt$logFC[m]
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      anno <- getAnnotation(grset)
      stat <- tt$t
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, anno$pos), stat = stat,
                           diff = tt$diff, ind.fdr = tt$adj.P.Val, is.sig = tt$adj.P.Val < fdr)
      names(annotated) <- rownames(tt)
    }, variability = {
      RSanno <- getAnnotation(grset)
      wholevar <- var(object)
      weights <- apply(object, 1, var)
      weights <- weights/mean(weights)
      annotated <- GRanges(as.character(RSanno$chr), IRanges(RSanno$pos, RSanno$pos), stat = weights,
                           diff = rep(0, nrow(object)), ind.fdr = rep(0, nrow(object)), is.sig = weights > quantile(weights, 0.95))
      names(annotated) <- rownames(object)
    }, ANOVA = {
      message("You are annotating in ANOVA mode: consider making the value of fdr quite small, e.g. 0.001")
      stopifnot(is.matrix(design))
      fit <- lmFit(object, design, ...)
      fit <- eBayes(fit)
      sqrtFs <- sqrt(fit$F)
      sqrtfdrs <- p.adjust(fit$F.p.value, method = "BH")
      nsig <- sum(sqrtfdrs < fdr)
      if (nsig == 0) {
        message("Your design returned no individually significant probes for ANOVA. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your design returned", nsig, 
                      "individually significant probes for ANOVA. We recommend the default setting of pcutoff in dmrcate(). Large numbers (e.g. > 100000) may warrant a smaller value of the argument passed to fdr"))
      }
      anno <- getAnnotation(grset)
      stat <- sqrtFs
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, anno$pos), stat = stat,
                           diff = 0, ind.fdr = sqrtfdrs, is.sig = sqrtfdrs < fdr)
      names(annotated) <-  rownames(object)
    }, diffVar = {
      stopifnot(is.matrix(design))
      if (!contrasts) {
        stopifnot(colnames(design)[1] == "(Intercept)")
      } else {
        stopifnot(!is.null(cont.matrix))
      }
      fitvar <- varFit(object, design = design, coef=varFitcoef)
      if (contrasts) {
        fitvar <- contrasts.varFit(fitvar, cont.matrix)
      }
      tt <- topVar(fitvar, coef = topVarcoef, number = nrow(object))
      nsig <- sum(tt$Adj.P.Value < fdr)
      if (nsig == 0) {
        message("Your contrast returned no individually significant probes. Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() to return DVMRs, but be warned there is an increased risk of Type I errors.")
      }
      if (nsig > 0 & nsig <= 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DVMRs, but be warned that doing this increases the risk of Type I errors."))
      }
      if (nsig > 100) {
        message(paste("Your contrast returned", nsig, 
                      "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))
      }
      m <- match(rownames(object), rownames(tt))
      tt <- tt[m, ]
      anno <- getAnnotation(grset)
      stat <- tt$t
      annotated <- GRanges(as.character(anno$chr), IRanges(anno$pos, anno$pos), stat = stat,
                           diff = 0, ind.fdr = tt$Adj.P.Value, is.sig = tt$Adj.P.Value < fdr)
      names(annotated) <-  rownames(tt)  
    })
    annotated <- sort(annotated)
    return(new("CpGannotated", ranges=annotated))
    
  }
  if (datatype == "sequencing") {
    stop("Sequencing mode is deprecated for cpg.annotate(). Please use sequencing.annotate().")
  }
  else {
    message("Error: datatype must be one of 'array' or 'sequencing'")
  }
}
