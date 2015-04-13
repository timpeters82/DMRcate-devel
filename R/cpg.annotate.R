cpg.annotate <- function (object, annotation = c(array = "IlluminaHumanMethylation450k", 
                                                 annotation = "ilmn12.hg19"), 
                          analysis.type = c("differential", "variability"), design, contrasts = FALSE, 
                          cont.matrix = NULL, coef, ...) 
{
  stopifnot(is.matrix(object))
  analysis.type <- match.arg(analysis.type)
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
    tt <- topTable(fit, coef=coef, number = nrow(object))
    nsig <- sum(tt$adj.P.Val < 0.05)
    if(nsig==0){message("Your contrast returned no individually significant probes. Set pcutoff manually in dmrcate() to return DMRs, but be warned there is an increased risk of Type I errors.")}
    if(nsig >0 & nsig <= 100) {message(paste("Your contrast returned", nsig, "individually significant probes; a small but real effect. Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases the risk of Type I errors."))}
    if(nsig > 100){message(paste("Your contrast returned", nsig, "individually significant probes. We recommend the default setting of pcutoff in dmrcate()."))}
    betafit <- lmFit(ilogit2(object), design, ...)
    if (contrasts) {
      betafit <- contrasts.fit(betafit, cont.matrix)
    }
    betafit <- eBayes(betafit)
    betatt <- topTable(betafit, coef=coef, number = nrow(object))
    m <- match(rownames(tt), rownames(betatt))
    tt$betafc <- betatt$logFC[m]
    m <- match(rownames(tt), rownames(object))
    object <- object[m, ]
    RSobject <- RatioSet(object, annotation = annotation)
    RSanno <- getAnnotation(RSobject)
    weights <- tt$t
    annotated <- data.frame(ID = rownames(object), weights = weights, 
                            CHR = RSanno$chr, pos = RSanno$pos, gene = RSanno$UCSC_RefGene_Name, 
                            group = RSanno$UCSC_RefGene_Group, betafc = tt$betafc, indfdr = tt$adj.P.Val)
  }, variability = {
    RSobject <- RatioSet(object, annotation = annotation)
    RSanno <- getAnnotation(RSobject)
    wholevar <- var(object)
    weights <- apply(object, 1, var)
    weights <- weights/mean(weights)
    annotated <- data.frame(ID = rownames(object), weights = weights, CHR = RSanno$chr, pos = RSanno$pos, gene = RSanno$UCSC_RefGene_Name, 
                            group = RSanno$UCSC_RefGene_Group, betafc = rep(0, nrow(object)), indfdr = rep(0, nrow(object)))
    
  })
  annotated <- annotated[order(annotated$CHR, annotated$pos), 
                         ]
  class(annotated) <- "annot"
  return(annotated)
}