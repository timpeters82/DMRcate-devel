cox.annotate <- function(object, 
                         surv, 
                         annotation = c(array="IlluminaHumanMethylation450k",
                                        annotation="ilmn12.hg19"), ...) {

  library(simulatorZ)
  if (!is.matrix(object)) object <- getM(object)
  stopifnot(!all( c("time", "status") %in% colnames(surv)))
  adj.p <- p.adjust(rowCoxTests(object, surv, ...)$p.value, "fdr")
  nsig <- sum(adj.p < 0.05)

  if (nsig == 0) {
    message("Your contrast returned no individually significant probes. ",
            "Set pcutoff manually in dmrcate() to return DMRs, ",
            "but be warned there is an increased risk of Type I errors.")
  } else if (nsig > 0 & nsig <= 100) {
    message("Your contrast returned ", nsig, 
            " individually significant probes; a small but real effect. ",
            "Consider manually setting the value(s) of pcutoff to return ",
            "more DMRs, but be warned that doing this increases the risk ",
            "of Type I errors.")
  } else if (nsig > 100) {
    message("Your contrast returned ", nsig, 
            " individually significant probes. ",
            "We recommend the default setting of pcutoff in dmrcate().")
  }
 
  betafit <- rowCoxTests(ilogit2(object), surv, ...)
  object$indfdr <- p.adjust(betafit$p.value, method="fdr")
  object$betafc <- scale(betafit$coef) ## normalized Cox effect
  RSobject <- RatioSet(object, annotation = annotation)
  RSanno <- getAnnotation(RSobject)
  weights <- sqrt(qchisq(betafit.p, 1, lower=F)) ## t ~ sqrt(x2)
  annotated <- data.frame(ID = rownames(object), 
                          weights = weights, 
                          CHR = RSanno$chr, 
                          pos = RSanno$pos, 
                          gene = RSanno$UCSC_RefGene_Name, 
                          group = RSanno$UCSC_RefGene_Group, 
                          betafc = tt$betafc, 
                          indfdr = tt$adj.P.Val)
  
  annotated <- annotated[order(annotated$CHR, annotated$pos), ]
  class(annotated) <- "annot"
  return(annotated)
}
