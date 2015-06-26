#' omnibus DMR finder starting from a m-values and a design/Surv matrix
#' 
#' @param x           matrix of M-values or SummarizedExperiment-derived object
#' @param design      design matrix (for limma) or Surv object/matrix (for Cox)
#' @param contrasts   apply contrasts to e.g. paired samples? (only for limma)
#' @param cont.matrix if contrasts=T, supply a contrast matrix here
#' @param coef        for which column of the design matrix shall we get DMRs?
#' @param pcutoff     p-value cutoff; if not specified, step from 10**-1:10**-8
#' @param betacutoff  DMRs must have at least this great a maxbetaFC difference
#' @param p.adjust.method   how to control the FDR, default is limma-style
#' @param what        (not used yet) whether to use limma or Cox PH to tag DMRs
#'
getDMRs <- function(x, design=NULL, contrasts=F, cont.matrix=NULL, 
                    coef=2, pcutoff=NULL, betacutoff=0.1, 
                    p.adjust.method="limma", what=c("limma", "cox"), ...) {

  ## support for Cox PH is coming:
  what <- tolower(match.arg(what))

  if (what == "cox") 
    stop("Cox regression DMRs are not (yet) supported.")
  if (what == "cox" && p.adjust.method == "limma")
    stop("Cannot use limma-style p-value correction for Cox regression DRMs.")

  if (is.null(pcutoff)) pcutoff <- 10 ** (-1 * seq(1, 8)) ## nice for bigWigs

  if (is(x, "SummarizedExperiment") || is(x, "RangedSummarizedExperiment")) {
    x <- prepMvals(x) ## would rather internalize its rowRanges within dmrcate
  }
  if (is.null(design)) {
    if ((is(x,"SummarizedExperiment") || is(x,"RangedSummarizedExperiment")) &&
        ("design" %in% names(exptData(x)) || "design" %in% names(metadata(x)))){
      if (is(x, "SummarizedExperiment")) design <- exptData(x)$design
      if (is(x, "RangedSummarizedExperiment")) design <- metadata(x)$design
    } else {
      stop("You need a design matrix (perhaps metadata(x)$design) to call DMRs")
    }
  }
  message("Annotating individual CpGs...")
  DMRannot <- switch(what, 
                     limma=cpg.annotate(x, design=design, coef=coef),
                     cox=cox.annotate(x, Surv(x$OS, x$OSevent)))

  message("Demarcating significant regions...")
  dmrcate(DMRannot, pcutoff=pcutoff, betacutoff=betacutoff, 
          p.adjust.method=p.adjust.method, ...)
}

getVMRs <- function(x, pcutoff=0.1, ...) 
{
  if (is(x, "SummarizedExperiment") || is(x, "RangedSummarizedExperiment")) {
    x <- prepMvals(x)
  }
  VMRannot <- cpg.annotate(x, analysis.type="variability", pcutoff=pcutoff)
  dmrcate(VMRannot, ...)
}

getDMRsAndVMRs <- function(x, design=NULL, contrasts=F, 
                           cont.matrix=NULL, coef=2, ...){
  mvals.x <- prepMvals(x)
  list(DMRs=getDMRs(mvals.x, design, coef, ...), VMRs=getVMRs(mvals.x, ...))
}
