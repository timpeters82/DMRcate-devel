dmrcate <- function(object, lambda=1000, C=2, 
                    p.adjust.method=c("limma","qvalue","fdr","BH"), 
                    pcutoff=0.05, consec=FALSE, conseclambda=10, 
                    betacutoff=NULL, parallel=FALSE) {

  ## {{{ ensure things are amenable to demarcation! 
  stopifnot(lambda >= 1)
  stopifnot(C >= 0.2)
  stopifnot(is(object, "annot") | ## for debugging, to avoid refitting:
            (is(object, "data.frame") && "indfdr" %in% names(object)))
  stopifnot(all(pcutoff >= 0) && all(pcutoff <= 1))
  if(consec & is.null(conseclambda)) {
    stop("Consecutive probe bandwidth must be specified")
  }
  ## }}}

  ## "limma" makes more sense as a p.adjust method than a cutoff
  p.adjust.method <- match.arg(p.adjust.method)

  ## support multiple p-value cutoffs to allow regional backoff...
  if (length(pcutoff) > 1) names(pcutoff) <- as.character(pcutoff)

  ## reshape the annotation results into a data.frame with appropriate columns
  object <- with(object, data.frame(ID, weights, CHR, pos, gene, group, 
                                    betafc, indfdr, stringsAsFactors=FALSE))
 
  ## parallel or serial by chromosome 
  message("Fitting by chromosome...")
  object$raw <- fitByChrom(object, parallel=parallel, 
                           lambda=lambda, consec=consec, 
                           conseclambda=conseclambda)

  ## support multiple p-value cutoffs
  if (p.adjust.method == "limma") {
    ## use limma method, but on q-values to avoid penalty
    object$qval <- qvalue(object$raw)$qvalues 
    if (sum(object$indfdr < max(pcutoff)) == 0) {
      stop("No significant regions")
    } else { 
      pcutoff <- unlist(sapply(c(pcutoff), 
                               function(p) {
                                 pp <- sort(object$qval)[sum(object$indfdr < p)]
                                 if(length(pp) == 0) pp <- 0
                                 return(pp)
                               }))
    }
    ## in case we need it later, record the fdr, too
    object$fdr <- p.adjust(object$raw, method="fdr")
  } else if (p.adjust.method == "qvalue") {
    object$fdr <- qvalue(object$raw)$qvalues 
  } else {
    object$fdr <- p.adjust(object$raw, method=p.adjust.method)
  } 

  ## cut on max(pcutoff) if multi-pcutoff
  if (all(object$fdr >= max(pcutoff))) {
    txt <- paste("No signficant regions found. Try increasing the value of",
                 "'pcutoff' in 'dmrcate' and/or 'cpg.annotate'.")
    stop(paste(strwrap(txt, exdent=2), collapse="\n"))
  } else {
    message("Identifying significant probes...")
    sigprobes <- object[object$fdr <= max(pcutoff), , drop = FALSE]
  }

  message("Demarcating regions...")
  results <- dmrcateByChr(sigprobes, pcutoff, lambda, parallel)
  if (!is.null(betacutoff)) {
    results <- results[which(abs(results$maxbetafc) >= betacutoff), ]
    if (nrow(results) == 0) {
      stop("No regions with at least ", betacutoff, " maximum differences!")
    }
    if (length(pcutoff) > 1) {
      ## fixup after limma fdr
      for (p in names(pcutoff)) {
        results$pcutoff[which(results$pcutoff == pcutoff[p])] <- as.numeric(p)
      }
    }
  }
  message("Done!")

  output <- list(input=object, results=results, cutoff=pcutoff, 
                 p.adjust.method=p.adjust.method, lambda=lambda, C=C)
  class(output) <- "dmrcate.output"
  return(output)

}
