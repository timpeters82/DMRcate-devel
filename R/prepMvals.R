#' this is the recommended user-facing probe-name version
#' 
#' @import minfi
#' 
#' @param grset         A GenomicRatioSet or something very much like it 
#' @param cutoff        Maximum absolute M-value before truncation (10)
#' @param returnBetas   Return beta values instead of M-values (FALSE) 
#' 
#' @return              A matrix of tidied and truncated M- or beta-values
#'
prepMvals <- function(grset, cutoff=10, returnBetas=FALSE) { 

  CpGs <- grep("^cg", rownames(grset))
  if (is(grset, "RangedSummarizedExperiment")) {
    xx <- prepM(rmSNPandCH(logit2(assays(grset)$Beta)[CpGs, ]))
  } else { 
    xx <- prepM(rmSNPandCH(getM(grset)[CpGs, ]))
  }
  if (returnBetas == TRUE) return(ilogit2(xx))
  else return(xx)
        
}
