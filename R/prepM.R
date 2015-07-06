#' this is the core of prepMvals but omits sanity checking; useful if collapsed
#'
#' @import impute 
#' 
#' @param M       M-value matrix
#' @param cutoff  extremes at which to truncate M-values 
#' 
#' @return        a matrix 
#' 
prepM <- function(M, cutoff=10) {
 
  if (is(M, "SummarizedExperiment")) {
    if (is(M, "RangedSummarizedExperiment")) M <- logit2(assays(M)$Beta)
    if (is(M, "GenomicRatioSet")) M <- getM(M) 
    else M <- -1 * assays(M)$exprs ## HELP data 
  }

  ## get rid of +/-Inf cells
  M[ which(M > cutoff) ] <- cutoff
  M[ which(M < -1 * cutoff) ] <- -1 * cutoff

  ## impute NAs via k-NN
  if (any(is.na(M))) {
    message("Imputing NAs...")
    M <- impute.knn(M)$data
  }
  return(M)

}
