#' collapse variable, and perhaps overlapping, rows of a matrix-like structure 
#'
#' @param x     anything descended from a SummarizedExperiment
#' @param z     the listed matrix indices within DMRs (NOT THE ACTUAL DMRs!)
#' @param f     the name of the function to be used to summarize x over z 
#'
#' @return      a list of assays, having been transformed
#' 
collapseMatrix <- function(x, z, fn, impute=TRUE) {
  fn <- selectMethod(fn, class(x))
  if (any(is.na(x)) && impute == TRUE) x <- impute.knn(x)$data
  do.call(rbind, lapply(z, function(w) fn(x[w, , drop=F])))
}
