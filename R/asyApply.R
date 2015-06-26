#' apply a function to each of the assays in an SE, efficiently
#'
#' @param x     anything descended from a SummarizedExperiment
#' @param fn    the actual function to run on each assay matrix
#'
#' @return      a list of assays having been transformed
#' 
asyApply <- function(x, fn, ..., parallel=FALSE) { 
    if(parallel) {
        mclapply(assays(x, withDimnames=FALSE), fn, ...)
    } else { 
        lapply(assays(x, withDimnames=FALSE), fn, ...)
    }
}
