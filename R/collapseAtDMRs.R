#' Collapse a SummarizedExperiment-like object over a GRanges of DMRs, sensibly.
#' Handy for things like model testing, DMR plotting, GenometriCorr, lars, etc.
#' 
#' FIXME: let the user provide their own summarizer function
#'
#' @import minfi  
#' @import impute
#' @import GenomicRanges 
#' 
#' @param x       anything descended from a SummarizedExperiment
#' @param y       differentially whatevered regions (DMRs) as a GRanges
#' @param how     how to summarize the assay values (default: median)
#' @param impute  should NA entries in input matrices be imputed? (yes)
#'
#' @return        an object with same colData but with new rowRanges & assays
#' 
collapseAtDMRs <- function(x, y, 
                           how=c("median", "mean", "sum", "max", "min"),
                           impute=TRUE) {

    ## check arguments
    stopifnot(is(x, "SummarizedExperiment") || 
              is(x, "RangedSummarizedExperiment"))
    stopifnot(class(x) %in% c("SummarizedExperiment","GenomicRatioSet"))
    if (class(y) == "dmrcate.output") y <- extractRanges(y)
    stopifnot(is(y, "GenomicRanges")) ## must be a GRanges
    if (is.null(names(y))) {
      if ("name" %in% names(mcols(y))) {
        names(y) <- make.unique(y$name)
      } else { 
        names(y) <- paste0("region", seq_len(length(y)))
      }
    }
    if (is(x, "GenomicRatioSet")) x <- x[ grep("^cg", rownames(x)), ] 

    ## find and index the runs
    yy <- subsetByOverlaps(y, x)
    xx <- subsetByOverlaps(x, yy)
    hitz <- findOverlaps(xx, yy)
    byDMR <- as(lapply(split(as(hitz, "data.frame"), subjectHits(hitz)),
                       function(z) z[, "queryHits"]), "List")
    names(byDMR) <- names(yy)

    ## obtain the summarizer
    how <- match.arg(how)
    fnBy <- switch(how,
                   "median"="colMedians",
                   "mean"="colMeans",
                   "sum"="colSums",
                   "max"="colMaxs",
                   "min"="colMins")

    ## loop through the assays as found in the SE 
    summarizeAsy <- function(xxy) collapseMatrix(xxy, byDMR, fnBy, impute=TRUE)
    summarizedAssays <- asyApply(xx, summarizeAsy)

    ## the following is a horrible mess 
    if (class(x) == "GenomicRatioSet") {
      
      ## shenanigans required to reassemble a GenomicRatioSet after collapsing
      if (!"Beta" %in% names(summarizedAssays)) summarizedAssays$Beta <- NULL
      if (!"M" %in% names(summarizedAssays)) summarizedAssays$M <- NULL

      ## now reassemble the GRSet
      res <- GenomicRatioSet(gr=yy,
                             Beta=summarizedAssays$Beta,
                             M=summarizedAssays$M, 
                             pData=colData(x),
                             annotation=annotation(x),
                             preprocessMethod=preprocessMethod(x))
    } else if (class(x) == "RangedSummarizedExperiment") {
      res <- SummarizedExperiment(assays=summarizedAssays,
                                  metadata=metadata(x),
                                  colData=colData(x), 
                                  rowRanges=yy)
    } else if (class(x) == "SummarizedExperiment") {
      res <- SummarizedExperiment(assays=summarizedAssays,
                                  exptData=exptData(x),
                                  colData=colData(x), 
                                  rowData=yy)
    } else {
      stop("Don't know how to summarize a ", class(x))
    }
    rownames(res) <- paste0(seqnames(yy), ":", start(yy), "-", end(yy))
    return(res)

}
