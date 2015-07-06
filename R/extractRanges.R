#' extract GRanges corresponding to DMRs from the results of a dmrcate run
#' 
#' @param   dmrcated  an object of class dmrcate.output
#' @param   bySign    boolean, shall the DMRs be extracted by sign?  (FALSE)
#' @param   atCutoff  a p-value cutoff desired for extraction, or NULL (NULL)
#' 
extractRanges <- function(dmrcated, bySign=FALSE, atCutoff=NULL, ...) {

  if (is(dmrcated, "GRanges") || is(dmrcated, "GRangesList")) {
    return(dmrcated)
  } else if (!is(dmrcated, "dmrcate.output")) {
    stop("Argument is not of class dmrcate.output; cannot process.")
  } else {
    dmrc <- DataFrame(cbind(extractCoords(dmrcated$results$hg19coord), 
                            dmrcated$results))
    dmrc$score <- dmrcated$results$maxbetafc 
    if (!is.null(atCutoff)) dmrc <- dmrc[dmrc$pcutoff == atCutoff, ]
    gr <- makeGRangesFromDataFrame(dmrc, keep.extra.columns=TRUE)
    genome(gr) <- "hg19" ## nothing else is supported, yet :-/
    if (bySign == FALSE) {
      # GRanges
      return(gr)
    } else if (bySign == TRUE) {
      # GRangesList
      grl <- split(gr, sign(score(gr)))
      names(grl) <- sub("1", "hyper", sub("-1", "hypo", names(grl)))
      return(grl)
    }
  }

}
