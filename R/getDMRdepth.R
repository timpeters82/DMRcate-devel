#' @param output DMRcate output (of class "dmrcate.output", as it happens)
#' @param stepSize multiplier for scores. Typically if using 10**-y, leave at 1 
#'
getDMRdepth <- function(output, stepSize=1, ...) { 

  stopifnot(class(output) == "dmrcate.output")

  bySign <- extractRanges(output, bySign=TRUE)

  binDepths <- function(x) {
    bins <- disjoin(x)
    bins$score <- countOverlaps(bins, x)
    if(median(sign(score(x)) < 0)) bins$score <- -1 * bins$score
    return(bins)
  }

  ## merge hyper and hypo, should never overlap (!)
  data(seqinfo.hg19)
  depths <- GRangesList(lapply(bySign, binDepths))
  ## add seqinfo so it can be dumped to a bigWig file
  seqinfo(depths) <- seqinfo.hg19[seqlevels(depths)] 

  fixOverlaps <- function(overlaps, segs) { # {{{
    segs$score <- rep(0, length(segs))
    for (i in seq_along(segs)) {
      segs[i]$score <- sum(subsetByOverlaps(overlaps, segs[i])$score)
    }
    return(segs[score(segs) != 0])
  } # }}}

  ## patch up any overlapping DMRs 
  overlapping <- NULL
  if(length(depths) > 1) {
    overlapping <- GRangesList(hyper=subsetByOverlaps(depths$hyper,depths$hypo),
                               hypo=subsetByOverlaps(depths$hypo,depths$hyper))
  } ## only makes sense if we have both kinds!
  if (length(unlist(overlapping)) > 0) {
    message("Overlapping DMRs with opposite signs found, fixing...")
    excludeDisputed <- function(x, y) x[-queryHits(findOverlaps(x, y))]
    hyper <- excludeDisputed(depths$hyper, depths$hypo)
    hypo <- excludeDisputed(depths$hypo, depths$hyper)
    depths$hyper <- hyper
    depths$hypo <- hypo 
    overlaps <- sort(unlist(overlapping))
    overlaps <- split(overlaps, seqnames(overlaps))
    overlaps <- overlaps[lapply(overlaps, length) > 0]
    segs <- lapply(overlaps, disjoin)
    fixed <- unlist(GRangesList(mapply(fixOverlaps, overlaps, segs)))
    depths <- sort(c(unlist(depths), fixed))
  } else { 
    depths <- sort(unlist(depths))
  }
  depths$score <- depths$score * stepSize
  return(depths)

}
