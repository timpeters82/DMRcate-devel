getSharedSegments <- function(..., minDepth=NULL) {
  if (is.null(minDepth)) minDepth <- length(...)
  if (length(...) < 2) stop('Only one set of ranges provided!')
  grl <- GRangesList(as(..., 'list'))
  merged <- unlist(grl)
  disjoint <- disjoin(merged)
  disjoint$depth <- countOverlaps(disjoint, merged)
  return(reduce(disjoint[ disjoint$depth >= minDepth]))
}

