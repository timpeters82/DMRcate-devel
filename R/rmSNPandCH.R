rmSNPandCH <- function(object, dist=2, mafcut=0.05, and=TRUE, rmcrosshyb=TRUE, rmXY=FALSE)
{
  stopifnot(is.matrix(object))
  dist <- as.integer(dist)
  stopifnot(0 <= mafcut & mafcut <= 1)
  
  env <- new.env(parent=emptyenv())
  data(dmrcatedata, envir=env)
  len0 <- vapply(env$snpsall$Distance, length, integer(1),
                 USE.NAMES=FALSE)
  len1 <- vapply(env$snpsall$MinorAlleleFrequency, length,
                 integer(1), USE.NAMES=FALSE)
  keep <- len0 == len1
  len0 <- len0[keep]
  env$snpsall <- env$snpsall[keep,]
  
  dist0 <- as.integer(unlist(env$snpsall$Distance, use.names=FALSE))
  distrange <- range(dist0)
  stopifnot(dist >= min(distrange) && dist <= max(distrange))
  
  test0 <- (dist0 >= -1) & (dist0 <= dist)
  test1 <-
    unlist(env$snpsall$MinorAlleleFrequency, use.names=FALSE) > mafcut
  test <- if (and) (test0 & test1) else (test0 | test1)
  ## 'any' by group
  ntrue <- cumsum(test)[cumsum(len0)]
  badidxs <- ntrue - c(0, head(ntrue, -1)) != 0
  
  badprobes <- rownames(env$snpsall)[badidxs]
  if(rmcrosshyb){
    badprobes <- union(badprobes, as.character(env$crosshyb))
  }
  
  if(rmXY){
    badprobes <- union(badprobes, env$XY.probes)
  }
  
  object[!(rownames(object) %in% badprobes),]
}

