rmSNPandCH <- function(object, dist=2, mafcut=0.05, and=TRUE, rmcrosshyb=TRUE, rmXY=FALSE)
{
  stopifnot(is.matrix(object))
  dist <- as.integer(dist)
  stopifnot(0 <= mafcut & mafcut <= 1)
  
  eh <- ExperimentHub()
  snpsall <- eh[["EH3130"]]
  
  dist0 <- as.integer(snpsall$distances)
  distrange <- range(dist0)
  stopifnot(dist >= min(distrange) && dist <= max(distrange))
  
  test0 <- (dist0 >= -1) & (dist0 <= dist)
  test1 <-  snpsall$mafs > mafcut
  test <- if (and) (test0 & test1) else (test0 | test1)
  
  badprobes <- snpsall$probe[test]
  if(rmcrosshyb){
    crosshyb <- eh[["EH3129"]]
    badprobes <- union(badprobes, as.character(crosshyb))
  }
  
  if(rmXY){
    XY.probes <- eh[['EH3131']]
    badprobes <- union(badprobes, XY.probes)
  }
  
  object[!(rownames(object) %in% badprobes),]
}

