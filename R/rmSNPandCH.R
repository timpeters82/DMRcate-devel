rmSNPandCH <- function(object, dist=2, mafcut=0.05, and=TRUE, rmcrosshyb=TRUE, rmXY=FALSE)
{
  stopifnot(is.matrix(object))
  dist <- as.integer(dist)
  stopifnot(0 <= mafcut & mafcut <= 1)
  
  env <- new.env(parent=emptyenv())
  data(dmrcatedata, envir=env)
  
  dist0 <- as.integer(env$snpsall$distances)
  distrange <- range(dist0)
  stopifnot(dist >= min(distrange) && dist <= max(distrange))
  
  test0 <- (dist0 >= -1) & (dist0 <= dist)
  test1 <-  env$snpsall$mafs > mafcut
  test <- if (and) (test0 & test1) else (test0 | test1)
  
  badprobes <- snpsall$probe[test]
  if(rmcrosshyb){
    badprobes <- union(badprobes, as.character(env$crosshyb))
  }
  
  if(rmXY){
    badprobes <- union(badprobes, env$XY.probes)
  }
  
  object[!(rownames(object) %in% badprobes),]
}

