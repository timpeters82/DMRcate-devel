fitByChrom <- function(x, parallel=FALSE, lambda=1000, C=2,  
                       consec=FALSE, conseclambda=10) {

  stopifnot(all(c("weights","CHR","ID") %in% colnames(x)))
  x$ID <- as.character(x$ID)
  if(parallel) {
    require("parallel") ## FIXME: finesse this
    pvals <- do.call(rbind, 
                     mclapply(split(x, x$CHR), 
                              fitChrom, 
                              lambda=lambda, 
                              consec=consec, 
                              conseclambda=conseclambda))
  } else { 
    pvals <- do.call(rbind, 
                     lapply(split(x, x$CHR), 
                            fitChrom, 
                            lambda=lambda, 
                            consec=consec, 
                            conseclambda=conseclambda))
  }

  ## ugly hack for combining pvals
  pvals$ID <- as.character(pvals$ID)
  rownames(pvals) <- pvals$ID
  if(!identical(pvals$ID, x$ID)) {
    pvals <- pvals[x$ID, ] 
  }
  return(pvals$pvals)

}
