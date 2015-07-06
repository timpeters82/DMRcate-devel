fitChrom <- function(xx, lambda=1000, C=2, consec=FALSE, conseclambda=10) {
  message("Fitting ", unique(xx$CHR), "...", sep="")
  if (consec) {
    stopifnot(conseclambda >= 1)
    lambda <- conseclambda
    pos <- seq_along(xx$pos)
  } else {
    pos <- xx$pos
  }
  sigma <- lambda / C
  lag <- lambda
  beta <- xx$weights
  df <- 1
  X2 <- beta ^ 2
  pvals <- KernelTest(pos = pos, X2 = X2, lambda = sigma, df = df)
  res <- data.frame(pvals = pvals, ID = xx$ID)
  rownames(res) <- res$ID
  return(res)
}
