fitParallel <- function(chr, object, consec, conseclambda, lambda, C){
  message(paste("Fitting ", chr, "...", sep = ""))
  chrIndex <- object$CHR %in% chr
  chromosome <- object[chrIndex,]
  if (consec) 
  {
    stopifnot(conseclambda >= 1)
    lambda <- conseclambda
    pos <- seq_along(chromosome$pos)
  }
  else 
  {
    pos <- chromosome$pos
  }
  sigma = lambda/C
  lag = lambda
  beta <- chromosome$weights
  df <- 1
  X2 <- beta^2
  pvalue <- KernelTest(pos = pos, X2 = X2, lambda = sigma, df = df)
  chromosome$raw <- pvalue
  chromosome
}