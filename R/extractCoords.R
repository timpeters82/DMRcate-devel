extractCoords <-
  function(coords.M)
  {
    M <- length(coords.M)
    A <- 3 # Should be three components
    split <- "[:-]"
    l <- strsplit(coords.M, split = split)
    stopifnot(all(sapply(l, length) == A))
    X.AM <- matrix(unlist(l), A, M)
    df.MA <- data.frame(t(X.AM), stringsAsFactors = FALSE)
    colnames(df.MA) <- c('chrom', 'chromStart', 'chromEnd')
    df.MA
  }

