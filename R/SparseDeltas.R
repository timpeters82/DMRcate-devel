SparseDeltas <-
  function(x, maxdiff)
  {
    n <- length(x)
    i <- seq(n)
    xout <- x + maxdiff
    j <- floor(approx(x = x, y = i, xout = xout, rule = 2)$y)
    
    s <- (j > i)
    i <- i[s]
    j <- j[s]
    
    d <- j - i
    ii <- rep(i, d)
    jj <- ii + unlist(lapply(d, seq))
    
    dx <- x[jj] - x[ii]
    triples <- cbind(i = ii, j = jj, dx = dx)

    triples
  }

