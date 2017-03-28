Segment <-
  function(jump.k)
  {
    K <- 1 + length(jump.k)
    kbreak.a <- which(jump.k)
    kstart.A <- c(1, 1 + kbreak.a)
    kend.A <- c(kbreak.a, K)
    ksegments.A2 <-
      cbind(start = kstart.A,
            end = kend.A
      )
    ksegments.A2
  }
