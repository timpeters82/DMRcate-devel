KernelTest <-
  function(pos,    # x values
           X2,     # Chi-square values, each with 'df' degrees of freedom
           lambda,     # Smoothing bandwidth
           df = 1  # 'nu'
           )
  {
    # Local kernel sums
    j <- KernelSums(x = pos, y = X2, scale = lambda)
    sky <- j[,'sky']
    sk <- j[,'sk']
    skk <- j[,'skk']
    
    # Satterthwaite approximation
    a <- df*skk/sk
    b <- sk*sk/skk
    
    outX2 <- sky/a # Approx. chi-square with b degrees of freedom
    
    P <- pchisq(outX2, df = b, lower.tail = FALSE)
    
    P
  }
