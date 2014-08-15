KernelSums <-
  function(x, 
           y,
           scale,
           K 
           )
  {
    if(missing(K))       
    {
      K <- function(dx) exp(-dx*dx/2) ## Gaussian
      attr(K, "support") <- 5
      #attr(K, "support") <- 1
    }
    support <- attr(K, "support")
    
    stopifnot(all(diff(x) >= 0))
    
    maxdiff <- support*scale
    triples <- SparseDeltas(x, maxdiff = maxdiff)
    #View(triples)
    
    i <- triples[,'i']
    j <- triples[,'j']
    dx <- triples[,'dx']
    
    k <- K(dx/scale)
      
    levels <- seq(length(x))
    ifac <- factor(i, levels = levels)
    jfac <- factor(j, levels = levels)
    
    itab <- table(ifac)
    jtab <- table(jfac)
    
    izero <- which(itab == 0)
    jzero <- which(jtab == 0)
    
    K0 <- K(0)
    
    SUM <-
      function(iterms, jterms)
      {
        S1 <- tapply(iterms, ifac, sum)
        S1[izero] <- 0
        S2 <- tapply(jterms, jfac, sum)
        S2[jzero] <- 0
        #browser()
        #tmp <- cbind(S1, S2); View(tmp)        
        S1 + S2
      }

    #browser()
    sky <- K0*y + SUM(k*y[j], k*y[i])
    sk <- K0 + SUM(k, k)
    k2 <- k*k
    skk <- K0*K0 + SUM(k2, k2)
    
    cbind(sky, sk, skk)
  }

