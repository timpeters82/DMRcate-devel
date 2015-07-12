dmrcate <-
  function(object, 
           lambda = 1000,
           C=2,
           p.adjust.method = "BH", 
           pcutoff = "limma", 
           consec = FALSE, 
           conseclambda = 10, 
           betacutoff = NULL
           ) 
  {
    
    ## Checks
    stopifnot(is(object, "annot"))
    stopifnot(lambda >= 1)
    stopifnot(pcutoff=="limma" | (0 <= pcutoff & pcutoff <= 1))
    stopifnot(C >= 0.2)
    if (consec & is.null(conseclambda)) 
    {
      stop("Consecutive probe bandwidth must be specified")
    }
    
    ## 'object' as data.frame
    object <- 
      data.frame(ID = object$ID, 
                 weights = object$weights, 
                 CHR = object$CHR, 
                 pos = object$pos, 
                 gene = object$gene, 
                 group = object$group, 
                 betafc = object$betafc,
                 indfdr = object$indfdr
      )    
    
    ## Loop over chromosomes
    chr.unique <- unique(c(as.character(object$CHR)))
    for (chr in chr.unique) 
    {
      message("Fitting ", chr, "...", sep = "")
      chrIndex <- object$CHR %in% chr      
      if (consec) 
      {
        stopifnot(conseclambda >= 1)
        lambda <- conseclambda
        pos <- seq_along(object$pos[chrIndex])
      }
      else 
      {
        pos <- object$pos[chrIndex]
      }
      sigma = lambda/C
      lag = lambda
      beta <- object$weights[chrIndex]
      df <- 1
      X2 <- beta^2
      pvalue <- KernelTest(pos = pos, X2 = X2, lambda = sigma, df = df)
      object$raw[chrIndex] <- pvalue
    }
    
    ## FDR correction
    object$fdr <- p.adjust(object$raw, method = p.adjust.method)
    if (pcutoff=="limma"){
      nsig <- sum(object$indfdr < 0.05)
      pcutoff <- sort(object$fdr)[nsig]
    }
    sigprobes <- object[object$fdr <= pcutoff, , drop = FALSE]
    
    if (nrow(sigprobes) == 0) 
    {
      txt <- "No significant regions found. Try increasing the value of\n 'pcutoff' in dmrcate()."
      stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
    }
    message("Demarcating regions...")
    ## Demarcate regions
    # Order: By chromosome, then position
    chr <- as.character(sigprobes$CHR)
    pos <- sigprobes$pos
    o <- order(chr, pos)
    sigprobes <- sigprobes[o,]
    # Identify region starts
    chr <- as.character(sigprobes$CHR)
    pos <- sigprobes$pos
    n <- nrow(sigprobes)    
    stopifnot(n >= 2)
    s <- seq(n-1)
    step.chr <- (chr[s] != chr[s+1])
    step.pos <- (pos[s+1] - pos[s] > lag)
    step <- (step.pos | step.chr)
    step <- c(TRUE, step)
    # Region factor
    region <- as.factor(cumsum(step))
    R <- length(levels(region))
    # Summaries
    no.probes <- c(table(region))
    maxbetafc <- c(tapply(sigprobes$betafc, region,
                          function (x) x[abs(x)==max(abs(x))]))
    maxbetafc <- sapply(maxbetafc, max) # Debug for ties
    minfdr <- c(tapply(sigprobes$fdr, region, min))
    meanfdr <- c(tapply(sigprobes$fdr, region, mean))
    # Coords
    U <- function(x) { u <- unique(x); stopifnot(length(u) == 1); u }
    chr <- tapply(chr, region, U)
    start <- tapply(pos, region, min)
    end <- tapply(pos, region, max)
    fmt <- "%s:%1d-%1d"
    hg19coord <- sprintf(fmt, chr, start, end)
    # Groups
    P <- function(x) 
    { 
      l <- strsplit(x, ";")
      u <- unique(unlist(l))
      paste(u, collapse = ",") 
    }
    group <- c(tapply(as.character(sigprobes$group), region, P))
    gene_assoc <- c(tapply(as.character(sigprobes$gene), region, P))
    results <- 
      data.frame(gene_assoc = gene_assoc,
                 group = group,
                 hg19coord = hg19coord,
                 no.probes = no.probes,
                 minpval = minfdr,
                 meanpval = meanfdr,
                 maxbetafc = maxbetafc,
                 row.names = seq(R),
                 stringsAsFactors = FALSE
      )
    results <- results[order(minfdr, -no.probes),, drop = FALSE]
    
    ## Wrap up
    message("Done!")
    output <- NULL
    output$input <- object
    output$results <- results
    output$cutoff <- pcutoff
    
    class(output) <- 'dmrcate.output'
    return(output)
}