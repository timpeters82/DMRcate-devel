dmrcate <-
  function(object, 
           lambda = 1000,
           C=NULL,
           p.adjust.method = "BH", 
           pcutoff = "fdr", 
           consec = FALSE, 
           conseclambda = 10, 
           betacutoff = NULL,
           min.cpgs=2,
           mc.cores=1
           ) 
  {
    
    ## Checks
    stopifnot(is(object, "annot"))
    stopifnot(lambda >= 1)
    stopifnot(pcutoff=="fdr" | (0 <= pcutoff & pcutoff <= 1))
    stopifnot(C >= 0.2)
    if (consec & is.null(conseclambda)) 
    {
      stop("Consecutive CpG bandwidth must be specified")
    }
    
    ## 'object' as data.frame
    object <- data.frame(ID=object$ID, 
                         weights=abs(object$stat), 
                         CHR=as.character(object$CHR),
                         pos=object$pos, 
                         betafc=object$betafc, 
                         indfdr=object$indfdr)
    if(is.null(C)){
      if(nrow(object) < 485513){
          C=2} else {
          C=50
        }
    }
    ## Loop over chromosomes
    lag = lambda
    chr.unique <- unique(c(as.character(object$CHR)))
    fitted <- mclapply(chr.unique, fitParallel, object=object, consec=consec,
                       conseclambda=conseclambda, lambda=lambda, C=C, mc.cores=mc.cores)
    object <- rbind.fill(fitted)
    ## FDR correction
    object$fdr <- p.adjust(object$raw, method = p.adjust.method)
    if (pcutoff=="fdr"){
      nsig <- sum(object$indfdr < 0.05)
      pcutoff <- sort(object$fdr)[nsig]
    }
    
    object$sig <- object$fdr <= pcutoff
    
    
    if (nrow(object) == 0) 
    {
      txt <- "No signficant regions found. Try increasing the value of\n    'pcutoff' in dmrcate() and/or 'fdr' in cpg.annotate()."
      stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
    }
    message("Demarcating regions...")
    ## Demarcate regions
    # Order: By chromosome, then position
    chr <- as.character(object$CHR)
    pos <- object$pos
    o <- order(chr, pos)
    object <- object[o,]
    # Identify region starts
    chr <- as.character(object$CHR)
    pos <- object$pos
    stopifnot(sum(object$sig) >= 2)
    bychr <- function(i){
      this.chr <- object[object$CHR %in% i,]
      n <- nrow(this.chr)    
      s <- seq(n-1)
      w <- which(this.chr$sig)
      step.pos <- rep(NA, length(w)-1)
      step.pos <- (pos[w][-1] - pos[w][-length(w)] > lag)
      step.pos <- c(TRUE, step.pos)
      step.dmr <- rep(NA, n)
      step.dmr[w] <- step.pos
      transitions <- paste(step.dmr[-length(step.dmr)], step.dmr[-1])
      gaps <- grep("NA FALSE", transitions)
      for (j in 1:length(gaps)){
        idx <- gaps[j]
        while(is.na(step.dmr[idx])){
          step.dmr[idx] <- FALSE
          idx <- idx - 1
        }
      }
      this.chr$step.dmr <- step.dmr
      this.chr
    }
    
    object <- lapply(unique(object$CHR), bychr)
    object <- rbind.fill(object)
    sigprobes <- object[!is.na(object$step.dmr),]
    
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
    no.cpgs <- c(table(region))
    maxbetafc <- c(tapply(sigprobes$betafc, region,
                          function (x) x[abs(x)==max(abs(x))]))
    maxbetafc <- sapply(maxbetafc, max) # Debug for ties
    minfdr <- c(tapply(sigprobes$fdr, region, min))
    Stouffer <- c(tapply(sigprobes$indfdr, region, function (x) pnorm(sum(qnorm(x)) / sqrt(length(x)))))
    meanbetafc <- c(tapply(sigprobes$betafc, region, mean))
    
    # Coords
    U <- function(x) { u <- unique(x); stopifnot(length(u) == 1); u }
    chr <- tapply(chr, region, U)
    start <- tapply(pos, region, min)
    end <- tapply(pos, region, max)
    fmt <- "%s:%1d-%1d"
    coord <- sprintf(fmt, chr, start, end)
    # Groups
    P <- function(x) 
    { 
      l <- strsplit(x, ";")
      u <- unique(unlist(l))
      paste(u, collapse = ",") 
    }
    results <- 
      data.frame(coord = coord,
                 no.cpgs = no.cpgs,
                 minfdr = minfdr,
                 Stouffer = Stouffer,
                 maxbetafc = maxbetafc,
                 meanbetafc = meanbetafc,
                 row.names = seq(R),
                 stringsAsFactors = FALSE
      )
    results <- results[order(Stouffer, -no.cpgs),, drop = FALSE]
    results <- results[results$no.cpgs >= min.cpgs,]
    ## Wrap up
    message("Done!")
    output <- NULL
    output$input <- object
    output$results <- results
    output$cutoff <- pcutoff
    class(output) <- "dmrcate.output"
    output
}