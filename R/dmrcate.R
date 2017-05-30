dmrcate <-
  function (object,
            lambda = 1000,
            C = NULL,
            p.adjust.method = "BH",
            pcutoff = "fdr",
            consec = FALSE,
            conseclambda = 10,
            betacutoff = NULL,
            min.cpgs = 2,
            mc.cores = 1
  )
  {
    ## Arguments
    stopifnot(is(object, "annot"))
    stopifnot(lambda >= 1)
    stopifnot(pcutoff == "fdr" | (0 <= pcutoff & pcutoff <= 1))
    stopifnot(C >= 0.2)
    if (consec & is.null(conseclambda)) {
      stop("Consecutive CpG bandwidth must be specified")
    }

    ## Modified 'object'
    object <-
      data.frame(ID = object$ID,
                 weights = abs(object$stat),
                 CHR = as.character(object$CHR),
                 pos = object$pos,
                 betafc = object$betafc,
                 indfdr = object$indfdr,
                 is.sig = object$is.sig
                 )
    # Order by position
    object <- object[order(object$CHR, object$pos),]

    # Automatic bandwidth specification
    if (is.null(C) & !consec) {
      if (nrow(object) < 900000) {
        C = 2
      }
      else {
        C = 50
      }
    }

    ## Handle 'consec' case
    if (consec)
      {
      lambda = conseclambda
      message(paste("Consecutive mode specified, lambda is now set at", conseclambda, "consecutive CpGs."))
      if (is.null(C)){
        stop("Error: argument C must be specified (in CpG sites) for consecutive mode.")
      }
      object$realcoordforconsec <- object$pos
      object$pos <- unlist(sapply(as.numeric(table(object$CHR)), function (x) 1:x))
    }

    ## Kernel (chi-squared) test via 'fitParallel'
    lag = lambda
    chr.unique <- unique(c(as.character(object$CHR)))
    fitted <-
      mclapply(chr.unique,
               fitParallel,
               object = object,
               consec = consec,
               conseclambda = conseclambda,
               lambda = lambda,
               C = C,
               mc.cores = mc.cores
      )
    object <- rbind.fill(fitted)

    ## FDR stuff
    object$fdr <- p.adjust(object$raw, method = p.adjust.method)
    if (pcutoff == "fdr")
      {
      nsig <- sum(object$is.sig)
      if (nsig == 0)
      {
        txt <- "The FDR you specified in cpg.annotate() returned no significant CpGs, hence there are no DMRs.\n    Try specifying a value of 'pcutoff' in dmrcate() and/or increasing 'fdr' in cpg.annotate()."
        stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
      }
      pcutoff <- sort(object$fdr)[nsig]
    }
    object$sig <- (object$fdr <= pcutoff)
    if (nrow(object) == 0)
    {
      txt <- "No signficant regions found. Try increasing the value of\n    'pcutoff' in dmrcate() and/or 'fdr' in cpg.annotate()."
      stop(paste(strwrap(txt, exdent = 2), collapse = "\n"))
    }

    ## Segmentation
    message("Demarcating regions...")
    # Define jump.k
    # K = number of significant CpGs
    # k = K - 1
    chr.N <- as.character(object$CHR)
    pos.N <- object$pos
    sig.N <- object$sig
    N <- length(sig.N)
    n.K <- which(sig.N)
    K <- length(n.K)
    stopifnot(K >= 2)
    pos.K <- pos.N[n.K]
    chr.K <- chr.N[n.K]
    jump_chr.k <- (chr.K[-1] != chr.K[-K])
    jump_pos.k <- (diff(pos.K) > lag)
    jump.k <- (jump_chr.k | jump_pos.k)
    # Segment using jump.k
    ksegments.A2 <- Segment(jump.k)
    A <- nrow(ksegments.A2)
    # Extract start/end indices
    kstart.A <- ksegments.A2[,"start"]
    kend.A <- ksegments.A2[,"end"]

    ## Regionwise stats
    # Positions
    realpos.K <- pos.K
    if(consec)
    {
      realpos.N <- object$realcoordforconsec
      realpos.K <- realpos.N[n.K]
    }

    # Per-DMR: Coordinates
    start.A <- realpos.K[kstart.A]
    end.A <- realpos.K[kend.A]
    chr.A <- chr.K[kstart.A]
    stopifnot(all(chr.K[kend.A] == chr.A))
    fmt <- "%s:%1d-%1d"
    coord.A <- sprintf(fmt, chr.A, start.A, end.A)

    # Region factor
    nstart.A <- n.K[kstart.A]
    nend.A <- n.K[kend.A]
    width.A <- nend.A + 1 - nstart.A
    a.Z <- rep(seq(A), width.A) # a.Z
    fn <-
      function(a)
        seq(from = nstart.A[a], to = nend.A[a])
    l.listA <- lapply(seq(A), fn)
    n.Z <- unlist(l.listA)
    region.N <- rep(NA_integer_, N)
    region.N[n.Z] <- a.Z
    levels <- seq(A)
    region.N <- factor(region.N, levels = levels)

    # Per-DMR: Number of CpGs
    no_cpg.A <- c(table(region.N))
    # Function to do regionwise summaries
    REGIONSTAT <-
      function(field,
               fn
               )
      {
        x.N <- object[[field]]
        x.R <- tapply(x.N, region.N, fn)
        c(x.R)
      }
    # results <- region-wise stats
    fn_Stouffer <- function(x) pnorm(sum(qnorm(x))/sqrt(length(x)))
    fn_max <- function(x) x[which.max(abs(x))]
    results <-
      data.frame(
        coord = coord.A,
        no.cpgs = no_cpg.A,
        minfdr = REGIONSTAT("fdr", min),
        Stouffer = REGIONSTAT("indfdr", fn_Stouffer),
        maxbetafc = REGIONSTAT("betafc", fn_max),
        meanbetafc = REGIONSTAT("betafc", mean),
        row.names = seq(A),
        stringsAsFactors = FALSE
      )

    # Order and filter DMRs
    
    keep <- (results$no.cpgs >= min.cpgs)
    results <- results[keep, ]
    if (!is.null(betacutoff))
    {
      keep <- (abs(results$meanbetafc) >= betacutoff)
      results <- results[keep,]
    }
    o <- order(results$Stouffer, -results$no.cpgs)
    results <- results[o,]
    message("Done!")


    ## Output list
    output <- NULL
    output$input <- object
    output$results <- results
    output$cutoff <- pcutoff
    class(output) <- "dmrcate.output"
    output
  }
