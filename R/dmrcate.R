dmrcate <-
  function (object,
            lambda = 1000,
            C = NULL,
            pcutoff = "fdr",
            consec = FALSE,
            conseclambda = 10,
            betacutoff = NULL,
            min.cpgs = 2
)
  {
    ## Arguments
    stopifnot(is(object, "CpGannotated"))
    stopifnot(lambda >= 1)
    stopifnot(pcutoff == "fdr" | (0 <= pcutoff & pcutoff <= 1))
    stopifnot(C >= 0.2)
    if (consec & is.null(conseclambda)) {
      stop("Consecutive CpG bandwidth must be specified")
    }

    ## Modified 'object'
    object <-
      data.frame(ID = names(object@ranges),
                 weights = abs(object@ranges$stat),
                 CHR = seqnames(object@ranges),
                 pos = start(object@ranges),
                 rawpval = object@ranges$rawpval,
                 diff = object@ranges$diff,
                 indfdr = object@ranges$ind.fdr,
                 is.sig = object@ranges$is.sig
                 )
    # Order by position
    object <- object[order(object$CHR, object$pos),]
    # Remove any chromosome with exactly 1 probe
    if (any(table(object$CHR)==1)) {
      torm <- names(which(table(object$CHR)==1))
      object <- object[!object$CHR %in% torm,]
    }
    
    # Automatic bandwidth specification
    if (is.null(C) & !consec) {
      C = 2
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
         lapply(chr.unique,
               fitParallel,
               object = object,
               consec = consec,
               conseclambda = conseclambda,
               lambda = lambda,
               C = C
    )
    object <- rbind.fill(fitted)
    object <- object[!is.na(object$is.sig),]

    ## FDR stuff
    object$fdr <- p.adjust(object$raw, method = "BH")
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
    fn_HMpval <- function (x) 1/mean(1/x)
    fn_Fisher <- function (x) pchisq((sum(log(x))*-2), df=length(x)*2, lower.tail=FALSE)
    fn_max <- function(x) x[which.max(abs(x))]
    results <-
      data.frame(
        coord = coord.A,
        no.cpgs = no_cpg.A,
        min_smoothed_fdr = REGIONSTAT("fdr", min),
        Stouffer = REGIONSTAT("rawpval", fn_Stouffer),
        HMpval = REGIONSTAT("rawpval", fn_HMpval),
        Fisher = REGIONSTAT("rawpval", fn_Fisher),
        maxdiff = REGIONSTAT("diff", fn_max),
        meandiff = REGIONSTAT("diff", mean),
        row.names = seq(A),
        stringsAsFactors = FALSE
      )
    #Correct DMR-wise "significances"
    results$Stouffer <- p.adjust(results$Stouffer, method = "fdr")
    results$Fisher <- p.adjust(results$Fisher, method = "fdr")
    results$HMFDR <- p.adjust(results$HMpval, method = "fdr")
    # Order and filter DMRs
    
    keep <- (results$no.cpgs >= min.cpgs)
    results <- results[keep, ]
    if(!(is.null(betacutoff))){
      message("Warning: betacutoff only meaningful for Illumina array data or WGBS results from DSS::DMLtest().")
      keep <- (abs(results$meandiff) > betacutoff)
      results <- results[keep, ]
    }
    o <- order(results$min_smoothed_fdr, -results$no.cpgs)
    results <- results[o,]
    message("Done!")
    return(new("DMResults", coord=results$coord, no.cpgs=results$no.cpgs, min_smoothed_fdr=results$min_smoothed_fdr,
        Stouffer=results$Stouffer, HMFDR=results$HMFDR, Fisher=results$Fisher, maxdiff=results$maxdiff, meandiff=results$meandiff))
  }
