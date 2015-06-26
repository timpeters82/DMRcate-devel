dmrcateByChr <- function(sigprobes, pcutoff, lambda, parallel=FALSE) { 

  ## might be easier using SE seqnames
  byChr <- split(sigprobes, sigprobes$CHR) 
  if (parallel == TRUE) {
    do.call(rbind, mclapply(byChr, dmrcateChr, pcutoff=pcutoff, lambda=lambda))
  } else { 
    do.call(rbind, lapply(byChr, dmrcateChr, pcutoff=pcutoff, lambda=lambda))
  }

} # }}}
