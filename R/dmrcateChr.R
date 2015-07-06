dmrcateChr <- function(byChr, pcutoff, lambda) { 

  ## this IS a per-chromosome function...
  chromosome <- unique(byChr$CHR)
  if (length(chromosome) != 1) return(NULL)

  ## here is where we loop through repeatedly iff length(pcutoff) > 1
  if (sum(pcutoff > 0) > 1) {
    do.call(rbind, 
            lapply(pcutoff[which(pcutoff > 0)] , 
                   function(p) dmrcateChr(byChr[byChr$fdr <= p, ], p, lambda)))
  } else { 
  
    ## "summits" are special
    if (nrow(byChr) < 2) {
      res <- with(byChr, 
                  data.frame(gene_assoc=gene,
                             group=group,
                             hg19coord=paste0(chromosome, ":", pos, "-", pos),
                             no.probes=1,
                             minpval=fdr,
                             meanpval=fdr,
                             maxbetafc=betafc,
                             pcutoff = rep(pcutoff, length(group)), 
                             stringsAsFactors = FALSE))

    } else { 

      message("Demarcating at ", pcutoff, " on ", chromosome, 
              " with lambda ", lambda, "...")
      sigprobes <- byChr[order(byChr$pos),]
      separateRegions <- c(TRUE, diff(sigprobes$pos) > lambda)
      region <- as.factor(cumsum(separateRegions))
      numRegions <- length(levels(region))

      # Summaries
      no.probes <- c(table(region))
      maxbetafc <- c(tapply(sigprobes$betafc, region,
                            function (x) x[abs(x) == max(abs(x))]))
      maxbetafc <- sapply(maxbetafc, max) # Debug for ties
      minfdr <- c(tapply(sigprobes$fdr, region, min))
      meanfdr <- c(tapply(sigprobes$fdr, region, mean))

      # Get genomic coordinates
      U <- function(x) { 
        u <- unique(x); 
        stopifnot(length(u) == 1)
        return(u)
      }

      # This could be turned into a GRanges right here...
      start <- tapply(sigprobes$pos, region, min)
      end <- tapply(sigprobes$pos, region, max)
      fmt <- "%s:%1d-%1d"
      hg19coord <- sprintf(fmt, rep(chromosome, numRegions), start, end)

      # Groups
      P <- function(x) paste(unique(unlist(strsplit(x, ";"))), collapse = ",")
      group <- c(tapply(as.character(sigprobes$group), region, P))
      gene_assoc <- c(tapply(as.character(sigprobes$gene), region, P))
      res <-
        data.frame(gene_assoc = gene_assoc,
                   group = group,
                   hg19coord = hg19coord,
                   no.probes = no.probes,
                   minpval = minfdr,
                   meanpval = meanfdr,
                   maxbetafc = maxbetafc,
                   pcutoff = rep(pcutoff, length(group)), 
                   stringsAsFactors = FALSE)
      res <- res[order(minfdr, -no.probes),, drop = FALSE]
    }
    return(res)
  }
}
