getSimpleMethHeatMap <- function(x, k=100, scale="none", asSNPs=F, binary=F,
                                 ColSideColors=NULL, rotate=F, distm="jaccard",
                                 ...){

  ## for DNA methylation data, this works well
  ward <- function(x) hclust(x, method="ward")

  if(is(x, "SummarizedExperiment") || is(x, "RangedSummarizedExperiment")){ #{{{
    x <- keepSeqlevels(x, paste0("chr", 1:22))
    x <- x[grep("^rs", rownames(x), invert=TRUE), ] 
    x <- x[grep("^ch[1234567890]", rownames(x), invert=TRUE), ] 
    if(all(grepl("^cg", rownames(x)))) {
      xxx <- rmSNPandCH(getBeta(x), mafcut=0.01)
      sds <- rowSds(xxx, na.rm=TRUE)
      names(sds) <- rownames(xxx)
      rowordering <- names(sds)[order(sds, decreasing=TRUE)]
      xx <- xxx[head(rowordering, k), ] 
    } else { 
      sds <- rowSds(getBeta(x), na.rm=TRUE)
      names(sds) <- rownames(x)
      rowordering <- names(sds)[order(sds, decreasing=TRUE)]
      xx <- getBeta(x)[head(rowordering, k), ]
    }
    if("indicator" %in% names(colData(x)) && is.null(ColSideColors)) {
      ColSideColors <- x$indicator
    } # }}}
  } else { # {{{
    sds <- rowSds(data.matrix(x))
    rowordering <- order(sds, decreasing=TRUE)
    xx <- x[head(rowordering, k), ]
  } # }}}
  if(asSNPs != TRUE) { # {{{
    jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", 
                              "yellow", "#FF7F00", "red", "#7F0000")) # }}}
  } else { ## {{{ as SNPs 
    xx <- round(xx * 2)
    jet <- colorRampPalette(c("blue","yellow","red"))
  } # }}}

  if (rotate == TRUE && !is.null(ColSideColors)) { # {{{
    hfun <- function(X) 
      heatmap(X, col=jet(64), scale=scale, hclustfun=ward,
              RowSideColors=ColSideColors, ...) # }}}
  } else if (rotate != TRUE && !is.null(ColSideColors)) { # {{{
    hfun <- function(X)
      heatmap(X, col=jet(64), scale=scale, hclustfun=ward,
              ColSideColors=ColSideColors, ...) # }}}
  } else { # {{{
     hfun <- function(X, ...) 
      heatmap(X, col=jet(64), scale=scale, hclustfun=ward,
              ...)
  } # }}}

  X <- xx 
  if (rotate) X <- t(X)
  hfun(X, ...)

}
