## the crux of these manipulations: make an array directly comparable to RNAseq
setAs("ExpressionSet", "SummarizedExperiment", 
      function(from) eSetToSE(from)) 
setMethod("exprs","SummarizedExperiment", 
          function(object) assays(object)$exprs)
setMethod("pData", "SummarizedExperiment", 
          function(object) colData(object))
setMethod("fData", "SummarizedExperiment", 
          function(object) rowData(object))
setMethod("featureNames", "SummarizedExperiment", 
          function(object) rownames(object))
setMethod("sampleNames", "SummarizedExperiment", 
          function(object) colnames(object))

## keep track of metadata 
setAs("MIAME", "SimpleList",
  function(from) { 
    to = list()
    for(i in slotNames(from)) if(i != '.__classVersion__') to[[i]]=slot(from, i)
    return(SimpleList(to))
  }
)

## this will eventually require a better approach; add Illumina expr support 
## although to be honest, it works pretty damned well for Affy and HELP arrays
eSetToSE <- function(from) {
  
  ## short circuit the search for Affy hgu133plus2
  ## since it is a common case for my work 
  if (annotation(from) == 'GPL570') {
    ## Affy hgu133plus2, the generic leukemia microarray...
    ## I could add others but really, I just don't care. 
    annotation(from) <- 'hgu133plus2'
  }
  if (grepl('^GPL', annotation(from))) {
    message('This coercion expects to find at least an Entrez ID via fData()')
    if ( all( c('CHR','MAPINFO') %in% fvarLabels(from)) ) {
      fdat <- fData(from)[, c('ID','CHR','MAPINFO','MAPINFO') ]
      names(fdat) <- c('name','chrom','chromStart','chromEnd')
      fdat[, 3:4] <- apply(fdat[, 3:4], 2, as.numeric)
      row.dat <- makeGRangesFromDataFrame(fdat, keep.extra.columns=TRUE)
    } else if (all(c('CHROMOSOME','RANGE_START','RANGE_STOP')  ## e.g. HELP
                   %in% fvarLabels(from))) { 
      fdat <- fData(from)[, c('ID','CHROMOSOME','RANGE_START','RANGE_STOP') ]
      names(fdat) <- c('name','chrom','chromStart','chromEnd')
      fdat[, 3:4] <- apply(fdat[, 3:4], 2, as.numeric)
      row.dat <- makeGRangesFromDataFrame(fdat)
      seqlevelsStyle(row.dat) <- 'UCSC' 
      chrs <- c(extractSeqlevelsByGroup("Homo_sapiens", "UCSC", group="auto"),
                extractSeqlevelsByGroup("Homo_sapiens", "UCSC", group="sex"))
      row.dat <- keepSeqlevels(row.dat, chrs)
      if (unique(fData(from)$BUILD) == 'June_2004') {
        message('This looks like hg17 data, attempting to lift over to hg19...')
        data(hg17ToHg19) ## included for just this purpose...
        row.dat$probe <- names(row.dat)
        row.dat <- GenomicRanges::unlist(liftOver(row.dat, hg17ToHg19))
        names(row.dat) <- mcols(row.dat)[,1]
        mcols(row.dat) <- NULL
        ## liftOver is imperfect, so use a GRL...
        row.dat <- split(row.dat, names(row.dat))
      }
    } else {
      stop(paste("Don't know how to handle platform", annotation(from)))
    }
  } else {
    chip <- annotation(from)
    db <- paste(chip, "db", sep = ".")
    library(db, character.only=TRUE)
    fdat <- select(get(db), 
                   keys = featureNames(from), 
                   columns = c("CHRLOC","CHRLOCEND","SYMBOL"))
    fdat$strand <- with(fdat, ifelse(sign(as.numeric(CHRLOC))=='-1', '-', '+'))
    names(fdat) <- sub('CHRLOCCHR', 'chrom', names(fdat))
    names(fdat) <- sub('CHRLOCEND', 'chromEnd', names(fdat))
    names(fdat) <- sub('CHRLOC', 'chromStart', names(fdat))
    names(fdat) <- sub('SYMBOL', 'geneSymbol', names(fdat))
    names(fdat) <- sub('PROBEID', 'probe', names(fdat))
    cols <-  c('chrom','chromStart','chromEnd','strand','geneSymbol','probe')
    chrs <- c(extractSeqlevelsByGroup("Homo_sapiens", "NCBI", group="auto"),
              extractSeqlevelsByGroup("Homo_sapiens", "NCBI", group="sex"))
    fdat <- fdat[ which(fdat$chrom %in% chrs), cols ] 
    fdat$chromStart <- abs(fdat$chromStart)
    fdat$chromEnd <- abs(fdat$chromEnd)
    row.dat <- makeGRangesFromDataFrame(fdat, keep.extra.columns=TRUE)
    names(row.dat) <- row.dat$geneSymbol
    seqlevelsStyle(row.dat) <- 'UCSC'
    row.dat <- split(row.dat, row.dat$probe)

  }
  asy.dat <- SimpleList()
  asy.dat$exprs = assayDataElement(from, 'exprs')[names(row.dat), ]
  SummarizedExperiment(assays=asy.dat,
                       rowData=row.dat,
                       colData=as(pData(from), 'DataFrame'),
                       metadata=as(experimentData(from), 'SimpleList'))
}
