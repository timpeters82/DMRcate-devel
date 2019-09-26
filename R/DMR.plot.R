DMR.plot <- function (ranges, dmr, CpGs, what = c("Beta", "M"), arraytype = c("EPIC", 
                                                                  "450K"), phen.col, genome = c("hg19", "hg38", "mm10"), ...) 
{
  env <- new.env(parent = emptyenv())
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  genome <- match.arg(genome)
  stopifnot(class(CpGs) %in% c("matrix", "BSseq", "GenomicRatioSet"))
  stopifnot(dmr %in% 1:length(ranges))
  group <- unique(names(phen.col))
  if (is(CpGs, "matrix") | is(CpGs, "GenomicRatioSet")) {
    if (is(CpGs, "matrix")) {
      if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
      if (arraytype == "EPIC") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
    }
    else {
      grset <- CpGs
    }
    CpGs <- getBeta(grset)
    RSanno <- getAnnotation(grset)
    RSanno <- RSanno[order(RSanno$chr, RSanno$pos), ]
    CpGs <- CpGs[rownames(RSanno), ]
    cpgs.ranges <- GRanges(RSanno$chr, IRanges(RSanno$pos, 
                                               RSanno$pos))
    values(cpgs.ranges) <- CpGs
    isbsseq <- FALSE
  }
  else {
    if (any(width(CpGs) > 1)) {
      stop("Error: all ranges in the BSseq object must be single nucleotides with width 1.")
    }
    if (is.null(rownames(colData(CpGs)))) {
      stop("Error: BSseq object must be annotated with colData with sample IDs as rownames of the data.frame.")
    }
    stopifnot(ncol(CpGs) == length(phen.col))
    cpgs.ranges <- CpGs
    isbsseq <- TRUE
  }
  ranges$ID <- paste0("DMR_", 1:length(ranges))
  ranges.reduce <- reduce(ranges + 5000)
  dmrs.inplot <- ranges[ranges %over% ranges.reduce[subjectHits(findOverlaps(ranges[dmr], 
                                                                             ranges.reduce))]]
  ranges.inplot <- ranges.reduce[ranges.reduce %over% dmrs.inplot]
  cpgs.ranges <- subsetByOverlaps(cpgs.ranges, ranges.inplot)
  if (isbsseq) {
    methRatios <- GRanges(seqnames(cpgs.ranges), ranges(cpgs.ranges), 
                          mcols = as.matrix(getCoverage(cpgs.ranges, type = "M"))/as.matrix(getCoverage(cpgs.ranges, 
                                                                                                        type = "Cov")))
  }
  else {
    methRatios <- cpgs.ranges
  }
  values(methRatios) <- as.matrix(values(methRatios))
  colnames(values(methRatios)) <- gsub("mcols.", "", colnames(values(methRatios)))
  dt.group <- lapply(unique(names(phen.col)), function(i) DataTrack(methRatios[, 
                                                                               names(phen.col) %in% i], name = i, background.title = phen.col[i], 
                                                                    type = "heatmap", showSampleNames = TRUE, ylim = c(0, 
                                                                                                                       1), genome = genome, gradient = c("blue", "white", 
                                                                                                                                                         "red")))
  dt.group <- c(dt.group, list(DataTrack(methRatios, groups = names(phen.col), 
                                         type = "smooth", aggregateGroups = TRUE, 
                                         aggregation = function (x) mean(x, na.rm=TRUE),
                                         col = phen.col[sort(group)], 
                                         ylim = c(0, 1), name = "Smoothed\n group means", na.rm=TRUE)))
  switch(genome, hg19 = {
    data(hg19.grt, envir = env);
    grt = env$hg19.grt
  }, hg38 = {
    data(hg38.grt, envir = env);
    grt = env$hg38.grt
  }, mm10 = {
    data(mm10.grt, envir = env);
    grt = env$mm10.grt
  })
  chromosome(grt) <- as.character(seqnames(methRatios)[1])
  extras <- list(AnnotationTrack(dmrs.inplot, name = "DMRs", 
                                 showFeatureId = TRUE, col = NULL, fill = "purple", id = dmrs.inplot$ID, 
                                 fontcolor = "black"))
  values(cpgs.ranges) <- NULL
  basetracks <- list(IdeogramTrack(genome = genome, chromosome = as.character(seqnames(ranges.inplot))), 
                     GenomeAxisTrack(), grt, AnnotationTrack(GRanges(seqnames(cpgs.ranges), 
                                                                     ranges(cpgs.ranges)), name = "CpGs", fill = "green", 
                                                             col = NULL, stacking = "dense"))
  suppressWarnings(plotTracks(c(basetracks, extras, dt.group), 
                              from = start(ranges.inplot), to = end(ranges.inplot), 
                              ...))
}