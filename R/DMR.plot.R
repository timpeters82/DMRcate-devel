DMR.plot <- function (ranges, dmr, CpGs, what = c("Beta", "M"), arraytype = c("EPICv2", "EPICv1", 
                                                                  "450K"), phen.col, genome = c("hg19", "hg38", "mm10"), labels = names(ranges), 
          group.means = FALSE, extra.ranges = NULL, extra.title = names(extra.ranges)) 
{
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  genome <- match.arg(genome)
  stopifnot(class(CpGs)[1] %in% c("matrix", "BSseq", "GenomicRatioSet"))
  if(arraytype=="EPICv2" & genome=="hg19"){
    stop("Error: genome must be hg38 for EPICv2 data.")
  }
  stopifnot(dmr %in% 1:length(ranges))
  IDs <- unique(names(phen.col))
  if (is(CpGs, "matrix") | is(CpGs, "GenomicRatioSet")) {
    if (is(CpGs, "matrix")) {
      if (arraytype == "450K") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
      if (arraytype == "EPICv1") {
        grset <- makeGenomicRatioSetFromMatrix(CpGs, 
                                               array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19", 
                                               mergeManifest = TRUE, what = what)
      }
      if (arraytype == "EPICv2") {
        if (!is(CpGs, "matrix")){
          stop("Error, CpGs argument must be a matrix for EPICv2")
        }
      }
    }
    else {
      grset <- CpGs
    }
    
    if (arraytype == "EPICv2") {
      message("Loading EPICv2 manifest...")
      ah <- AnnotationHub()
      EPICv2manifest <- ah[["AH116484"]]
      EPICv2manifest <- EPICv2manifest[EPICv2manifest$CHR!="chr0",]
      RSanno <- GRanges(EPICv2manifest$CHR, IRanges(EPICv2manifest$MAPINFO, EPICv2manifest$MAPINFO))
      names(RSanno) <- rownames(EPICv2manifest)
      RSanno <- sort(sortSeqlevels(RSanno))
      RSanno <- RSanno[names(RSanno) %in% rownames(CpGs)]
      CpGs <- CpGs[names(RSanno),]
      cpgs.ranges <- RSanno
      
    } else {
          CpGs <- getBeta(grset)
          RSanno <- getAnnotation(grset)
          RSanno <- RSanno[order(RSanno$chr, RSanno$pos), ]
          CpGs <- CpGs[rownames(RSanno), ]
          cpgs.ranges <- GRanges(RSanno$chr, IRanges(RSanno$pos, 
                                                     RSanno$pos))
          
    }
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
    isbsseq <- TRUE
  }
  ranges$ID <- rep("", length(ranges))
  ranges.reduce <- reduce(ranges + 5000)
  dmrs.inplot <- ranges[queryHits(findOverlaps(ranges, ranges.reduce[subjectHits(findOverlaps(ranges[dmr], 
                                                                             ranges.reduce))]))]
  ranges.inplot <- ranges.reduce[queryHits(findOverlaps(ranges.reduce, dmrs.inplot))]
  if (is(CpGs, "matrix")) {
    cpgs.ranges <- subsetByOverlaps(cpgs.ranges, ranges.inplot)
  }
  else {
    cpgs.ranges <- subsetByOverlaps(CpGs, ranges.inplot)
  }
  genome(cpgs.ranges) <- genome
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
                                                                                                                       1), genome = genome, gradient = colorRampPalette(c("black", 
                                                                                                                                                                          "orange"))(20)))
  if (!group.means) {
    smoothedmeans <- DataTrack(methRatios, groups = 1:ncol(values(methRatios)), 
                               type = "smooth", col = phen.col, ylim = c(0, 1), 
                               name = "Smoothed methylation", na.rm = TRUE, legend = FALSE)
    suppressMessages(setPar(smoothedmeans, "groupAnnotation", 
                            "feature"))
    dt.group <- c(dt.group, list(smoothedmeans))
  }
  else {
    smoothedmeans <- DataTrack(methRatios, groups = names(phen.col), 
                               type = "smooth", aggregateGroups = TRUE, aggregation = function(x) mean(x, 
                                                                                                       na.rm = TRUE), col = phen.col[sort(unique(names(phen.col)))], 
                               ylim = c(0, 1), name = "Smoothed\n group means", 
                               na.rm = TRUE)
    dt.group <- c(dt.group, list(smoothedmeans))
  }
  suppressWarnings(switch(genome, hg19 = {
    ensembl <- useEnsembl(host = "https://grch37.ensembl.org", 
                          biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    grt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = as.character(seqnames(methRatios)[1]), 
                                  start = min(start(ranges.inplot)), end = max(start(ranges.inplot)), 
                                  name = "ENSEMBL", biomart = ensembl)
  }, hg38 = {
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    grt <- BiomartGeneRegionTrack(genome = "hg38", chromosome = as.character(seqnames(methRatios)[1]), 
                                  start = min(start(ranges.inplot)), end = max(start(ranges.inplot)), 
                                  name = "ENSEMBL", biomart = ensembl)
  }, mm10 = {
    ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    grt <- BiomartGeneRegionTrack(genome = "mm10", chromosome = as.character(seqnames(methRatios)[1]), 
                                  start = min(start(ranges.inplot)), end = max(start(ranges.inplot)), 
                                  name = "ENSEMBL", biomart = ensembl)
  }))
  suppressMessages(grt <- setPar(grt, "collapseTranscripts", 
                                 "meta"))
  suppressMessages(grt <- setPar(grt, "exonAnnotation", "symbol"))
  suppressMessages(grt <- setPar(grt, "fontcolor.item", "black"))
  suppressMessages(grt <- setPar(grt, "cex", 0.6))
  suppressMessages(grt <- setPar(grt, "rotation.item", 45))
  suppressMessages(grt <- setPar(grt, "rotation.title", 0))
  cpgs.track <- AnnotationTrack(GRanges(seqnames(cpgs.ranges), 
                                        ranges(cpgs.ranges)), name = "CpGs", fill = "green", 
                                stacking = "dense", rotation.title = 0)
  suppressMessages(cpgs.track <- setPar(cpgs.track, "lty", 
                                        "blank"))
  if (all(is.null(names(dmrs.inplot)))) {
    names(dmrs.inplot) <- rep("DMR", length(dmrs.inplot))
  }
  dmrs.track <- AnnotationTrack(dmrs.inplot, name = "DMRs", 
                                showFeatureId = TRUE, fill = "purple", id = names(dmrs.inplot), 
                                fontcolor = "white", rotation.title = 0)
  if (!is.null(extra.ranges)) {
    extra.ranges <- extra.ranges[extra.ranges %over% dmrs.inplot]
    extras.track <- AnnotationTrack(extra.ranges, showFeatureId = TRUE, 
                                    name = extra.title, fill = "pink", id = names(extra.ranges), 
                                    rotation.title = 0)
    basetracks <- list(IdeogramTrack(genome = genome, chromosome = as.character(seqnames(ranges.inplot))), 
                       GenomeAxisTrack(), grt, cpgs.track, dmrs.track, 
                       extras.track)
  }
  else {
    basetracks <- list(IdeogramTrack(genome = genome, chromosome = as.character(seqnames(ranges.inplot))), 
                       GenomeAxisTrack(), grt, cpgs.track, dmrs.track)
  }
  trackstoplot <- c(basetracks, dt.group)
  suppressWarnings(plotTracks(trackstoplot, from = min(start(ranges.inplot)), 
                              to = max(end(ranges.inplot))))
}
