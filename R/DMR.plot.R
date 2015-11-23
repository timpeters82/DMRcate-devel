DMR.plot <- function(ranges, dmr, CpGs, phen.col,
         genome=c("hg19", "hg38", "mm10"), array.annotation=c(array="IlluminaHumanMethylation450k",
                                                                      annotation="ilmn12.hg19"), 
         samps=NULL, ...)
{
  stopifnot(class(CpGs) %in% c("matrix", "GRanges"))
  stopifnot(dmr %in% 1:length(ranges))
  if(is.null(samps)){samps=1:length(phen.col)}
  group <- unique(names(phen.col))
  if(is.matrix(CpGs)){
    RSobject <- RatioSet(CpGs, annotation=array.annotation)
    RSanno <- getAnnotation(RSobject)
    RSanno <- RSanno[order(RSanno$chr, RSanno$pos),]
    CpGs <- CpGs[rownames(RSanno),]
    colnames(CpGs) <- paste(colnames(CpGs), ".C", sep='')
    cov <- matrix(1, nrow(CpGs), ncol(CpGs), dimnames = list(rownames(CpGs), sub(".C", ".cov", colnames(CpGs))))
    cpgs.ranges <- GRanges(RSanno$chr, IRanges(RSanno$pos, RSanno$pos))
    dummy <- matrix(0, nrow=nrow(CpGs), ncol=2*ncol(CpGs))
    dummy[,seq(1, 2*ncol(CpGs), 2)] <- CpGs
    dummy[,seq(2, 2*ncol(CpGs), 2)] <- cov
    colnames(dummy)[seq(2, 2*ncol(CpGs), 2)] <- colnames(cov)
    colnames(dummy)[seq(1, 2*ncol(CpGs), 2)] <- colnames(CpGs)
    values(cpgs.ranges) <- dummy 
  } else {
    stopifnot(length(colnames(values(CpGs)))/2 == length(phen.col))
    if(!all(gsub(".*\\.", "", colnames(values(CpGs))) == rep(c("C", "cov"), length(phen.col)))){
        stop("Error: Column names of values(ranges) might not be in the correct format. Must be c('<sample1>.C', '<sample1>.cov', '<sample2>.C', '<sample2>.cov'...) and so on for all samples.")
      }
    cpgs.ranges <- CpGs
  }
  ranges$ID <- paste0("DMR_", 1:length(ranges))
  ranges.reduce <- reduce(ranges+5000)
  dmrs.inplot <- ranges[ranges %over% ranges.reduce[subjectHits(findOverlaps(ranges[dmr], ranges.reduce))]]
  ranges.inplot <- ranges.reduce[ranges.reduce %over% dmrs.inplot]
  cpgs.ranges <- subsetByOverlaps(cpgs.ranges, ranges.inplot)
  methRatios <- as.data.frame(values(cpgs.ranges)[,grep("C$", colnames(values(cpgs.ranges)))])/
                as.data.frame(values(cpgs.ranges)[,grep("cov$", colnames(values(cpgs.ranges)))])
  methRatios <- GRanges(cpgs.ranges, mcols=methRatios)
  mcols(methRatios) <- mcols(methRatios)[samps]
  names(mcols(methRatios)) <- gsub("mcols.", "", gsub("*.C", "", names(mcols(methRatios))))
  phen.col <- phen.col[samps]
  
 
  dt.group <- lapply(unique(names(phen.col)), function(i)
        DataTrack(methRatios[,names(phen.col) %in% i], name=i, background.title=phen.col[i],
                  type="heatmap", showSampleNames=TRUE, ylim=c(0, 1), genome=genome,
                  gradient=c("blue", "white", "red")))
 
  dt.group <- c(dt.group, list(DataTrack(methRatios, groups=names(phen.col), type="b",
                                             aggregateGroups=TRUE, col=phen.col[sort(group)], ylim=c(0, 1), name="Group means")))
     
  switch(genome, 
         hg19={tx=tx.hg19},
         hg38={tx=tx.hg38},
         mm10={tx=tx.mm10}
  )
  extras <- list(AnnotationTrack(dmrs.inplot, name="DMRs", showFeatureId=TRUE, col=NULL, fill="purple", id=dmrs.inplot$ID, 
                                 fontcolor="black"))
  extras <- endoapply(extras, function(x) {
    chromosome(x) <- as.character(seqnames(methRatios[dmr]))
    x
  })
  values(cpgs.ranges) <- NULL
    basetracks <- list(IdeogramTrack(genome = "hg19", chromosome = as.character(seqnames(ranges.inplot))),
                       GenomeAxisTrack(),
                       GeneRegionTrack(subsetByOverlaps(tx, ranges.inplot), name = "Gene", showId=TRUE, geneSymbol=TRUE, symbol = subsetByOverlaps(tx, ranges.inplot)$gene_name, col=NULL, fill="lightblue", transcriptAnnotation = "symbol", shape="arrow"),
                       AnnotationTrack(cpgs.ranges, name="CpGs", fill="green", col=NULL, stacking="dense"))
    plotTracks(c(basetracks, extras, dt.group), from=start(ranges.inplot), to=end(ranges.inplot), ...)
  }
  
  
