getCollapsedBetas <- function (annot, ranges=NULL) 
{
  if(!is(annot, "CpGannotated")){
    stop("Error: annot is not a CpGannotated object. Please create one with cpg.annotate() or sequencing.annotate()")
  }
  if(!is.null(ranges)){
    if(!is(ranges, "GenomicRanges")){
      stop("Error: ranges is not a GRanges object")
    }
  }
  betas <- annot@betas
  if(!is.null(ranges)){
    idx <- annot@ranges %over% ranges
    if(sum(idx)==0){
      stop("Error: Ranges given do not overlap any CpGs in the given CpGannotated object.")
    }
    betas <- betas[idx,]
  }
  betas
}