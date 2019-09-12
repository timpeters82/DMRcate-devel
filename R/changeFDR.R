changeFDR <- function (annot, FDR) 
{
  if(!is(annot, "CpGannotated")){
    stop("Error: annot is not a CpGannotated object. Please create one with cpg.annotate() or sequencing.annotate()")
  }
  if(FDR <=0 | FDR >=1){
    stop("Error: please enter an appropriate FDR value > 0 or < 1.")
  }
  annot@is.sig <- annot@ind.fdr < FDR
  cat(paste0("Threshold is now set at FDR=", FDR, ", resulting in ", 
                 sum(annot@is.sig), " significantly differential CpGs."))
  annot
  
}