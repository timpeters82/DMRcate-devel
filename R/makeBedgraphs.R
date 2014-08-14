makeBedgraphs <- function(dmrcoutput, betas,
                          annotation=c(array="IlluminaHumanMethylation450k",
                                       annotation="ilmn12.hg19"),
                          samps=NULL)
{
  stopifnot(is(dmrcoutput, 'dmrcate.output'))
  stopifnot(is.matrix(betas))
  stopifnot(all(samps %in% 1:ncol(betas)))
  if (is.null(samps)) {
    samps <- 1:ncol(betas)
  }
  sigprobes <- dmrcoutput$input[dmrcoutput$input$fdr < dmrcoutput$cutoff,]
  betas <- betas[as.character(sigprobes$ID),]
  RSobject <- RatioSet(betas, annotation=annotation)
  RSanno <- getAnnotation(RSobject)
  probes <- matrix(nrow=nrow(betas), ncol=3)
  ## Taking strandedness into account
  probes[RSanno$strand == "-",] <-
    cbind(as.character(sigprobes$CHR)[RSanno$strand == "-"],
          as.integer(as.character(sigprobes$pos))[RSanno$strand == "-"],
          as.integer(as.character(sigprobes$pos))[RSanno$strand == "-"] + 1)
  probes[RSanno$strand == "+",] <-
    cbind(as.character(sigprobes$CHR)[RSanno$strand == "+"],
          as.integer(as.character(sigprobes$pos))[RSanno$strand == "+"] - 1,
          as.integer(as.character(sigprobes$pos))[RSanno$strand == "+"])
  for(i in samps) {
    write.table(cbind(probes, betas[,i]),
                paste(colnames(betas)[i], ".bedGraph", sep=''), sep='\t',
                row.names=FALSE, col.names=FALSE, quote=FALSE)
  } 
  message(length(samps), " bedGraphs written to ", getwd(), ".", sep='')
}