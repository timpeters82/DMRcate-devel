DMR.plot <- function(dmrcoutput, dmr, betas, phen.col,
                     annotation=c(array="IlluminaHumanMethylation450k",
                                  annotation="ilmn12.hg19"),
                     samps=NULL, toscale=FALSE, plotmedians=FALSE, ...)
{

  stopifnot(is(dmrcoutput, 'dmrcate.output')) 
  stopifnot(is.matrix(betas))
  stopifnot((length(dmr) == 1) && (dmr %in% 1:nrow(dmrcoutput$results)))
  stopifnot(all(samps %in% 1:ncol(betas)))
  stopifnot(ncol(betas) == length(phen.col))
  if(dmrcoutput$results$no.probes[dmr] < 2){
    stop("Region must have 2 or more CpGs")
  }
  coords <- dmrcoutput$results$hg19coord[dmr]
  chr <- sub(":.*", "", coords)
  bookends <- sub(".*:", "", coords)
  startcpg <- as.integer(sub("-.*", "", bookends))
  stopcpg <- as.integer(sub(".*-", "", bookends))
  RSobject <- RatioSet(betas, annotation=annotation)
  RSanno <- getAnnotation(RSobject)
  cpgs <- rownames(RSanno)[RSanno$chr %in% chr &
                             RSanno$pos >= startcpg & RSanno$pos <= stopcpg]
  cpgs <- cpgs[order(RSanno[cpgs,"pos"])]
  if(is.null(samps)){samps <- 1:ncol(betas)}
  betas <- betas[as.character(cpgs), samps]
  m <- match(as.character(cpgs), rownames(RSanno))
  clusta <- data.frame(
    gene=RSanno$UCSC_RefGene_Name[m],
    group=RSanno$UCSC_RefGene_Group[m], pos=RSanno$pos[m])
  if(!toscale) {
    plot(1:nrow(clusta), , type='n', xlab=paste(chr, "consecutive probes"),
         ylab="beta values", ylim=c(0,1.15), ...)
    for(i in 1:nrow(clusta)) {
      points(rep(i, ncol(betas)), betas[i,], col=phen.col, ...)
    } 
    
    if(plotmedians){
      medians <- matrix(0, nrow(clusta), ncol=length(unique(phen.col)))
      colnames(medians) <- unique(phen.col)
      for(j in colnames(medians)) {
        for (i in 1:nrow(clusta)){
          medians[i,j] <- median(betas[i,phen.col==j])
        }
        lines(1:nrow(clusta), medians[,j], col=j, lwd=2)
      }
    }
      
    abline(1, 0, col='gray')
    abline(0.5, 0, col='gray')
    abline(0, 0, col='gray')
    empty <- rep(-1, nrow(clusta))
    empty[grep("5'UTR", clusta$group)] <- 1
    points(1:length(empty), empty, adj=0.6, col='blue4', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("3'UTR", clusta$group)] <- 1.01
    points(1:length(empty), empty, adj=0.8, col='cyan', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("Body", clusta$group)] <- 1.02
    points(1:length(empty), empty, col='firebrick1', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS1500", clusta$group)] <- 1.03
    points(1:length(empty), empty, adj=0.2, col='chartreuse', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS200", clusta$group)] <- 1.04
    points(1:length(empty), empty, adj=0.4, col='forestgreen', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("1stExon", clusta$group)] <- 1.05
    points(1:length(empty), empty, adj=1, col='magenta', pch=15)
    empty <- rep('', nrow(clusta))
    geneidxs <- cumsum(rle(as.character(clusta$gene))$lengths)
    geneidxs <- ceiling(c(geneidxs[1]/2, geneidxs[2:length(geneidxs)] -
                            diff(geneidxs)/2))
    empty[geneidxs] <- as.character(clusta$gene[geneidxs])
    genevector <- empty
    text(1:nrow(clusta), 1.03, genevector, pos=3, offset=1)
    }
  if (toscale) {
    plot(clusta$pos, 1:nrow(clusta), type='n',
         xlab=paste(chr, "hg19 coords"), ylab="beta values",
         ylim=c(0,1.15), ...)
    for(i in 1:nrow(clusta)){
      points(rep(clusta$pos[i], ncol(betas)), betas[i,], col=phen.col, ...)
    } 
    
    if(plotmedians){
      medians <- matrix(0, nrow(clusta), ncol=length(unique(phen.col)))
      colnames(medians) <- unique(phen.col)
      for(j in colnames(medians)) {
        for (i in 1:nrow(clusta)){
          medians[i,j] <- median(betas[i,phen.col==j])
        }
        lines(clusta$pos, medians[,j], col=j, lwd=2)
      }
    }
    
    abline(1, 0, col='gray')
    abline(0.5, 0, col='gray')
    abline(0, 0, col='gray')
    empty <- rep(-1, nrow(clusta))
    empty[grep("5'UTR", clusta$group)] <- 1
    points(clusta$pos, empty, adj=0.6, col='blue4', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("3'UTR", clusta$group)] <- 1.01
    points(clusta$pos, empty, adj=0.8, col='cyan', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("Body", clusta$group)] <- 1.02
    points(clusta$pos, empty, col='firebrick1', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS1500", clusta$group)] <- 1.03
    points(clusta$pos, empty, adj=0.2, col='chartreuse', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("TSS200", clusta$group)] <- 1.04
    points(clusta$pos, empty, adj=0.4, col='forestgreen', pch=15)
    empty <- rep(-1, nrow(clusta))
    empty[grep("1stExon", clusta$group)] <- 1.05
    points(clusta$pos, empty, adj=1, col='magenta', pch=15)
    empty <- rep('', nrow(clusta))
    geneidxs <- cumsum(rle(as.character(clusta$gene))$lengths)
    geneidxs <- ceiling(c(geneidxs[1]/2, geneidxs[2:length(geneidxs)] -
                            diff(geneidxs)/2))
    empty[geneidxs] <- as.character(clusta$gene[geneidxs])
    genevector <- empty
    text(clusta$pos, 1.03, genevector, pos=3, offset=1)
  }
}
