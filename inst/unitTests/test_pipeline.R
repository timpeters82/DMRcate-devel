testpipeline <- function(){
  data(dmrcatedata)
  myMs <- logit2(myBetas)
  # Test rmSNPandCH()
  #########################################
  checkException(rmSNPandCH(object=list()))
  checkException(rmSNPandCH(myMs, dist=-8))
  checkException(rmSNPandCH(myMs, dist=60))
  checkException(rmSNPandCH(myMs, mafcut=-0.05))
  checkException(rmSNPandCH(myMs, mafcut=1.05))
  checkTrue(nrow(rmSNPandCH(myMs, and=TRUE)) != nrow(rmSNPandCH(myMs, and=FALSE)))
  checkTrue(nrow(rmSNPandCH(myMs, rmcrosshyb=TRUE)) != nrow(rmSNPandCH(myMs, rmcrosshyb=FALSE)))
  checkTrue(nrow(rmSNPandCH(myMs, rmcrosshyb=TRUE)) < nrow(myMs))
  ###############################################################
  
  myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
  patient <- factor(sub("-.*", "", colnames(myMs)))
  type <- factor(sub(".*-", "", colnames(myMs)))
  design <- model.matrix(~patient + type) 
  nointerceptdesign <- design
  nointerceptdesign <- nointerceptdesign[,-1]
    
  #Test annotate()
  ##################################################################################
  checkException(cpg.annotate(datatype="array", object=list(), what="M", arraytype = "450K", analysis.type="differential", design=design, 
                          coef=39))
  checkException(cpg.annotate(datatype="array", myMs.noSNPs, what="M", arraytype = "550K", analysis.type="differential", design=design, 
                          coef=39))
  checkException(cpg.annotate(datatype="array", myMs.noSNPs, what="M", arraytype = "450K", analysis.type="deferential", design=design, 
                          coef=39))
  checkException(cpg.annotate(datatype="array", myMs.noSNPs, what="M", arraytype = "450K", analysis.type="differential", design=design,
                          coef=40))
  checkException(cpg.annotate(datatype="array", myMs.noSNPs, what="M", arraytype = "450K", analysis.type="differential", design=nointerceptdesign,
                          coef=39))
    
  #checkEquals(nrow(myMs.noSNPs), length(cpg.annotate(datatype="array", object=myMs.noSNPs, what="M", arraytype = "450K", analysis.type="differential", design=design, 
  #                                               coef=39)$ID))
  #checkEquals(class(cpg.annotate(datatype="array", myMs.noSNPs, what="M", arraytype = "450K", analysis.type="differential", design=design, 
  #                           coef=39)), "annot")
  ###################################################################################
  
  myannotation <- cpg.annotate("array", myMs.noSNPs, what="M", arraytype = "450K", analysis.type="differential", design=design, 
                           coef=39, pcutoff=0.01)
  
  #Test dmrcate()
  #########################################################################
  checkException(dmrcate(object=myannotation$ID))
  checkException(dmrcate(myannotation, lambda=-100))
  checkException(dmrcate(myannotation, p.adjust.method="BLEH"))
  checkException(dmrcate(myannotation, pcutoff=-0.01))
  checkException(dmrcate(myannotation, C=0.1))
  checkException(dmrcate(myannotation, pcutoff=1.01))
  checkException(dmrcate(myannotation, consec=TRUE, conseclambda=-1))
  checkIdentical(attributes(dmrcate(myannotation))$names, c("input", "results", "cutoff"))
  #########################################################################
  
  dmrcoutput <- dmrcate(myannotation, lambda=1000)
  
  results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
  groups <- c(Tumour="magenta", Normal="forestgreen")
  cols <- groups[as.character(type)]
  
  
  #Test DMR.plot
  ###############################################################################
  checkException(DMR.plot(dmrcoutput=dmrcoutputwrong, dmr=1, CpGs=myBetas, what="M", arraytype = "450K", phen.col=cols, genome="hg19"))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1.5, CpGs=myBetas, what="M", arraytype = "450K", phen.col=cols, genome="hg19"))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, CpGs=factor(c(1, 2)), what="M", arraytype = "450K", phen.col=cols, genome="hg19"))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, CpGs=myBetas, what="M", arraytype = "450K", phen.col=cols, genome="hg19", samps=c(1, 79)))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, CpGs=myBetas, what="M", arraytype = "450K", phen.col=cols, annotation=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg20")))
  badphen.col <- c(rep("orange", 39), rep("blue", 38))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, CpGs=myBetas, what="M", arraytype = "450K", genome="hg19", phen.col=badphen.col))
  ##############################################################################
}
