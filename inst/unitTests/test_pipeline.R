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
  c.matrix <- makeContrasts(patient2675 - patient2679, levels=nointerceptdesign)
  
  #Test annotate()
  ##################################################################################
  checkException(cpg.annotate(object=list(), analysis.type="differential", design=design, 
                          coef=39))
  checkException(cpg.annotate(myMs.noSNPs, annotation=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg20"), analysis.type="differential", design=design, 
                          coef=39))
  checkException(cpg.annotate(myMs.noSNPs, analysis.type="deferential", design=design, 
                          coef=39))
  checkException(cpg.annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=40))
  checkException(cpg.annotate(myMs.noSNPs, analysis.type="differential", design=nointerceptdesign,
                          coef=39))
  checkException(cpg.annotate(myMs.noSNPs, analysis.type="differential", design=nointerceptdesign,
                              contrasts=TRUE, cont.matrix=c.matrix, coef=39, pcutoff=1.01))
  
  checkEquals(nrow(myMs.noSNPs), length(cpg.annotate(myMs.noSNPs, analysis.type="differential", design=design, 
                                                 coef=39)$ID))
  checkEquals(class(cpg.annotate(myMs.noSNPs, analysis.type="differential", design=design, 
                             coef=39)), "annot")
  ###################################################################################
  
  myannotation <- cpg.annotate(myMs.noSNPs, analysis.type="differential", design=design, 
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
  
  #Test makeBedgraphs()
  ################################################################################
  dmrcoutputwrong <- dmrcoutput
  class(dmrcoutputwrong) <- "list"
  checkException(makeBedgraphs(dmrcoutput=dmrcoutputwrong, betas=myBetas, samps=c(1, 39)))
  checkException(makeBedgraphs(dmrcoutput=dmrcoutput, betas=list(), samps=c(1,39)))
  checkException(makeBedgraphs(dmrcoutput=dmrcoutput, betas=myBetas, annotation=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg20")))
  checkException(makeBedgraphs(dmrcoutput=dmrcoutput, betas=myBetas, samps=c(1,79)))
  ###############################################################################
  
  phen.col <- c(rep("orange", 38), rep("blue", 38))
  
  #Test DMR.plot
  ###############################################################################
  checkException(DMR.plot(dmrcoutput=dmrcoutputwrong, dmr=1, betas=myBetas, phen.col=phen.col))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1.5, betas=myBetas, phen.col=phen.col))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=list(), phen.col=phen.col))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=myBetas, phen.col=phen.col, samps=c(1, 79)))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=myBetas, phen.col=phen.col, annotation=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg20")))
  badphen.col <- c(rep("orange", 39), rep("blue", 38))
  checkException(DMR.plot(dmrcoutput=dmrcoutput, dmr=1, betas=myBetas, phen.col=badphen.col))
  ##############################################################################
}
