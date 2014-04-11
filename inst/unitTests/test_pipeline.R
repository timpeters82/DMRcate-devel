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
  
  #Test annotate()
  ##################################################################################
  checkException(annotate(object=list(), analysis.type="differential", design=design, 
                          coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01))
  checkException(annotate(myMs.noSNPs, annotation=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg20"), analysis.type="differential", design=design, 
                          coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="deferential", design=design, 
                          coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=40, diff.metric="FC", paired=TRUE, pcutoff=0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=39, diff.metric="FC", paired=FALSE, pcutoff=0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=39, diff.metric="FC", paired=TRUE, pcutoff=-0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=39, diff.metric="FC", paired=TRUE, pcutoff=1.01))
  checkException(annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=39, diff.metric="FC", paired=TRUE, betacutoff=0.2, pcutoff=0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=39, diff.metric="FC", paired=TRUE, betacutoff=-0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="differential", design=design,
                          coef=39, diff.metric="FC", paired=TRUE, betacutoff=1.01))
  checkException(annotate(myMs.noSNPs, analysis.type="variability", quantcut=-0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="variability", quantcut=1.01))
  checkException(annotate(myMs.noSNPs, analysis.type="hypermethylation", cut=-0.01))
  checkException(annotate(myMs.noSNPs, analysis.type="hypermethylation", cut=1.01))
  checkEquals(nrow(myMs.noSNPs), length(annotate(myMs.noSNPs, analysis.type="differential", design=design, 
                                                 coef=39, diff.metric="FC", paired=TRUE, pcutoff=1)$ID))
  checkEquals(class(annotate(myMs.noSNPs, analysis.type="differential", design=design, 
                             coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01)), "annot")
  ###################################################################################
  
  myannotation <- annotate(myMs.noSNPs, analysis.type="differential", design=design, 
                           coef=39, diff.metric="FC", paired=TRUE, pcutoff=0.01)
  
  #Test dmrcate()
  #########################################################################
  checkException(dmrcate(object=myannotation$ID))
  checkException(dmrcate(myannotation, bw=-100))
  checkException(dmrcate(myannotation, p.adjust.method="BLEH"))
  checkException(dmrcate(myannotation, pcutoff=-0.01))
  checkException(dmrcate(myannotation, pcutoff=1.01))
  checkException(dmrcate(myannotation, consec=TRUE, consecbw=-1))
  checkIdentical(attributes(dmrcate(myannotation))$names, c("input", "results", "cutoff"))
  #########################################################################
  
  dmrcoutput <- dmrcate(myannotation, bw=1000)
  
  #Test makeBedgraphs()
  ################################################################################
  dmrcoutputwrong <- dmrcoutput
  attributes(dmrcoutputwrong)$names <- c("input", "results", "pcutoff")
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
