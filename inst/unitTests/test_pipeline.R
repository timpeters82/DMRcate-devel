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
    
  checkEquals(nrow(myMs.noSNPs), length(cpg.annotate(datatype="array", object=myMs.noSNPs, what="M", arraytype = "450K", analysis.type="differential", design=design, 
                                                 coef=39)$ID))
  checkEquals(class(cpg.annotate(datatype="array", myMs.noSNPs, what="M", arraytype = "450K", analysis.type="differential", design=design, 
                             coef=39)), "annot")
  ###################################################################################
 
}
