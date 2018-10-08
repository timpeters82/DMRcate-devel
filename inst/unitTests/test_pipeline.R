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

}
