#' Generate model predictions based on the posterior
#'
#' @param runObj an object created by runJAGSmodel
#' @param predType the value to be predicted
#' @param at the x locations at which predictions should be made
#' @return a 3 dimensional array (sims, population, at value)
#' @examples
#' predictVal(m7,"smoltAtDat",seq(0,10000,by=100))
#' @export
predictVal <- function(runObj,predType="smoltAtTrap",at=NULL){
  x <- getPostDraws(runObj)
  bdat <- runObj$bdat
  nn <- length(xx)
  pVal <- array(NA,dim=c(numSims,bdat$Npops,nn))
  for(i in 1:bdat$Npops){
    pOut <- invLogit(x$pOutPopL[,i])
    for(j in 1:nn){
      predParr <- xx[j]/(1/x$prod[,i] + xx[j]/x$cap[,i])
      wsSurvMig <- invLogit(x$winterSurvMigPopL[,i]+x$winterSurvSlope[,i]*predParr/bdat$habVar[i])
      wsSurvRes <- invLogit(x$winterSurvResPopL[,i]+x$winterSurvSlope[,i]*predParr/bdat$habVar[i])
      if(predType=="parr"){
        pVal[,i,j] <- predParr
      }else if(predType=="parrAtTrap"){
        pVal[,i,j] <- predParr * pOut
      }else if(predType=="smoltAtTrap"){
        pVal[,i,j] <- predParr * (1-pOut) * wsSurvRes
      }else if(predType=="smoltAtDam"){
        pVal[,i,j] <- predParr * (1-pOut) * wsSurvRes +
          predParr * pOut * wsSurvMig
      }else if(predType=="winterSurvMig"){
        pVal[,i,j] <- wsSurvMig
      }else if(predType=="winterSurvRes"){
        pVal[,i,j] <- wsSurvRes
      }else if(predType=="winterSurvAgg"){
        pVal[,i,j] <- wsSurvMig * pOut + wsSurvRes * (1-pOut)
      }else{
        stop("ERROR: unrecognized name")
      }
    }
  }
  names(pVal) <- predType
  pVal
}
