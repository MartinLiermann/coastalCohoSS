#' create default model priors
#'
#' @param modParams a list with any priors that you would like to override.
#' @return a list of priors.
#' @examples
#' priors <- createDefaultPriors()
#' @export
createDefaultPriorsOld <- function(modParams=NULL){
  tmpList <- list(
    # priors
    prodMuPrior = c(mu=100, tau = 0.001),
    prodSDPrior = c(mu=100, tau = 0.001),
    capSlopePrior = c(mu=1, tau = 0.001),
    oceanSurvMuPrior = c(mu=0, tau = 0.001),
    sexRatioMuPrior = c(mu=0, tau = 0.001)
  )
  if(!is.null(modParams)){
    for(param in names(modParams)){
      if(param %in% names(tmpList)){
        tmpList[param] <- modParams[param]
      }else{
        stop(paste("Error: ",param," in modParams is not a valid parameter.",sep=""))
      }
    }
  }
  tmpList
}
