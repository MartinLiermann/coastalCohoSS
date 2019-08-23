# a function to to run the JAGS model
#' @export
runJAGSmodel <- function(bMod, dat, priors, calcInits, MCMCsims=10000, projection=FALSE){

  require(R2jags)

  saveList <- c("prod","cap","logProdMu","logProdSD","capSlope","logCapSD",
                "SRresidSD","yearEffect","yearEffectTau",
                "smolt","escapement","spawnersWild")

  if("smoltObs" %in% names(dat$jagsDat)) saveList <- c(saveList,"oceanSurv","oceanSurvPopL","yearEffectOS","yearEffectTauOS")
    
  # write the model to the local directory
  writeLines(bMod,con="bMod.txt") # writeLine seems to work best (cat, and write gave truncated results)

  mm <- MCMCsims/1000 # thin to 1000

  # see the generic run model for an example of how to run in parallel
  jagsDat <- c(dat$jagsDat,priors)

  m1 <- jags(data=jagsDat, inits=calcInits, parameters.to.save=saveList, model.file="bMod.txt",
                n.chains=3, n.iter=1100*mm, n.burnin=100*mm, n.thin=mm,
                DIC=TRUE, digits=5)

  list(JAGSout=m1, dat=dat, priors=priors, calcInits=calcInits, MCMCsims=MCMCsims, projection=projection)
}
