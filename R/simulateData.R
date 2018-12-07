#' Simulate data based on the fitted model
#'
#' @param runObj an object created by runJAGSmodel
#' @param years the number of years to project forward in time
#' @param numSims the number of forward simulations
#' @param obsDat should observation data be generated?
#' @return predictions
#' @examples
#' simulate(runObj,years=10,numSims=100)
#' @export
## -- NOT WORKING -- IN DEVELOPMENT --
## notes:
##  - can simulate based on multiple draws from the posterior or a single draw (or list that we construct).
##  - parameter should be how many times we simulate for each draw / parameter vect.
##  - create a separate sub function that generates data for testing the model and other posterior predictive checks.
##  - takes last state as an initial values (again, this is just another parameter in the parameter vector)
##  - need to decide whether to simulate from all levels of the hierarchical process.
##      probably not the population level, but yes for the year level.
##      otherwise we would not be including all of the year to year variability in the projections.
##      If there's correlation in the year to year variabiliy for the different parameters, we will need to
##      include it to make the projections realistic. For example, resident winter survival and pOut?
##      use a bivariate normal distribution for the prior? This means estimating another parameter (right?).
##      currently:
##        pOutPopL[pop] ~ dnorm(pOutMu, pOutTau)
##        winterSurvResPopL[pop] ~ dnorm(winterSurvResMu, winterSurvResTau)



simulateData <- function(runObj, years=25, numSims=10, obsDat=FALSE){
  x <- getPostDraws(runObj)
  bdat <- runObj$dat$jagsDat

  # create necessary vectors
  muParr <- numeric(bdat$N)
  parr <- numeric(bdat$N)
  migrantParr <- numeric(bdat$N)
  residentParr <- numeric(bdat$N)
  migrantSmoltAtDam <- numeric(bdat$N)
  residentSmoltAtDam <- numeric(bdat$N)
  smoltAtDam <- numeric(bdat$N)

  prod <- apply(x$prod,2,median)
  cap <- apply(x$cap,2,median)
  pOut <- apply(x$pOut,2,median)
  winterSurvMig <- apply(x$winterSurv,2,median)
  winterSurvRes <- apply(x$winterSurv,2,median)
  moveSurv <- apply(x$moveSurvMig,2,median)
  SRresidTau <- apply(x$SRresidTau,2,median)

  for(i in 1:bdat$N){
    ### process model ###
    # parr
    muParr[i] <- bdat$spawnersObs[i]/(1/prod[bdat$stock[i]] + spawnersObs[i]/cap[bdat$stock[i]])
    parr[i] <- rlnorm(log(muParr[i]), SRresidTau[stock[i]]^(-1/2))
    #exp(yearEffect[year[i]])
    # parr after fall migration
    migrantParr[i] <- parr[i]*pOut[i]
    residentParr[i] <- parr[i] - migrantParr[i]
    # smolt
    residentSmolt[i] <- residentParr[i]*winterSurv[i]
    # smolt at dam.
    #     moveSurvMig and moveSurvRes are the survivals incurred while moving from the winter
    #     rearing grounds to the dam for fish that overwinter below and above the trap respectively
    migrantSmoltAtDam[i] <- migrantParr[i] * winterSurv[i] * moveSurvMig[i]
    residentSmoltAtDam[i] <- residentParr[i] * winterSurv[i] * moveSurvRes[i]
    smoltAtDam[i] <- migrantSmoltAtDam[i] + residentSmoltAtDam[i]

    ### observation model ###
    # trap data
    parrOutObs[i] <- rlnorm(log(migrantParr[i]),bdat$parrOutObsTau[i]^(-1/2))
    smoltOutObs[i] <- rlnorm(log(residentSmolt[i]),bdat$smoltOutObsTau[i]^(-1/2))
    # trap to dam survival for fall (parr) and spring (smolt) outmigrants
    LparrSurvObs[i] <- rnorm(logit(migrantSmoltAtDam[i]/migrantParr[i]),bdat$LparrSurvObsTau[i]^(-1/2))
    LsmoltSurvObs[i] <- rnorm(logit(residentSmoltAtDam[i]/residentSmolt[i]),bdat$LsmoltSurvObsTau[i]^(-1/2))
    # summer parr (efishing) to dam survival (assumes no tagging to fall migration mortality)
    LsummerParrToSmoltSurvObs[i] <- rnorm(logit(smoltAtDam[i]/parr[i]),
                                         bdat$LsummerParrToSmoltSurvObsTau[i]^(-1/2))
    # winter parr (efishing) to dam survival (assumes no fall migration to tagging mortality)
    #LwinterParrToSmoltSurvObs[i] ~ dnorm(logit(residentSmoltAtDam[i]/residentParr[i]),
    #                                     LwinterParrToSmoltSurvObsTau[i])
  }
}
