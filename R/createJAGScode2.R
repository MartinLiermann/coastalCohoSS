#' creates the JAGS code for the model
#'
#' @examples
#' bmod <- createJAGScode()
#' @export
createJAGScode2 <- function(){

  # Create a text string with the JAGS model specification
  mod1 <-   paste("
model
{
  ###################################
  ########### PROCESS MODEL #########
  ###################################

  for(i in 1:N){ # iterate over all years and populations

    smolt[i] ~ dlnorm(log(muSmolt[i]), SRresidTau[stock[i]])
    muSmolt[i] <- spawnersFem[i]/(1/prod[stock[i]]^3 + spawnersFem[i]^3/cap[stock[i]]^3)^(1/3) *
                 exp(yearEffect[year[i]])

    escapement[i] <- smolt[i] * oceanSurv[i]
    escapementFem[i] <- escapement[i] * sexRatio[i]

    logit(oceanSurv[i]) <- oceanSurvL[i] + yearEffectOS[year[i]]
    oceanSurvL[i] ~ dnorm(oceanSurvPopL[stock[i]], oceanSurvPopTau[stock[i]])

    logit(sexRatio[i]) <- sexRatioL[i]
    sexRatioL[i] ~ dnorm(sexRatioPopL[stock[i]], sexRatioPopTau[stock[i]])
  }

  # offset escapement to produce spawners
  # spawners are in calendar years and escapement are in brood years
  for(pop in 1:Npops){
    for(i in (escStartInd[pop]+3):escStopInd[pop]){
      spawnersWild[i] <- escapement[i-3]
      spawnersFemWild[i] <- spawnersWild[i] * sexRatio[i]
      spawnersFem[i] <- spawnersFemWild[i] / (1-pHOS[i])
    }
    # fill in the missing years with a vague prior
    for(i in escStartInd[pop]:(escStartInd[pop]+2)){
      spawnersWild[i]  ~ dlnorm(0,0.0001)
      spawnersFemWild[i] <- spawnersWild[i] * sexRatio[i]
      spawnersFem[i] <- spawnersFemWild[i] / (1-pHOS[i])
    }
  }

  ### population specific priors ###
  for(pop in 1:Npops){
    prod[pop] ~ dlnorm(logProdMu, logProdTau)
    cap[pop] ~ dlnorm(log(capSlope*habVar[pop]),logCapTau)

    oceanSurvPopL[pop] ~ dnorm(oceanSurvMu, oceanSurvTau)
    sexRatioPopL[pop] ~ dnorm(sexRatioMu, sexRatioTau)

    SRresidTau[pop] ~ dgamma(0.001,0.001)
    SRresidSD[pop] <- 1.0/sqrt(SRresidTau[pop])

    oceanSurvPopTau[pop] ~ dgamma(0.001,0.001)
    sexRatioPopTau[pop] ~ dgamma(0.001,0.001)
  }

  ### Hyper-priors (i.e. priors describing the distributions of parameters that vary by population)
  # productivity
  logProdMu ~ dnorm(prodMuPrior[1],prodMuPrior[2])
  logProdSD ~ dnorm(prodSDPrior[1],prodSDPrior[2]) # dt(0,1,1)T(0,) # half cauchy with var=tau=sd=1
  logProdTau <- 1.0/(logProdSD*logProdSD)

  # capacity
  capSlope ~ dlnorm(capSlopePrior[1],capSlopePrior[2])
  logCapSD ~ dt(0,1,1)T(0,)  # half cauchy with var=tau=sd=1 # dunif(0,10) #
  logCapTau <- 1.0/(logCapSD*logCapSD)

  # ocean survival
  oceanSurvMu ~ dnorm(oceanSurvMuPrior[1],oceanSurvMuPrior[2])
  oceanSurvSD ~ dt(0,1,1)T(0,)  # half cauchy with var=tau=sd=1
  oceanSurvTau <- pow(oceanSurvSD,-2)

  # sex ratio
  sexRatioMu ~ dnorm(sexRatioMuPrior[1],sexRatioMuPrior[2])
  sexRatioSD ~ dt(0,1,1)T(0,)  # half cauchy with var=tau=sd=1
  sexRatioTau <- pow(sexRatioSD,-2)

  ### year effects (SR resids), smolt residuals and ocean survival common to all populations
  for(k in 1:Nyears){
    yearEffectTmp[k] ~ dnorm(0,yearEffectTau)
    yearEffect[k] <- yearEffectTmp[k]-mean(yearEffectTmp)
    yearEffectTmpOS[k] ~ dnorm(0,yearEffectTauOS)
    yearEffectOS[k] <- yearEffectTmpOS[k]-mean(yearEffectTmpOS)
  }
  yearEffectTau ~ dgamma(0.001,0.001)
  yearEffectTauOS ~ dgamma(0.001,0.001)

  #######################################
  ########### OBSERVATION MODEL #########
  #######################################

  # smolt data
  for(i in 1:Nsmolt){ # smolt trap count (spring)
    smoltObs[i] ~ dlnorm(log(smolt[smoltInd[i]]),smoltObsTau)
  }
  smoltObsTau ~ dgamma(0.001,0.001)

  # escapement data (assume approximate CV of 0.15)
  for(i in 1:Nesc){
    escapementObs[i] ~ dlnorm(log(spawnersWild[escInd[i]]),1/(0.15^2))
  }

  # sex ratio (assume a logit normal sigma 0.38 based on binom(p=0.5, N=30))
  for(i in 1:NsexRatio){
    sexRatioObsL[i] ~ dnorm(logit(sexRatio[sexRatioInd[i]]),1/(0.38^2))
  }

 }
",sep="")

  mod1
}

