#' create function that generates initial values for the Grand Ronde chinook model
#'
#' @param bdat the model data.
#' @return a function which generates a list of initial values.
#' @examples
#' initFunc <- createInitValsFunc(bdat)
#'
#' @export
#
####################################################
# WHAT I NEED TO DO NEXT
####################################################
createInitValsFunc <- function(dat, includePlots=FALSE,includeFry=FALSE){
  bdat <- dat$jagsDat
  sites <- dat$siteNames
  popNames <- tapply(names(bdat$stock),list(bdat$stock),unique)
  function(){
    logit <- function(x) log(x/(1-x))
    invLogit <- function(x) 1/(1+exp(-x))
    
    smDat <- "smoltObs" %in% names(bdat)

    ##### expand the observation data to include all combinations
    #####   of years and populations so it matches the parameters
    edat <- expandAllData(bdat)

    #####  spawners to trap model (prepare values) #####
    spawners <- edat$escapementObs
    Rtemp <- spawners/(1-edat$HR)/edat$oceanSurv
    smolt <- rep(NA,length(spawners))
    if(smDat){
      smolt <- edat$smoltObs
    } else {
      for(i in 1:length(sites)){
        smolt[bdat$stock==i] <- c(Rtemp[bdat$stock==i][-(1:3)],rep(NA,3))
      }
    }
    if(includeFry){
      lwinterSurv <- rep(-5,bdat$Npops)
    }
    habVar <- bdat$habVar

    ##### estimate the SR params #####
    srF <- function(p,ss) ss/(1/p[1]^3 + ss^3/p[2]^3)^(1/3)
    fitBH <- function(ss,rr){
      ssq <- function(pp,ss,rr){
        p <- exp(pp)
        sum((rr-srF(p,ss))^2)
      }
      bFit <- nlm(ssq,p=log(c(median(rr/ss),quantile(rr,prob=0.8))),ss,rr)
      bFit
    }
    habDat <- habVar[edat$stock]
    escDat <- spawners/(1-edat$pHOS)
    escDat <- pmax(1,escDat)
    srDat <- data.frame(S=escDat/habDat, R=smolt/habDat, stock=edat$stock)
    srDat <- srDat[!is.na(srDat$S) & !is.na(srDat$R),] # remove NAs
    srParams <- array(NA,dim=c(bdat$Npops,2))
    
    if(includePlots){
      windows(10,10)
      if(smDat){
        par(mfrow=c(3,2))
      }else {
        par(mfrow=c(5,5))
      }
        
      xx <- 0:round(max(srDat$S))
    }
    fryProd <- rep(NA,bdat$Npops)
    if(includeFry){
      fObs <- rep(NA,bdat$N)
      fObs[bdat$fryInd] <- bdat$fryObs
    }
    for(st in 1:bdat$Npops){
      sInd <- srDat$stock==st
      xLim <- c(0,max(srDat$S[sInd]))
      yLim <- c(0,max(srDat$R[sInd]))
      r1 <- fitBH(srDat$S[sInd],srDat$R[sInd])
      srParams[st,] <- exp(r1$estimate)
      if(includePlots){
        plot(srDat$S[sInd],srDat$R[sInd],pch=16,col="black",xlab="Spawners",ylab="Smolt",bty="l",main=popNames[st],xlim=xLim,ylim=yLim)
        lines(xx,srF(exp(r1$estimate),ss=xx),lwd=2,col=rgb(0.2,0.8,0.2,0.5))
        grid()
      }
      if(includeFry){
        fryProd[st] <- mean(fObs[bdat$stock==st]/bdat$escapementObs[bdat$stock==st],na.rm=TRUE)
        fry <- rep(NA,bdat$N)
        fry[bdat$fryInd] <- bdat$fryObs
      }
      
    }

    # now calculate the model specific parameters
    compTau <- function(x) 1/var(x,na.rm=T)
    prod <- srParams[,1]
    cap <- srParams[,2]
    capSlope <- median(cap/bdat$habVar)
    escObs <- edat$escapementObs
    esc <- rep(NA,length(escObs))
    for(i in 1:bdat$Npops){
      pInd <- which(edat$stock==i & edat$year %in% (1:3))
      esc[pInd] <- escObs[pInd]
    }

    pInd <- bdat$stock
    if(smDat) pMax <- 800 else pMax <- 20
    initList <- list(
      prod = pmin(prod,pMax),
      cap = cap,
      capSlope = capSlope,
      yearEffectTmp = rep(0,bdat$Nyears),
      yearEffectTmpOS = rep(0,bdat$Nyears)
    )
    if(includeFry) initList <- c(initList,
      list(
        lwinterSurv = rep(-4,bdat$Npops),
        fry = fry,
        fryProd = fryProd))
    
    initList
  }
}
