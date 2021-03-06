#' create a data list for the JAGS models
#'
#' @param dataDir the directory of the file together.csv.
#' @param includeSmolt if TRUE, smolt data will be included
#'
#' @return a list with the data necessary to run the JAGS model.
#' @examples
#' dat <- createJAGSprojectionData("data dir goes here")
#'
#' @export
createJAGSdata <- function(dataDir="",dataType="smolt",includeSmolt=TRUE,habType="area",femaleEsc=FALSE,includeFry=FALSE){
  #dataDir <- "Y:/sept2018_sept2019/oregonCoho/data"

  library(dplyr)
  library(tidyr)
  
  if(dataType=="OCN"){
    fdat <- read.csv(file=paste(dataDir,"fishDataOCN.csv",sep="/"),stringsAsFactors=FALSE)
    hdat <- read.csv(file=paste(dataDir,"habitatDataOCN.csv",sep="/"),stringsAsFactors=FALSE)
    hrdat <- read.csv(file=paste(dataDir,"harvestRatesOCN.csv",sep="/"),stringsAsFactors=FALSE)
    osdat <- read.csv(file=paste(dataDir,"oceanSurvivalOCN.csv",sep="/"),stringsAsFactors=FALSE)
    
    sites <- unique(fdat$Site)
    stockNum <- 1:length(sites)
    names(stockNum) <- sites
    stock <- stockNum[fdat$Site]
    yrRange <- range(fdat$broodYear)
    year <- fdat$broodYear-yrRange[1]+1
    
    escapementObs <- fdat$escapementObs
    escInd <- which(!is.na(escapementObs))
    
    # create escStartInd and escStopInd for summing age specific escapement
    # starts at year minYear+5 and ends at year maxYear+3 (since we only need age 4&5)
    escStartInd <- rep(-1,length(sites))
    escStopInd <- rep(-1,length(sites))
    for(i in 1:length(sites)){
      popNam <- sites[i]
      escStartInd[i] <- which(fdat$broodYear==min(fdat$broodYear[fdat$Site==popNam]) & fdat$Site==popNam)
      escStopInd[i] <- which(fdat$broodYear==max(fdat$broodYear[fdat$Site==popNam]) & fdat$Site==popNam)
    }
    
    # harvest rates
    HRdat <- fdat %>% left_join(hrdat,by="broodYear") 
    
    # the habitat variable
    habVar <- hdat[match(sites,hdat$Site),habType]
    
    # ocean survival (may just be a bunch of 1's)
    OSdat <- fdat %>% left_join(osdat,by="broodYear") 
    
    pHOS <- fdat$pHOS

    jdat <- list(
      stock = stock,
      year = year,
      escapementObs = escapementObs[escInd],
      escInd = escInd,
      pHOS = pHOS,
      HR = HRdat$HR,
      oceanSurv = OSdat$OS,
      habVar = habVar,
      escStartInd = escStartInd,
      escStopInd = escStopInd,
      N = length(fdat$Site),
      Npops = length(sites),
      Nyears = yrRange[2]-yrRange[1]+1,
      Nesc = length(escInd)
    )
  } else if(dataType=="smolt") {
    fdat <- read.csv(file=paste(dataDir,"fishData.csv",sep="/"),stringsAsFactors=FALSE)
    hdat <- read.csv(file=paste(dataDir,"habitatData.csv",sep="/"),stringsAsFactors=FALSE)
    hrdat <- read.csv(file=paste(dataDir,"harvestRatesOCN.csv",sep="/"),stringsAsFactors=FALSE)
    
    # remove leading years without female spawner estimates
    sites <- unique(fdat$Site)
    rmRows <- c()
    for(site in sites){
      yr <- min(fdat$BroodYear[fdat$Site==site])
      while(is.na(fdat$FemaleParentsWild[fdat$BroodYear==yr & fdat$Site==site])){
        rmRows <- c(rmRows,which(fdat$BroodYear==yr & fdat$Site==site))
        yr <- yr+1
      }
    }
    fdat <- fdat[-rmRows,]
    
    # add harvest data
    fdat <- fdat %>% left_join(hrdat,by=c("BroodYear"="broodYear"))

    # calculate the estimated proportion of hatchery origin fish 
    #  this is based on females because we did not have the male hatchery fish numbers
    pHOS <- fdat$FemaleParentsHatchery/(fdat$FemaleParentsHatchery + fdat$FemaleParentsWild)
    
    # calculate the total escapement and observed sex ratio
    # because escapementObs includes natural and hatchery origin adjust for 
    # hatchery fish
    escapementObs <- (fdat$FemaleParentsWild + fdat$MaleParentsWild) / (1-pHOS)
    sexRatioObs <- fdat$FemaleParentsWild / escapementObs

    # fill in missing sex ratios using within population average (males are missing for initial years)
    for(site in sites){
      sexRatioObs[is.na(sexRatioObs) & fdat$site==site] <- mean(sexRatioObs[fdat$site==site],na.rm=TRUE)
    }
    
    # use the filled in sex ratios to fill in total escapement 
    naInd <- is.na(escapementObs)
    escapementObs[naInd] <- (fdat$FemaleParentsWild[naInd] / sexRatioObs[naInd])/(1-pHOS[naInd]) 

    # if femaleESc = TRUE then just use females for escapement
    if(femaleEsc){
      escapementObs <- fdat$FemaleParentsWild
    }
    
    # smolt observations
    smoltObs <- fdat$Smolts
    fryObs <- fdat$FryOutmigrant
    
    # create indices for the observed data (in case there are missing values).
    #  Here there should not be missing values.
    escInd <- which(!is.na(escapementObs))
    sexRatioInd <- which(!is.na(sexRatioObs))
    smoltInd <- which(!is.na(smoltObs))
    fryInd <- which(!is.na(fryObs))
    
    # create year index for year effect
    yrRange <- range(fdat$BroodYear)
  
    # create stock and year indices
    stockNum <- 1:length(sites)
    names(stockNum) <- sites
    stock <- stockNum[fdat$Site]
    year <- fdat$BroodYear-yrRange[1]+1
  
    # the habitat variable
    habVar <- hdat[match(sites,hdat$Site),habType]
  
    # create escStartInd and escStopInd for summing age specific escapement
    # starts at year minYear+5 and ends at year maxYear+3 (since we only need age 4&5)
    escStartInd <- rep(-1,length(sites))
    escStopInd <- rep(-1,length(sites))
    for(i in 1:length(sites)){
      popNam <- sites[i]
      escStartInd[i] <- which(fdat$BroodYear==min(fdat$BroodYear[fdat$Site==popNam]) & fdat$Site==popNam)
      escStopInd[i] <- which(fdat$BroodYear==max(fdat$BroodYear[fdat$Site==popNam]) & fdat$Site==popNam)
    }
  
    jdat <- list(
      stock = stock,
      year = year,
      escapementObs = escapementObs[escInd],
      escInd = escInd,
      pHOS = pHOS,
      HR = fdat$HR,
      habVar = habVar,
      escStartInd = escStartInd,
      escStopInd = escStopInd,
      N = length(fdat$Site),
      Npops = length(sites),
      Nyears = yrRange[2]-yrRange[1]+1,
      Nesc = length(escInd)
    )
    if(includeSmolt){
      jdat$smoltObs <- smoltObs[smoltInd]
      jdat$smoltInd <- smoltInd
      jdat$Nsmolt <- length(smoltInd)
    }else{
      jdat$oceanSurv <- rep(1,length(pHOS))
    }
    if(includeFry){
      jdat$fryObs <- fryObs[fryInd]
      jdat$fryInd <- fryInd
      jdat$Nfry <- length(fryInd)
    }
  } else {
    stop("Error: dataType must be smolt or OCN")
  }
  # for now just fill in the missing pHOS and HR values.
  # we could also model this using a logit-normal hierarchical model
  # logit(pHOS[i]) ~ norm(pHOSmuPop[stock[i]],pHOStau[stock[i]])
  # then we could include an observation model
  jdat$pHOS <- fillGapsKern(jdat,"pHOS")
  jdat$HR <- fillGapsKern(jdat,"HR")

  list(jagsDat=jdat, fdat=fdat, hdat=hdat, siteNames=sites)
}
