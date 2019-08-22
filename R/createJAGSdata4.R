#' create a data list for the JAGS models
#'
#' @param dataDir the directory of the file together.csv.
#'
#' @return a list with the data necessary to run the JAGS model.
#' @examples
#' dat <- createJAGSprojectionData("data dir goes here")
#'
#' @export
createJAGSdata4 <- function(dataDir=""){
  #dataDir <- "Y:/sept2018_sept2019/oregonCoho/data"
  fdat <- read.csv(file=paste(dataDir,"fishData.csv",sep="/"),stringsAsFactors=FALSE)
  fdat2 <- read.csv(file=paste(dataDir,"allPops.csv",sep="/"),stringsAsFactors=FALSE)
  hdat <- read.csv(file=paste(dataDir,"habitatData.csv",sep="/"),stringsAsFactors=FALSE)
  aPops <- names(fdat2)[3:23]
  ndat <- fdat2 %>% gather(key="site",value="escapement",aPops)
  
  #### stopped here need to keep on going! ####
  
  
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
  
  # calculate the total escapement and observed sex ratio
  escapementObs <- fdat$FemaleParentsWild + fdat$MaleParentsWild
  sexRatioObs <- fdat$FemaleParentsWild / escapementObs
  
  # fill in missing sex ratios using within population average (males are missing for initial years)
  for(site in sites){
    sexRatioObs[is.na(sexRatioObs) & fdat$site==site] <- mean(sexRatioObs[fdat$site==site],na.rm=TRUE)
  }
  
  # use the filled in sex ratios to fill in total escapement
  naInd <- is.na(escapementObs)
  escapementObs[naInd] <- fdat$FemaleParentsWild[naInd] / sexRatioObs[naInd]
  
  # smolt observations
  smoltObs <- fdat$Smolts
  
  # create indices for the observed data (in case there are missing values).
  #  Here there should not be missing values.
  escInd <- which(!is.na(escapementObs))
  sexRatioInd <- which(!is.na(sexRatioObs))
  smoltInd <- which(!is.na(smoltObs))
  
  # calculate the estimated proportion of hatchery origin fish 
  #  this is based on females because we did not have the male hatchery fish numbers
  pHOS <- fdat$FemaleParentsHatchery/(fdat$FemaleParentsHatchery + fdat$FemaleParentsWild)

  # create year index for year effect
  yrRange <- range(fdat$BroodYear)

  # create stock and year indices
  stockNum <- 1:length(sites)
  names(stockNum) <- sites
  stock <- stockNum[fdat$Site]
  year <- fdat$BroodYear-yrRange[1]+1

  # the habitat variable
  habVar <- hdat$area[match(sites,hdat$Site)]

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
    sexRatioObsL = logit(sexRatioObs[sexRatioInd]),
    sexRatioInd = sexRatioInd,
    smoltObs = smoltObs[smoltInd],
    smoltInd = smoltInd,
    pHOS = pHOS,
    habVar = habVar,
    escStartInd = escStartInd,
    escStopInd = escStopInd,
    N = length(fdat$Site),
    Npops = length(sites),
    Nyears = yrRange[2]-yrRange[1]+1,
    Nesc = length(escInd),
    NsexRatio = length(sexRatioInd),
    Nsmolt = length(smoltInd)
  )
  # for now just fill in the missing pHOS values.
  # we could also model this using a logit-normal hierarchical model
  # logit(pHOS[i]) ~ norm(pHOSmuPop[stock[i]],pHOStau[stock[i]])
  # then we could include an observation model
  jdat$pHOS <- fillGapsKern(jdat,"pHOS")

  list(jagsDat=jdat, fdat=fdat, hdat=hdat, siteNames=sites)
}
