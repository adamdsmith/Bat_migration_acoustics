# Load nightly weather/radar data file
load("./Data/weatherVars.rda") # Created by nightly_weather_ASOS1.R
load("./Data/nightSky.rda") # Created by nighly_weather_ASOS5.R
load("./Data/tempResids.rda") # Created by temperature_GAMM.R
load("./Data/NEXRAD_reflectivity.rda") # Created by radar_assessment.R

# Load packages
toLoad = c("plyr")
instant_pkgs(toLoad); rm(toLoad)

# Set years of interest
years <- as.character(2009:2012)

# Creating the final combined data set
allWeather <- vector(mode="list", length=length(years))
names(allWeather) <- years
for (i in 1:length(years)) {
  allWeather[[i]] <- join(weatherVars[[i]], tempResids[[i]], by = "night")
  allWeather[[i]] <- join(allWeather[[i]], nightSky[[i]], by = "night")
  allWeather[[i]] <- join(allWeather[[i]], reflectivity[[i]], by = "night")
}

# Second, creating a dataset containing the single observations ~ 30 min prior to 
# sunset, which will then be averaged over the 4-5 stations to
# create a regional variable

# Get ASOS1 and ASOS5 data sets
load("./Data/ASOS1.rda") # Created by ASOS1_read.R and modified by nightly_weather_ASOS1.R
load("./Data/ASOS5.rda") # Created by nightly_weather_ASOS5.R
preSet <- list(ASOS=ASOS, ASOS5=ASOS5)

for (i in 1:length(years)) { 
  
  # Modifying ASOS-1
  preSet[["ASOS"]][[i]] <- subset(preSet[["ASOS"]][[i]], sinceset < -0.945 & sinceset > -0.955)
  
  # Calculating average of pertinent variables from up to 5 regional stations (all ASOS with 50 km) 
  # from observations recorded 30 minutes before sunset (~ 1 hour before civil sunset) for each night,
  # and calculating three difference variables (i.e., change in temperature, atmospheric pressure, 
  # and relative humidity in the past 24 hours)
  
  # First, create object containing circular mean direction to add subsequently...
  # Second, calculate means of linear variables
  # Third, join them
  tempwdir <- as.data.frame(cast(melt(subset(preSet[["ASOS"]][[i]], wdir > 0), id = "night", measure = "wdir", na.rm = T),
                                 night ~ variable, wdirfun))
  preSet[["ASOS"]][[i]] <- as.data.frame(cast(melt(preSet[["ASOS"]][[i]], id = "night",
                                                   measure = c("time", "sinceset", "wsp", "cTempC", "rh", "cWsp", "cMb"), na.rm=T), night ~ variable, mean))
  preSet[["ASOS"]][[i]] <- join(preSet[["ASOS"]][[i]], tempwdir)
  names(preSet[["ASOS"]][[i]]) <- c("night", "time", "sinceset", "uncent_wsp", "setTemp", "setRh", "setWsp", "setMb", "setWdir")
  
  # Fourth, calculate new variables from observations averaged across stations
  # Calculating differences variables (i.e., 24-hour change in temperature and  
  # atmospheric pressure, and relative humidity)
  for (j in 1:nrow(preSet[["ASOS"]][[i]])) {
    if (j == 1) {
      preSet[["ASOS"]][[i]]$setdTemp[j] <- NA
      preSet[["ASOS"]][[i]]$setdMb[j] <- NA
      preSet[["ASOS"]][[i]]$setdRh[j] <- NA
    } else {
      diff <- as.duration(preSet[["ASOS"]][[i]]$night[j] - preSet[["ASOS"]][[i]]$night[j-1])
      if (diff == duration(1, "days")) {
        preSet[["ASOS"]][[i]]$setdTemp[j] <- preSet[["ASOS"]][[i]]$setTemp[j] - preSet[["ASOS"]][[i]]$setTemp[j-1]
        preSet[["ASOS"]][[i]]$setdMb[j] <- preSet[["ASOS"]][[i]]$setMb[j] - preSet[["ASOS"]][[i]]$setMb[j-1]
        preSet[["ASOS"]][[i]]$setdRh[j] <- preSet[["ASOS"]][[i]]$setRh[j] - preSet[["ASOS"]][[i]]$setRh[j-1]
      } else {
        preSet[["ASOS"]][[i]]$setdTemp[j] <- NA
        preSet[["ASOS"]][[i]]$setdMb[j] <- NA
        preSet[["ASOS"]][[i]]$setdRh[j] <- NA
      }
    }
  }
  # set reference direction for wind profit calculation below
  ref_dir <- 135 
  
  # Wind direction (for posterity; not used in analysis) and wind profit
  preSet[["ASOS"]][[i]] <- within(preSet[["ASOS"]][[i]], {
    setWdir <- round(ifelse(setWdir < 0, 360 + setWdir, setWdir), 0)
    setWp12 <- round(ifelse(setWdir >= 180,
                            12 - sqrt(12^2 + uncent_wsp^2 - 2 * 12 * uncent_wsp * cos(rad(setWdir - 180 - ref_dir))),
                            12 - sqrt(12^2 + uncent_wsp^2 - 2 * 12 * uncent_wsp * cos(rad(setWdir + 180 - ref_dir)))), 2)
    time <- sinceset <- NULL
  })
  
  
  # Modifying ASOS-5
  preSet[["ASOS5"]][[i]] <- subset(preSet[["ASOS5"]][[i]], sinceset < -0.92)
  
  # Calculating average sky cover and visibility for pre-set observations
  preSet[["ASOS5"]][[i]] <- as.data.frame(cast(melt(preSet[["ASOS5"]][[i]], id = "night", measure = c("skycov", "vis"), na.rm=T),
                                               night ~ variable, mean))
  names(preSet[["ASOS5"]][[i]]) <- c("night", "setSkycov", "setVis")
}

# Third, creating a dataset containing single observations of atmospheric pressure
# ~ 6.5 hr prior to sunset, which will then be averaged over the 4-5 stations to
# create a regional variable indicating the 6-hour pressure tendency prior to sunset

set6h <- ASOS
for (i in 1:length(years)) { 
  
  # Modifying ASOS-1
  set6h[[i]] <- subset(set6h[[i]][, c("station", "night", "cMb", "sinceset")], sinceset < -6.945 & sinceset >= -6.955)
  
  set6h[[i]] <- as.data.frame(cast(melt(set6h[[i]], id = "night",
                                        measure = c("cMb"), na.rm=T), night ~ variable, mean))
  names(set6h[[i]]) <- c("night", "set6Mb") 
}

# Fourth, final join of all weather variables, and creation of final few variables
for (i in 1:length(years)) {
  allWeather[[i]] <- join(allWeather[[i]], preSet[["ASOS"]][[i]], by = "night")
  allWeather[[i]] <- join(allWeather[[i]], preSet[["ASOS5"]][[i]], by = "night")
  allWeather[[i]] <- join(allWeather[[i]], set6h[[i]], by = "night")
  allWeather[[i]] <- within(allWeather[[i]], {
    setd6Mb <- round(setMb - set6Mb, 2)
    set6Mb <- uncent_wsp <- NULL
  })
}

save(allWeather, file = "./Data/allWeather.rda")