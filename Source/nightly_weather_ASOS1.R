# Loading required packages
toLoad <- c("plotrix", "circular", "plyr", "lubridate", "reshape", "AED", 
            "PerformanceAnalytics", "ggplot2")
instant_pkgs(toLoad); rm(toLoad)

# Defining any necessary functions
# 1) Function for centering wind speed, temperature and pressure among stations
center_wsp = function(x) {
  x.obs <- x[!is.na(x)]
  return(x - mean(x.obs))
}

# 2) Function for calculating mean wind directions
mean_wdir <- function(wdir) {
  meandir <- mean.circular(as.circular(wdir, units = "degrees"))
  return(meandir)
}

# Loading data files into R
# ASOS 1-minute data (to be summarized below)
load("./Data/ASOS1.rda")

# Centering a few variables to account for systematic site differences
# Specifically, stations likely differ systematically in pressure,
# temperature, and wind speed...
tempASOS <- arrange(ldply(ASOS), station, night, sinceset)
# First, centering temperature
cTempC <- vector()
tempTemp <- with(tempASOS, tapply(tempC, station, center_wsp))
for (i in 1:length(tempTemp)) {
  cTempC <- append(cTempC, tempTemp[[i]])
}
# Second, centering station pressure
cMb <- vector()
tempMb <- with(tempASOS, tapply(mb, station, center_wsp))
for (i in 1:length(tempMb)) {
  cMb <- append(cMb, tempMb[[i]])
}
# Finally, centering wind speed
cWsp <- vector()
tempWsp <- with(tempASOS, tapply(wsp, station, center_wsp))
for (i in 1:length(tempWsp)) {
  cWsp <- append(cWsp, tempWsp[[i]])
}
# Consolidating
tempASOS <- data.frame(tempASOS, cTempC, cMb, cWsp)

# Returning to list form
years <- c("2009", "2010", "2011", "2012") # Change as necessary
ASOS <- vector(mode="list", length=length(years))
names(ASOS) <- years
for (i in 1:length(years)) {
  ASOS[[i]] <- arrange(subset(tempASOS, year(night) == years[i]), night, sinceset, station)
}
save(ASOS, file = "./Data/final_ASOS1.rda")

# Calculating divisor for proportion of night with precipitation variable.
# The reason a divisor is needed is twofold:
#   1) night length (and thus, # hours in a night) changes seasonally
#   2) data records can be incomplete in a given night (e.g., observations 
#      missing for one or more hours)
# Thus, the divisor is calculated on a nightly basis as the number of hours
# in which at least a single precipitation observation is recorded

divisors <- vector(mode="list", length=length(years))
names(divisors) <- years
for (i in 1:length(years)) {
  tempPrecip <- subset(ASOS[[i]], precip >= 0)
  divisors[[i]] <- as.data.frame(as.table(ifelse(table(tempPrecip$night, tempPrecip$hrset) > 0, 1, 0)))
  divisors[[i]] <- cast(melt(divisors[[i]], id = "Var1", measure = "Freq", na.rm=T), Var1 ~ variable, sum)
  names(divisors[[i]]) <- c("night", "divisor")
  divisors[[i]]$night <- ymd(divisors[[i]]$night)
}

# Calculating hourly sum of precipitation (across all stations) for each night
hrPrecipSum <- vector(mode="list", length=length(years))
names(hrPrecipSum) <- years
hrsWprecip <- vector(mode="list", length=length(years))
names(hrsWprecip) <- years
for (i in 1:length(years)) {
  hrPrecipSum[[i]] <- cast(melt(ASOS[[i]], id = c("night", "hrset"), measure = "precip", na.rm=T),
                          night ~ hrset ~ variable, sum)
  
  # Converting array (which is constructed as a contingency table) into a data frame for 
  # additional manipulation
  hrPrecipSum[[i]] <- arrange(as.data.frame(as.table(hrPrecipSum[[i]])), night, hrset)
  hrPrecipSum[[i]] <- within(hrPrecipSum[[i]], {
    precipSum <- Freq
    variable <- Freq <- NULL
    precipYN <- ifelse(precipSum > 0, 1, 0) # indicator of whether any precip was measured
  })
  hrsWprecip[[i]] <- cast(melt(hrPrecipSum[[i]], id = "night", measure = "precipYN", na.rm=T), night ~ variable, sum)
  hrsWprecip[[i]]$night <- ymd(hrsWprecip[[i]]$night)
  hrsWprecip[[i]] <- join(hrsWprecip[[i]], divisors[[i]], by = "night")
  hrsWprecip[[i]] <- within(hrsWprecip[[i]], {
    propPrecip <- round(precipYN / divisor, 2)
    precipYN <- divisor <- NULL
  })
}

### NIGHTLY AVERAGES

# Calculating nightly averages of pertinent variables (across all stations) 
# NOTE: recall that a couple of variables have been previously centered to make
# them more comparable among sites (i.e., temperature, wind speed, temperature)

nightASOS <- vector(mode="list", length=length(years))
names(nightASOS) <- years  
for (i in 1:length(years)) {
  nightASOS[[i]] <- subset(ASOS[[i]], hrset != -1)
  tempwdir <- as.data.frame(cast(melt(subset(nightASOS[[i]], wdir > 0), id = "night", measure = "wdir", na.rm = T),
                                night ~ variable, mean_wdir))
  
  nightASOS[[i]] <- as.data.frame(cast(melt(nightASOS[[i]], id = "night",
                                           measure = c("cTempC", "rh", "cWsp", "cMb", "wprof12"), na.rm=T),
                                      night ~ variable, mean))
  nightASOS[[i]] <- join(nightASOS[[i]], tempwdir)
  names(nightASOS[[i]]) <- c("night", "nightTemp", "nightRh", "nightWsp", "nightMb", "nightWp12", "nightWdir")
  
  # Calculating nightly average wind direction... 
  nightASOS[[i]] <- within(nightASOS[[i]], {
    nightWdir <- round(ifelse(nightWdir < 0, 360 + nightWdir, nightWdir), 0)
  })
  
  # Calculating differences variables (i.e., 24-hour change in temperature and  
  # atmospheric pressure, and relative humidity)
  for (j in 1:nrow(nightASOS[[i]])) {
    if (j == 1) {
      nightASOS[[i]]$nightdTemp[j] <- NA
      nightASOS[[i]]$nightdMb[j] <- NA
      nightASOS[[i]]$nightdRh[j] <- NA
    } else {
      diff <- as.duration(nightASOS[[i]]$night[j] - nightASOS[[i]]$night[j-1])
      if (diff == duration(1, "days")) {
        nightASOS[[i]]$nightdTemp[j] <- nightASOS[[i]]$nightTemp[j] - nightASOS[[i]]$nightTemp[j-1]
        nightASOS[[i]]$nightdMb[j] <- nightASOS[[i]]$nightMb[j] - nightASOS[[i]]$nightMb[j-1]
        nightASOS[[i]]$nightdRh[j] <- nightASOS[[i]]$nightRh[j] - nightASOS[[i]]$nightRh[j-1]
      } else {
        nightASOS[[i]]$nightdTemp[j] <- NA
        nightASOS[[i]]$nightdMb[j] <- NA
        nightASOS[[i]]$nightdRh[j] <- NA
      }
    }
  }
}

# Consolidating all variables into a single list
weatherVars <- vector(mode="list", length=length(years))
names(weatherVars) <- years
for (i in 1:length(years)) {
  nightN <- count(subset(ASOS[[i]], hrset != -1), "night"); names(nightN) <- c("night", "nightN")
  weatherVars[[i]] <- nightN #
  weatherVars[[i]] <- join(weatherVars[[i]], nightASOS[[i]], by = "night")
  weatherVars[[i]] <- join(weatherVars[[i]], hrsWprecip[[i]], by = "night")
}

save(weatherVars, file = "./Data/weatherVars.rda")
