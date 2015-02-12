# Load necessary packages
toLoad <- c("plyr", "ggplot2", "lubridate")
instant_pkgs(toLoad); rm(toLoad)

# How many weather samples are nightly averages based on? 
sample_summary <- ldply(allWeather)
sample_summary <- within(sample_summary, {
  year <- as.factor(year(night))
  dd <- decimal_date(night)-year(night)
})

#Testing VIF of nightly average weather variables
testNightVIF <- within(sample_summary, {
  .id <- night <- nightN <- dbZ <- Eta <- logEta <- tempresid <- year <- NULL
  nightWdir <- Z <- VCP <- precip150 <- setTemp <- setRh <- 
	setWsp <- setMb <- setWdir <- setdMb <- setdTemp <- setdRh <- setWp12 <- 
	setSkycov <- setVis <- pmode <- setd6Mb <- NULL
})
cat("Nightly average weather variables")
print(corvif(testNightVIF))
pairs(na.omit(sample_summary[,c(3:7, 9:12, 14:16, 35)]), lower.panel=panel.cor)

#Testing VIF of sunset weather variables
testSunsetVIF <- within(sample_summary, {
  .id <- night <- nightN <- dbZ <- Eta <- logEta <- tempresid <- year <- NULL
  # uncomment next 3 lines to evaluate only preSunset variables
  setdMb <- setWdir <- VCP <- nightTemp <- nightRh <- nightWsp <- nightMb <- 
	nightWp12 <- nightWdir <- nightdTemp <- nightdRh <- 
	nightdMb <- propPrecip <- skycov <- vis <- NULL
})
cat("Pre-sunset weather variables")
print(corvif(testSunsetVIF))
pairs(na.omit(sample_summary[,c(18, 22:26, 28, 30:35)]), lower.panel=panel.cor)