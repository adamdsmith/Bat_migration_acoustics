# Load necessary packages
toLoad = c("reshape", "lubridate", "plyr")
instant_pkgs(toLoad); rm(toLoad)

# Combine with weather data
years <- as.character(2010:2012)
allBats <- vector(mode="list", length=length(years))
names(allBats) <- years

for (i in 1:length(years)) {
  allBats[[i]] <- subset(bats, year(night) == years[i])
  allBats[[i]] <- as.data.frame(cast(melt(allBats[[i]], id = c("site", "night"), 
                                          measure = c("passes", "hif", "lof"), na.rm=T),
                                    site + night ~ variable, sum))
  
  # Combining weather data
  year <- years[i]
  allBats[[i]] <- arrange(join(allBats[[i]], allWeather[[year]], by = "night"), night, site)
}

#  Finalize bats and weather data set                      
batsWeather<- ldply(allBats)[, c(2:6, 10:17, 19:21, 23, 27, 30:39)]
batsWeather<- within(batsWeather, {
  doy <- yday(night)
  year <- factor(year(night))
  det_sets <- ifelse(year == 2010, 1, 0)
  siteyear <- as.factor(as.character(interaction(site, year)))
})
# Truncating to set 31 Oct as last recording night of each year
batsWeather <- arrange(subset(batsWeather, month(night) < 11), siteyear, doy)
save(batsWeather, file = "./Data/batsWeather.rda")
