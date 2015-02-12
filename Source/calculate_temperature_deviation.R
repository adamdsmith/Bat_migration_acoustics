# Loading required packages
toLoad <- c("lubridate", "mFilter", "mgcv", "nlme", "plyr", "ggplot2", "reshape")
instant_pkgs(toLoad); rm(toLoad)

years <- c("2009", "2010", "2011", "2012")
# Using daily maximum temperature data from Providence, RI (KPVD) to generate a GAM model 
# from which nightly temperature residuals will be calculated.  

# We evaluated temperature residuals in two ways:
# (1) based on the deviation from the 30-year average for a given night (actually decimal  
#     date, to accommodate leap years); and
# (2) based on a "moving" generalized additive model (GAM) that terminates on the date
#     of interest; the GAM is based on the previous 365 days of temperature observations 
#     from the GSOD for the nearest airport with an ASOS station (KPVD in this case)

# Using a moving GAM allows the model to accommodate any time trend (i.e., warming or cooling
# trends) over the past year while keeping each night at the same location in the dataset 
# (i.e., the terminal datum).  This also permits the prediction of the coming night's 
# temperature departure from expectation over the previous year

# Loading daily maximum temperatures for 2008 - 2012 (2 Aug 2008 - 3 Dec 2012)
temps <- read.csv("./Data/dailytemps.csv", header=T)
temps <- within(temps, {
  night <- ymd(as.character(date))
  tempC <- round(5 * (tempf - 32) / 9, 1) # convert to Celsius
  doy <- as.POSIXlt(night)$yday + 1
  year <- year(night)
  dec_date <- round(ifelse(is.nan(decimal_date(night) - year), 0, 
                          decimal_date(night) - year), 3)
  tempf <- date <- NULL
})
temps$time <- 1:nrow(temps)

# Loading daily maximum temperatures for 1982 - 2012
tempsall <- read.csv("./Data/30_yr_hitemp_history.csv", header=T)
tempsall <- within(tempsall, {
  night <- ymd(as.character(date))
  tempC <- round(5 * (hitemp - 32) / 9, 1) # convert to Celsius
  doy <- as.POSIXlt(night)$yday + 1
  year <- year(night)
  dec_date <- round(ifelse(is.nan(decimal_date(night) - year), 0, 
                          decimal_date(night) - year), 3)
  hitemp <- date <- NULL
})

# First, looking at the data
qplot(dec_date, tempC, data=temps)
qplot(dec_date, tempC, data=tempsall)

# In GAM based on historic data, we ignore (potentially important)
# time trends to calculate a 30-year average based only on decimal date,
# which accomodates leap years

# Fit basic model to all historic data; no correlation yet
gam1a <- gamm(tempC ~ s(dec_date, bs="cc"), data=tempsall)
summary(gam1a$gam)

## look at autocorrelation in residuals:
acf(resid(gam1a$lme, type = "normalized"))
## ...wait... look at the plot, only then do...
pacf(resid(gam1a$lme, type = "normalized"))
## Sharp cut-off in PACF - AR(3) might be best for historic data

## Fitting AR1s
gam2a <- gamm(tempC ~ s(dec_date, bs="cc", k=20), data=tempsall,
             correlation = corAR1(form = ~1|year))

## Fitting AR2
gam3a <- gamm(tempC ~ s(dec_date, bs="cc", k=20), data=tempsall,
             correlation = corARMA(form = ~1|year, p = 2)) 

## Fitting AR3
gam4a <- gamm(tempC ~ s(dec_date, bs="cc", k=20), data=tempsall,
             correlation = corARMA(form = ~1|year, p = 3))

## Comparing models
## AR3 is preferred model, but is slower
anova(gam1a$lme, gam2a$lme, gam3a$lme, gam4a$lme)

## Diagnostic plots of AR3 model 
## Heavier tails (i.e., more large residuals) than expected,
## but otherwise okay
with(tempsall, tsDiagGamm(gam4a, timevar = doy, observed = tempC))

## Plotting AR3 model, with superimposed residuals
## Largest residuals seem to occur in Jan/Feb, so
## even better implications for our application
plot(gam4a$gam, residuals = TRUE, pch = 19, cex = 0.75)

## Getting unique decimal date values to be encountered
uniq_dd <- data.frame(dec_date = unique(tempsall$dec_date))

## Getting predicted temperatures from AR3 model
avg_temp <- predict.gam(gam4a$gam, newdata=uniq_dd)

## Linking to dates of interest, and calculating
## deviation from 30-year average (i.e., "normal") high temperature
avg_temps <- data.frame(uniq_dd, avg_temp)
temps <- join(temps, avg_temps, by = "dec_date")
temps$tempDev <- temps$tempC - temps$avg_temp

## Adding temperature residual calculated from GAM on last 365 days
## (see below for procedure by which these were calculated)
tempresids = dget(paste(setsave, "tempResids", sep=""))
for (i in 1:length(years)) {
  tempresids[[i]] = join(tempresids[[i]], temps[, c(5, 8)], by = "night")
}
save(tempresids, file = "./Data/tempResids.rda")

## Calculating daily "running residuals", which fits a GAM using the maximum daily temperatures
## from the previous year to calculate the residual temperature for a given day

## Note: changed spline type of seasonal trend (i.e., doy) since data so no inclusion of 
## spline for time trend is necessary (i.e., data are no longer confined to be cyclical)

## wrapper for function if needed for public use
## ifelse("mgcv" %in% rownames(installed.packages()), require(mgcv), install.packages("mgcv"))
## ifelse("lubridate" %in% rownames(installed.packages()), require(lubridate), install.packages("lubridate"))
pb <- winProgressBar(title = "progress bar", min = 0, max = nrow(temps), width = 300)
temps$tempresid <- 0
startyear <- 2009 #user controlled
for (i in 1:nrow(temps)) {
  if (year(temps[i,"night"]) < startyear) {temps$tempresid[i] = NA} else {
    if (month(temps[i,"night"]) < 8 | month(temps[i,"night"]) > 11) {temps$tempresid[i] = NA} else {
      #Find date 365 prior to current
      begin <- which(temps$night == (temps$night[i] - days(364))) 
      # Create new subset of data including only previous year
      current_temps <- temps[begin:i,] 
      #cat("Number of records used in GAM:", nrow(current_temps), "\n")
      current_gam <- gamm(tempC ~ s(doy, k=20), data=current_temps, correlation = corARMA(p = 2)) 
      temps$tempresid[i] <- tail(resid(current_gam$gam), 1)  
    }
  }
  setWinProgressBar(pb, i, title=paste(round(i/nrow(temps)*100, 0), "% done"))
}
Sys.sleep(5)
close(pb)

ggplot(subset(temps, year(night) > 2008), aes(doy, tempresid)) + geom_point() + 
  geom_line() + facet_grid(year~.) + xlim(210,340) 
            
## Creating list of temperature residuals for further use in analysis
tempResids <- vector(mode="list", length=length(years))
names(tempResids) <- years
for (i in 1:length(years)) {
  tempResids[[i]] <- within(subset(temps, year == years[i] & month(night) >= 8 & month(night) <= 11), {
    tempC <- time <- year <- doy <- NULL
    }) 
}  
     
save(tempResids, file = "./Data/tempResids.rda")
