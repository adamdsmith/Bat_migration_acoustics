# Loading necessary packages
toLoad <- c("plyr", "glmmADMB", "splines", "mgcv", "dsm", "grid")
#toLoad <- c("glmmADMB", "mgcv", "plyr", "reshape", "RCurl", "nlme", "ggplot2", "party")
instant_pkgs(toLoad); rm(toLoad)

# Define the nights to use and filter bats and weather data
# Identify nights with at least three active detectors;  
# the number of active detectors will be used as an offset, so add it
# 2010 (n = 23), 2011 (n = 54), 2012 (n = 85)
useNights <- subset(ddply(batsWeather, .(night), "nrow"), nrow >= 3)
batsWeather <- batsWeather[batsWeather$night %in% useNights$night, ]
batsWeather <- join(batsWeather, useNights)

# Create bat data set with average nightly weather conditions
nightBats <- within(batsWeather, {ldetectors = log(nrow)})[, c(1:8, 10:16, 31:32, 34)]

# Sum passes by night; these counts will be used for nightly and pre-sunset analyses
nightlyPasses <- ddply(nightBats, .(night), summarize,
                       totalPasses = sum(passes),
                       hiPasses = sum(hif),
                       loPasses = sum(lof))

nightBats <- join(nightlyPasses, unique(nightBats[,  c(2, 6:18)]))

# Center and scale (1 SD) nightly weather; retain doy (as variable "time") for use in corARMA structure
scaledNightBats <- data.frame(time = nightBats$doy, nightBats[, c(1:4, 15, 17)], 
                              apply(nightBats[,c(5:14, 16)], 2, scale))
scaledNightBats <- arrange(scaledNightBats, year, time)

# Create bat data set with pre-sunset weather conditions
sunsetBats <- within(batsWeather, {ldetectors = log(nrow)})[, c(1:5, 14, 17:20, 22, 24:28, 31:32, 34)]
sunsetBats <- join(nightlyPasses, unique(sunsetBats[, c(2, 6:19)]))

# Center and scale (1 SD) pre-sunset weather; retain doy (as variable "time") for use in corARMA structure
scaledSunsetBats <- data.frame(time = sunsetBats$doy, sunsetBats[, c(1:4, 16, 18)], 
                               apply(sunsetBats[,c(5:6, 8:15, 17)], 2, scale),
                               precip150 = round(scale(sunsetBats[,7], scale=FALSE), 7))
scaledSunsetBats <- within(scaledSunsetBats, {
  precip150 <- as.factor(precip150)
})
scaledSunsetBats <- arrange(scaledSunsetBats, year, time)

################################################################################################
# Objective 1 - Associations between nightly atmospheric conditions and nightly bat activity
################################################################################################
# Hi frequency bats
# Get approximation of theta for use in GAMM
HF_nightGLMM <- glmmadmb(hiPasses ~ ns(doy, df = 3) + nightWsp + nightMb + 
                           nightWp12 + nightdTemp + nightdMb + nightdRh + 
                           propPrecip + tempDev + vis + tempDev:doy + 
                           nightdTemp:doy + offset(ldetectors) + (1|year), 
                         data = scaledNightBats,
                         family="nbinom")
HF_nightTheta <- HF_nightGLMM$alpha

# Fit the GAMM
HF_nightGAMM <- gamm(hiPasses ~ s(doy) + tempDev + nightdTemp + doy:tempDev + doy:nightdTemp + 
                       nightWsp + nightMb + nightWp12 + nightdRh + nightdMb + propPrecip + 
                       vis + offset(ldetectors),
                     random=list(year=~1),
                     data = scaledNightBats, 
                     family=negbin(HF_nightTheta),
                     correlation=corARMA(form = ~time|year, p = 1),
                     niterPQL = 100)

# Generating desired output
summary(HF_nightGAMM$gam)
arrange(data.frame(Variable = row.names(summary(HF_nightGAMM$gam)$p.table),
                   round(summary(HF_nightGAMM$gam)$p.table, 4)), -abs(Estimate))

# Evaluate GAMM fit
par(mfrow=c(2,1))
acf(resid(HF_nightGAMM$lme, type="normalized"))
pacf(resid(HF_nightGAMM$lme, type="normalized"))
par(mfrow=c(1,1))
rqgam.check(HF_nightGAMM$gam)

# GAMM - Low frequency bats
# Get approximation of theta
LF_nightGLMM <- glmmadmb(loPasses ~ ns(doy, df=3) + nightWsp + nightMb +
                           nightWp12 + nightdTemp + nightdMb + nightdRh + 
                           propPrecip + tempDev + vis + tempDev:doy + 
                           nightdTemp:doy + offset(ldetectors) + (1|year),
                         data = scaledNightBats,
                         family="nbinom")
LF_nightTheta <- LF_nightGLMM$alpha

# Fit the GAMM
# Smooth term for doy indicated linear relationship, so refit with GLMM
LF_nightGAMM <- gamm(loPasses ~ doy + tempDev + nightdTemp + doy:tempDev + 
                       doy:nightdTemp + nightWsp + nightMb + nightWp12 + 
                       nightdRh + nightdMb + propPrecip + vis + offset(ldetectors),
                     random=list(year=~1),
                     data = scaledNightBats, 
                     family=negbin(LF_nightTheta),
                     correlation=corARMA(form = ~time|year, p = 1),
                     niterPQL = 100)

# Generating desired output
summary(LF_nightGAMM$gam)
arrange(data.frame(Variable = row.names(summary(LF_nightGAMM$gam)$p.table),
                   round(summary(LF_nightGAMM$gam)$p.table, 4)), -abs(Estimate))

# Step 4; evaluate GAMM fit
par(mfrow=c(2,1))
acf(resid(LF_nightGAMM$lme, type="normalized"))
pacf(resid(LF_nightGAMM$lme, type="normalized"))
par(mfrow=c(1,1))
rqgam.check(LF_nightGAMM$gam)


################################################################################################
# Objective 2 - Predict nightly bat activity as a function of pre-set atmospheric conditions
################################################################################################

# GAMM - High frequency bat activity

# Actually reduces to a GLMM, estimated via the nlme package
# with glmmPQL (to retain autoregressive structure)
# We fit with gamm (mgcv) nonetheless, to better visualize the model fit

# Get approximation of theta
# This model has been updated to reflect the final variable evaluation (i.e., after variable reduction
# using penalized spline terms); it initially included all variables and interactions of interest, 
# including a natural cubic spline for day of year
HF_sunsetGLMM <- glmmadmb(hiPasses ~ doy + setWsp + setWp12 +
                            setd6Mb + tempDev + tempDev:doy +
                            offset(ldetectors) + (1|year),
                          data = scaledSunsetBats,
                          family="nbinom")

HF_sunsetTheta <- HF_sunsetGLMM$alpha

# Fit GAMM
# Since spline for doy and its interaction with tempDev indicated linear was best fit, 
# we refit with GLMM (although still using gamm in mgcv for continuity) to provide a
# comparison of variable importance
HF_sunsetGAMM <- gamm(hiPasses ~ doy + setWsp + setWp12 + 
                        setd6Mb + tempDev + tempDev:doy +
                        offset(ldetectors),
                      random = list(year = ~1), 
                      data = scaledSunsetBats, 
                      family=negbin(HF_sunsetTheta), 
                      correlation=corARMA(form = ~time|year, p=1),
                      niterPQL=100)
# Next line prints the table of smoothing parameters (for getting EDF values)
# It was used to evaluate whether a given term was penalized out of the model
# All other variables retained a linear form; precip150 was not evaluated in this manner
# summary(HF_sunsetGAMM$gam)$s.table

# Generating desired output
arrange(data.frame(Variable = row.names(summary(HF_sunsetGAMM$gam)$p.table),
                   round(summary(HF_sunsetGAMM$gam)$p.table, 4)), -abs(Estimate))

# Evaluate GAMM fit
par(mfrow=c(2,1))
acf(resid(HF_sunsetGAMM$lme, type="normalized"))
pacf(resid(HF_sunsetGAMM$lme, type="normalized"))
par(mfrow=c(1,1))
rqgam.check(HF_sunsetGAMM$gam)

# Leave-one-out cross validate predictions on a nightly basis
# 162 unique nights of data collection
# Leave out one night, refit regional mode, and then use it to predict activity for
# that night
HF_sunsetCV = data.frame()
nCV = nrow(scaledSunsetBats)
for (i in 1:nCV) {
  cat("\nProcessing record", i, "of", nCV)
  tmpDat = scaledSunsetBats[-i, ]
  newDat = scaledSunsetBats[i, ]
  tmpGAMM <- gamm(hiPasses ~ doy + setWsp + setWp12 + setd6Mb + tempDev + 
                    tempDev:doy + offset(ldetectors),
                  random = list(year = ~1),
                  data = tmpDat, 
                  family=negbin(HF_sunsetTheta), # New version from updated HF_sunsetGLMM above
                  correlation=corARMA(form = ~time|year, p=1),
                  niterPQL=100, verbosePQL=FALSE)
  HF_sunsetCV = rbind(HF_sunsetCV, 
                      data.frame(newDat, 
                                 CVfit = predict.gam(tmpGAMM$gam, newdata = newDat, type="response"),
                                 CV_lp = predict.gam(tmpGAMM$gam, newdata = newDat),
                                 Fullfit = predict.gam(HF_sunsetGAMM$gam, newdata = newDat, type="response")))
}

# Place observed passes and their associated CV predictions into quartiles for comparison
HF_sunsetCV = within(HF_sunsetCV, {
  passRate = hiPasses/exp(ldetectors)
  CVrate = CVfit/exp(ldetectors)
  rateClass = as.factor(ifelse(passRate <= quantile(passRate, 0.25), 1, 
                               ifelse(passRate <= quantile(passRate, 0.5), 2,
                                      ifelse(passRate <= quantile(passRate, 0.75), 3, 4))))
  CVrateClass = as.factor(ifelse(CVrate <= quantile(CVrate, 0.25), 1, 
                                 ifelse(CVrate <= quantile(CVrate, 0.5), 2,
                                        ifelse(CVrate <= quantile(CVrate, 0.75), 3, 4))))
})

# GAMM - Low frequency bat activity

# Get approximation of theta
# This model has been updated to reflect the final variable evaluation (i.e., after variable reduction
# using penalized spline terms); it initially included all variables and interactions of interest, 
# including a third-order polynomial for day of year
LF_sunsetGLMM <- glmmadmb(loPasses ~ doy + setdRh + tempDev + setdTemp +
                            setWp12 + tempDev:doy + setdTemp:doy +
                            setVis + offset(ldetectors), 
                          data = scaledSunsetBats, 
                          family="nbinom")
LF_sunsetTheta <- LF_sunsetGLMM$alpha

# Fit the GAMM (reduces to GLMM)
# Since spline for doy and its interaction with tempDev indicated linear was best fit, 
# we refit with GLMM (although still using gamm in mgcv for continuity) to provide a
# comparison of variable importance
LF_sunsetGAMM <- gamm(loPasses ~ doy + setdRh + tempDev + setdTemp +
                        setWp12 + tempDev:doy + setdTemp:doy +
                        setVis + offset(ldetectors),
                      random = list(year = ~1), 
                      data = scaledSunsetBats, 
                      family=negbin(LF_sunsetTheta), 
                      correlation=corARMA(form = ~time|year, p=1),
                      niterPQL=100)

# Next line prints the table of smoothing parameters (for getting EDF values)
# It was used to evaluate whether a given term was penalized out of the model
# All other variables retained a linear form; precip150 was not evaluated in this manner
# summary(LF_sunsetGAMM$gam)$s.table

# Generating desired output
arrange(data.frame(Variable = row.names(summary(LF_sunsetGAMM$gam)$p.table),
                   round(summary(LF_sunsetGAMM$gam)$p.table, 4)), -abs(Estimate))

# Evaluate GLMM fit
# Evaluate GAMM fit
par(mfrow=c(2,1))
acf(resid(LF_sunsetGAMM$lme, type="normalized"))
pacf(resid(LF_sunsetGAMM$lme, type="normalized"))
par(mfrow=c(1,1))
rqgam.check(LF_sunsetGAMM$gam)

# Leave-one-out cross validate predictions on a nightly basis
# 162 unique nights of data collection
# Leave out one night, refit regional mode, and then use it to predict activity for
# that night
LF_sunsetCV = data.frame()
nCV = nrow(scaledSunsetBats)
for (i in 1:nCV) {
  cat("\nProcessing record", i, "of", nCV)
  tmpDat = scaledSunsetBats[-i, ]
  newDat = scaledSunsetBats[i, ]
  tmpGAMM <- gamm(loPasses ~ doy + setdRh + tempDev + setdTemp +
                    setWp12 + tempDev:doy + setdTemp:doy +
                    setVis + offset(ldetectors),
                  random = list(year = ~1),
                  data = tmpDat, 
                  family=negbin(LF_sunsetTheta), # New version from updated LF_sunsetGLMM above
                  correlation=corARMA(form = ~time|year, p=1),
                  niterPQL=100, verbosePQL=FALSE)
  LF_sunsetCV = rbind(LF_sunsetCV, 
                      data.frame(newDat, 
                                 CVfit = predict.gam(tmpGAMM$gam, newdata = newDat, type="response"),
                                 CV_lp = predict.gam(tmpGAMM$gam, newdata = newDat),
                                 Fullfit = predict.gam(LF_sunsetGAMM$gam, newdata = newDat, type="response")))
}

# Place observed passes and their associated CV predictions into quartiles for comparison
LF_sunsetCV = within(LF_sunsetCV, {
  passRate = loPasses/exp(ldetectors)
  CVrate = CVfit/exp(ldetectors)
  rateClass = as.factor(ifelse(passRate <= quantile(passRate, 0.25), 1, 
                               ifelse(passRate <= quantile(passRate, 0.5), 2,
                                      ifelse(passRate <= quantile(passRate, 0.75), 3, 4))))
  CVrateClass = as.factor(ifelse(CVrate <= quantile(CVrate, 0.25), 1, 
                                 ifelse(CVrate <= quantile(CVrate, 0.5), 2,
                                        ifelse(CVrate <= quantile(CVrate, 0.75), 3, 4))))
})

# Save relevant data
save(HF_nightGLMM, HF_nightGAMM, HF_sunsetGLMM, HF_sunsetGAMM, HF_sunsetCV, 
     LF_nightGLMM, LF_nightGAMM, LF_sunsetGLMM, LF_sunsetGAMM, LF_sunsetCV, 
     nightBats, sunsetBats, scaledNightBats, scaledSunsetBats,
     file = "./Output/bat_activity_GAMMs.rda")