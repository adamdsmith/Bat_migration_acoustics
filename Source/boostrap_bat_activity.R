# Load bat and weather data
load("./Data/batsWeather.rda")


# Estimate GLMM based on the full, site-specific data set, but omit a few variables that
# proved unimportant for simplicity 
fullGLMM <- glmmadmb(passes ~ ns(doy, df = 3) + site + year + tempDev + nightdTemp +  
                       nightWp12 + nightMb + nightWsp + vis, 
                     data = batsWeather,
                     family="nbinom")
# Use GAMM to estimate autoregressive structure for simulation
fullGAMM <- gamm(passes ~ s(doy) + site + year + tempDev + nightdTemp +  
                       nightWp12 + nightMb + nightWsp + vis, 
                     data = batsWeather,
                 family=negbin(fullGLMM$alpha),
                 correlation=corARMA(form = ~doy|siteyear, p = 1),
                 niterPQL = 100)

# This takes a while!  
# Nsim > 1000 to accomodate those iterations that do not converge
Nsim = 1070
simCoefs = data.frame()

# View GLMM Output
summary(fullGLMM)
# Coefficients:
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       6.76761    0.42482   15.93  < 2e-16 ***
# ns(doy, df = 3)1 -0.82650    0.21006   -3.93  8.3e-05 ***
# ns(doy, df = 3)2 -3.84810    0.51442   -7.48  7.4e-14 ***
# ns(doy, df = 3)3 -2.75278    0.20222  -13.61  < 2e-16 ***
# siteKURZ         -1.16511    0.20461   -5.69  1.2e-08 ***
# siteNINI         -0.37410    0.19156   -1.95   0.0508 .  
# siteSACH         -1.80250    0.16772  -10.75  < 2e-16 ***
# siteTRS2         -1.65501    0.17184   -9.63  < 2e-16 ***
# siteTRST          0.06620    0.16474    0.40   0.6878    
# siteWASH         -2.92119    0.20985  -13.92  < 2e-16 ***
# year2011         -2.38624    0.14978  -15.93  < 2e-16 ***
# year2012         -1.56912    0.18969   -8.27  < 2e-16 ***
# tempDev           0.05855    0.01576    3.72   0.0002 ***
# nightdTemp        0.03561    0.01811    1.97   0.0492 *  
# nightWp12         0.29093    0.03497    8.32  < 2e-16 ***
# nightMb           0.02844    0.00865    3.29   0.0010 ** 
# nightWsp         -0.10912    0.03849   -2.84   0.0046 ** 
# vis               0.08385    0.03365    2.49   0.0127 *  
  
set.seed(62)

for (loop in 1:Nsim) {
  
  cat("\nIteration:", loop, "\n")
  
  # Generating 14 samples of 91 days
  # 14 site/year combinations; I did not preserve clustering
  # i.e., I did not bootstrap a given 91 day sample only from the matching site/year data,
  # but rather from the whole data set.
  
  # Bootstrapping weather data 
  # Bootstrapped collectively to preserve relationships among variables
  weather <- data.frame()
  for (i in 1:14) {
    rows <- sample(1:775, 91)  ###CHECK THIS NUMBER###
    # This bootstrapping preserves the relationships among these four variables,
    # nightWp12, nightWsp, tempDev, and vis, particularly the atypical
    # relationship between nightWsp & nightWp12
    temporary <- batsWeather[rows,c("tempDev", "nightdTemp", "nightWp12", "nightMb", "nightWsp", "vis")]
    weather <- rbind(weather, temporary)
  }

  # Adding day of year (doy) and its natural splines (degree = 3)
  doy <- data.frame(rep(215:305, 14), ns(rep(215:305, 14), df = 3))
  names(doy) <- c("doy", "doy1", "doy2", "doy3")

  # Adding factor design variables
  
  # Site; KTLP as reference class
  sitestring <- as.factor(sort(rep(substr(unique(batsWeather$siteyear), 1, 4), 91)))
  site <- data.frame(site = sitestring, model.matrix(~ sitestring)[,2:7])
  names(site) <- c("site", as.character(sort(unique(batsWeather$site))[2:7]))

  #Year; 2010 as reference class
  yearstring <- vector()
    for (i in 1:14) {
      temp <- rep(substr(unique(batsWeather$siteyear), 6, 9)[i], 91)
      yearstring <- c(yearstring, temp)
    }
  suppressWarnings(year <- data.frame(year = yearstring, model.matrix(~ yearstring)[,2:3]))
  names(year) <- c("year", "year11", "year12")
  
  # Adding site and year combination factor
  siteyear <- interaction(site$site, year$year)
  
  # Adding AR(2) structure
  rho <- 0.3111548 # estimated from GAMM structured like fullGLMM with AR1 correlation structure
  ar1.sim <- as.vector(replicate(14, arima.sim(model=list(ar=rho), n=91)))

  # Consolidating
  xdata <- cbind(weather, doy, site, year, siteyear, ar1.sim)

  # Calculating linear predictor
  # Parameter estimates from comparable GLMM
  GLMMcoefs <- coef(fullGLMM)
  
  #Intercept 
  B0 <- coef(fullGLMM)["(Intercept)"]

  # Day of year natural splines
  B1 <- GLMMcoefs["ns(doy, df = 3)1"]
  B1a <- GLMMcoefs["ns(doy, df = 3)2"]
  B1b <- GLMMcoefs["ns(doy, df = 3)3"]

  # Site
  B2 <-GLMMcoefs["siteKURZ"]
  B3 <- GLMMcoefs["siteNINI"]
  B4 <- GLMMcoefs["siteSACH"]
  B5 <- GLMMcoefs["siteTRS2"]
  B6 <- GLMMcoefs["siteTRST"]
  B7 <- GLMMcoefs["siteWASH"]
  
  # Year
  B8 <- GLMMcoefs["year2011"]
  B9 <- GLMMcoefs["year2012"]

  # Weather variables
  B10 <- GLMMcoefs["tempDev"]
  B11 <- GLMMcoefs["nightdTemp"]
  B12 <- GLMMcoefs["nightWp12"]
  B13 <- GLMMcoefs["nightMb"]
  B14 <- GLMMcoefs["nightWsp"]
  B15 <- GLMMcoefs["vis"]

  # Linear predictor
  eta <- with(xdata, B0 + B1*doy1 + B1a*doy2 + B1b*doy3 + B2*KURZ + B3*NINI + B4*SACH + 
                B5*TRS2 + B6*TRST + B7*WASH + B8*year11 + B9*year12 + B10*tempDev +
                B11*nightdTemp + B12*nightWp12 + B13*nightMb + B14*nightWsp + B15*vis +
                ar1.sim)

  mu <- rnegbin(1274, mu=exp(eta), theta=fullGLMM$alpha)

  batsim <- cbind(passes = mu, xdata)

  #Estimate theta
  estTheta <- tryCatch.W.E(glmmadmb(passes ~ doy1 + doy2 + doy3 + tempDev +
                                       nightdTemp + nightWp12 + nightMb + nightWsp + 
                                       vis + (1|siteyear),
                       family = "nbinom", data = batsim))
  if(!is.null(estTheta$warning)) {next}
  theta <- estTheta$value$alpha
  
  #Constructing four versions of the model: (1) full simulated data set with known theta; (2) full
  #data set with estimated theta; (3) data set clipped to match missing pattern of actual data with 
  #known theta; and (4) clipped data set and estimated theta
  
  #(1) Full data set; known theta
  testFull <- tryCatch.W.E(gamm(passes ~ doy1 + doy2 + doy3 + tempDev + nightdTemp + 
                                  nightWp12 + nightMb + nightWsp + vis, 
                                random = list(siteyear = ~1), 
                                correlation = corARMA(form = ~doy|siteyear, p=1),
                                niterPQL = 100, family=negbin(fullGLMM$alpha), data=batsim))
  
  #(2) Full data set; estimated theta
  testFull2 <- tryCatch.W.E(gamm(passes ~ doy1 + doy2 + doy3 + tempDev + nightdTemp + 
                                   nightWp12 + nightMb + nightWsp + vis, 
                                 random = list(siteyear = ~1), 
                                 correlation = corARMA(form = ~doy|siteyear, p=1),
                                 niterPQL = 100, family=negbin(theta), data=batsim))
    
  # Now, clipping out the sites, years, and doys actually recorded
  clipBats <- batsim[which(with(batsim, interaction(site, year, doy)) %in% 
                        with(batsWeather, interaction(site, year, doy))),]

  #(3) Clipped data set; known theta
  testClip <- tryCatch.W.E(gamm(passes ~ doy1 + doy2 + doy3 + tempDev + nightdTemp + 
                                  nightWp12 + nightMb + nightWsp + vis, 
                                random = list(siteyear = ~1), 
                                correlation = corARMA(form = ~doy|siteyear, p=1),
                                niterPQL = 100, family=negbin(fullGLMM$alpha), data=clipBats))
    
  #(4) Clipped data set; estimated theta
  testClip2 <- tryCatch.W.E(gamm(passes ~ doy1 + doy2 + doy3 + tempDev + nightdTemp + 
                                   nightWp12 + nightMb + nightWsp + vis, 
                                 random = list(siteyear = ~1), 
                                 correlation = corARMA(form = ~doy|siteyear, p=1),
                                 niterPQL = 100, family=negbin(theta), data=clipBats))
  
  # Checking whether any of the four fitted gamms failed
  # If so, proceeding to next iteration
  # If not, storing coefficients
  if(sum(is.error(testFull$value), is.error(testFull2$value), 
         is.error(testClip$value), is.error(testClip2$value)) >= 1) {next}
  if(!is.null(c(testFull$warning, testFull2$warning, testClip$warning, testClip2$warning))) {next}
  simCoefs <- rbind(simCoefs, 
                   data.frame(simNum = loop, theta = fullGLMM$alpha, dataset = "full", 
                              t(coef(testFull$value$gam))),
                   data.frame(simNum = loop, theta = theta, dataset = "full", 
                              t(coef(testFull2$value$gam))), 
                   data.frame(simNum = loop, theta = fullGLMM$alpha, dataset = "clipped", 
                              t(coef(testClip$value$gam))),
                   data.frame(simNum = loop, theta = theta, dataset = "clipped", 
                              t(coef(testClip2$value$gam))))
  
}

# Keep only first 1000 successful bootstraps
simCoefs <- arrange(simCoefs[1:4000, ], dataset, theta)

# Melt simulation results for visualization
simSummary <- melt(subset(simCoefs, dataset == "clipped" & theta != fullGLMM$alpha)[,c(2,5:13)])

save(fullGLMM, GLMMcoefs, simCoefs, simSummary, file = "./Output/bat_bootstrap.rda")
