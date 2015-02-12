# Load bat and weather data
load("./Data/batsWeather.rda")

# Load raw bat pass data
load("./Data/bats.rda")

# Load nightly bat composition data
load("./Data/batPhenology.rda")

# Load necessary packages
toLoad <- c("ggplot2", "plyr", "reshape", "lubridate", "splines", "glmmADMB", "mgcv")
instant_pkgs(toLoad); rm(toLoad)

# Setting theme
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.minor = element_blank(),
             panel.grid.major= element_blank(),
             panel.border = element_blank(),
             panel.background= element_blank(),
             axis.line = element_line(color = 'black'),
             legend.position = 'top')


# Center and scale (1 SD) nightly weather; retain doy (as variable "time") for use in corARMA structure
scaledfullBats <- data.frame(time = batsWeather$doy, batsWeather[, c(1, 3, 15, 17, 29, 31)],
                             apply(batsWeather[,c(6:8, 10:14, 16, 32)], 2, scale),
                             det_sets = factor(batsWeather$det_sets))
scaledfullBats <- arrange(scaledfullBats, year, time)

# Estimate GAMM based on the full, site-specific data set to calculate the adjustment for 
# change in detection setting after 2010.  
fullGLMM <- glmmadmb(passes ~ ns(doy, df = 3) + site + tempDev + nightdTemp + doy:tempDev + 
                          doy:nightdTemp + nightWsp + nightMb + nightWp12 + nightdRh + 
                          nightdMb + propPrecip + vis + det_sets + (1|year), 
                         data = scaledfullBats,
                         family="nbinom")
fullTheta <- fullGLMM$alpha

fullGAMM <- gamm(passes ~ s(doy) + site + tempDev + nightdTemp + doy:tempDev + doy:nightdTemp + 
                   nightWsp + nightMb + nightWp12 + nightdRh + nightdMb + propPrecip + 
                   vis + det_sets,
                 random=list(year=~1),
                 data = scaledfullBats, 
                 family=negbin(fullTheta),
                 correlation=corARMA(form = ~time|siteyear, p = 1),
                 niterPQL = 100)

# Create data set with bat phenolgy data
# For a given day of the year, calculate total passes for a given species/classification,
# sum them across years, and standardize by the total detector nights, also summed
# over all years
nSites <- ddply(bats, .(night, site), summarize, passes = sum(passes))
activeNights <- ddply(nSites, .(night), "nrow")
activeNights <- within(activeNights, {
  doy <- yday(night)
})

batPhenology <- join(batPhenology, activeNights, by="night")
batPhenology <- within(batPhenology, {
  # Estimated adjustment to bat rates to account for more sensitive microphones
  # in 2010; based on nightly GAMM combined for all bat frequencies (allnightgam)
  # Adjustment was very similar for the two frequencies, so we used the combined 
  # estimate
  doy <- yday(night)
  detect_adj <- ifelse(year(night) == 2010, exp(2.187797), 1)
})

# Applying detector setting adjustment to put counts on approximately the same scale
# Specifically, divided 2010 species group counts by adjustment factor
batPhenology[,2:9] <- apply(batPhenology[,2:9], 2, function(x) round(x/batPhenology$detect_adj))

# Summing counts for a given species/classification for each day of the year,
# and dividing by the number of active detectors
adj_batPhenology <- as.data.frame(cast(melt(batPhenology, id = "doy", measure = names(batPhenology)[2:9], na.rm=T), 
                                  doy ~ variable, sum))
adj_batPhenology <- data.frame(adj_batPhenology, numdetects = ddply(batsWeather, .(yday(night)), "nrow")[, 2])
adj_batPhenology[, 2:9] <- apply(adj_batPhenology[, 2:9], 2, function(x) round(x / adj_batPhenology$numdetects))
adj_batPhenology <- melt(adj_batPhenology[, c("EPFU", "LABO", "LACI", "LANO", "MYSP", "PESU", "doy")], id = "doy")
adj_batPhenology$group <- ifelse(adj_batPhenology$variable %in% c("LABO", "PESU", "MYSP"), "HF", "LF")

### FIGURE A1 ###
phenPlots <- vector(mode = "list", length = 2); names(phenPlots) <- c("HF", "LF")
for (i in c("HF", "LF")) {
  index <- which(c("HF", "LF") == i)
  tmpDat <- subset(adj_batPhenology, group == i)
  
  phenPlots[[i]] <- ggplot(tmpDat, aes(x = doy, y = value, linetype = variable)) + 
    xlab("Day of year") + 
    ylab("# of passes / detector / night") + 
    geom_line() +
    annotate("text", x = max(tmpDat$doy), y = max(tmpDat$value), label = LETTERS[index], size = 10) + 
    scale_x_continuous(breaks=c(215, 243, 273, 304), labels=c("3 Aug", "31 Aug", "30 Sep", "31 Oct")) +
    scale_linetype("") +
    theme(legend.justification=c(0.5,1), legend.position=c(0.5, 1),
          legend.direction = "horizontal")
  
}
png(file="./Output/figureA1.png", width = 6.5, height = 8, units = "in", res = 600)
multiplot(plotlist = phenPlots, layout = matrix(1:2, ncol=1))
dev.off()
