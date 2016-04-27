##### FIGURE 3 ######
theme_set(theme_classic(base_size = 13))

# Variables to plot, in desired plotting order
plotVars <- c("doy", "nightWp12", "tempDev", "nightdTemp", "nightdMb")
# Full variable names of variables to be plotted, in order of appearance
varNames <- c("Day of year", "Wind profit (m/s)", 
              expression(paste("Temperature deviation ("^"o", "C)", sep="")),
              expression(paste(Delta, " temperature ("^"o", "C)", sep="")),
              expression(paste(Delta, " pressure (mb)")))

# List to hold ggplot objects
figure3 <- vector(mode = "list", length = length(plotVars)); names(figure3) <- plotVars

# GAMM models to use
gamList <- list(HF_nightGAMM, LF_nightGAMM)

# String of variables used in the models, excluding offset
modelVars <- c("tempDev", "nightdTemp", "nightWsp", "nightMb", "nightWp12", "nightdRh", 
               "nightdMb", "propPrecip", "vis", "doy")

for (var in plotVars) {
  
  # Specify new data for prediction and intervals
  # It contains raw and scaled values of variable of interest over its range
  # in the data; all other variables will be held at their mean (i.e., scaled value of 0)
  
  # Have to accommodate interactions for a couple of variables
  interx <- var %in% c("nightdTemp", "tempDev")
  if (interx) {
    doymean <- mean(nightBats[, "doy"]); doysd = sd(nightBats[, "doy"])
    newDat <- expand.grid(var_scaled = unique(scaledNightBats[, var]), 
                          doy = (seq(230, 290, 30) - doymean) / doysd,
                          ldetectors = log(3),
                          var1 = 0, var2 = 0, var3 = 0, var4 = 0,
                          var5 = 0, var6 = 0, var7 = 0, var8 = 0)
    newDat$var_raw <- rep(unique(nightBats[, var]), 3)
    names(newDat) <- c(var, "doy", "ldetectors", 
                       modelVars[-which(modelVars == var | modelVars == "doy")], 
                       paste0(var, "_raw"))
  } else {
    newDat <- expand.grid(var_scaled = unique(scaledNightBats[, var]), 
                          ldetectors = log(3),
                          var1 = 0, var2 = 0, var3 = 0, var4 = 0,
                          var5 = 0, var6 = 0, var7 = 0, var8 = 0,
                          var9 = 0)
    newDat$var_raw <- unique(nightBats[, var])
    names(newDat) <- c(var, "ldetectors", modelVars[-which(modelVars == var)], paste0(var, "_raw"))
  }
  
  # Cycle through each model for this variable
  plotDat <- data.frame()
  pResid <- data.frame()
  
  for (g in 1:length(gamList)) {
    model <- gamList[[g]]$gam
    
    # Predicting over the range of current variable, all others at means
    out <- predict.gam(model, newdata = newDat, se.fit = TRUE)
    # Centering Additive predictor relative to overall mean
    out$fit <- out$fit - coef(model)["(Intercept)"] - log(3)
    lcl <- with(out, fit - se.fit)
    ucl <- with(out, fit + se.fit)
    grp <- ifelse(g == 1, "hiF", "lowF")
    tmpDat <- data.frame(grp = grp,
                         var = newDat[, paste0(var, "_raw")],
                         doy = newDat[, "doy"], 
                         fit = out$fit, lcl = lcl, ucl = ucl)
    fvTerms <- predict(model, type = "terms")
    attr(fvTerms, "dimnames")[[2]] <- gsub("s\\(|\\)", "", attr(fvTerms, "dimnames")[[2]])
    tmpResid <- data.frame(grp = grp,
                           var = nightBats[, var],
                           fit = fvTerms[, var] + residuals(model, type = "working"))
    plotDat <- rbind(plotDat, tmpDat)
    pResid <- rbind(pResid, tmpResid)
  }
  
  # Create plot
  
  # Get some parameters for labels
  ymin <- floor(ifelse(interx, min(plotDat$fit), min(pResid$fit)))
  ymax <- ceiling(ifelse(interx, max(plotDat$fit), max(pResid$fit)))
  y_range <- ymax - ymin
  xmin <- min(plotDat$var); xmax <- max(plotDat$var); x_range <- xmax - xmin
  
  if (interx) {
    
    # Get rid of unimportant interactions
    if (var == "nightdTemp") {
      plotDat <- subset(plotDat, grp == "hiF" | grp == "lowF" & doy == (260 - doymean) / doysd)
    } else {
      plotDat <- subset(plotDat, grp == "lowF" | grp == "hiF" & doy == (260 - doymean) / doysd)
    }
    
    p <- ggplot(plotDat, aes(x = var, y = fit, group = interaction(grp, as.factor(doy)))) + 
      #      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3) + 
      geom_line(aes(colour = interaction(grp, as.factor(doy)), 
                    linetype = grp), size = 1.25) +
      scale_colour_manual("", values = c("gray70", rep("gray35", 2), "black")) +
      annotate("segment", x = xmin + 0.65*x_range, xend = xmin + 0.75*x_range, size = 1.25, 
               y = ymin + 0.225*y_range, yend = ymin + 0.225*y_range, color = "gray70") +
      annotate("segment", x = xmin + 0.65*x_range, xend = xmin + 0.75*x_range, size = 1.25,  
               y = ymin + 0.15*y_range, yend = ymin + 0.15*y_range, color = "gray35") +
      annotate("segment", x = xmin + 0.65*x_range, xend = xmin + 0.75*x_range,  size = 1.25, 
               y = ymin + 0.075*y_range, yend = ymin + 0.075*y_range, color = "black") +
      annotate("text", x = xmin + 0.775*x_range, y = ymin + 0.225*y_range, hjust=0, 
               label = "Early", size = 4) + 
      annotate("text", x = xmin + 0.775*x_range, y = ymin + 0.15*y_range,  hjust=0, 
               label = "Middle", size = 4) + 
      annotate("text", x = xmin + 0.775*x_range, y = ymin + 0.075*y_range,  hjust=0, 
               label = "Late", size = 4)   
    #facet_grid(grp ~ .)

  } else {
    p <- ggplot(plotDat, aes(x = var, y = fit, group = grp)) + 
        geom_ribbon(aes(ymin = lcl, ymax = ucl, group = grp), alpha = 0.3) + 
        geom_line(aes(linetype = grp), size = 1.25) +
        geom_point(data = pResid, aes(fill = grp), shape = 21, size = 0.75) + 
        scale_fill_manual(values = c("black", NA))
    }
  
  # Tidy the plot
  p <- p + 
    geom_rug(sides = "b") + 
    scale_linetype_manual("", labels = c("High frequency", "Low frequency"),
                          values = c("solid", "dashed")) +
    xlab(varNames[which(plotVars == var)]) +
    theme(legend.position = "none", plot.margin = unit(c(0.1, 0.22, 0.1, 0), "cm"))
  #  theme(legend.justification=c(0.5,1), legend.position=c(0.5, 1),
  #        legend.key = element_blank(), legend.direction = "horizontal",
  #        legend.key.width = unit(0.05, units = "npc"))
  
  # Label panels
  index <- which(plotVars == var)
  p <- p + annotate("text", x = min(plotDat$var), y = ymax,
                    label = letters[index], hjust = 0, vjust=0.6,
                    size = 6)
  
  # Axis manipulation for multipanel plot
  p <-  p +
    if(index %in% c(1, 4)){
      scale_y_continuous("Additive predictor", limits = c(ymin, ymax), 
                         breaks = seq(ymin, ymax, 1))
    } else {
      scale_y_continuous("",  limits = c(ymin, ymax), 
                         breaks = seq(ymin, ymax, 1))
    }
  
  if(var == "doy") {
    p <- p + scale_x_continuous(breaks=c(227,258,288),
                                labels=c("15 Aug", "15 Sep", "15 Oct"))
  }
  
  figure3[[var]] <- p + 
      # Workaround for current ggplot2 bug
      theme(axis.line.x = element_line(),
            axis.line.y = element_line())
  
}


#tiff(file = "./Output/figure3.tif", width = 8, height = 4.8, units = "in", 
#     compression = "lzw", res = 1000)
multiplot(plotlist = figure3, layout = matrix(1:6, ncol=3, byrow=T))
#dev.off()

#setEPS()
#postscript("./Output/figure3.eps", width = 10, height = 6)
#multiplot(plotlist = figure3, layout = matrix(1:6, ncol=3, byrow=T))
#dev.off()
