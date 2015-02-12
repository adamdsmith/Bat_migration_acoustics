##### FIGURE 3 ######

# Variables to plot, in desired plotting order
plotVars <- c("doy", "setWp12", "tempDev", "setdTemp", "setd6Mb", "setVis")
# Full variable names of variables to be plotted, in order of appearance
varNames <- c("Day of year", "Wind profit (m/s)", 
              expression(paste("Temperature deviation ("^"o", "C)", sep="")),
              expression(paste(Delta, " temperature ("^"o", "C)", sep="")),
              expression(paste(Delta, " pressure (mb)")),
              "Visibility (mi)")

# List to hold ggplot objects
figure3 <- vector(mode = "list", length = length(plotVars)); names(figure3) <- plotVars

# GAMM models to use
gamList <- list(HF_sunsetGAMM, LF_sunsetGAMM)

# String of variables used in the models, excluding offset
modelVars <- c("tempDev", "setdTemp", "setWsp", "setMb", "setWp12", "setdRh", 
               "setd6Mb", "setVis", "doy")

for (var in plotVars) {
  
  plotDat <- data.frame()
  
  # Specify new data for prediction and intervals
  # It contains raw and scaled values of variable of interest over its range
  # in the data; all other variables will be held at their mean (i.e., scaled value of 0)
  
  # Have to accommodate interactions for a couple of variables
  interx <- var %in% "tempDev"
  if (interx) {
    doymean = mean(sunsetBats[, "doy"]); doysd = sd(sunsetBats[, "doy"])
    newDat <- expand.grid(var_scaled = unique(scaledSunsetBats[, var]), 
                          doy = (seq(230, 290, 30) - doymean) / doysd,
                          ldetectors = log(3),
                          var1 = 0, var2 = 0, var3 = 0, var4 = 0,
                          var5 = 0, var6 = 0, var7 = 0)
    newDat$var_raw <- rep(unique(sunsetBats[, var]), 3)
    names(newDat) <- c(var, "doy", "ldetectors", 
                       modelVars[-which(modelVars == var | modelVars == "doy")], 
                       paste0(var, "_raw"))
  } else {
    newDat <- expand.grid(var_scaled = unique(scaledSunsetBats[, var]), 
                          ldetectors = log(3),
                          var1 = 0, var2 = 0, var3 = 0, var4 = 0,
                          var5 = 0, var6 = 0, var7 = 0, var8 = 0)
    newDat$var_raw <- unique(sunsetBats[, var])
    names(newDat) <- c(var, "ldetectors", modelVars[-which(modelVars == var)], paste0(var, "_raw"))
  }
  
  # Cycle through each model for this variable
  for (g in 1:length(gamList)) {
    model <- gamList[[g]]$gam
    
    # Predicting over the range of current variable, all others at means
    out <- predict.gam(model, newdata = newDat, se.fit = TRUE)
    # Centering linear predictor relative to overall mean
    out$fit <- out$fit - coef(model)["(Intercept)"] - log(3)
    lcl <- with(out, fit - se.fit)
    ucl <- with(out, fit + se.fit)
    tmpDat <- data.frame(grp = ifelse(g == 1, "hiF", "lowF"),
                         var = newDat[, paste0(var, "_raw")],
                         doy = newDat[, "doy"], 
                         fit = out$fit, lcl = lcl, ucl = ucl)
    
    plotDat <- rbind(plotDat, tmpDat)
  }
  
  # Create plot
  
  # Get some parameters for labels
  ymin <- floor(ifelse(interx, min(plotDat$fit), min(plotDat$lcl)))
  ymax <- ceiling(ifelse(interx, max(plotDat$fit), max(plotDat$ucl)))
  yrange <- ymax - ymin
  xmin <- min(plotDat$var); xmax <- max(plotDat$var); xrange <- xmax - xmin
  
  if (interx) {
    
    p <- ggplot(plotDat, aes(x = var, y = fit, group = interaction(grp, as.factor(doy)))) + 
      #      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3) + 
      geom_line(aes(colour = interaction(grp, as.factor(doy)), 
                    linetype = grp), size = 1.25) +
      scale_colour_manual("", values = c(rep("gray70", 2), rep("gray35", 2), rep("black", 2))) +
      annotate("segment", x = xmin + 0.65*xrange, xend = xmin + 0.75*xrange, size = 1.25, 
               y = ymin + 0.225*yrange, yend = ymin + 0.225*yrange, color = "gray70") +
      annotate("segment", x = xmin + 0.65*xrange, xend = xmin + 0.75*xrange, size = 1.25,  
               y = ymin + 0.15*yrange, yend = ymin + 0.15*yrange, color = "gray35") +
      annotate("segment", x = xmin + 0.65*xrange, xend = xmin + 0.75*xrange,  size = 1.25, 
               y = ymin + 0.075*yrange, yend = ymin + 0.075*yrange, color = "black") +
      annotate("text", x = xmin + 0.775*xrange, y = ymin + 0.225*yrange, hjust=0, 
               label = "Early", size = 4) + 
      annotate("text", x = xmin + 0.775*xrange, y = ymin + 0.15*yrange,  hjust=0, 
               label = "Middle", size = 4) + 
      annotate("text", x = xmin + 0.775*xrange, y = ymin + 0.075*yrange,  hjust=0, 
               label = "Late", size = 4)   
    #facet_grid(grp ~ .)
    
  } else {
    p <- ggplot(plotDat, aes(x = var, y = fit, group = grp)) + 
      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.3) + 
      geom_line(aes(linetype = grp), size = 1.25)
  }
  
  # Tidy the plot
  p <- p + 
    geom_rug(sides = "b") + 
    scale_linetype_manual("", labels = c("High frequency", "Low frequency"),
                          values = c("solid", "dashed")) +
    xlab(varNames[which(plotVars == var)]) +
    theme(legend.position = "none", plot.margin = unit(rep(0.25, 4), "cm"))
  #  theme(legend.justification=c(0.5,1), legend.position=c(0.5, 1),
  #        legend.key = element_blank(), legend.direction = "horizontal",
  #        legend.key.width = unit(0.05, units = "npc"))
  
  # Label panels
  index <- which(plotVars == var)
  p <- p + annotate("text", x = max(plotDat$var), y = ymax,
                    label = LETTERS[index], hjust = 1, vjust=0.6,
                    size = 8)
  
  # Axis manipulation for multipanel plot
  p <-  p +
    if(index %in% c(1, 3, 5)){
      scale_y_continuous("Linear predictor", limits = c(ymin, ymax), 
                         breaks = seq(ymin, ymax, 1))
    } else {
      scale_y_continuous("",  limits = c(ymin, ymax), 
                         breaks = seq(ymin, ymax, 1))
    }
  
  if(var == "doy") {
    p <- p + scale_x_continuous(breaks=c(227,258,288),
                                labels=c("15 Aug", "15 Sep", "15 Oct"))
  }
  
  figure3[[var]] <- p
  
}


#png(file = "./Output/figure3.png", width = 6.5, height = 9, units = "in", res = 900)
multiplot(plotlist = figure3, layout = matrix(1:6, ncol=2, byrow=T))
#dev.off()
