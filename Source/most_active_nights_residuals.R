# Figure to examine how well most active nights were modelled

par(mfcol = c(2, 1), mar = c(4, 4.5, 0.5, 0.5), mgp = c(2.75, 1, 0), cex.axis = 1.75, cex.lab = 2)

gamList <- list(hinightgam, lownightgam)
for (i in 1:length(gamList)) {
  model <- gamList[[i]]
  resi <- resid(model$lme, type = "normalized")
  fits <- fitted(model$lme)
  obs <- model$gam$model[1][, 1]
  cutoff <- quantile(obs, probs = 0.95)
  temp <- data.frame(resi, fits, obs)
  with(temp, plot(resi ~ fits, ylab = "Normalized residual", bty = "l", 
                  pch = ifelse(obs < cutoff, 1, 19), xlab = ifelse(i == 2, "Fitted value (linear predictor)", ""), 
                  yaxt="n"))
  axis(2, las=1)
  coords = par("usr")
  xrange = coords[2] - coords[1]
  text(coords[1] + 0.95*xrange, coords[4], label = LETTERS[i], pos = 1, cex = 2.5) 
  legend( x="topleft", bty = "n", 
          legend="Most active 5% of nights",
          pch=19)
}
layout(1)
