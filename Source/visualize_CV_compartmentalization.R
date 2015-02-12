# Scatterplot to visualize compartmentalization of CV predictions
#png(file = "./Output/GAMM_predicted_vs_observed.png", width = 9, height = 6.5, units = "in", res = 600)
layout(matrix(1:2, nrow=2))
CVlist <- list(HF_sunsetCV, LF_sunsetCV)
labels <- c("(HF)", "(LF)")

par(mar = c(3.5, 4.5, 0.5, 0.5), mgp = c(2.25, 0.75, 0), cex = 0.9, cex.axis = 1.25, cex.lab = 1.25)

for (i in 1:length(CVlist)) {
  model <- CVlist[[i]]
  label <- LETTERS[i]
  with(model, plot(CVrate, passRate, xlab = paste("Predicted bat pass rate", labels[i]),
                   yaxt = "n", ylab="", bty = "l"))
  axis(2, las=1)
  mtext(paste("Observed bat pass rate", labels[i]), 2, line = 3, cex = 1.25)
  abline(v = quantile(model$CVrate, c(0.25, 0.5, 0.75)), col = "red", 
         lty = c("dotted", "dashed", "solid"))
  # Adding quartile lines for predictions
  abline(h = quantile(model$passRate, c(0.25, 0.5, 0.75)), col = "blue", 
         lty = c("dotted", "dashed", "solid"))
  text(x = max(model$CVrate), y = max(model$passRate),label = label, cex = 1.75)
}
#dev.off()