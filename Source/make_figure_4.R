theme_set(theme_bw(base_size = 18))
theme_update(panel.grid.minor = element_blank(),
             panel.grid.major= element_blank(),
             panel.border = element_blank(),
             panel.background= element_blank(),
             axis.line = element_line(color = 'black'),
             legend.position = 'none')

### FIGURE 4 ###

# Getting confusion matrices in table form
HFcm <- with(HF_sunsetCV, table(rateClass, CVrateClass))
LFcm <- with(LF_sunsetCV, table(rateClass, CVrateClass))

# Converting raw numbers to percentages
HFcmPerc <- HFcm
LFcmPerc <- LFcm
for (i in 1:4) {
  HFcmPerc[, i] <- HFcmPerc[, i] / colSums(HFcmPerc)[i]
  LFcmPerc[, i] <- LFcmPerc[, i] / colSums(LFcmPerc)[i]
}

plotCVdat <- data.frame(rbind(melt(HFcmPerc), melt(LFcmPerc)), 
                        grp = rep(c("HF", "LF"), each=16),
                        label = rep(letters[1:2], each=16))
plotCVdat <- within(plotCVdat, {
  obs <- factor(rateClass)
  CVrateClass <- factor(CVrateClass)
})

# Some trickery to improve the placement of grouped bars
# Obsolete; retained for posterity's sake
#plotCVdat$barpl <- c(rep(1, 4), rep(3.25, 4), rep(5.5, 4), rep(7.75, 4),
#                     rep(2, 4), rep(4.25, 4), rep(6.5, 4), rep(8.75, 4))

labels <- c("Low\n(\u2264 25%)", "Low/Med\n(26 - 50%)", "Med/High\n(51 - 75%)", 
                     "High\n(> 75%)")
p <- ggplot(plotCVdat, aes(x = CVrateClass, y = value, fill = obs)) +
  geom_bar(position = "dodge", stat="identity") +
  geom_bar(position = "dodge", stat="identity", color = "black", show_guide=FALSE) +
  scale_fill_manual("Observed bat activity", 
                    values = rev(c("#cccccc", "#969696", "#636363", "#252525")),
                    labels = c("Low", "Low/Med", "Med/High", "High")) + 
  scale_y_continuous("Proportion within activity class", limits = c(0, 0.7), expand = c(0,0), 
                     breaks = seq(0, 0.7, 0.1), labels = insert_minor(c(0, 0.7), 0.2, 0.1)) +
  scale_x_discrete("Predicted bat activity (percentiles)", labels = labels) + 
  #  scale_x_discrete("Predicted bat activity (percentiles)", labels = rep(c("HF", "LF"), 4),
  #                   limits = sort(unique(plotCVdat$barpl))) + 
  #  annotate("text", x = 1.5, y = -0.1, vjust = 0.75, label = labels[1]) +
  #  annotate("text", x = 3.75, y = -0.1, vjust = 0.75, label = labels[2]) +
  #  annotate("text", x = 6, y = -0.1, vjust = 0.75, label = labels[3]) +
  #  annotate("text", x = 8.25, y = -0.1, vjust = 0.75, label = labels[4]) +
  facet_grid(grp ~ .) +
  theme(axis.ticks.x = element_blank(), 
        axis.title.x = element_text(vjust = 0),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        legend.key = element_rect(color = "black"),
        legend.justification=c(0.5,1), legend.position=c(0.5, 1),
        legend.direction = "horizontal",
        panel.margin = unit(1, "lines")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))

# Label panels
xlab <- xrange(p, 1)[1] + diff(xrange(p, 1) * 0.01)
p <- p + geom_text(aes(x = xlab, y = 0.7, label = label),
                   hjust = 0, vjust=1, size = 9)
  
tiff(file = "./Output/figure4.tif", width = 5.75, height = 5.5, units = "in", res = 1000)
print(p)
dev.off()

