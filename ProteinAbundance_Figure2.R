# Ensure that the "ProteinAbundance_PartTwo_Analysis"
# has been run, and all required dataframes are properly
# loaded into R


# Plotting all the points into one scatterplot

x <- molcell[,1:22]

# Find median value for protein copy number for each ORF
x$Median <- apply(x[,2:22], 1, median, na.rm = TRUE)

# Order ORFs based on protein copy number
x <- x[order(x$Median),]

# Assign a value of 1 to nrow(DF) for each ORF based on protein copy number
x$ORD <- seq(1:nrow(x))

# Melt the data to have it in a format compatible for GGPLOT
xMELT <- melt(x, id = c("ORF", "ORD", "Median"))

# The result is a DF in long data format, where all the ORFs are laid out, and
# a value is a assigned to each ORF for protein molecules per cell if a data
# set has one reported. Even if an ORF has no value, it is included in the DF
# and represented as NA

kellys_palette <- c('grey', '#222222', '#F3C300', '#875692', '#F38400', 
                    '#A1CAF1',
                    '#BE0032', '#C2B280', '#848482', '#008856', '#E68FAC', 
                    '#0067A5', 
                    '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', 
                    '#882D17',
                    '#068481', '#654522', '#E25822', '#002366')


p.allORF <- ggplot(xMELT, aes(x = ORD, y = value))

dev.new(height = 13, width = 16)

p.allORF + geom_point(aes(colour = factor(variable), alpha = 0.5), size = 3) +
  scale_y_log10() +
  scale_colour_manual(values = kellys_palette) +
  theme(panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(colour = "grey", fill = NA, size = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())





# Plotting each independent dataset overlaying the "all_scatterplot".
# You need to specify which dataset you want to be coloured in blue,
# and this is specified after "variable" in the first line for the
# subset function

melt_subset <- subset(xMELT, variable == "PENG")
dev.new(height = 5, width = 5.7)
p.allORF + geom_point(aes(colour = "grey"), size = 1.5) +
           geom_point(data = melt_subset, colour = "#4070CD", size = 1.5) +
           scale_y_log10() +
           scale_colour_manual(values = kellys_palette) +
           guides(fill = FALSE) +
           theme(legend.position="none") +
           theme(panel.background = element_rect(fill = "#FFFFFF"),
                 panel.border = element_rect(colour = "grey", fill = NA, size = 2),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())
dev.print(png, file = "PENG.png", height = 500, width = 570)
dev.off()




# making the violin plot

# Figure 2
x <- molcell
x_bxplot <- x[,2:22]
x_bxplot[x_bxplot < 1] = NA

# violin plot
library(vioplot)
xMelt <- melt(x_bxplot)
p <- ggplot(xMelt, aes(x = variable, y = value))

p + geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75)) +
    scale_y_log10() +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())






# Making the comparison with TAP and low-scale studies
lowscale <- read.csv("smallscale.csv", na.strings = c("", "NA"), stringsAsFactors = F)

unified <- subset(molcell, select = c("ORF", "Median", "GHA"))

unif_merge <- merge(unified, lowscale, by.x = "ORF", by.y = "Systematic.Name",
                   all.x = T, all.y = T)
colnames(unif_merge) <- c("ORF", "UnifiedMedian", "MPC_Tap",
                          "StandardName", 
                          "Method", "MPC_smallscale", "MPC_cal", "XNAME")






# Plot against low-scale small-scale
plot(log(MPC_smallscale) ~ log(UnifiedMedian), data = unif_merge,
     pch = 16, cex = 1, col = "grey", xlim = c(0,14), ylim = c(0, 14))
abline(lm(log(MPC_smallscale) ~ log(UnifiedMedian), data = unif_merge),
       col = "firebrick", lwd = 2)


plot(log(MPC_Tap) ~ log(UnifiedMedian), data = unif_merge,
     pch = 16, cex = 1, col = "grey", xlim = c(0,14), ylim = c(0, 14))
abline(lm(log(MPC_Tap) ~ log(UnifiedMedian), data = unif_merge),
       col = "firebrick", lwd = 2)





# Now let's plot the log2 ratios between calibration and small scale / tap

unif_merge <- merge(unified, lowscale, by.x = "ORF", by.y = "Systematic.Name",
                   all.x = T, all.y = T)
colnames(unif_merge) <- c("ORF", "UnifiedMedian", "MPC_Tap",
                          "StandardName", 
                          "Method", "MPC_smallscale", "MPC_cal", "XNAME")
unif_merge$MPC_smallscale <- as.numeric(unif_merge$MPC_smallscale)
unif_merge$MPC_Tap <- as.numeric(unif_merge$MPC_Tap)

lowscale <- subset(unif_merge, select = c("ORF", "MPC_Tap",
                                          "MPC_smallscale"))

cal_set <- molcell[,1:6]
cal_set$cal_mean <- apply(molcell[,2:6], 1, mean, na.rm = T)

cal_merge <- merge(cal_set, lowscale, by.x = "ORF", by.y = "ORF",
                   all.x = T, all.y = T)
cal_merge$log2SS <- log2(cal_merge$cal_mean/cal_merge$MPC_smallscale)
cal_merge$log2TAP <- log2(cal_merge$cal_mean/cal_merge$MPC_Tap)

log2TAP_df <- subset(cal_merge, select = c("ORF", "log2TAP"))
log2TAP_df <- log2TAP_df[complete.cases(log2TAP_df$log2TAP),]
log2TAP_df <- log2TAP_df[order(log2TAP_df$log2TAP),]
log2TAP_df$ord <- seq(1:nrow(log2TAP_df))

log2SS <- subset(cal_merge, select = c("ORF", "log2SS"))
log2SS <- log2SS[complete.cases(log2SS$log2SS),]
log2SS <- log2SS[order(log2SS$log2SS),]
log2SS$ord <- seq(1:nrow(log2SS))

dev.new(width = 5, height = 5)
plot(log2TAP ~ ord, data = log2TAP_df, pch = 16, cex = 1, col = "grey")
abline(v = median(log2TAP_df$ord), lty = 3)

dev.new(width = 5, height = 5)
plot(log2SS ~ ord, data = log2SS, pch = 16, cex = 1, col = "grey")
abline(v = median(log2SS$ord), lty = 3)


