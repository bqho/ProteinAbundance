# Ensure that the "ProteinAbundance_PartTwo_Analysis"
# has been run, and all required dataframes are properly
# loaded into R


corPair <- combined_raw[,2:ncol(combined_raw)]
corPair <- cor(log(corPair), use = "pairwise")
corPair[lower.tri(corPair, diag = TRUE)] <- NA
corPair_Melt <- melt(corPair, id = rownames(corPair))



# Visualize and obtain heatmap colours corresponding to correlation value
p <- ggplot(corPair_Melt, aes(x = Var1, y = Var2))

p + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "#005b96", high = "#EDED4C",
                      breaks = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
                                 0.8, 0.85, 0.9, 0.95, 1.0),
                      limits = c(0.5,1), guide = "legend") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# Also get the scatterplot matrix
scatterplotMatrix(log(combined_raw[,2:ncol(combined_raw)]),
                  smoother = FALSE,
                  transform = FALSE,
                  cex = 0.5,
                  pch = 16,
                  lwd = 1,
                  col = "DarkGrey",
                  legend.plot = TRUE)
