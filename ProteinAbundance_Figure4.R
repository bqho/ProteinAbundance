# Ensure that the "ProteinAbundance_PartTwo_Analysis"
# has been run, and all required dataframes are properly
# loaded into R


# Making the violin plot for just the total proteome [median of all data]
x <- molcell

x_median <- subset(x, select = c("ORF", "Median"))
xmedian_melt <- melt(x_median)

dev.new(height = 8, width = 3)
p_median <- ggplot(xmedian_melt, aes(x = variable, y = value))

p_median + geom_violin(scale = "width", draw_quantiles = seq(1:9)/10) +
    scale_y_log10() +
    geom_hline(yintercept = 1000) +
    geom_hline(yintercept = 10000) +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())


# Side note, SAFE analysis was performed in CytoScape as outlined in
# the materials and methods, using the median abundance measurements