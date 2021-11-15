# Ensure that the "ProteinAbundance_PartTwo_Analysis"
# has been run, and all required dataframes are properly
# loaded into R


tagComp_sub <- molcell[,1:22]
tagComp_sub$logMS <- log(apply(tagComp_sub[,2:12], 1, median, na.rm = TRUE))
tagComp_sub$logGFP <- log(apply(tagComp_sub[,c(13:17,19:21)], 1, median, na.rm = TRUE))

tagComp_sub$MS <- apply(tagComp_sub[,2:12], 1, median, na.rm = TRUE)
tagComp_sub$GFP <- apply(tagComp_sub[,c(13:17,19:21)], 1, median, na.rm = TRUE)

# Make a separate data frame from the means
tagComp_sub <- subset(tagComp_sub, select = c("ORF", "logMS", "logGFP",
                                              "MS", "GFP"))



# Merge with the TAP-tag data, for comparison downstream of this analysis
tagvar <- merge(tagComp_sub, subset(molcell, 
                               select = c("ORF", "GHA")), by = "ORF")
tagvar$logTAP <- log(tagvar$GHA)
tagvar$logFC <- log10(tagvar$MS/tagvar$GFP)
colnames(tagvar) <- c("ORF", "logMS", "logGFP", "MS", "GFP", "TAP", "logTAP", "logFC")
tagvar <- subset(tagvar, logMS > 0)
tagvar <- subset(tagvar, GFP > 0)




# Obtain the linear regression model for log transformed MS and GFP mean values
tagvar_fit <- lm(logMS ~ logGFP, data = tagvar) # log-transformed




# you can get Cook Distance with simple function, then merge the statistics 
# with the tagvar dataframe, which is where the statistics was obtained. This
# is where we determine how different tagged proteins are with MS quantified
# proteins

cooksD <- data.frame(cooks.distance(tagvar_fit))
cooksD <- data.frame(rownames = as.numeric(rownames(cooksD)), cooksD)
StudRes <- data.frame(studres(tagvar_fit)) # getting the stud residuals as well
StudRes <- data.frame(rownames = as.numeric(rownames(StudRes)), StudRes)
tagvar <- data.frame(rownames = as.numeric(rownames(tagvar)), tagvar)

tagvar_stats <- Reduce(function(x,y) merge(x,y, all=TRUE),
                       list(tagvar, cooksD, StudRes))

colnames(tagvar_stats) <- c("rownames", "ORF", "logMS", "logGFP",
                            "MS", "GFP",
                            "TAP", "logTAP", "logFC",
                            "cookDistance", "studentRes")



tagvar_outliers <- subset(tagvar_stats, studentRes > 2 | studentRes < -2)
tagvar_outliers <- tagvar_outliers[!duplicated(tagvar_outliers$ORF),]

# genes that do not change
unch <- tagvar_stats[!(tagvar_stats$ORF %in% tagvar_outliers$ORF),]

remove(tagvar, StudRes, cooksD) # clean up environment





# Generate the final table

# Apply gene names to the systematic ORF names
ORF_df <- subset(SGD_ORFcomp, select = c("ORF", "Gene"))
colnames(ORF_df) <- c("ORF", "Gene")
newData <- merge(tagvar_outliers, ORF_df, by.x = "ORF", 
                 by.y = "ORF", all.x = T)

# Then we want to change all NA's (those without gene names) back to their
# systematic names
newData[is.na(newData$Gene),"Gene"] <- as.character(newData[is.na(newData$Gene),"ORF"])


# We are also going to want to combine this with the Yof tag stuff

Yof_sub <- subset(molcell, select = c("ORF", "YOF"))
newData2 <- merge(newData, Yof_sub, all.x = TRUE)

# write.csv(newData2, file = "tagvar_median20170729.csv")









# data visualization

# Plotting the bubble plot
dev.new(width = 5, height = 5)
plot(hatvalues(tagvar_fit), rstudent(tagvar_fit), type = "n",
     main = "Influence Plot for TAG Variation Data",
     xlab = "HAT values",
     ylab = "Studentized Residuals")
cook <- sqrt(cooks.distance(tagvar_fit))
points(hatvalues(tagvar_fit), rstudent(tagvar_fit), cex = 10*cook/max(cook))
abline(h = c(-2,0,2), lty = 2)

# Also plot MSmean vs GFPmean, but highlight those that are outliers
p <- ggplot(data = tagvar_stats, aes(x = GFP, y = MS))

dev.new(width = 5, height = 4.5)
p + geom_point(colour = "grey") +
    geom_point(data = tagvar_outliers, colour = "#FF1654") +
    scale_x_log10() +
    scale_y_log10() +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

# Show the log-fold changes for just the proteins that are identified
# as outliers

# First order the dataframe
tagvar_out_order <- tagvar_outliers[order(tagvar_outliers$MS),]
tagvar_out_order$ORD <- seq(1:nrow(tagvar_out_order))

dev.new(width = 8, height = 3)
p <- ggplot(data = tagvar_out_order, aes(x = ORD, y = logFC))
p + geom_bar(stat = "identity") +
    theme(panel.background = element_rect(fill = "#FFFFFF"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())





# We could also show the above scatterplot as a heatmap instead, which might
# help us cluster which proteins have reduced abundance when tagged
tagvar_outlier_heatmap <- subset(tagvar_outliers,
                                 select = c("logMS", "logGFP", "logTAP"))
rownames(tagvar_outlier_heatmap) <- tagvar_outliers$ORF

# Order the heatmap, so that we are looking at increasing mass spec abundance
tagvar_outlier_heatmap <- tagvar_outlier_heatmap[order(tagvar_outlier_heatmap[,1]),]
tagvar_outlier_heatmap <- data.matrix(tagvar_outlier_heatmap)


# creates a own color palette, black to yellow
my_palette <- colorRampPalette(c("#000000", "yellow"))(n = 20)

dev.new(width = 3, height = 10)
heatmap.2(tagvar_outlier_heatmap,
          col = my_palette,
          Colv = FALSE,
          Rowv = FALSE,
          dendrogram = "none",
          density.info = "none",
          trace = "none",
          na.rm = TRUE)


GFP_comp <- tagvar_stats[complete.cases(tagvar_stats$logGFP),]
TAP_comp_wtinGFP <- GFP_comp[complete.cases(GFP_comp$logTAP),]


# looking at the greatest change (fold-change) of the unaffected

affected <- tagvar_out_order$ORF

unaff <- molcell[!(molcell$ORF %in% affected),]


unaff$MS <- apply(unaff[,2:12], 1, median, na.rm = T)
unaff$GFP <- apply(unaff[,c(13:17,19:21)], 1, median, na.rm = T)


# and the following are the fold changes

FC <- unaff$MS/unaff$GFP