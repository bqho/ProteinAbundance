# Ensure that the "ProteinAbundance_PartTwo_Analysis"
# has been run, and all required dataframes are properly
# loaded into R


# Binning of the data does not require the molecules/cell data, but just the
# mean (or median) molecules per cell, and the coefficient of variation

prebin.df <- subset(molcell, select = c("Median", "Variance", 
                                        "CoefficientVariation"))
prebin.df <- subset(prebin.df, CoefficientVariation > 0)
prebin.df <- prebin.df[complete.cases(prebin.df),]

# Assign bin values to the data frame, which is ordered by increasing mean
num.bins <- 10 # enter the number of bins you want
v <- seq(from = 1, to = nrow(prebin.df), 
         by = (nrow(prebin.df)-1)/num.bins)
prebin.df$bin <- findInterval(1:nrow(prebin.df), v, all.inside = TRUE)


bin.data <- function(x, num.bins){
    postbin = 0
    bin <- seq(1:num.bins)
    
    cbind.fill <- function(...){
        nm <- list(...) 
        nm <- lapply(nm, as.matrix)
        n <- max(sapply(nm, nrow)) 
        do.call(cbind, lapply(nm, function (y) 
            rbind(y, matrix(, n-nrow(y), ncol(y))))) 
    }
    
    for(i in 1:length(bin)){
        temp.subset <- subset(x, bin == i)
        postbin = cbind.fill(postbin, temp.subset$CoefficientVariation)
    }
    
    return(postbin[,2:ncol(postbin)])
}

# This postbin.df is the final DF that puts all the coefficient of variation in
# one data frame
postbin.df <- data.frame(bin.data(prebin.df, num.bins))
colnames(postbin.df) <- c("bin1", "bin2", "bin3", "bin4", "bin5", 
                          "bin6", "bin7", "bin8", "bin9", "bin10")



#######------------------------------------------------------ Data Visualization
# We will want to plot a boxplot with the binned data as a beeswarm

dev.new(width = 5, height = 5)
beeswarm(postbin.df, method = c("swarm"), corral = c("wrap"), pch = 16,
         cex = 0.9, col = rgb(0.14, 0.14, 0.14, 0.4))
bxplot(postbin.df, add = TRUE, col = "#FF1654")

dev.off()

# Clean up global environment, remove variables

remove(mean.negrm, median.negrm, cv.negrm, var.negrm)
remove(prebin.df)
remove(final_ALL)
remove(num.bins, v)

finalDF <- merge(molcell, SGD_ORFcomp, all.x = T, all.y = T)
# write.csv(finalDF, file = "molcell_filter.csv")



