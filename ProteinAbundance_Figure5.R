# Ensure that the "ProteinAbundance_PartTwo_Analysis"
# has been run, and all required dataframes are properly
# loaded into R


# NOTES: All excel/CSV files containing the raw data from the original studies
# were placed in a single folder for ease of access. Each CSV corresponding
# to the study investigated is named in a specific manner such that copying all
# of PART 1 of this script would load all dataframes quickly for analysis. To
# have this feature, all CSVs must be in the same folder (remember to set the
# proper working directory), and name all CSVs exactly as how it is written
# below (case sensitive)



# Comparing with mRNA microarrays, and RNA-seq, and ribosome translation data
# First, let's load the data for microarray

Causton_mRNA <- read.csv("Causton_Conditions.csv", na.strings = c("", "NA"),
                         stringsAsFactors = FALSE)
Causton_Keep <- colnames(Causton_mRNA[,grep("Call", colnames(Causton_mRNA),
                                           invert = TRUE)])
Causton_mRNA <- subset(Causton_mRNA, select = Causton_Keep)

x <- subset(Causton_mRNA, 
            select=colnames(Causton_mRNA[,grep(".0.", colnames(Causton_mRNA),
                                               fixed = TRUE)]))
cor(x, use = "pairwise")

# It seems all the data correlate well with one another within the same study.
# For their mRNA abundance measurements, we will just use the Heat = 0 timepoint

Causton_mRNA <- read.csv("Causton_Conditions.csv", na.strings = c("", "NA"),
                         stringsAsFactors = FALSE)
Causton_mRNA <- subset(Causton_mRNA, select = c("ORF", "Heat.0...A."))
colnames(Causton_mRNA) <- c("ORF", "Causton2001")




# Loading Fritz Roth Data
Roth_microarray <- read.csv("FRoth_mRNAAbund_ver2.csv", na.strings = c("", "NA"),
							stringsAsFactors = FALSE)
Roth_microarray <- subset(Roth_microarray, select = c("ORF", "FY4a", 
													  "thresholdFY4"))

# get rid of mRNA below detection limit
mRNA_below <- subset(Roth_microarray, thresholdFY4 == 1)
mRNA_confident <- Roth_microarray[!Roth_microarray$ORF %in% mRNA_below$ORF,]

# And it is the mRNA_confident that we care about




# Loading Lipson Data
Lipson_microarray <- read.csv("lipson_2009.csv", na.strings = c("", "NA"),
							  stringsAsFactors = FALSE)
Lipson_microarray <- subset(Lipson_microarray, 
							select = c("ORF", "Avg", "MA", "Nagal."))

# All the microarray data that we will use is now loaded into the workspace
# and all we need to do now is merge the data

microarray_merge <- Reduce(function(x,y) merge(x, y, all = TRUE),
                     	   list(Causton_mRNA, mRNA_confident, 
                     	   	    Lipson_microarray))

colnames(microarray_merge) <- c("ORF", "Causton_MA", "Roth_MA", "RothThreshold",
								"Lipson_RNAseq", "Lipson_MA", "Nagal")

# The correlation between these data sets themselves aren't that great. I 
# wonder if taking the mean of all these data sets for our comparison is in 
# fact a good idea. We should probably normalize the data first, but I think
# it might be best to analyze our protein abundance data with each individual
# data set. 

microarray <- subset(microarray_merge, select = c("ORF", "Causton_MA",
												  "Roth_MA", "Lipson_MA"))

molcell_prot <- subset(molcell, select = c("ORF", "Median"))

molcell_merge <- merge(microarray, molcell_prot, all.x = TRUE, all.y = TRUE)



# Now we can plot and find correlation values

dev.new(width = 5, height = 5)
Causton_fit <- lm(log(Causton_MA) ~ log(Median), data = molcell_merge)
plot(log(Causton_MA) ~ log(Median), data = molcell_merge,
	 pch = 16, col = "grey")
abline(Causton_fit, lty = 6, col = "firebrick", lwd = 2)

cor(log(molcell_merge$Median), log(molcell_merge$Causton_MA), 
	use = "complete.obs")



dev.new(width = 5, height = 5)
Roth_fit <- lm(log(Roth_MA) ~ log(Median), data = molcell_merge)
plot(log(Roth_MA) ~ log(Median), data = molcell_merge,
	 pch = 16, col = "grey")
abline(Roth_fit, lty = 6, col = "firebrick", lwd = 2)

cor(log(molcell_merge$Median), log(molcell_merge$Roth_MA), 
	use = "complete.obs")



dev.new(width = 5, height = 5)
Lipson_fit <- lm(log(Lipson_MA) ~ log(Median), data = molcell_merge)
plot(log(Lipson_MA) ~ log(Median), data = molcell_merge,
	 pch = 16, col = "grey")
abline(Lipson_fit, lty = 6, col = "firebrick", lwd = 2)

cor(log(molcell_merge$Median), log(molcell_merge$Lipson_MA), 
	use = "complete.obs")













# Now repeat for RNA - seq data ___________________________________________________

Yassour_RS <- read.csv("Yassour_2009.csv", na.strings = c("", "NA", "Inf", 0.00),
					   stringsAsFactors = FALSE)
Yassour_RS <- subset(Yassour_RS, select = c("GENE", "Name", "Status", 
											"YPD0.1..rpkb.", "YPD0.2..rpkb.",
											"YPD15.1..rpkb."))
colnames(Yassour_RS) <- c("ORF", "Gene", "Status", "YPD0.1", "YPD0.2", 
						  "YPD15.1")

# They all correlate highly with one another, so just use any YPD RNA-seq data


Lipson_RS <- read.csv("lipson_2009.csv", na.strings = c("", "NA", "0.00"),
					  stringsAsFactors = FALSE)
Lipson_RS <- subset(Lipson_RS, select = c("ORF", "Avg", "Nagal."))
Lipson_RS <- Lipson_RS[Lipson_RS$Avg > 0,]

# This data frame actually includes two RNA-seq studies; Lipson and Nagal. et al

RNAseq <- merge(Yassour_RS, Lipson_RS, all.x = TRUE, all.y = TRUE)

molcell_merge <- merge(molcell_prot, RNAseq, all.x = TRUE, all.y = TRUE)


dev.new(width = 5, height = 5)
Yassour_fit <- lm(log(YPD0.1) ~ log(Median), data = molcell_merge)
plot(log(YPD0.1) ~ log(Median), data = molcell_merge,
	 pch = 16, col = "grey")
abline(Yassour_fit, lty = 6, col = "firebrick", lwd = 2)

cor(log(molcell_merge$Median), log(molcell_merge$YPD0.1), 
	use = "complete.obs")


dev.new(width = 5, height = 5)
LipsonRS_fit <- lm(log(Avg) ~ log(Median), data = molcell_merge)
plot(log(Avg) ~ log(Median), data = molcell_merge,
	 pch = 16, col = "grey")
abline(LipsonRS_fit, lty = 6, col = "firebrick", lwd = 2)

cor(log(molcell_merge$Median), log(molcell_merge$Avg), 
	use = "complete.obs")


dev.new(width = 5, height = 5)
Nagal_fit <- lm(log(Nagal.) ~ log(Median), data = molcell_merge)
plot(log(Nagal.) ~ log(Median), data = molcell_merge,
	 pch = 16, col = "grey")
abline(Nagal_fit, lty = 6, col = "firebrick", lwd = 2)

cor(log(molcell_merge$Median), log(molcell_merge$Nagal.), 
	use = "complete.obs")













# Comparisons with ribosomal profiling_________________________________________

# Ensure working directory is set to folder containing all the ribosomal
# profiling datasets we need

#####
Albert <- read.csv("Albert_GSM1335348.csv", na.strings = c("", "NA"),
				   stringsAsFactors = FALSE)
colnames(Albert) <- c("ORF", "Albert_RPKM")
Albert$Albert_RPKM <- as.numeric(Albert$Albert_RPKM)


#####
Brar <- read.csv("Brar_GSE34082.csv", na.strings = c("", "NA"),
				 stringsAsFactors = FALSE)
Brar <- subset(Brar, select = c("description",
								"gb15.vegetative.exponential.rib.fp.RPKM"))
colnames(Brar) <- c("ORF", "Brar_RPKM")

Brar$ORF <- lapply(Brar$ORF, function(x) strsplit(x," ")[[1]][1])
Brar$ORF <- as.character(Brar$ORF)


#####
Ingolia1 <- read.csv("Ingolia_RichQuant1.csv", na.strings = c("", "NA"),
					 stringsAsFactors = FALSE)
Ingolia1 <- subset(Ingolia1, select = c("yorf", "norm"))
colnames(Ingolia1) <- c("ORF", "Ingolia1")

Ingolia2 <- read.csv("Ingolia_RichQuant2.csv", na.strings = c("", "NA"),
					 stringsAsFactors = FALSE)
Ingolia2 <- subset(Ingolia2, select = c("yorf", "norm"))
colnames(Ingolia2) <- c("ORF", "Ingolia2")

Ingolia <- merge(Ingolia1, Ingolia2, all = T)
Ingolia$Ingolia_RPKM <- apply(Ingolia[,2:3], 1, mean, na.rm = T)
Ingolia <- subset(Ingolia, select = c("ORF" ,"Ingolia_RPKM"))


#####
Pop <- read.csv("Pop_GSE63789.csv", na.strings = c("", "NA"),
				stringsAsFactors = FALSE)
Pop <- subset(Pop, select = c("X.Name", "Sum.FP"))
colnames(Pop) <- c("ORF", "Pop_RPKM")
Pop$Pop_RPKM <- as.numeric(Pop$Pop_RPKM)


#####
Weinberg <- read.csv("Weinberg_RiboZeroRPKM.csv", na.strings = c("", "NA"),
					 stringsAsFactors = FALSE)
colnames(Weinberg) <- c("ORF", "Weinberg_RPKM")



# Remove unnecessary dataframes
remove(Ingolia1, Ingolia2)


# Merge the various ribosomal profiling data into one dataframe
ribosome_footprint <- Reduce(function(x,y) merge(x, y, all = TRUE),
                    		 list(Albert, Brar, Ingolia, Pop, Weinberg))


# Now merge with the protein abundance dataset
molcell_sub <- subset(molcell, select = c("ORF", "Median"))
molcell_fp <- merge(molcell_sub, ribosome_footprint, all.x = TRUE)
ln_molcell_fp <- log(molcell_fp[,2:ncol(molcell_fp)])

ln_molcell_fp[mapply(is.infinite, ln_molcell_fp)] <- NA

cor(ln_molcell_fp, use = "pairwise.complete.obs")


# Data visualization
dev.new(width = 5, height = 5)
Ingolia_fit <- lm(Ingolia_RPKM ~ Median, data = ln_molcell_fp)
plot(Ingolia_RPKM ~ Median, data = ln_molcell_fp, pch = 16, col = "grey")
abline(Ingolia_fit, lty = 6, col = "firebrick", lwd = 2)


dev.new(width = 5, height = 5)
Albert_fit <- lm(Albert_RPKM ~ Median, data = ln_molcell_fp)
plot(Albert_RPKM ~ Median, data = ln_molcell_fp, pch = 16, col = "grey")
abline(Albert_fit, lty = 6, col = "firebrick", lwd = 2)


dev.new(width = 5, height = 5)
Brar_fit <- lm(Brar_RPKM ~ Median, data = ln_molcell_fp)
plot(Brar_RPKM ~ Median, data = ln_molcell_fp, pch = 16, col = "grey")
abline(Brar_fit, lty = 6, col = "firebrick", lwd = 2)


dev.new(width = 5, height = 5)
Pop_fit <- lm(Pop_RPKM ~ Median, data = ln_molcell_fp)
plot(Pop_RPKM ~ Median, data = ln_molcell_fp, pch = 16, col = "grey")
abline(Pop_fit, lty = 6, col = "firebrick", lwd = 2)


dev.new(width = 5, height = 5)
Weinberg_fit <- lm(Weinberg_RPKM ~ Median, data = ln_molcell_fp)
plot(Weinberg_RPKM ~ Median, data = ln_molcell_fp, pch = 16, col = "grey")
abline(Weinberg_fit, lty = 6, col = "firebrick", lwd = 2)











# 3D PLOTTING_________________________________________________________________

# First, get the molcell_merge dataframe from RNA comparisons scripts.
# Also retrieve the molcell_fp dataframe from ribosome profiling script

molcell_merge_rna <- subset(molcell_merge, 
							select = c("ORF", "Median", "YPD0.1",
									   "Avg", "Nagal."))

colnames(molcell_merge_rna) <- c("ORF", "Median", "Yassour", "Lipson",
								 "Nagal")

molcell_fp_v2 <- subset(molcell_fp, select = c("ORF", "Albert_RPKM",
											   "Brar_RPKM", "Ingolia_RPKM",
											   "Pop_RPKM", "Weinberg_RPKM"))

molcell_RNA_FP <- merge(molcell_fp_v2, molcell_merge_rna, all.x = TRUE)


# Mode-shift normalize the RNA data sets
get.mode <- function(x, y){
  x.mode = NULL
  
  hist.info <- hist(log10(x[,y]), breaks = 50)
  hist.mids <- hist.info$mids
  x.mode <- hist.mids[which(hist.info$counts == max(hist.info$counts))]
  
  return(10^x.mode)
}

Albert.mode <- get.mode(molcell_RNA_FP, 2)[1]
Brar.mode <- get.mode(molcell_RNA_FP, 3)
Ingolia.mode <- get.mode(molcell_RNA_FP, 4)
Pop.mode <- get.mode(molcell_RNA_FP, 5)
Weinberg.mode <- get.mode(molcell_RNA_FP, 6)
Median.mode <- get.mode(molcell_RNA_FP, 7)
Yassour.mode <- get.mode(molcell_RNA_FP, 8)
Lipson.mode <- get.mode(molcell_RNA_FP, 9)
Nagal.mode <- get.mode(molcell_RNA_FP, 10)


MS_RNA <- subset(molcell_RNA_FP, select = c("ORF", "Median"))
MS_RNA$Albert <- sapply(molcell_RNA_FP[,2], function(x) x*(100/Albert.mode))
MS_RNA$Brar <- sapply(molcell_RNA_FP[,3], function(x) x*(100/Brar.mode))
MS_RNA$Ingolia <- sapply(molcell_RNA_FP[,4], function(x) x*(100/Ingolia.mode))
MS_RNA$Pop <- sapply(molcell_RNA_FP[,5], function(x) x*(100/Pop.mode))
MS_RNA$Weinberg <- sapply(molcell_RNA_FP[,6], function(x) x*(100/Weinberg.mode))
MS_RNA$Mean2 <- sapply(molcell_RNA_FP[,7], function(x) x*(100/Median.mode))
MS_RNA$Yassour <- sapply(molcell_RNA_FP[,8], function(x) x*(100/Yassour.mode))
MS_RNA$Lipson <- sapply(molcell_RNA_FP[,9], function(x) x*(100/Lipson.mode))
MS_RNA$Nagal <- sapply(molcell_RNA_FP[,10], function(x) x*(100/Nagal.mode))

MS_RNA$RNAseq <- apply(log(MS_RNA[,9:11]), 1, mean, na.rm = TRUE)
MS_RNA$RiboFP <- apply(log(MS_RNA[,3:7]), 1, mean, na.rm = TRUE)

MS_RNA2 <- subset(MS_RNA, select = c("ORF", "Mean2", "RNAseq", "RiboFP", "Median"))

lnMS_RNA2 <- MS_RNA2
lnMS_RNA2$Mean2 <- log(lnMS_RNA2$Mean2)
lnMS_RNA2[mapply(is.infinite, lnMS_RNA2)] <- NA


# Get the graph that compares RNA or riboFP to median abundance
RNA_fit <- lm(RNAseq ~ log(Median), data = lnMS_RNA2)
plot(RNAseq ~ log(Median), data = lnMS_RNA2, pch = 16, col = "grey")
abline(RNA_fit, lty = 6, col = "firebrick", lwd = 2)


FP_fit <- lm(RiboFP ~ log(Median), data = lnMS_RNA2)
plot(RiboFP ~ log(Median), data = lnMS_RNA2, pch = 16, col = "grey")
abline(FP_fit, lty = 6, col = "firebrick", lwd = 2)


# Merge data to Gene Names
# keep all.x = T to retain all rows in your data
ORF_df <- SGD_ORFcomp
ORF_df <- subset(SGD_ORFcomp, select = c("ORF", "Gene"))
colnames(ORF_df) <- c("ORF", "Gene")
data <- lnMS_RNA2
newData <- merge(data, ORF_df, by.x = "ORF", by.y = "ORF", all.x = T)


# replace "NA's" with ORF name
newData[is.na(newData$Gene), "Gene"] <- newData[is.na(newData$Gene), "ORF"]

newData <- subset(newData, select = c("Gene", "RNAseq", "RiboFP", "Mean2"))


# mahalanobis distances (alpha = 0.001, chi-square = 16.2662 df = 3)
# but to find things that deviate from distribution, use alpha = 0.01
# |________> (alpha = 0.01, crit val = 11.34, df = 3)
lnMS_RNA2_mat <- data.matrix(newData[,2:4])
rownames(lnMS_RNA2_mat) <- newData$Gene

lnMS_RNA2_mat <- lnMS_RNA2_mat[complete.cases(lnMS_RNA2_mat),]

Sx <- cov(lnMS_RNA2_mat)
mal_dist <- mahalanobis(lnMS_RNA2_mat, center = colMeans(lnMS_RNA2_mat),
						cov = cov(lnMS_RNA2_mat))

outliers <- data.frame(mal_dist)
outliers <- subset(outliers, mal_dist > 12.838)

# y will represent the entire dataset, x is the outliers
x <- newData[newData$Gene %in% rownames(outliers),]
y <- newData


# The following code will generate the table summarizing the results of
# this analysis
outliers_df <- data.frame(Gene = rownames(outliers), outliers)

out_df <- newData[newData$Gene %in% rownames(outliers),]
out_df2 <- merge(out_df, outliers_df, by.x = "Gene", by.y = "Gene")


# keep this dataframe, because we will soon add the cluster below


#Heatmap with k-means
mapkm.getk <-function(expmat, kRange){
    set.seed(100)
    clustOb <- lapply(kRange, function(k) { kmeans(expmat, k) })
    score <- sapply(clustOb, function(ob) { ob$tot.withinss })
    #scoreAll <- sapply(clustOb, function(ob) { ob$totss })
    plot(kRange, score, main=paste0("Sum of squares w/in clusters for a k-means range of ",min(kRange)," to ",max(kRange)), xlab="k", ylab="Total w/in sum of squares", pch=20)
}

z_df <- x
set.seed(100)
kclust <- kmeans(z_df[,2:4], 5, iter.max = 100000, nstart = 1)
z_df2 <- z_df
z_df2$Cluster <- kclust$cluster
z_df2 <- z_df2[order(z_df2$Cluster),]

z_df2_mat <- data.matrix(z_df2[,2:5])
rownames(z_df2_mat) <- z_df2$Gene

# now we want ot merge z_df2, and out_df2
mal <- subset(out_df2, select = c("Gene", "mal_dist"))
outliers <- z_df2

final_outliers <- merge(mal, outliers, by.x = "Gene", by.y = "Gene")

# doing the heatmap on all clusters at once is too computationally
# expensive. So break it up, subset each by cluster

outliersmat <- z_df2_mat

cluster1 <- outliersmat[which(outliersmat[,4] == 1),]
cluster2 <- outliersmat[which(outliersmat[,4] == 2),]
cluster3 <- outliersmat[which(outliersmat[,4] == 3),]
cluster4 <- outliersmat[which(outliersmat[,4] == 4),]
cluster5 <- outliersmat[which(outliersmat[,4] == 5),]


colours = c(seq(0, 12, length = 20))
my_palette <- colorRampPalette(c("yellow", "#ff8040" , "#800040"))(n = 19)


hm <- heatmap.2(outliersmat,
		  		col = my_palette,
		  		breaks = colours,
		  		symm = F, symkey = F, scale = "none",
		  		Colv = F, Rowv = F,
		  		density.info = "none",
		  		trace = "none",
		  		dendrogram = "none")


write.csv(outliersmat, file = "outliers.csv")


## Load scatterplot3d
library(scatterplot3d)

scatterplot3d(x = y$Mean2, y = y$RNAseq, z = y$RiboFP,
			  xlim = c(0, 12), ylim = c(0, 12), pch = 16, color = "grey",
			  angle = 45, zlim = c(0, 12))
par(new = TRUE)
scatterplot3d(cluster1[,3], cluster1[,1], cluster1[,2],
			  xlim = c(0, 12), ylim = c(0, 12), pch = 16, color = "#7654A2",
			  angle = 45, zlim = c(0, 12))
par(new = TRUE)
scatterplot3d(cluster2[,3], cluster2[,1], cluster2[,2],
			  xlim = c(0, 12), ylim = c(0, 12), pch = 16, color = "#F15A29",
			  angle = 45, zlim = c(0, 12))
par(new = TRUE)
scatterplot3d(cluster3[,3], cluster3[,1], cluster3[,2],
			  xlim = c(0, 12), ylim = c(0, 12), pch = 16, color = "#3A52A4",
			  angle = 45, zlim = c(0, 12))
par(new = TRUE)
scatterplot3d(cluster4[,3], cluster4[,1], cluster4[,2],
			  xlim = c(0, 12), ylim = c(0, 12), pch = 16, color = "#00A79D",
			  angle = 45, zlim = c(0, 12))
par(new = TRUE)
scatterplot3d(cluster5[,3], cluster5[,1], cluster5[,2],
			  xlim = c(0, 12), ylim = c(0, 12), pch = 16, color = "#FBB040",
			  angle = 45, zlim = c(0, 12))


setEPS()
postscript("whatever.eps")
plot3d(y$Mean2, y$RNAseq, y$RiboFP, col = "grey",
	   size = 4, xlim = c(0, 10), ylim = c(0, 10), zlim = c(0, 10))
dev.off()










































outliers_final <- data.frame(ORF = rownames(outliersmat), outliersmat)



halflife <- read.csv(file.choose(), na.strings = c("", "NA"),
					 stringsAsFactors = FALSE)
halflife2 <- read.csv(file.choose(), na.strings = c("", "NA"),
					 stringsAsFactors = FALSE)
halflife_merge <- merge(outliers_final, halflife, by.x = "ORF",
						by.y = "Gene.Name", all.x = TRUE)
halflife_merge2 <- merge(outliers_final, halflife2, by.x = "ORF",
						by.y = "Gene.name", all.x = TRUE)


x <- halflife
x <- subset(x, select = c("ORF", "Gene.Name", "Corrected.Half.Life"))
y <- halflife2
y <- subset(y, select = c("ENSG", "t1.2..min."))
halflife_meas <- merge(x, y, by.x = "ORF",
					   by.y = "ENSG", all.x = T, all.y = T)
colnames(halflife_meas) <- c("ORF", "Gene", "A", "B")

halflife_meas$A <- as.numeric(halflife_meas$A)
halflife_meas$B <- as.numeric(halflife_meas$B)

























































# Try clustering them first, via heatmap and hierarchical clustering
colours = c(seq(2, 10, length = 20))
my_palette <- colorRampPalette(c("white", "red"))(n = 19)

lnMS_RNA2_mat <- data.matrix(newData[,2:4])
rownames(lnMS_RNA2_mat) <-newData$Gene


dev.new(height = 4, width = 4)

#plot(RiboFP ~ RNAseq, data = MS_RNA2,
#	 cex = 0.5, pch = 16, col = "dark grey")
#abline(lm(RiboFP ~ RNAseq, data = MS_RNA2), lty = 3)

hm <- heatmap.2(lnMS_RNA2_mat,
		  		col = my_palette,
		  		breaks = colours,
		  		symm = F, symkey = F, scale = "none",
		  		Colv = F,
		  		density.info = "none",
		  		trace = "none")


write.csv(lnMS_RNA2_mat[rev(hm$rowInd), hm$colInd],
		  file = "lnMS_RNA2_mat.csv")



















