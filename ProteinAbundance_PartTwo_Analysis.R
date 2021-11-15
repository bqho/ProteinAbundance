#____________________________TABLE OF CONTENTS__________________________________
#
# loading data for downstream analysis ................................. line 31
# merging data into unified dataframe ................................. line 389
# apply mode shift normalization ...................................... line 439
# remove autofluorescence ............................................. line 506
# calculate molecules per cell ........................................ line 537
# assemble molcell dataframe .......................................... line 669
#
#
#
#
#

# The following packages will be required:

library(stringr) 
library(ggplot2) 
library(reshape2) 
library(beeswarm)
library(car)
library(raster) 
#library(qpcR)
library(LSD)
library(psych)
library(gplots)
library(cluster)
library(MASS)
library(RColorBrewer)

#_______________________________________________________________________________

#     ###
#      ##
#      ##            LOADING THE DATA FOR DOWNSTREAM ANALYSIS
#      ##
#    ######

#_______________________________________________________________________________

# NOTES: All excel/CSV files containing the raw data from the original studies
# were placed in a single folder for ease of access. Each CSV corresponding
# to the study investigated is named in a specific manner such that copying all
# of PART 1 of this script would load all dataframes quickly for analysis. To
# have this feature, all CSVs must be in the same folder (remember to set the
# proper working directory), and name all CSVs exactly as how it is written
# below (case sensitive)

# Otherwise, load each dataset individually, and you can search through your
# files manually with file.choose()

# Chong Data (2015) - Raw GFP Values
Chong2015 <- read.csv("Chong2015_complete.csv", na.strings = c("", "NA"), 
                      stringsAsFactors = FALSE)
Chong2015$Mean <- apply(Chong2015[,3:5], 1, mean, na.rm = TRUE)
Chong2015 <- subset(Chong2015, select = c("ORF", "Mean"))
colnames(Chong2015) <- c("ORF", "Chong2015")




# Ghaemmaghami Data (2003) - Molecules per cell
Ghaemmaghami2003 <- read.csv("Ghaemmaghami2003_Nature.csv", 
							 stringsAsFactors = FALSE,
                             na.strings = c("", "NA"))
Ghaemmaghami2003 <- subset(Ghaemmaghami2003,
                           select = c("ORF", "Protein.Molecules.Cell"))
colnames(Ghaemmaghami2003) <- c("ORF", "Ghaemmaghami2003")




# Kulak Data (2014) - Molecules per cell
Kulak2014 <- read.csv("Kulak2014_NatMethods.csv", stringsAsFactors = FALSE,
                      na.strings = c("", "NA"))
Kulak2014 <- subset(Kulak2014, select = c("ORF", "Copy.number"))
colnames(Kulak2014) <- c("ORF", "Kulak2014")

# Kulak data has some rows with numerous ORFs because they grouped all the
# paralogs together. Let's split these up
para <- Kulak2014[grep(";", Kulak2014$ORF, invert = FALSE),]
s <- strsplit(as.character(para$ORF), ';')
para <- data.frame(ORF = unlist(s),
                   Kulak2014 = rep(para$Kulak, sapply(s, FUN=length)))
Kulak2014 <- Kulak2014[grep(";", Kulak2014$ORF, invert = TRUE),]
Kulak2014 <- rbind(Kulak2014, para)




# Newman Data (2006) - Raw GFP values
Newman2006 <- read.csv("Newman2006_Nature.csv", stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
Newman2006 <- subset(Newman2006, select = c("ORF", "Ave..YEPD"))
colnames(Newman2006) <- c("ORF", "Newman2006")




# Tkach Data (2012) - Raw GFP values
Tkach2012 <- read.csv("Tkach2012_NatCellBio.csv", stringsAsFactors = FALSE,
                      na.strings = c("", "NA"))
Tkach2012 <- subset(Tkach2012, select = c("Systematic.ORF", "Control"))
colnames(Tkach2012) <- c("ORF", "Tkach2012")




# Breker Data (2012) - Raw GFP values
Breker2013 <- read.csv("Breker2013_JCellBiol.csv", stringsAsFactors = FALSE,
                       na.strings = c("", "NA"))
Breker2013 <- subset(Breker2013, select = c("ORF", "Control.Median"))
colnames(Breker2013) <- c("ORF", "Breker2013")




# Peng Data (2007) - Molecules per cell
Lu2007 <- read.csv("Lu2007_NatBioTech.csv", stringsAsFactors = FALSE,
                     na.strings = c("", "NA"))
Lu2007 <- subset(Lu2007, select = c("YORF", "APEX"))
Lu2007$APEX <- as.numeric(Lu2007$APEX)
colnames(Lu2007) <- c("ORF", "Lu2007")




# Mazumder Data (2012) - Raw GFP values
# Mazumder has some duplicate values that we need to remove
Mazumder2013 <- read.csv("Samson2013_NAR.csv", stringsAsFactors = FALSE,
                         na.strings = c("", "NA"))
dup.mazumder <- Mazumder2013[which(duplicated(Mazumder2013$ORF)),]
Mazumder2013 <- Mazumder2013[!(Mazumder2013$ORF %in% dup.mazumder$ORF),]
remove(dup.mazumder)

Mazumder2013 <- subset(Mazumder2013, select = c("ORF", "Cont.mean.Int"))
colnames(Mazumder2013) <- c("ORF", "Mazumder2013")




# Denervaud Data (2013) - Raw GFP values
# ** requires a little bit of processing first, since it's not measured at 
# ** one time point, but rather at numerous time points
Denervaud2013 <- read.csv("Denervaud2013_PNAS.csv", stringsAsFactors = FALSE,
                          na.strings = c("", "NA"))
Denervaud.control <- Denervaud2013[,1:19]
Denervaud.control <- Denervaud.control[complete.cases(Denervaud.control),]
Denervaud.control$Median <- as.numeric(apply(Denervaud.control, 1, 
                                             function(x) median(x, 
                                                                na.rm = TRUE)))
Denervaud2013 <- subset(Denervaud.control, select = c("ORF", "Median"))
remove(Denervaud.control)
colnames(Denervaud2013) <- c("ORF", "Denervaud2013")




# Lawless Data
Lawless2016 <- read.csv("Lawless2016_MCP.csv", stringsAsFactors = FALSE,
                        na.strings = c("", "NA"))
Lawless2016 <- subset(Lawless2016, select = c("Protein", "Final.Quant"))
colnames(Lawless2016) <- c("ORF", "Lawless2016")

# Lawless data has some rows with numerous ORFs because they grouped all the
# paralogs together. Let's split these up
para <- Lawless2016[grep("_", Lawless2016$ORF, invert = FALSE),]
s <- strsplit(as.character(para$ORF), '_')
para <- data.frame(ORF = unlist(s),
                   Lawless2016 = rep(para$Lawless, sapply(s, FUN=length)))
Lawless2016 <- Lawless2016[grep("_", Lawless2016$ORF, invert = TRUE),]
Lawless2016 <- rbind(Lawless2016, para)
remove(para, s)




# DeGodoy Data
deGodoy2008 <- read.csv("deGodoy2008_Nature.csv", stringsAsFactors = FALSE,
                        na.strings = c("", "NA"))
deGodoy2008 <- subset(deGodoy2008, select = c("ORF.Name", "Intensity.L"))
colnames(deGodoy2008) <- c("ORF", "deGodoy2008")




# Peng 2012
Peng2012 <- read.csv("Peng2012_NatMethods.csv", stringsAsFactors = FALSE,
                     na.strings = c("", "NA", "n.d."))
Peng2012 <- subset(Peng2012, select = c("ORF", "Average.of.all.data"))
colnames(Peng2012) <- c("ORF", "Peng2012")
Peng2012$Peng2012 <- as.numeric(Peng2012$Peng2012)




# Nagaraj 2011
Nagaraj2011 <- read.csv("Nagaraj2012_MCP.csv", stringsAsFactors = FALSE,
                        na.strings = c("", "NA"))
Nagaraj2011 <- subset(Nagaraj2011, select = c("Majority.Protein.IDs",
                                              "Intensity", "Intensity.run01",
                                              "Intensity.run02",
                                              "Intensity.run03",
                                              "Intensity.run04",
                                              "Intensity.run05"))
Nagaraj2011$Mean <- apply(Nagaraj2011[,2:ncol(Nagaraj2011)], 1, mean, na.rm = T)
Nagaraj2011 <- subset(Nagaraj2011, select = c("Majority.Protein.IDs", "Mean"))
colnames(Nagaraj2011) <- c("ORF", "Nagaraj2011")

# And get rid of the grouped paralogs
para <- Nagaraj2011[grep(";", Nagaraj2011$ORF, invert = FALSE),]
s <- strsplit(as.character(para$ORF), ';')
para <- data.frame(ORF = unlist(s),
                   Nagaraj2011 = rep(para$Nagaraj, sapply(s, FUN=length)))
Nagaraj2011 <- Nagaraj2011[grep(";", Nagaraj2011$ORF, invert = TRUE),]
Nagaraj2011 <- rbind(Nagaraj2011, para)
remove(para, s)

# Some values reported as 0, need to remove for analysis
Nagaraj2011 <- Nagaraj2011[is.finite(log(Nagaraj2011$Nagaraj2011)),]




# Lee 2011
# This data set has 3 biological replicates. We will average them all
Lee2011 <- read.csv("Lee2011_MolSystBio.csv", stringsAsFactors = FALSE,
                    na.strings = c("", "NA"))
Lee2011_2 <- read.csv("Lee2011_MolSystBio_2.csv", stringsAsFactors = FALSE,
                    na.strings = c("", "NA"))
Lee2011_3 <- read.csv("Lee2011_MolSystBio_3.csv", stringsAsFactors = FALSE,
                    na.strings = c("", "NA"))

ORFonly <- function(x){
    return(strsplit(x, " ")[[1]][1])
}

Lee2011$ORF <- as.character(lapply(Lee2011$Protein.Description, ORFonly))
Lee2011_2$ORF <- as.character(lapply(Lee2011_2$Protein.Description, ORFonly))
Lee2011_3$ORF <- as.character(lapply(Lee2011_3$Protein.Description, ORFonly))
Lee2011 <- subset(Lee2011, select = c("ORF", "TMT.130.Intensity.0.min"))
colnames(Lee2011) <- c("ORF", "Lee1")
Lee2011_2 <- subset(Lee2011_2, select = c("ORF", "TMT.130.Intensity.0.min"))
colnames(Lee2011_2) <- c("ORF", "Lee2")
Lee2011_3 <- subset(Lee2011_3, select = c("ORF", "TMT.130.Intensity.0.min"))
colnames(Lee2011_3) <- c("ORF", "Lee3")
Lee2011 <- Reduce(function(x,y) merge(x,y, all=TRUE),
                  list(Lee2011, Lee2011_2, Lee2011_3))

# Then take the mean for Lee's data
Lee2011$Lee <- apply(Lee2011[,2:4], 1, mean, na.rm = TRUE)
Lee2011 <- subset(Lee2011, select = c("ORF", "Lee"))
colnames(Lee2011) <- c("ORF", "Lee2011")
Lee2011 <- Lee2011[is.finite(log(Lee2011$Lee)),]




# Thakur 2011
Thakur2011 <- read.csv("Thakur2011_MCP.csv", stringsAsFactors = F,
                         na.strings = c("NA", ""))

Thakur2011 <- subset(Thakur2011, select = c("Protein.IDs", "Intensity.Exp1",
                                      "Intensity.Exp2", "Intensity.Exp3"))
Thakur2011$Thakur <- apply(Thakur2011[,2:ncol(Thakur2011)], 1, mean, na.rm = TRUE)
Thakur2011 <- subset(Thakur2011, select = c("Protein.IDs", "Thakur"))
colnames(Thakur2011) <- c("ORF", "Thakur2011")
Thakur2011 <- subset(Thakur2011, Thakur2011 > 0)

# And get rid of the grouped paralogs
para <- Thakur2011[grep(";", Thakur2011$ORF, invert = FALSE),]
s <- strsplit(as.character(para$ORF), ';')
para <- data.frame(ORF = unlist(s),
                   Thakur2011 = rep(para$Thakur, sapply(s, FUN=length)))
Thakur2011 <- Thakur2011[grep(";", Thakur2011$ORF, invert = TRUE),]
Thakur2011 <- rbind(Thakur2011, para)
remove(para, s)




# Yofe 2016
Yofe2016 <- read.csv("Yofe2016_NatBiotech.csv", 
					 na.strings = c("", "N/D", "NA"),
					 stringsAsFactors = FALSE)
Yofe2016 <- subset(Yofe2016, select = c("ORF", "SWAT.GFP.intensity.median..a.u..",
										"Seamless.GFP.intensity.median..a.u.."))
colnames(Yofe2016) <- c("ORF", "YofeSWAT", "YofeSeam")
Yofe2016 <- Yofe2016[Yofe2016$YofeSWAT > 0,]
Yofe2016 <- Yofe2016[Yofe2016$YofeSeam > 0,]

YofeSeam <- subset(Yofe2016, select = c("ORF", "YofeSeam"))
colnames(YofeSeam) <- c("ORF", "Yofe2016")
YofeSeam <- YofeSeam[complete.cases(YofeSeam$ORF),]




# Paulo et al 2016
# Correction, Paulo is not needed for analysis! Since it is not abundance meas.
Paulo2016 <- read.csv("Paulo2016_Jprot.csv", na.strings = c("", "NA"),
					  stringsAsFactors = FALSE)
Paulo2016 <- subset(Paulo2016, select = c("ORF", "glucose"))
colnames(Paulo2016) <- c("ORF", "Paulo2016")




# Davidson et al 2011
Davidson2011 <- read.csv("davidson2011_mboc.csv", na.strings = c("", "NA"),
						 stringsAsFactors = FALSE)
Davidson2011 <- subset(Davidson2011, select = c("ID2", "AVG.Mean"))
colnames(Davidson2011) <- c("ORF", "Davidson2011")




# Lee 2007 Dataset from Yeast Journal
Lee2007 <- read.csv("Lee2007_Yeast.csv", na.strings = c("", "NA"),
					stringsAsFactors = FALSE)
Lee2007 <- subset(Lee2007, select = c("ORF", "noMMS_norm"))
colnames(Lee2007) <- c("ORF", "Lee2007")




# Webb 2013 Dataset from J Proteome Research Journal
Webb2013_1 <- read.csv("Webb2013_JProtRes_1.csv", na.strings = c("", "NA", "*"),
					   stringsAsFactors = FALSE)
Webb2013_1 <- subset(Webb2013_1, select = c("Locus", "EMPAI"))
Webb2013_1 <- Webb2013_1[complete.cases(Webb2013_1$Locus),]
colnames(Webb2013_1) <- c("ORF", "EMPAI1")


Webb2013_2 <- read.csv("Webb2013_JProtRes_2.csv", na.strings = c("", "NA", "*"),
					   stringsAsFactors = FALSE)
Webb2013_2 <- subset(Webb2013_2, select = c("Locus", "EMPAI"))
Webb2013_2 <- Webb2013_2[complete.cases(Webb2013_2$Locus),]
colnames(Webb2013_2) <- c("ORF", "EMPAI2")


Webb2013_3 <- read.csv("Webb2013_JProtRes_3.csv", na.strings = c("", "NA", "*"),
					   stringsAsFactors = FALSE)
Webb2013_3 <- subset(Webb2013_3, select = c("Locus", "EMPAI"))
Webb2013_3 <- Webb2013_3[complete.cases(Webb2013_3$Locus),]
colnames(Webb2013_3) <- c("ORF", "EMPAI3")


Webb2013_combined <- Reduce(function(x,y) merge(x, y, all = TRUE),
							list(Webb2013_1, Webb2013_2, Webb2013_3))
Webb2013_combined$Mean <- apply(Webb2013_combined[,2:4], 1, mean, na.rm = TRUE)

Webb2013 <- subset(Webb2013_combined, select = c("ORF", "Mean"))
colnames(Webb2013) <- c("ORF", "Webb2013")




# Lahtvee 2017 dataset from cell systems
Lahtvee2017 <- read.csv("Lahtvee2017_CellSys.csv", na.strings = c("", "NA"),
                        stringsAsFactors = FALSE)
Lahtvee2017 <- subset(Lahtvee2017, select = c("Gene.ID", "REF"))
Lahtvee2017$Lahtvee2017 <- 13*(Lahtvee2017$REF)
Lahtvee2017 <- subset(Lahtvee2017, select = c("Gene.ID", "Lahtvee2017"))
colnames(Lahtvee2017) <- c("ORF", "Lahtvee2017")




# Picotti 2013 dataset from Nature
Picotti2013 <- read.csv("Picotti2013_Nature.csv", na.strings = c("", "NA"),
                        stringsAsFactors = FALSE)
Picotti2013 <- subset(Picotti2013, select = c("Accession", "BY", "BY.2"))
Picotti2013$BY_t <- 2^Picotti2013$BY
Picotti2013$BY2_t <- 2^Picotti2013$BY.2
Picotti2013$Avg <- apply(Picotti2013[,4:5], 1, mean, na.rm = T)
Picotti2013 <- subset(Picotti2013, select = c("Accession", "Avg"))
colnames(Picotti2013) <- c("ORF", "Picotti2013")










#_______________________________________________________________________________

#    ####
#   ##  ##
#      ##      MERGE ALL DATASETS (BOTH RAW AND MPC) INTO ONE DATAFRAME           
#    ##
#   ######

#_______________________________________________________________________________

#######---------------------------------------------- Merge data subset proteome
SGD_ORF <- read.csv("SGD_ORFeome.csv", na.strings = c("", "NA"), 
                    stringsAsFactors = FALSE)
SGD_ORFcomp <- subset(SGD_ORF, 
                      Qualifier == "Verified" | Qualifier == "Uncharacterized")
colnames(SGD_ORFcomp) <- c("ORF", "Gene", "Qualifier")


combined_raw <- Reduce(function(x,y) merge(x,y, all=TRUE),
                   	   list(Kulak2014, Lu2007, Lawless2016, Peng2012, 
                            Lahtvee2017,
                   	   		deGodoy2008, Nagaraj2011, Lee2011, Thakur2011,
                   	   		Webb2013, Picotti2013, Tkach2012, Breker2013, 
                            Mazumder2013, Denervaud2013, Chong2015, Newman2006, 
                   	   		YofeSeam, Davidson2011, Lee2007, Ghaemmaghami2003))

combined_raw <- combined_raw[combined_raw$ORF %in% SGD_ORFcomp$ORF,]
combined_raw <- combined_raw[!duplicated(combined_raw$ORF),]

# just re-organize the data to the format we want for the figures

combined_raw <- subset(combined_raw, 
					   select = c("ORF", "Lu2007", "Peng2012", "Kulak2014",
								  "Lawless2016", "Lahtvee2017", "deGodoy2008", 
                                  "Lee2011", "Thakur2011", "Nagaraj2011", 
                                  "Picotti2013", "Webb2013",
								  "Tkach2012", "Breker2013", 
								  "Denervaud2013", "Mazumder2013",
								  "Chong2015", "Yofe2016", "Newman2006",
								  "Lee2007", "Davidson2011", 
								  "Ghaemmaghami2003"))

remove(Webb2013_1, Webb2013_2, Webb2013_3, Webb2013_combined, Yofe2016,
       Lee2011_2, Lee2011_3, Paulo2016)

#######---------------------------------------------- Merge data subset proteome




#######------------------------------------------------- Apply GetMode and Shift
deGodoy.mode <- get.mode(combined_raw, 7)
Lee.mode <- get.mode(combined_raw, 8)
Thakur.mode <- get.mode(combined_raw, 9)
Nagaraj.mode <- get.mode(combined_raw, 10)
Picotti.mode <- get.mode(combined_raw, 11)
Webb.mode <- get.mode(combined_raw, 12)
Tkach.mode <- get.mode(combined_raw, 13)
Breker.mode <- get.mode(combined_raw, 14)
Denervaud.mode <- get.mode(combined_raw, 15)
Mazumder.mode <- get.mode(combined_raw, 16)
Chong.mode <- get.mode(combined_raw, 17)
Yofe.mode <- get.mode(combined_raw, 18)
Newman.mode <- get.mode(combined_raw, 19)
Lee2007.mode <- get.mode(combined_raw, 20)[2]
Davidson.mode <- get.mode(combined_raw, 21)
dev.off()

# now we have all the "modes" for each data frame, and we can now shift all the
# data values so the modes are matching. Centre all modes to A.U. 100

ModeShift.DF <- data.frame(combined_raw[,1])
ModeShift.DF$deGodoy <- sapply(combined_raw[,7], 
							   function(x) x*(100/deGodoy.mode))
ModeShift.DF$Lee <- sapply(combined_raw[,8], 
						   function(x) x*(100/Lee.mode))
ModeShift.DF$Thakur <- sapply(combined_raw[,9], 
							  function(x) x*(100/Thakur.mode))
ModeShift.DF$Nagaraj <- sapply(combined_raw[,10], 
							   function(x) x*(100/Nagaraj.mode))
ModeShift.DF$Picotti <- sapply(combined_raw[,11], 
							   function(x) x*(100/Picotti.mode))
ModeShift.DF$Webb <- sapply(combined_raw[,12], 
							function(x) x*(100/Webb.mode))
ModeShift.DF$Tkach <- sapply(combined_raw[,13], 
							 function(x) x*(100/Tkach.mode))
ModeShift.DF$Breker <- sapply(combined_raw[,14], 
							  function(x) x*(100/Breker.mode))
ModeShift.DF$Denervaud <- sapply(combined_raw[,15],
								 function(x)x*(100/Denervaud.mode))
ModeShift.DF$Mazumder <- sapply(combined_raw[,16],
							    function(x) x*(100/Mazumder.mode))
ModeShift.DF$Chong <- sapply(combined_raw[,17],
						     function(x) x*(100/Chong.mode))
ModeShift.DF$Yofe <- sapply(combined_raw[,18], 
							function(x) x*(100/Yofe.mode))
ModeShift.DF$Newman <- sapply(combined_raw[,19], 
							  function(x) x*(100/Newman.mode))
ModeShift.DF$Lee2007 <- sapply(combined_raw[,20], 
							   function(x) x*(100/Lee2007.mode))
ModeShift.DF$Davidson <- sapply(combined_raw[,21], 
								function(x) x*(100/Davidson.mode))


colnames(ModeShift.DF) <- c("ORF", "deGodoyAU", "LeeAU", "ThakurAU","NagarajAU", 
                            "PicottiAU", "WebbAU", "TkachAU", "BrekerAU", 
                            "DenervaudAU", "MazumderAU", "ChongAU", "YofeAU",
                            "NewmanAU", "Lee2007AU", "DavidsonAU")

# And to document the data, write the CSV and save table
# write.csv(ModeShift.DF, file = "ModeShifted_GFPdatasets.csv")

#######------------------------------------------------- Apply GetMode and Shift




#________________________________________________________Remove autofluorescence

calibration_set <- subset(combined_raw, 
                          select = c("ORF", "Kulak2014", "Lu2007",
                                     "Lawless2016", "Peng2012", "Lahtvee2017"))

# identify which ORFs to remove
auto <- read.csv("AutoFluor.csv", stringsAsFactors = FALSE,
                 na.strings = c("", "NA"))

auto_Mode <- ModeShift.DF[ModeShift.DF$ORF %in% auto$ORF,]

max(auto_Mode$Chong, na.rm = T)
# so it looks like the autofluorescence values are at < 106.5494 AU

GFP_mat <- data.matrix(ModeShift.DF[8:16])
rownames(GFP_mat) <- ModeShift.DF$ORF
GFP_mat[GFP_mat < 106.5494] = NA

MS_AU <- ModeShift.DF[,1:7]

GFP_DF <- data.frame(ORF = rownames(GFP_mat),
                           GFP_mat)

Mode_final <- merge(MS_AU, GFP_DF, all.x = T, all.y = T)

#________________________________________________________Remove autofluorescence




#__________________________________________________ Calculate molecules per cell
# Merging the data, such that all data sets are assembled into one DF
compare.df <- merge(calibration_set, Mode_final, all.y = T)
compare.df <- compare.df[!duplicated(compare.df$ORF),]

compare.df$Calibration <- apply(compare.df[,2:6], 1, mean, na.rm = T)
compare.df$GFPAUmean <- apply(compare.df[,13:21], 1, mean, na.rm = T)
compare.df$MSAUmean <- apply(compare.df[,7:12], 1, mean, na.rm = T)
compare.df$ALLmean <- apply(compare.df[,7:21], 1, mean, na.rm = T)

dev.new(width = 4, height = 4)
heatscatter(log(compare.df$ALLmean), log(compare.df$Calibration),
            xlim = c(-2,12), ylim = c(0, 15))
abline(lm(log(Calibration) ~ log(ALLmean), data = compare.df), lwd = 1,
       lty = 2, col = "firebrick")

cor(log(compare.df$ALLmean), log(compare.df$Calibration), use = "complete.obs")
lm_mode <- lm(log(Calibration) ~ log(ALLmean), data = compare.df)
dev.off()

# Merging the data, such that all data sets are assembled into one DF
compare.df <- merge(calibration_set, Mode_final, all.y = T)

compare.df$Calibration <- apply(compare.df[,2:6], 1, mean, na.rm = T)


fit_DGD <- lm(log(compare.df[,22]) ~ log(compare.df[,7]),
              data = compare.df)
DGD <- fit_DGD$coefficients


fit_LEE2 <- lm(log(compare.df[,22]) ~ log(compare.df[,8]),
              data = compare.df)
LEE2 <- fit_LEE2$coefficients


fit_THAK <- lm(log(compare.df[,22]) ~ log(compare.df[,9]),
              data = compare.df)
THAK <- fit_THAK$coefficients


fit_NAG <- lm(log(compare.df[,22]) ~ log(compare.df[,10]),
              data = compare.df)
NAG <- fit_NAG$coefficients


fit_PIC <- lm(log(compare.df[,22]) ~ log(compare.df[,11]),
              data = compare.df)
PIC <- fit_PIC$coefficients


fit_WEBB <- lm(log(compare.df[,22]) ~ log(compare.df[,12]),
              data = compare.df)
WEBB <- fit_WEBB$coefficients


fit_TKAC <- lm(log(compare.df[,22]) ~ log(compare.df[,13]),
              data = compare.df)
TKAC <- fit_TKAC$coefficients


fit_BRE <- lm(log(compare.df[,22]) ~ log(compare.df[,14]),
              data = compare.df)
BRE <- fit_BRE$coefficients


fit_DEN <- lm(log(compare.df[,22]) ~ log(compare.df[,15]),
              data = compare.df)
DEN <- fit_DEN$coefficients


fit_MAZ <- lm(log(compare.df[,22]) ~ log(compare.df[,16]),
              data = compare.df)
MAZ <- fit_MAZ$coefficients


fit_CHO <- lm(log(compare.df[,22]) ~ log(compare.df[,17]),
              data = compare.df)
CHO <- fit_CHO$coefficients


fit_YOF <- lm(log(compare.df[,22]) ~ log(compare.df[,18]),
              data = compare.df)
YOF <- fit_YOF$coefficients


fit_NEW <- lm(log(compare.df[,22]) ~ log(compare.df[,19]),
              data = compare.df)
NEW <- fit_NEW$coefficients


fit_LEE <- lm(log(compare.df[,22]) ~ log(compare.df[,20]),
              data = compare.df)
LEE <- fit_LEE$coefficients


fit_DAV <- lm(log(compare.df[,22]) ~ log(compare.df[,21]),
              data = compare.df)
DAV <- fit_DAV$coefficients



# Generate new DF, preparing to calculate molecules/cell for the AU data
MolCell_AU <- data.frame(Mode_final$ORF)

MolCell_AU$deGodoy <- exp(1)^(DGD[2]*log(Mode_final[,2]) + DGD[1])
MolCell_AU$Lee <- exp(1)^(LEE2[2]*log(Mode_final[,3]) + LEE2[1])
MolCell_AU$Thakur <- exp(1)^(THAK[2]*log(Mode_final[,4]) + THAK[1])
MolCell_AU$Nagaraj <- exp(1)^(NAG[2]*log(Mode_final[,5]) + NAG[1])
MolCell_AU$Picotti <- exp(1)^(PIC[2]*log(Mode_final[,6]) + PIC[1])
MolCell_AU$Webb <- exp(1)^(WEBB[2]*log(Mode_final[,7]) + WEBB[1])
MolCell_AU$Tkach <- exp(1)^(TKAC[2]*log(Mode_final[,8]) + TKAC[1])
MolCell_AU$Breker <- exp(1)^(BRE[2]*log(Mode_final[,9]) + BRE[1])
MolCell_AU$Denervaud <- exp(1)^(DEN[2]*log(Mode_final[,10]) + DEN[1])
MolCell_AU$Mazumder <- exp(1)^(MAZ[2]*log(Mode_final[,11]) + MAZ[1])
MolCell_AU$Chong <- exp(1)^(CHO[2]*log(Mode_final[,12]) + CHO[1])
MolCell_AU$Yofe <- exp(1)^(YOF[2]*log(Mode_final[,13]) + YOF[1])
MolCell_AU$Newman <- exp(1)^(NEW[2]*log(Mode_final[,14]) + NEW[1])
MolCell_AU$Lee2007 <- exp(1)^(LEE[2]*log(Mode_final[,15]) + LEE[1])
MolCell_AU$Davidson <- exp(1)^(DAV[2]*log(Mode_final[,16]) + DAV[1])


colnames(MolCell_AU) <- c("ORF", "deGodoy", "Lee", "Thakur", "Nagaraj",
                           "Picotti", "Webb", "Tkach", "Breker", 
                           "Denervaud", "Mazumder", "Chong", "Yofe",
                           "Newman", "Lee2007", "Davidson")

#__________________________________________________ Calculate molecules per cell




#__________________________________________________ Assemble final mpc dataframe

merge1 <- merge(calibration_set, MolCell_AU, all.x = T, all.y = T)

Ghaemmaghami <- subset(combined_raw, select = c("ORF", "Ghaemmaghami2003"))

merge2 <- merge(merge1, Ghaemmaghami, all.x = T, all.y = T)

# Now we have the molcell dataframe. Let's see what it looks like, when
# we do the heatscatter

final_ALL <- merge2
colnames(final_ALL) <- c("ORF", "KUL", "LU", "LAW", "PENG", "LAHT", "DGD",
                         "LEE2", "THAK", "NAG", "PIC", "WEB", "TKA", "BRE",
                         "DEN", "MAZ", "CHO", "YOF", "NEW", "LEE", "DAV",
                         "GHA")


finalDF <- merge(final_ALL, SGD_ORFcomp, all.x = T, all.y = T)
# write.csv(finalDF, file = "final_molcell.csv")


final <- final_ALL
final$lnMSMean <- log(apply(final[,7:12], 1, mean, na.rm = T))
final$lnGFPMean <- log(apply(final[,13:21], 1, mean, na.rm = T))
final$lnCalibration <- log(apply(final[,2:6], 1, mean, na.rm = T))
final$lnALL <- log(apply(final[,7:21], 1, mean, na.rm = T))

dev.new(width = 4, height = 4)
heatscatter(final$lnALL, final$lnCalibration, cex = 0.1)
abline(lm(lnCalibration ~ lnALL, data = final), lwd = 1,
       lty = 2, col = "firebrick")

cor(final$lnCalibration, final$lnALL, use = "complete.obs")
dev.off()

#__________________________________________________ Assemble final mpc dataframe




#__________________________________________________________________ Calculations

final_ALL$Mean <- apply(final_ALL[,2:22], 1, mean.negrm)
final_ALL$Median <- apply(final_ALL[,2:22], 1, median.negrm)
final_ALL$Variance <- apply(final_ALL[,2:22], 1, var.negrm)
final_ALL$CoefficientVariation <- apply(final_ALL[,2:22], 1, 
                                        cv.negrm)
final_ALL$CV.complete <- apply(final_ALL[,2:22], 1, cv)

# With all the data now ready for further analysis and data visualization, let's
# order the dataframe by increasing protein abundance (either by mean or median)

molcell <- final_ALL # rename final_ALL, so that final_ALL is does not get 
                     # modified in any way. Our "master" table

molcell <- molcell[order(molcell$Median),] # Choose to order by mean or median
molcell$Order <- seq(1:nrow(molcell))


remove(auto, auto_Mode, Breker.mode, Chong.mode, Davidson.mode,
       deGodoy.mode, Denervaud.mode, Lee2007.mode, 
       Newman.mode, Picotti.mode, Thakur.mode, Tkach.mode, 
       Webb.mode, Yofe.mode, Mazumder.mode, Nagaraj.mode, 
       Lee.mode, merge1, merge2, GFP_mat, final, finalDF, GFP_DF,
       compare.df, MS_AU, MolCell_AU, ModeShift.DF, SGD_ORF)

#__________________________________________________________________ Calculations
