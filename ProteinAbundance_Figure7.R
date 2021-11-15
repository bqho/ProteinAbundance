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

#_______________________________________________________________________________

#      ####
#        ##
#        ##                 Load the datasets for conditional abundance analysis
#        ##
#        ##
#      ######

#_______________________________________________________________________________


# Chong data - HU and Rapamycin
chongSTRESS <- read.csv("Chong2015_complete.csv", na.strings = c(""),
                        stringsAsFactors = FALSE)
chongSTRESS <- subset(chongSTRESS, select = c("ORF", "WT1", "WT2", "WT3", 
                                              "HU160", "RAP700"))
chongSTRESS$UNTmean <- apply(chongSTRESS[,2:4], 1, mean, na.rm = T)
chongSTRESS <- subset(chongSTRESS, select = c("ORF", "UNTmean", "HU160",
                                              "RAP700"))
colnames(chongSTRESS) <- c("ORF", "chongUNT", "chongHU", "chongRAP700")


# Tkach data - MMS and HU
tkachSTRESS <- read.csv("Tkach2012_NatureCellBio.csv", na.strings = c("", "NA"),
                        stringsAsFactors = FALSE)
tkachSTRESS <- subset(tkachSTRESS, select = c("Systematic.ORF", 
                                              "Control", "HU", "MMS"))
colnames(tkachSTRESS) <- c("ORF", "tkachUNT", "tkachHU", "tkachMMS")
tkachSTRESS <- tkachSTRESS[!duplicated(tkachSTRESS$ORF),]


# Mazumder data - MMS
mazumderSTRESS <- read.csv("Samson2013_NAR.csv", na.strings = c("", "NA"),
                           stringsAsFactors = FALSE)
mazumderSTRESS <- subset(mazumderSTRESS, select = c("ORF", 
                                                    "Cont.mean.Int",
                                                    "MMS.mean.Int"))
colnames(mazumderSTRESS) <- c("ORF", "mazumderUNT", "mazumderMMS")
mazumderSTRESS <- mazumderSTRESS[!duplicated(mazumderSTRESS$ORF),]


# Breker data - DTT, H2O2, starvation
brekerSTRESS <- read.csv("Breker2013_JCellBiol.csv", na.strings = c("", "NA"),
                         stringsAsFactors = FALSE)
brekerSTRESS <- subset(brekerSTRESS, select = c("ORF", 
                                                "Control.Median",
                                                "DTT.Median", 
                                                "H2O2.Median", 
                                                "Starvation.Median"))
colnames(brekerSTRESS) <- c("ORF", "brekerUNT", "brekerDTT", "brekerH2O2", 
                            "brekerSTARVE")


# Denervaud data - MMS (low/high/Pulse), HU, UV, Sorbitol
denervaudMMS_High <- read.csv("Denervaud_HighMMS.csv", na.strings = c("", "NA"),
                        stringsAsFactors = FALSE)
denervaudMMS_High <- subset(denervaudMMS_High, select = c("ORF", 
                                                          "X.0.5", "X2.1666667"))
denervaudMMS_High$ORF <- gsub("'", "", denervaudMMS_High$ORF)


# Davidson Data - Stationary vs Exponential Phase
davidson_EXP <- read.csv("davidson2011_EXP.csv", stringsAsFactors = FALSE,
                         na.strings = c("", "NA"))
davidson_EXP <- subset(davidson_EXP, select = c("ID2", "AVG.Mean"))
colnames(davidson_EXP) <- c("ORF", "EXP")
davidson_SP <- read.csv("davidson2011_SP.csv", stringsAsFactors = FALSE,
                        na.strings = c("", "NA"))
davidson_SP <- subset(davidson_SP, select = c("ORF.Name", "AVG.mean"))
colnames(davidson_SP) <- c("ORF", "SP")
davidsonSTRESS <- merge(davidson_EXP, davidson_SP, all.x = TRUE, all.y = TRUE)


# Lee Data - MMS
leeSTRESS <- read.csv("Lee2007_Yeast.csv", stringsAsFactors = FALSE,
                   na.strings = c("", "NA", "<NA>"))
leeSTRESS <- subset(leeSTRESS, select = c("ORF", "noMMS_norm", "MMS_norm"))
colnames(leeSTRESS) <- c("ORF", "leeUNT", "leeMMS")








#_______________________________________________________________________________

#      ######
#          ##
#         ###                      Convert arbitrary units to molecules per cell
#      ####
#      ##
#      ######

#_______________________________________________________________________________

# Make a data frame with all stresses, merging all the data from each condtional
# study, examining protein abundance

comb_rawStress <- Reduce(function(x,y) merge(x,y, all=TRUE),
                       list(chongSTRESS, tkachSTRESS, leeSTRESS, mazumderSTRESS,
                            brekerSTRESS, denervaudMMS_High, davidsonSTRESS))
comb_rawStress <- comb_rawStress[comb_rawStress$ORF %in% SGD_ORFcomp$ORF,]




chongUNT_mode <- get.mode(comb_rawStress, 2)
chongHU_mode <- get.mode(comb_rawStress, 3)
chongRAP_mode <- get.mode(comb_rawStress, 4)
tkachUNT_mode <- get.mode(comb_rawStress, 5)
tkachHU_mode <-get.mode(comb_rawStress, 6)
tkachMMS_mode <- get.mode(comb_rawStress, 7)
leeUNT_mode <- get.mode(comb_rawStress, 8)[2]
leeMMS_mode <- get.mode(comb_rawStress, 9)
mazumderUNT_mode <- get.mode(comb_rawStress, 10)
mazumderMMS_mode <- get.mode(comb_rawStress, 11)
brekerUNT_mode <- get.mode(comb_rawStress, 12)
brekerDTT_mode <- get.mode(comb_rawStress, 13)
brekerH2O2_mode <- get.mode(comb_rawStress, 14)
brekerSTARVE_mode <- get.mode(comb_rawStress, 15)
denervaudMMSHighUNT_mode <- get.mode(comb_rawStress, 16)
denervaudMMSHigh_mode <- get.mode(comb_rawStress, 17)
davidsonUNT_mode <- get.mode(comb_rawStress, 18)
davidsonSP_mode <- get.mode(comb_rawStress, 19)

# Assemble into one "stress data" matrix

ms_stress <- data.frame(ORF = comb_rawStress$ORF)
ms_stress$choUNT <- sapply(comb_rawStress[,2], 
                           function(x) x*(100/chongUNT_mode))
ms_stress$choHU <- sapply(comb_rawStress[,3], 
                           function(x) x*(100/chongHU_mode))
ms_stress$choRAP <- sapply(comb_rawStress[,4], 
                           function(x) x*(100/chongRAP_mode))

ms_stress$tkaUNT <- sapply(comb_rawStress[,5], 
                           function(x) x*(100/tkachUNT_mode))
ms_stress$tkaHU <- sapply(comb_rawStress[,6], 
                           function(x) x*(100/tkachHU_mode))
ms_stress$tkaMMS <- sapply(comb_rawStress[,7], 
                           function(x) x*(100/tkachMMS_mode))

ms_stress$leeUNT <- sapply(comb_rawStress[,8], 
                           function(x) x*(100/leeUNT_mode))
ms_stress$leeMMS <- sapply(comb_rawStress[,9], 
                           function(x) x*(100/leeMMS_mode))

ms_stress$mazUNT <- sapply(comb_rawStress[,10], 
                           function(x) x*(100/mazumderUNT_mode))
ms_stress$mazMMS <- sapply(comb_rawStress[,11], 
                           function(x) x*(100/mazumderMMS_mode))

ms_stress$breUNT <- sapply(comb_rawStress[,12], 
                           function(x) x*(100/brekerUNT_mode))
ms_stress$breDTT <- sapply(comb_rawStress[,13], 
                           function(x) x*(100/brekerDTT_mode))
ms_stress$breH2O2 <- sapply(comb_rawStress[,14], 
                           function(x) x*(100/brekerH2O2_mode))
ms_stress$breSTARVE <- sapply(comb_rawStress[,15], 
                           function(x) x*(100/brekerSTARVE_mode))

ms_stress$denUNT <- sapply(comb_rawStress[,16], 
                           function(x) x*(100/denervaudMMSHighUNT_mode))
ms_stress$denMMS <- sapply(comb_rawStress[,17], 
                           function(x) x*(100/denervaudMMSHigh_mode))

ms_stress$davUNT <- sapply(comb_rawStress[,18], 
                           function(x) x*(100/davidsonUNT_mode))
ms_stress$davSP <- sapply(comb_rawStress[,19], 
                           function(x) x*(100/davidsonSP_mode))


# And now remove autofluorescence proteins (106.5494)
stress_mat <- data.matrix(ms_stress[,2:19])
rownames(stress_mat) <- ms_stress$ORF

stress_mat[stress_mat < 106.5494] = NA

stress_df <- data.frame(ORF = rownames(stress_mat), stress_mat)

# convert to molecules per cell
mpc_stress <- data.frame(ORF = stress_df$ORF)
mpc_stress$choUNT <- exp(1)^(CHO[2]*log(stress_df[,2]) + CHO[1])
mpc_stress$choHU <- exp(1)^(CHO[2]*log(stress_df[,3]) + CHO[1])
mpc_stress$choRAP <- exp(1)^(CHO[2]*log(stress_df[,4]) + CHO[1])

mpc_stress$tkaUNT <- exp(1)^(TKAC[2]*log(stress_df[,5]) + TKAC[1])
mpc_stress$tkaHU <- exp(1)^(TKAC[2]*log(stress_df[,6]) + TKAC[1])
mpc_stress$tkaMMS <- exp(1)^(TKAC[2]*log(stress_df[,7]) + TKAC[1])

mpc_stress$leeUNT <- exp(1)^(LEE[2]*log(stress_df[,8]) + LEE[1])
mpc_stress$leeMMS <- exp(1)^(LEE[2]*log(stress_df[,9]) + LEE[1])

mpc_stress$mazUNT <- exp(1)^(MAZ[2]*log(stress_df[,10]) + MAZ[1])
mpc_stress$mazMMS <- exp(1)^(MAZ[2]*log(stress_df[,11]) + MAZ[1])

mpc_stress$breUNT <- exp(1)^(BRE[2]*log(stress_df[,12]) + BRE[1])
mpc_stress$breDTT <- exp(1)^(BRE[2]*log(stress_df[,13]) + BRE[1])
mpc_stress$breH2O2 <- exp(1)^(BRE[2]*log(stress_df[,14]) + BRE[1])
mpc_stress$breSTARVE <- exp(1)^(BRE[2]*log(stress_df[,15]) + BRE[1])

mpc_stress$denUNT <- exp(1)^(DEN[2]*log(stress_df[,16]) + DEN[1])
mpc_stress$denMMS <- exp(1)^(DEN[2]*log(stress_df[,17]) + DEN[1])

mpc_stress$davUNT <- exp(1)^(DAV[2]*log(stress_df[,18]) + DAV[1])
mpc_stress$davSP <- exp(1)^(DAV[2]*log(stress_df[,19]) + DAV[1])

combinedSTRESS <- mpc_stress

# Rearrange for the order of the table

combinedSTRESS <- subset(combinedSTRESS, 
                         select = c("ORF", "breUNT", "breDTT",
                                    "breH2O2", "breSTARVE",
                                    "choUNT", "choHU", "choRAP",
                                    "davUNT", "davSP", 
                                    "denUNT", "denMMS",
                                    "leeUNT", "leeMMS", "mazUNT",
                                    "mazMMS", "tkaUNT", "tkaHU",
                                    "tkaMMS"))

# then merge with gene names
final_combined <- merge(SGD_ORFcomp, combinedSTRESS, by.x = "ORF",
                        by.y = "ORF", all.x = T)
final_combined <- final_combined[!duplicated(final_combined$ORF),]

# write.csv(final_combined, file = "stress_mpc.csv")


# Now determine the fold-changes

Ratio <- subset(final_combined, select = c("ORF", "Gene"))
Ratio$bDTT <- (final_combined$breDTT/final_combined$breUNT)
Ratio$bH2O2 <- (final_combined$breH2O2/final_combined$breUNT)
Ratio$bSTARVE <- (final_combined$breSTARVE/final_combined$breUNT)

Ratio$cHU <- (final_combined$choHU/final_combined$choUNT)
Ratio$cRAP <- (final_combined$choRAP/final_combined$choUNT)

Ratio$dSP <- (final_combined$davSP/final_combined$davUNT)

Ratio$dMMS <- (final_combined$denMMS/final_combined$denUNT)

Ratio$mMMS <- (final_combined$mazMMS/final_combined$mazUNT)

Ratio$lMMS <-  (final_combined$leeMMS/final_combined$leeUNT)

Ratio$tHU <- (final_combined$tkaHU/final_combined$tkaUNT)
Ratio$tMMS <- (final_combined$tkaMMS/final_combined$tkaUNT)



# this is just to confirm that these molecules per cell are
# equivalent to the previously measured observations

# test <- read.csv(file.choose(), na.strings = c("", "NA"),
#                 stringsAsFactors = FALSE)

# x <- merge(final_combined, test, by.x = "ORF", by.y = "ORF",
#      all.x = T, all.y = T)



# We will determine whether there are protein changes based on
# a simple standard deviation analysis

final_combined <- subset(final_combined,
                         select = c("ORF", "Gene", "breUNT",
                                    "choUNT", "davUNT", "denUNT",
                                    "mazUNT", "leeUNT", "tkaUNT",
                                    "breDTT", "breH2O2", "breSTARVE",
                                    "choHU", "choRAP", "davSP", "denMMS",
                                    "mazMMS", "leeMMS", "tkaHU", "tkaMMS"))

final_combined$Mean <- apply(log(final_combined[,3:9]), 1, mean, na.rm = T)
final_combined$SD <- apply(log(final_combined[,3:9]), 1, sd, na.rm = T)
final_combined$UP <- final_combined$Mean + 2*final_combined$SD
final_combined$DOWN <- final_combined$Mean - 2*final_combined$SD



final_combined$uDTT <- log(final_combined$breDTT) - final_combined$UP
final_combined$uH2O2 <- log(final_combined$breH2O2) - final_combined$UP
final_combined$uSTARVE <- log(final_combined$breSTARVE) - final_combined$UP

final_combined$uHU <- log(final_combined$choHU) - final_combined$UP
final_combined$uRAP <- log(final_combined$choRAP) - final_combined$UP

final_combined$uSP <- log(final_combined$davSP) - final_combined$UP

final_combined$udMMS <- log(final_combined$denMMS) - final_combined$UP

final_combined$umMMS <- log(final_combined$mazMMS) - final_combined$UP

final_combined$ulMMS <- log(final_combined$leeMMS) - final_combined$UP

final_combined$utHU <- log(final_combined$tkaHU) - final_combined$UP
final_combined$utMMS <- log(final_combined$tkaMMS) - final_combined$UP



final_combined$dDTT <- log(final_combined$breDTT) - final_combined$DOWN
final_combined$dH2O2 <- log(final_combined$breH2O2) - final_combined$DOWN
final_combined$dSTARVE <- log(final_combined$breSTARVE) - final_combined$DOWN

final_combined$dHU <- log(final_combined$choHU) - final_combined$DOWN
final_combined$dRAP <- log(final_combined$choRAP) - final_combined$DOWN

final_combined$dSP <- log(final_combined$davSP) - final_combined$DOWN

final_combined$ddMMS <- log(final_combined$denMMS) - final_combined$DOWN

final_combined$dmMMS <- log(final_combined$mazMMS) - final_combined$DOWN

final_combined$dlMMS <- log(final_combined$leeMMS) - final_combined$DOWN

final_combined$dtHU <- log(final_combined$tkaHU) - final_combined$DOWN
final_combined$dtMMS <- log(final_combined$tkaMMS) - final_combined$DOWN


UP_mat <- final_combined[,c(25:35)]
rownames(UP_mat) <- final_combined$ORF
UP_mat[UP_mat < 0] = 0
UP_mat[UP_mat > 0] = 1

DOWN_mat <- final_combined[,c(36:46)]
rownames(DOWN_mat) <- final_combined$ORF
DOWN_mat[DOWN_mat > 0] = 0
DOWN_mat[DOWN_mat < 0] = 1

all_mat <- data.matrix(final_combined[,10:20])
rownames(all_mat) <- final_combined$ORF

UP <- UP_mat*log(all_mat)
DOWN <- DOWN_mat*log(all_mat)

allchanges <- (UP_mat + DOWN_mat)*log(all_mat)


# write.csv(UP_mat, file = "up.csv")
# write.csv(DOWN_mat, file = "down.csv")


UP_df <- data.frame(ORF = rownames(UP), UP)
DOWN_df <- data.frame(ORF = rownames(DOWN), DOWN)
ALL_df <- data.frame(ORF = rownames(all_mat), all_mat)

final_combined_sub <- subset(final_combined, 
                             select = c("ORF", "Mean", "UP", "DOWN"))
final_combined_sub <- final_combined_sub[order(final_combined_sub$Mean),]
final_combined_sub$ord <- as.numeric(seq(1:nrow(final_combined_sub)))

merge_final <- merge(final_combined_sub, UP_df, all.x = T, all.y = T)
merge_final <- merge(merge_final, DOWN_df, all.x = T, all.y = T)

UP_sub <- subset(merge_final, select = c("ORF", "ord", "uDTT", "uH2O2",
                                         "uSTARVE", "uHU", "uRAP", "uSP",
                                         "udMMS", "umMMS", "ulMMS", "utHU", "utMMS"))
UP_melt <- melt(UP_sub, id = c("ORF", "ord"))

DOWN_sub <- subset(merge_final, select = c("ORF", "ord", "dDTT", "dH2O2",
                                         "dSTARVE", "dHU", "dRAP", "dSP",
                                         "ddMMS", "dmMMS", "dlMMS", "dtHU", "dtMMS"))
DOWN_melt <- melt(DOWN_sub, id = c("ORF", "ord"))

ALL_sub <- subset(merge_final, select = c("ORF", "ord", "UP", "DOWN"))



p1 <- ggplot(ALL_sub, aes(x = ord))
p1 + geom_errorbar(data = ALL_sub, aes(ymin = UP, ymax = DOWN),
                   colour = "grey") +
     geom_point(data = UP_melt, aes(y = value), colour = "#e2725b", size = 1) +
     geom_point(data = DOWN_melt, aes(y = value), colour = "#8787ff", size = 1) +
     scale_x_continuous(limits = c(1550,1600))


# Find out which changes also have a two fold change, at least

ratio_mat <- Ratio[,3:ncol(Ratio)]
rownames(ratio_mat) <- Ratio$ORF

FC_change <- (UP_mat + DOWN_mat)*ratio_mat

FC_up <- UP_mat*ratio_mat
FC_down <- DOWN_mat*ratio_mat

FC_change[FC_change < 0.01] = NA
FC_change[FC_change > 2] = 100
FC_change[FC_change < 0.5] = 100
FC_change[FC_change < 100] = NA