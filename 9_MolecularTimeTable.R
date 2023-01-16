# FILE DESCRIPTION: 9_MolecularTimeTable.R -----------------------------------
#
# DESCRIPTION : implementation of Ueda et al 2004, 
# Molecular-timetable methods for detection of body time and rhythm disorders 
# from single-time-point genome-wide expression profiles
#
# USAGE: 
#
#  OPTIONS:  none
#  REQUIREMENTS:  data.table, lubridate, stringr
#  BUGS: --
#  NOTES:  ---
#  AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
#  COMPANY:  EPFL, Lausanne, Switzerland
#  VERSION:  1
#  CREATED:  10.04.2018
#  REVISION: 10.04.2018

setwd('~/Desktop/BitBucket/AroundTheClock_Aug_Oct2017/')
setwd('~/Desktop/AroundTheClock_Aug_Oct2017/')
source('0_common_functions.R')

# INPUTS ----------------------------------------------------------------------
set.seed(123)

INFO <- readRDS('Rds/sampleInfo.Rds')
INFO <- as.data.table(INFO)
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'
cyclingGenes <- readRDS('Rds/cyclGenesList.Rds')
# list all tissues
allTissues <- unique(INFO$Tissue)
tissueColor <- c('brown3', 'darkgoldenrod1', 'darkgreen', 'darkorchid3')
names(tissueColor) <- allTissues

# counts: in this case, I put dgrps and w- into one table for normalization and
# batch correction, so there wouldn't be batch effect bacause of it
# prepare info table
whiteDGRP <- readRDS(paste0(saveRDSdir, 'whiteDGRP.Rds'))
# get only white- and profiled DGRPs
whiteDgrpProf <- lapply(whiteDGRP, function(x) x[, !grepl('_n_', colnames(x))])
# many methods require matrix with samples in rows and genes in columns, here 
# you go 
whiteDgrpProfMatr <- lapply(whiteDgrpProf, function(x) as.data.frame(t(x)))

tissueOutliers <- list(FB = c('FB43_Oct_2', 'FB30_Oct_2',
                              'BRB2_FB_ZT11_29655'),
                       GUT = c('BRB5_Gut_ZT11_29655', 'BRB5_Gut_ZT23_29655',
                               'BRB5_Gut_ZT11_25205', 'BRB4_Gut_ZT9_28211', 
                               'GUT31_Jan_6', 'BRB4_Gut_ZT21_29655'),
                       BRAIN = c('BRB8_Brain_ZT1_29655', 'BRAIN13_Oct_3',
                                 'BRAIN37_Oct_3', 'BRAIN1_Oct_3', 
                                 'BRAIN28_Oct_1', 'BRB7_Brain_ZT21_28211'),
                       MT = c())

# DGRP samples: "normal", which biological time we're interested in
dgrpNorm <- lapply(whiteDGRP, function(x) x[, grepl('_n_', colnames(x))])
# many methods require matrix with samples in rows and genes in columns, here 
# you go 
dgrpNormMatr <- lapply(dgrpNorm, function(x) as.data.frame(t(x)))

# WHITE & DGRP: READ-IN COUNTS, NORMALIZE, BATCH CORRECT ----------------------
tissue <- 'GUT'
outliers <- tissueOutliers[[tissue]]
colorForTissue <- tissueColor[tissue]
INFOTissue <- INFO[Tissue == tissue]
setkey(INFOTissue, rightGT_name)
cyclingGenesTiss <- cyclingGenes[[tissue]]
cyclingGenesTiss <- cyclingGenesTiss[cyclingGenesTiss$ADJ.P < 0.05, ]

whiteDgrpProf <- whiteDgrpProfMatr[[tissue]]
whiteDgrpProfnoOut <- whiteDgrpProf[!rownames(whiteDgrpProf) %in% outliers, ]
dgrpNorm <- dgrpNormMatr[[tissue]]
# time for both sets
whiteDgrpProfnoOutTime <- INFOTissue[rownames(whiteDgrpProfnoOut)]$Time 
dgrpNormTime <- INFOTissue[rownames(dgrpNorm)]$Time

# Read-in Ueda genes ----------------------------------------------------------
uedaTimeIndGenes <- fread('inputs_for_Rscripts/Ueda_Table3.csv', header = T)
length(intersect(rownames(cyclingGenesTiss), uedaTimeIndGenes$Symbol))

# Time-indicating genes & their molecular peak time----------------------------
# To select time-indicating genes whose expression exhibits circadian 
# rhythmicity with high amplitude, the expression profile of each gene was 
# analyzed through two filters, one for circadian rhythmicity and the other for
# high amplitude.In paper: "we first calculated the correlation over time 
# between 12-point time courses under LD (or DD) conditions and cosine curves 
# of defined periods and phases generate bundles of cosine waves". They had 4
# days of observation: 2 in LD and 2 in DD. Our set up differs from that, 
# however, I'll try to apply their strategy

# LD
LD <- whiteDgrpProfnoOutTime[whiteDgrpProfnoOutTime <= 24]
whiteDgrpProfnoOutLD <- whiteDgrpProfnoOut[whiteDgrpProfnoOutTime <= 24, ]
corCosLD <- corCosBundleMolPeakTime(LD, whiteDgrpProfnoOutLD)
stDevLD <- apply(whiteDgrpProfnoOutLD, 2, function(x) sd(x)/mean(x))

# DD
DD <- whiteDgrpProfnoOutTime[whiteDgrpProfnoOutTime > 24]
whiteDgrpProfnoOutDD <- whiteDgrpProfnoOut[whiteDgrpProfnoOutTime > 24, ]
corCosDD <- corCosBundleMolPeakTime(DD, whiteDgrpProfnoOutDD)
stDevDD <- apply(whiteDgrpProfnoOutDD, 2, function(x) sd(x)/mean(x))

# together: LD + DD
corCosLDDD <- corCosBundleMolPeakTime(whiteDgrpProfnoOutTime,
                                      whiteDgrpProfnoOut)
stDevLDDD <- apply(whiteDgrpProfnoOut, 2, function(x) sd(x)/mean(x))

# time indicating genes. In the paper: 1) best correlation values in LD and 
# DD conditions were both above the cutoff correlation value of 0.8, 
# 2) coefficients of variation in LD and DD conditions were both above the
# cutoff value of 0.15 (high amplitude). But for us it seems to be too strict,
# because it gives only 2-3 genes. Since we're going to fit cos curves into
# molecular peak time (MPT) of TIGs afterwards, we need as many genes as 
# possible
corCutOff <- 0.5
stDevCutOff <- 0.15
TIG <- colnames(whiteDgrpProfnoOut)[which(corCosLD$maxCor >= corCutOff &
                                          corCosDD$maxCor >= corCutOff &
                                          stDevLD >= stDevCutOff &
                                          stDevDD >= stDevCutOff)]
TIG <- colnames(whiteDgrpProfnoOut)[which(corCosLDDD$maxCor >= corCutOff &
                                          stDevLDDD >= stDevCutOff)]
#TIG <- rownames(cyclingGenesTiss[1:50, ])

# Plot time profiles with expression of time indicating genes -----------------
# create DF with the expression values of TIGs in every sample + TIG peak time
TIGexpr <- t(whiteDgrpProfnoOut[, TIG])
TIGexpr <- melt(TIGexpr)
colnames(TIGexpr) <- c('TIG', 'Sample', 'Expr')
TIGexpr$TIG <- as.character(TIGexpr$TIG)
TIGexpr$Sample <- as.character(TIGexpr$Sample)
TIGexpr$Time <- factor(INFOTissue[TIGexpr$Sample, Time] %% 24,
                       levels = sort(unique(INFOTissue[TIGexpr$Sample, Time] 
                                            %% 24)))
TIGexpr$MPT <- corCosLDDD[TIGexpr$TIG, 'MPT']
TIGexpr <- as.data.table(TIGexpr)

# plot expression profiles of TIG for every time point where we have reference
# samples
MPTallSamp <- ggplot(TIGexpr, 
                     aes(x = MPT, y = Expr, 
                         color = factor(MPT, levels = sort(unique(MPT))))) +
              geom_point(shape = 20, cex = 2) + mashaGgplot2Theme +
              theme(legend.position = "none")
facet_multiple(plot = MPTallSamp, facets = 'Time', ncol = 4, nrow = 3)

# plot expression profiles of TIG for reference sample
MPTbySamp <- ggplot(TIGexpr, 
                    aes(x = MPT, y = Expr, 
                        color = factor(MPT, levels = sort(unique(MPT))))) +
             geom_point(shape = 20, cex = 2) + mashaGgplot2Theme +
             theme(legend.position = "none")
facet_multiple(plot = MPTbySamp, facets = 'Sample', ncol = 4, nrow = 3)

# Get physiological time for white- and profiled DGRPs ------------------------
physioTime <- sapply(rownames(whiteDgrpProfnoOut), 
                     function(x) calcBodyTime(corCosLDDD[TIG, ]$MPT, 
                                              unlist(whiteDgrpProfnoOut[x,
                                                                        TIG])))
molTTrealPredTimeDF <- cbind(rownames(whiteDgrpProfnoOut), tissue,
                             'MolTT', 
                             INFOTissue[rownames(whiteDgrpProfnoOut)]$Time,
                             physioTime)

# Get physiological time for profiled once DGRPs ------------------------------
# best parameters
corCutOff <- 0.5
stDevCutOff <- 0.15
timeDiffCutOff <- 4
for (tissue in allTissues) {
  outliers <- "" #tissueOutliers[[tissue]]
  colorForTissue <- tissueColor[tissue]
  INFOTissue <- INFO[Tissue == tissue]
  setkey(INFOTissue, rightGT_name)
  
  whiteDgrpProf <- whiteDgrpProfMatr[[tissue]]
  whiteDgrpProfnoOut <- whiteDgrpProf[!rownames(whiteDgrpProf) %in% outliers, ]
  dgrpNorm <- dgrpNormMatr[[tissue]]
  # time for both sets
  whiteDgrpProfnoOutTime <- INFOTissue[rownames(whiteDgrpProfnoOut)]$Time 
  dgrpNormTime <- INFOTissue[rownames(dgrpNorm)]$Time
  
  corCosLDDD <- corCosBundleMolPeakTime(whiteDgrpProfnoOutTime,
                                        whiteDgrpProfnoOut)
  stDevLDDD <- apply(whiteDgrpProfnoOut, 2, function(x) sd(x) / mean(x))
  
  TIG <- colnames(whiteDgrpProfnoOut)[which(corCosLDDD$maxCor >= 0.5 &
                                            stDevLDDD >= 0.15)]
  
  physioTime <- sapply(rownames(dgrpNorm), 
                       function(x) calcBodyTime(corCosLDDD[TIG, ]$MPT, 
                                                unlist(dgrpNorm[x, TIG])))
  molTTrealPredTimeDF <- data.frame(Sample = rownames(dgrpNorm), 
                                    timeObs = INFOTissue[rownames(dgrpNorm)]$Time,
                                    timePred = physioTime)
  molTTrealPredTimeDF$timeError <- apply(molTTrealPredTimeDF, 1, 
                                         function(x) diffInTime(as.numeric(x[2]), 
                                                                as.numeric(x[3])))
  
  plot(molTTrealPredTimeDF[, 2:3], pch = 20, bty = 'n',
       col = ifelse(abs(molTTrealPredTimeDF$timeError) > timeDiffCutOff, 
                    tissueColor[tissue], 'black'), cex = 2, 
       xlim = c(0, 25), ylim = c(0, 25), cex.lab = 2, cex.axis = 2,
       xlab = 'Time of harversting, h', ylab = 'Physiological time, h',
       main = paste0('Physiological and harversting time, ', tissue))
  abline(1, 1, lty = 2)
  text(molTTrealPredTimeDF[abs(molTTrealPredTimeDF$timeError) > 
                               timeDiffCutOff, 2:3],
       labels = gsub('.*_n_', '', 
                     molTTrealPredTimeDF$Sample[abs(molTTrealPredTimeDF$timeError) > 
                                         timeDiffCutOff]),
       pos = 3)
  grid()
  
  saveRDS(molTTrealPredTimeDF, paste0(saveRDSdir, tissue, 
                                      '_physioPredTimeDF_MTT_withOutliers_CV.Rds'))
}