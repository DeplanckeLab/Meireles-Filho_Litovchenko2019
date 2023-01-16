# FILE DESCRIPTION: 10_PhysioTimePredMTTonDGRP.R ------------------------------
#
# DESCRIPTION : In the script 10_PhysioTimePredModelsEval.R we identified that
#               the best method to predict physiological time is LD_DD_05,
#               molecular time table. So in this script we're going to predict
#               our DGRP samples.
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
#  CREATED:  
#  REVISION: 

setwd('~/Desktop/AroundTheClock_Aug_Oct2017/')
source('0_common_functions.R')
library(neuralnet)

# INPUTS ----------------------------------------------------------------------
set.seed(123)

INFO <- readRDS('Rds/sampleInfo.Rds')
INFO <- as.data.table(INFO)
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'
cyclingGenes <- readRDS('Rds/cyclGenesList.Rds')

# bloomington to ral 
bloomRal <- 'inputs_for_Rscripts/wolb_and_inv_DGRP2.txt'
bloomRal <- fread(bloomRal, header = T)
bloomRal[, line := gsub('DGRP-', 'line_', line)]
bloomRal[, line := gsub('line_0', 'line_', line)]
setnames(bloomRal, 'bloom', 'ID')
bloomRal[, ID := as.character(ID)]
setkey(bloomRal, ID)

# counts: in this case, I put dgrps and w- into one table for normalization and
# batch correction, so there wouldn't be batch effect bacause of it
# prepare info table
whiteDGRP <- readRDS(paste0(saveRDSdir, 'whiteDGRP.Rds'))
# get only white- and profiled DGRPs
whiteDgrpProf <- lapply(whiteDGRP, function(x) x[, !grepl('_n_', colnames(x))])
# according to Harbison et al 2019 (Behavior Genetics), line 25205 has a period
# of 18h and 20h in males and females respectfully, if estimated with MESA. So
# I exclude it from the training set
whiteDgrpProf <- lapply(whiteDgrpProf, 
                        function(x) x[, !grepl('_25205$', colnames(x))])
# many methods require matrix with samples in rows and genes in columns, here 
# you go 
whiteDgrpProf <- lapply(whiteDgrpProf, function(x) as.data.frame(t(x)))

# counts for DGRPs only to predict physiological time. I also add 25205 here
# for the reasons described above
dgrpOneTimePoint <- lapply(whiteDGRP, 
                           function(x) x[, grepl('_n_|_25205$', colnames(x))])
# many methods require matrix with samples in rows and genes in columns, here 
# you go 
dgrpOneTimePoint <- lapply(dgrpOneTimePoint, function(x) as.data.frame(t(x)))

# list all tissues
allTissues <- sort(unique(INFO$Tissue))
tissueColor <- c("#AB505E", "#05839C", "#0D6759", "#626499")
names(tissueColor) <- allTissues

# parameters for calculating time indicating genes
TIGname <- 'TIG_LDDD_05'
outliersName <- 'Without'
corCut <- 0.5
stDevCut <- 0.2

# Calculate time indicating genes ---------------------------------------------
# time indicating genes for all tissues
TIGs <- list()
# molecular peak time
MPTs <- list()
# for every tissue separetly 
for (tissue in allTissues) {
  message(paste('Started', tissue, 'at', Sys.time()))
  
  INFOTissue <- INFO[Tissue == tissue]
  setkey(INFOTissue, rightGT_name)
  
  # expression in this tissue for w- and profiled dgrps
  whiteDgrpProfTiss <- whiteDgrpProf[[tissue]]
  whiteDgrpProfTissTime <- INFOTissue[rownames(whiteDgrpProfTiss)]$Time
  
  # order cycling genes
  cyclingGenesTiss <- cyclingGenes[[tissue]]
  cyclingGenesTiss <- cyclingGenesTiss[cyclingGenesTiss$ADJ.P < 0.05, ]
  cyclingGenesTiss <- cyclingGenesTiss[order(-cyclingGenesTiss$AMP, 
                                             cyclingGenesTiss$ADJ.P), ]
  
  ## LD
  #LD <- whiteDgrpProfTissTime[whiteDgrpProfTissTime <= 24]
  #whiteDgrpProfTissLD <- whiteDgrpProfTiss[whiteDgrpProfTissTime <= 24, ]
  #corCosLD <- corCosBundleMolPeakTime(LD, whiteDgrpProfTissLD)
  #stDevLD <- apply(whiteDgrpProfTissLD, 2, function(x) sd(x) / mean(x))
  ## DD
  #DD <- whiteDgrpProfTissTime[whiteDgrpProfTissTime > 24]
  #whiteDgrpProfTissDD <- whiteDgrpProfTiss[whiteDgrpProfTissTime > 24, ]
  #corCosDD <- corCosBundleMolPeakTime(DD, whiteDgrpProfTissDD)
  #stDevDD <- apply(whiteDgrpProfTissDD, 2, function(x) sd(x) / mean(x))
  #TIG_LD_DD_05 <- colnames(whiteDgrpProfTiss)[which(corCosLD$maxCor >= corCut &
  #                                                  corCosDD$maxCor >= corCut &
  #                                                  stDevLD >= stDevCut &
  #                                                  stDevDD >= stDevCut)]
  ## add to the final list
  #TIGs[[length(TIGs) + 1]] <- TIG_LD_DD_05
  
  ## molecular peak times for time indicating genes
  #MPTs1Tiss <- sapply(TIG_LD_DD_05, function(x) mean(corCosLD[x, 'MPT'],
  #                                                   corCosDD[x, 'MPT']))
  #MPTs[[length(MPTs) + 1]] <- MPTs1Tiss
  
  corCosLDDD <- corCosBundleMolPeakTime(whiteDgrpProfTissTime,
                                        whiteDgrpProfTiss)
  stDevLDDD <- apply(whiteDgrpProfTiss, 2, function(x) sd(x) / mean(x))
  TIG_LDDD_05 <- colnames(whiteDgrpProfTiss)[which(corCosLDDD$maxCor >= corCut &
                                                   stDevLDDD >= stDevCut)]
  
  # add to the final list
  TIGs[[length(TIGs) + 1]] <- TIG_LDDD_05
  
  ## molecular peak times for time indicating genes
  MPTs1Tiss <- sapply(TIG_LDDD_05, function(x) corCosLDDD[x, 'MPT'])
  MPTs[[length(MPTs) + 1]] <- MPTs1Tiss
}
names(TIGs) <- allTissues
names(MPTs) <- allTissues

# Predict physiological time for DGRP samples ---------------------------------
dgrpOneTP_PredTimeDF <- data.frame()
# for every tissue separetly 
for (tissue in allTissues) {
  INFOTissue <- INFO[Tissue == tissue]
  setkey(INFOTissue, rightGT_name)
  
  # expression in this tissue for w- and profiled dgrps
  dgrpOneTimePointTiss <- dgrpOneTimePoint[[tissue]]
  dgrpOneTimePointTissTime <- INFOTissue[rownames(dgrpOneTimePointTiss)]$Time
  
  # TIGs and MPTs for this tissue
  tigs <- TIGs[[tissue]]
  # tigs %in% colnames(whiteDgrpProfnoOut) because cycling gene detection 
  # was on done only on w-, so then added DGRPs, certain genes might be
  # too low at expression
  tigs <- tigs[tigs %in% colnames(dgrpOneTimePointTiss)]
  mpts <- MPTs[[tissue]][tigs]
  
  # calculate physiological time
  physioTime <- sapply(rownames(dgrpOneTimePointTiss), 
                       function(x) calcBodyTime(mpts, 
                                                unlist(dgrpOneTimePointTiss[x,
                                                                            tigs])))
  # add to final data frame
  dgrpOneTP_PredTimeDF <- rbind(dgrpOneTP_PredTimeDF,
                                cbind(rownames(dgrpOneTimePointTiss), tissue,
                                     'MolTT', TIGname, outliersName, 
                                     dgrpOneTimePointTissTime,
                                     physioTime))
}
colnames(dgrpOneTP_PredTimeDF) <- c('Sample', 'Tissue', 'Method', 'TIGname', 
                                    'outliersName', "timeObs", "timePred")
dgrpOneTP_PredTimeDF <- as.data.table(dgrpOneTP_PredTimeDF)
dgrpOneTP_PredTimeDF[, Sample := as.character(Sample)]
dgrpOneTP_PredTimeDF[, Tissue := as.character(Tissue)]
dgrpOneTP_PredTimeDF[, timeObs := as.double(as.character(timeObs))]
dgrpOneTP_PredTimeDF[, timePred := as.double(as.character(timePred))]
timeErr <- sapply(1:nrow(dgrpOneTP_PredTimeDF),
                  function(x) diffInTime(dgrpOneTP_PredTimeDF$timeObs[x],
                                         dgrpOneTP_PredTimeDF$timePred[x]))
dgrpOneTP_PredTimeDF[, timeError := timeErr]
setkey(INFO, rightGT_name)
dgrpOneTP_PredTimeDF[, Time := INFO[Sample]$Time]
dgrpOneTP_PredTimeDF[, DGRP := INFO[Sample]$right_GT]
dgrpOneTP_PredTimeDF[, RAL := gsub('line_', '', bloomRal[DGRP]$line)]
dgrpOneTP_PredTimeDF[, RAL := factor(RAL, 
                                     levels = sort(as.integer(unique(RAL))))]
saveRDS(dgrpOneTP_PredTimeDF, 
        paste0(saveRDSdir, 'physioPredTimeDF_MTT_withoutOutliers_CV_no25205.Rds'))