# FILE DESCRIPTION: 26_PhysioTimePredModelsEval.R -----------------------------
#
# DESCRIPTION :
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
#  CREATED:  08.04.2018
#  REVISION: 08.04.2018

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

# list all tissues
allTissues <- sort(unique(INFO$Tissue))
tissueColor <- c("#AB505E", "#05839C", "#0D6759", "#626499")
names(tissueColor) <- allTissues

tissueOutliers <- list(With = list(FB = c('FB43_Oct_2', 'FB30_Oct_2', 
                                          'BRB2_FB_ZT11_29655'),
                                   GUT = c('BRB5_Gut_ZT11_29655', 
                                           'BRB5_Gut_ZT23_29655', 
                                           'BRB5_Gut_ZT11_25205', 
                                           'BRB4_Gut_ZT9_28211', 
                                           'GUT31_Jan_6', 
                                           'BRB4_Gut_ZT21_29655'),
                                   BRAIN = c('BRB8_Brain_ZT1_29655', 
                                             'BRAIN13_Oct_3', 
                                             'BRAIN37_Oct_3',
                                             'BRAIN1_Oct_3', 
                                             'BRAIN28_Oct_1', 
                                             'BRB7_Brain_ZT21_28211'),
                                   MT = c()),
                       Without = list(FB = c(), GUT = c(), BRAIN = c(), 
                                      MT = c()))
# LASSO, EVALUATE THE MODEL ---------------------------------------------------
lassoRealPredTimeDF <- c() # DF will hold real and predicted time
for (outliersName in names(tissueOutliers)) {
  message(paste('Started', outliersName, 'at', Sys.time()))
  for (tissue in allTissues) {
    message(paste('\tStarted', tissue, 'at', Sys.time()))
    outliers <- tissueOutliers[[outliersName]][[tissue]]
    colorForTissue <- tissueColor[tissue]
    INFOTissue <- INFO[Tissue == tissue]
    setkey(INFOTissue, rightGT_name)
    
    whiteDgrpProfnoOut <- whiteDgrpProf[[tissue]]
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[!rownames(whiteDgrpProfnoOut) %in%
                                               outliers, ]
    
    quadCoeffs <- list() # put coefficients for the first model (quadrant) here
    lassoCoeffs_1_12 <- list()
    lassoCoeffs_12_24 <- list()
    for (i in 1:nrow(whiteDgrpProfnoOut)) {
      message(paste0('LASSO', tissue, ': ', i))
      manualSel <- colnames(whiteDgrpProfnoOut)
      # divide into train and test dataset
      trainTestSets <- divideIntoTrainAndTest(whiteDgrpProfnoOut, INFOTissue,
                                              leaveOneOut = i, genesOI = manualSel)
      refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time
      refTrainTimeFact <- as.factor(ifelse(refTrainTime %% 24 > 12, 1, 0))
      
      # 1st STEP: TRAIN MODEL FOR DETERMINE, TO WHICH PART OF THE DAY SAMPLE
      # BELONGS
      lassoFit_DayTime <- lassoOnTrainSet(trainTestSets$train, refTrainTimeFact, 
                                          F, family = 'binomial')
      dayTimeCoeffs <- coef(lassoFit_DayTime)[, 1][coef(lassoFit_DayTime)[, 1] != 0]
      quadCoeffs[[length(quadCoeffs) + 1]] <- dayTimeCoeffs
      
      # 2nd STEP: TRAN MODELS FOR BEFORE AND AFTER 12
      # create day time specific train sets
      trainSet_1_12 <- trainTestSets$train[refTrainTime  %% 24 <= 12, ]
      refTrainTime_1_12 <- INFOTissue[rownames(trainSet_1_12)]$Time %% 24
      trainSet_12_24 <- trainTestSets$train[refTrainTime  %% 24 >= 12, ]
      refTrainTime_12_24 <- INFOTissue[rownames(trainSet_12_24)]$Time %% 24
      # fit 
      lassoFit_1_12 <- lassoOnTrainSet(trainSet_1_12, refTrainTime_1_12, F)
      coeffs_1_12 <- coef(lassoFit_1_12)[, 1][coef(lassoFit_1_12)[, 1] != 0]
      lassoCoeffs_1_12[[length(lassoCoeffs_1_12) + 1]] <- coeffs_1_12
      lassoFit_12_24 <- lassoOnTrainSet(trainSet_12_24, refTrainTime_12_24, F)
      coeffs_12_24 <- coef(lassoFit_12_24)[, 1][coef(lassoFit_12_24)[, 1] != 0]
      lassoCoeffs_12_24[[length(lassoCoeffs_12_24) + 1]] <- coeffs_12_24
      
      # 3rd STEP: predict
      # run prediction on test
      testPredict <- predict(lassoFit_DayTime, newx = trainTestSets$tes, 
                             s = "lambda.min", family = 'binomial')
      predTestTime <- as.factor(ifelse(testPredict[, 1] > 0, 1, 0))
      if (predTestTime == 0) {
        testPredTime <- predict(lassoFit_1_12, newx = trainTestSets$test, 
                                s = "lambda.min")
      } else {
        testPredTime <- predict(lassoFit_12_24, newx = trainTestSets$test, 
                                s = "lambda.min")
      }
      # record predicted and harversted time
      refTestTime <- INFOTissue[rownames(trainTestSets$test)]$Time %% 24
      lassoRealPredTimeDF <- rbind(lassoRealPredTimeDF,  
                                   c(rownames(whiteDgrpProfnoOut)[i],
                                     tissue, 'LASSO', '', outliersName,
                                     refTestTime, testPredTime[, 1]))
    }
  }
}
saveRDS(lassoRealPredTimeDF, 
        paste0(saveRDSdir, 'LASSO_realPredTimeDF_no25205.Rds'))

# ZEITZEIGER, EVALUATE THE MODEL ----------------------------------------------
optSumabsv <- 1.5
optSPC <- 2

zzRealPredTimeDF <- c()
for (outliersName in names(tissueOutliers)) {
  message(paste('Started', outliersName, 'at', Sys.time()))
  for (tissue in allTissues) {
    message(paste('\tStarted', tissue, 'at', Sys.time()))
    outliers <- tissueOutliers[[outliersName]][[tissue]]
    colorForTissue <- tissueColor[tissue]
    INFOTissue <- INFO[Tissue == tissue]
    setkey(INFOTissue, rightGT_name)
    
    whiteDgrpProfnoOut <- whiteDgrpProf[[tissue]]
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[!rownames(whiteDgrpProfnoOut) %in%
                                              outliers, ]
    
    for (i in 1:nrow(whiteDgrpProfnoOut)) {
      message(paste0('ZEITZEIGER ', tissue, ': ', i))
      manualSel <- colnames(whiteDgrpProfnoOut)
      # divide into train and test dataset
      trainTestSets <- divideIntoTrainAndTest(whiteDgrpProfnoOut, INFOTissue,
                                              leaveOneOut = i, genesOI = manualSel)
      refTrainTime <- (INFOTissue[rownames(trainTestSets$train)]$Time %% 24) / 24
      
      # fit a periodic smoothing spline to the behavior of each feature as a function
      # of time
      fitWhiteDgrpProf <- zeitzeigerFit(trainTestSets$train, refTrainTime)
      # fits to calculate sparse principal components (SPCs) for how the features
      # change over time
      spcWhiteDgrpProf <- zeitzeigerSpc(fitWhiteDgrpProf$xFitMean, 
                                        fitWhiteDgrpProf$xFitResid, 
                                        sumabsv = optSumabsv)
      # uses the training data and the SPCs to predict the corresponding time for
      # each test observation
      predFitTest <- zeitzeigerPredict(as.matrix(trainTestSets$train),
                                       refTrainTime,
                                       as.matrix(trainTestSets$test),
                                       spcWhiteDgrpProf, nSpc = optSPC)
      # record predicted and harversted time
      refTestTime <- INFOTissue[rownames(trainTestSets$test)]$Time %% 24
      zzRealPredTimeDF <- rbind(zzRealPredTimeDF,
                                c(rownames(whiteDgrpProfnoOut)[i], tissue, 
                                  'ZeitZeiger', '', outliersName,
                                  refTestTime, 24 * predFitTest$timePred))
    }
  }
}
saveRDS(zzRealPredTimeDF, 
        paste0(saveRDSdir, 'ZeitZeiger_realPredTimeDF_no25205.Rds'))

# MOLECULAR TIME TABLE, EVALUATE THE MODEL ------------------------------------
molTTrealPredTimeDF <- c()
MPTlistTotal <- list()
for (outliersName in names(tissueOutliers)) {
  message(paste('Started', outliersName, 'at', Sys.time()))
  for (tissue in allTissues) {
    message(paste('\tStarted', tissue, 'at', Sys.time()))
    outliers <- tissueOutliers[[outliersName]][[tissue]]
    colorForTissue <- tissueColor[tissue]
    INFOTissue <- INFO[Tissue == tissue]
    setkey(INFOTissue, rightGT_name)
    # order cycling genes
    cyclingGenesTiss <- cyclingGenes[[tissue]]
    cyclingGenesTiss <- cyclingGenesTiss[cyclingGenesTiss$ADJ.P < 0.05, ]
    cyclingGenesTiss <- cyclingGenesTiss[order(-cyclingGenesTiss$AMP, 
                                               cyclingGenesTiss$ADJ.P), ]
    whiteDgrpProfnoOut <- whiteDgrpProf[[tissue]]
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[!rownames(whiteDgrpProfnoOut) %in%
                                               outliers, ]
    whiteDgrpProfnoOutTime <- INFOTissue[rownames(whiteDgrpProfnoOut)]$Time 
    
    # LD
    LD <- whiteDgrpProfnoOutTime[whiteDgrpProfnoOutTime <= 24]
    whiteDgrpProfnoOutLD <- whiteDgrpProfnoOut[whiteDgrpProfnoOutTime <= 24, ]
    corCosLD <- corCosBundleMolPeakTime(LD, whiteDgrpProfnoOutLD)
    stDevLD <- apply(whiteDgrpProfnoOutLD, 2, function(x) sd(x) / mean(x))
    
    # DD
    DD <- whiteDgrpProfnoOutTime[whiteDgrpProfnoOutTime > 24]
    whiteDgrpProfnoOutDD <- whiteDgrpProfnoOut[whiteDgrpProfnoOutTime > 24, ]
    corCosDD <- corCosBundleMolPeakTime(DD, whiteDgrpProfnoOutDD)
    stDevDD <- apply(whiteDgrpProfnoOutDD, 2, function(x) sd(x) / mean(x))
    
    # together: LD + DD
    corCosLDDD <- corCosBundleMolPeakTime(whiteDgrpProfnoOutTime,
                                          whiteDgrpProfnoOut)
    stDevLDDD <- apply(whiteDgrpProfnoOut, 2, function(x) sd(x) / mean(x))
    
    # for this method, I'll try several sets of genes
    # on LD and DD separetly, cut off on correlation 0.5 and 0.8
    TIG_LD_DD_08 <- colnames(whiteDgrpProfnoOut)[which(corCosLD$maxCor >= 0.8 &
                                                         corCosDD$maxCor >= 0.8 &
                                                         stDevLD >= 0.20 &
                                                         stDevDD >= 0.20)]
    TIG_LD_DD_05 <- colnames(whiteDgrpProfnoOut)[which(corCosLD$maxCor >= 0.5 &
                                                         corCosDD$maxCor >= 0.5 &
                                                         stDevLD >= 0.20 &
                                                         stDevDD >= 0.20)]
    # on LD and DD together
    TIG_LDDD_08 <- colnames(whiteDgrpProfnoOut)[which(corCosLDDD$maxCor >= 0.8 &
                                                        stDevLDDD >= 0.20)]
    TIG_LDDD_05 <- colnames(whiteDgrpProfnoOut)[which(corCosLDDD$maxCor >= 0.5 &
                                                        stDevLDDD >= 0.20)]
    # on top 25/50/75/100/all cycling genes
    TIG_cycl25 <- rownames(cyclingGenesTiss[1:25, ])
    TIG_cycl50 <- rownames(cyclingGenesTiss[1:50, ])
    TIG_cycl75 <- rownames(cyclingGenesTiss[1:75, ])
    TIG_cycl100 <- rownames(cyclingGenesTiss[1:100, ])
    TIG_cycl <- rownames(cyclingGenesTiss)
    # On Ueda genes
    uedaTimeIndGenes <- fread('inputs_for_Rscripts/Ueda_Table3.csv', header = T)
    TIG_Ueda <- intersect(uedaTimeIndGenes$Symbol,
                          colnames(whiteDgrpProfnoOut))
    # On overlap between Ueda genes and our genes
    TIG_UedaUs <- intersect(rownames(cyclingGenesTiss),
                            uedaTimeIndGenes$Symbol)
    realPredTime <- c() # put real and predicted time here
    
    TIGlist <- list(TIG_LD_DD_08, TIG_LD_DD_05, TIG_LDDD_08, TIG_LDDD_05,
                    TIG_cycl25, TIG_cycl50, TIG_cycl75, TIG_cycl100, TIG_cycl,
                    TIG_Ueda, TIG_UedaUs)
    names(TIGlist) <- c("TIG_LD_DD_08", "TIG_LD_DD_05", "TIG_LDDD_08",
                        "TIG_LDDD_05", "TIG_cycl25", "TIG_cycl50",
                        "TIG_cycl75", "TIG_cycl100", "TIG_cycl",
                        "TIG_Ueda", "TIG_UedaUs")
    MPTlist <- list(TIG_LD_DD_08 = sapply(TIG_LD_DD_08, 
                                          function(x) mean(corCosLD[x, 'MPT'],
                                                           corCosDD[x, 'MPT'])),
                    TIG_LD_DD_05 = sapply(TIG_LD_DD_05, 
                                          function(x) mean(corCosLD[x, 'MPT'],
                                                           corCosDD[x, 'MPT'])))
    MPTlist <- c(MPTlist, lapply(TIGlist[-c(1, 2)], 
                                 function(x) corCosLDDD[x, ]$MPT))
    
    toAdd <- lapply(1:length(MPTlist), 
                    function(x) data.table(Tissue = tolower(tissue),
                                           TIGtype = names(TIGlist)[x],
                                  TIG = TIGlist[[x]],
                                  MPT = MPTlist[[x]],
                                  out = outliersName))
    toAdd <- do.call(rbind, toAdd)
    MPTlistTotal <- rbind(MPTlistTotal, toAdd)
    
    for (TIGname in names(TIGlist)) {
      tigs <- TIGlist[[TIGname]]
      mpts <- MPTlist[[TIGname]][tigs %in% colnames(whiteDgrpProfnoOut)]
      tigs <- tigs[tigs %in% colnames(whiteDgrpProfnoOut)]
      # tigs %in% colnames(whiteDgrpProfnoOut) because cycling gene detection 
      # was on done only on w-, so then added DGRPs, certain genes might be
      # too low at expression
      physioTime <- sapply(rownames(whiteDgrpProfnoOut), 
                           function(x) calcBodyTime(mpts, 
                                                    unlist(whiteDgrpProfnoOut[x,
                                                                              tigs])))
      molTTrealPredTimeDF <- rbind(molTTrealPredTimeDF,
                                   cbind(rownames(whiteDgrpProfnoOut), tissue,
                                         'MolTT', TIGname, outliersName, 
                                         INFOTissue[rownames(whiteDgrpProfnoOut)]$Time,
                                         physioTime))
    }
  }
}
saveRDS(MPTlistTotal, paste0(saveRDSdir, 'MolecularTimeTable_MPT_no25205.Rds'))
saveRDS(molTTrealPredTimeDF, 
        paste0(saveRDSdir, 'MolecularTimeTable_realPredTimeDF_cv_no25205.Rds'))

# Directional statistics ------------------------------------------------------
directStatsRealPredTimeDF <- c() # DF will hold real and predicted time
for (outliersName in names(tissueOutliers)) {
  message(paste('Started', outliersName, 'at', Sys.time()))
  for (tissue in allTissues) {
    message(paste('\tStarted', tissue, 'at', Sys.time()))
    outliers <- tissueOutliers[[outliersName]][[tissue]]
    colorForTissue <- tissueColor[tissue]
    INFOTissue <- INFO[Tissue == tissue]
    setkey(INFOTissue, rightGT_name)
    
    # get DGRP and proviled samples for tissue
    whiteDgrpProfnoOut <- whiteDgrpProf[[tissue]]
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[!rownames(whiteDgrpProfnoOut) %in%
                                               outliers, ]
    
    # leave one out approach
    for (i in 1:nrow(whiteDgrpProfnoOut)) {
      message(paste0('Directional statistics ', tissue, ': ', outliersName, ' ', 
                     i))
      # all genes here, because no selection is made
      manualSel <- colnames(whiteDgrpProfnoOut)
      # divide into train and test dataset
      trainTestSets <- divideIntoTrainAndTest(whiteDgrpProfnoOut, INFOTissue,
                                              leaveOneOut = i, 
                                              genesOI = manualSel)
      refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time
      testTime <- INFOTissue[rownames(trainTestSets$test)]$Time
      
      modelPredict <- fitDirectionalModel(trainTestSets$train, refTrainTime, 
                                          trainTestSets$test, testTime)
      modelPredict <- c(rownames(trainTestSets$test), tissue, 'DirectStat',
                        'DirectStat', outliersName, modelPredict)
      directStatsRealPredTimeDF <- rbind(directStatsRealPredTimeDF, 
                                         modelPredict)
    }
  }
}
saveRDS(directStatsRealPredTimeDF, 
        paste0(saveRDSdir, 'DirectStats_realPredTimeDF_no25205.Rds'))

# NEURAL NET ------------------------------------------------------------------
neuralNetPredTimeDF <- c()
for (outliersName in names(tissueOutliers)) {
  for (tissue in allTissues) {
    outliers <- tissueOutliers[[outliersName]][[tissue]]
    colorForTissue <- tissueColor[tissue]
    INFOTissue <- INFO[Tissue == tissue]
    setkey(INFOTissue, rightGT_name)
    
    cyclGenesTiss <- cyclingGenes[[tissue]]
    cyclGenesTiss <- cyclGenesTiss[cyclGenesTiss$ADJ.P < 0.05, ]
    cyclGenesTiss <- cyclGenesTiss[order(-cyclGenesTiss$AMP, 
                                         cyclGenesTiss$ADJ.P), ]
    
    whiteDgrpProfnoOut <- whiteDgrpProf[[tissue]]
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[, colnames(whiteDgrpProfnoOut) %in% 
                                               rownames(cyclGenesTiss)]
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[!rownames(whiteDgrpProfnoOut) %in%
                                               outliers, ]
    neurons <- list(n12_6_3 = c(12, 6, 3), n12_6 = c(12, 6), n12 = 12)
    for (neuronName in names(neurons)) {
      hiddenNeurons <- neurons[[neuronName]]
      for (i in 1:nrow(whiteDgrpProfnoOut)) {
        message(paste0('NEURAL NET ', tissue, ': ', i))
        # divide into train and test dataset
        trainTestSets <- divideIntoTrainAndTest(whiteDgrpProfnoOut, INFOTissue,
                                                leaveOneOut = i)
        refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time %% 24
        
        trainSet <- cbind(refTrainTime, trainTestSets$train)
        testSet <- trainTestSets$test
        colnames(trainSet)[1] <- 'timeOfSample'
        colnames(trainSet)[2:ncol(trainSet)] <- paste0('GENE_', 2:ncol(trainSet))
        colnames(testSet)[2:ncol(testSet)] <- paste0('GENE_', 2:ncol(testSet))
        
        predictorGenes <- colnames(trainSet)[-1]
        modelFormula <- as.formula(paste("timeOfSample ~", 
                                         paste(predictorGenes, collapse = " + ")))
        neuralNet <- NULL
        preditNN <- NULL
        while (is.null(neuralNet) | is.null(preditNN)) {
          tryCatch({neuralNet <- neuralnet(modelFormula, data = trainSet,
                                           hidden = hiddenNeurons,
                                           linear.output = T)
          preditNN <- compute(neuralNet, testSet)}, 
          warning = function(w) {print('WARNING! trying again')
            neuralNet = NULL},
          error = function(e) {print('ERROR! trying again')
            neuralNet = NULL})
        }
        neuralNetPredTimeDF <- rbind(neuralNetPredTimeDF,
                                     cbind(rownames(whiteDgrpProfnoOut)[i], tissue,
                                           'NeuralNet', neuronName, outliersName, 
                                           INFOTissue[rownames(whiteDgrpProfnoOut)[i]]$Time,
                                           preditNN$net.result[, 1]))
        neuralNet <- NULL
        preditNN <- NULL
      }
    }
  }
}

saveRDS(neuralNetPredTimeDF, 
        paste0(saveRDSdir, 'NeuralNet_realPredTimeDF_no25205.Rds'))

# Compare all the models by time diff and std dev -----------------------------
fileSuff <- '_realPredTimeDF_no25205.Rds'
lassoRealPredTimeDF <- readRDS(paste0(saveRDSdir, 'LASSO', fileSuff))
molTTrealPredTimeDF <- readRDS(paste0(saveRDSdir, 'MolecularTimeTable',
                                      '_realPredTimeDF_cv_no25205.Rds'))
zzRealPredTimeDF <- readRDS(paste0(saveRDSdir, 'ZeitZeiger', fileSuff))
neuralNetPredTimeDF <- readRDS(paste0(saveRDSdir, 'NeuralNet', fileSuff))
directStatsRealPredTimeDF <- readRDS(paste0(saveRDSdir, 'DirectStats', fileSuff))

# merge all methods together, reformat a bit
realPredTimeDF <- rbind(lassoRealPredTimeDF, molTTrealPredTimeDF, 
                        zzRealPredTimeDF, neuralNetPredTimeDF, 
                        directStatsRealPredTimeDF)
realPredTimeDF <- as.data.table(realPredTimeDF)
colnames(realPredTimeDF) <- c('Sample', 'Tissue', 'Method', 'Submethod',
                              'OutliersStatus', 'HarvTime', 'PredTime')
realPredTimeDF[, HarvTime := as.numeric(HarvTime) %% 24]
realPredTimeDF[, PredTime := as.numeric(PredTime)]
realPredTimeDF[, TimeDiff := apply(realPredTimeDF, 1, 
                                   function(x) diffInTime(as.numeric(x['HarvTime']),
                                                          as.numeric(x['PredTime'])))]
realPredTimeDF[Method == 'LASSO' | Method == 'ZeitZeiger', 
               Submethod := Method]
realPredTimeDF[Method == 'MolTT', Submethod := gsub('TIG_', '', Submethod)]

# calculate standart deviation and mean
realPredTimeSDmean <- realPredTimeDF[, .(Mean = mean(TimeDiff), 
                                         StDev = sd(TimeDiff)), 
                                     by = .(OutliersStatus, Submethod, Method)]
realPredTimeSDmean <- realPredTimeSDmean[order(Mean, StDev), ]

message('Top 3 best methods')
realPredTimeSDmean[OutliersStatus == 'Without'][1:3, ]

# plot - all tissues together
ggplot(realPredTimeDF[OutliersStatus == 'Without', ], 
       aes(x = Submethod, y = TimeDiff, fill = Method)) + 
  geom_boxplot() + scale_fill_brewer(palette = "Greys") +
  facet_wrap(~ Method, scales = "free_x", ncol = 5) +
  ylab('Difference between predicted and harversting time') +
  mashaGgplot2Theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position="none")

# plot - tissue by tissue
ggplot(realPredTimeDF[OutliersStatus == 'Without', ],
       aes(x = Submethod, y = TimeDiff, fill = Tissue)) + 
  geom_boxplot() + scale_fill_manual(values = tissueColor) +
  facet_grid(Tissue ~ Method, scales = "free_x") +
  ylab('Difference between predicted and harversting time') +
  mashaGgplot2Theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")