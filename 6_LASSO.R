# FILE DESCRIPTION ------------------------------------------------------------
#
# FILE: 7_LASSO.R
#
# DESCRIPTION : 
# Performs modelling of the physiological time based on gene expression with 
# use of LASSO
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
#  CREATED:  22.01.2018
#  REVISION: 22.01.2018
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

tissueOutliers <- list(FB = c('FB43_Oct_2', 'FB30_Oct_2', 'BRB2_FB_ZT11_29655'),
                       GUT = c('BRB5_Gut_ZT11_29655', 'BRB5_Gut_ZT23_29655', 
                               'BRB5_Gut_ZT11_25205', 'BRB4_Gut_ZT9_28211', 
                               'GUT31_Jan_6', 'BRB4_Gut_ZT21_29655'),
                       BRAIN = c('BRB8_Brain_ZT1_29655', 'BRAIN13_Oct_3', 
                                 'BRAIN37_Oct_3', 'BRAIN1_Oct_3', 
                                 'BRAIN28_Oct_1', 'BRB7_Brain_ZT21_28211'))
tissueOutliers <- list(FB = c('BRB2_FB_ZT11_29655', 'FB43_Oct_2'),
                       GUT = c('BRB5_Gut_ZT23_29655', 'BRB5_Gut_ZT11_25205', 
                               'BRB4_Gut_ZT9_28211', 'GUT6_Jan_4',
                               BRAIN = 'BRAIN1_Oct_3'))
tissueOutliers <- list(FB = c(), GUT = c(), BRAIN = c())

tissue <- 'FB'
outliers <- tissueOutliers[[tissue]]
colorForTissue <- tissueColor[tissue]
INFOTissue <- INFO[Tissue == tissue]
setkey(INFOTissue, rightGT_name)
cyclingGenesTiss <- cyclingGenes[[tissue]]

# WHITE & DGRP: READ-IN COUNTS, NORMALIZE, BATCH CORRECT ----------------------
# in this case, I put dgrps and w- into one table for normalization and 
# batch correction, so there wouldn't be batch effect bacause of it
# prepare info table
whiteDGRP <- readRDS(paste0(saveRDSdir, 'whiteDGRP.Rds'))[[tissue]]
# separate into the training set + test set : white and profiled DGRPs
# and set of interest - DGRPs profiled once
whiteDgrpProf <- whiteDGRP[, !grepl('_n_', colnames(whiteDGRP))]
newRownames <- colnames(whiteDgrpProf)
whiteDgrpProf <- as.data.frame(t(whiteDgrpProf))
rownames(whiteDgrpProf) <- newRownames
whiteDgrpProfnoOut <- whiteDgrpProf[!rownames(whiteDgrpProf) %in% outliers, ]

dgrpNorm <- whiteDGRP[, grepl('_n_', colnames(whiteDGRP))]
newRownames <- colnames(dgrpNorm)
dgrpNorm <- as.data.frame(t(dgrpNorm))
rownames(dgrpNorm) <- newRownames

# time for both sets
whiteDgrpProfTime <- INFOTissue[rownames(whiteDgrpProfnoOut)]$Time 
dgrpNormTime <- INFOTissue[rownames(whiteDgrpProfnoOut)]$Time

# LASSO, STEP 1: EVALUATE THE MODEL -------------------------------------------
# put difference in time here
realPredTimeDiff <- c()
# put coefficients for the first model (quadrant) here
quadCoeffs <- list()
lassoCoeffs_1_12 <- list()
lassoCoeffs_12_24 <- list()
for (i in 1:nrow(whiteDgrpProfnoOut)) {
  print(i)
  manualSel <- colnames(whiteDgrpProfnoOut)
  # divide into train and test dataset
  trainTestSets <- divideIntoTrainAndTest(whiteDgrpProfnoOut, INFOTissue,
                                          leaveOneOut = i, genesOI = manualSel)
  refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time
  refTrainTimeFact <- as.factor(ifelse(refTrainTime %% 24 > 12, 1, 0))
  
  # 1st STEP: TRAIN MODEL FOR DETERMINE, TO WHICH PART OF THE DAY SAMPLE
  # BELONGS
  # run training
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
  testPredict <- predict(lassoFit_DayTime, newx = trainTestSets$test, 
                         s = "lambda.min", family = 'binomial')
  predTestTime <- as.factor(ifelse(testPredict[, 1] > 0, 1, 0))
  if (predTestTime == 0) {
    testPredTime <- predict(lassoFit_1_12, newx = trainTestSets$test, 
                            s = "lambda.min")
  } else {
    testPredTime <- predict(lassoFit_12_24, newx = trainTestSets$test, 
                            s = "lambda.min")
  }
  
  # get mean square error
  refTestTime <- INFOTissue[rownames(trainTestSets$test)]$Time %% 24
  realPredTimeDiff <- c(realPredTimeDiff, 
                        diffInTime(testPredTime[, 1], refTestTime))
}
names(realPredTimeDiff) <- rownames(whiteDgrpProfnoOut)
png(paste0(savePlotsDir, tissue, '_Lasso7_predictRealTimeDiff.png'),
    width = 800, height = 800)
boxplot(realPredTimeDiff, 
        xlab = 'Difference between predicted and real time', 
        col = tissueColor[tissue], border = 'black', 
        names = 'Difference in time', horizontal = T, 
        main = paste('Difference between predicted and real time for LASSO 7',
                     'leave one out, no outliers,', tissue, sep = '\n'))
dev.off()
saveRDS(realPredTimeDiff, paste0(saveRDSdir, tissue, '_realPredTimeDiff_Lasso.Rds'))

# if barplot with all 3 tissues is needed
realPredTimeDiffAllTiss <- data.frame()
for (tiss in c('BRAIN', 'GUT', 'FB')) {
  realPredTimeDiffTiss <- readRDS(paste0(saveRDSdir, tiss, 
                                         '_realPredTimeDiff.Rds'))
  realPredTimeDiffAllTiss <- rbind(realPredTimeDiffAllTiss,
                                   data.frame(tissue = tiss,
                                              realPredTimeDiffTiss))
}
colnames(realPredTimeDiffAllTiss) <- c('tissue', 'timeDiff')
realPredTimeDiffAllTiss$timeDiff <- as.numeric(realPredTimeDiffAllTiss$timeDiff)
png(paste0(savePlotsDir, '_allTissLasso7_predictRealTimeDiff.png'),
    width = 800, height = 500)
boxplot(timeDiff ~ tissue, realPredTimeDiffAllTiss,
        xlab = 'Difference between predicted and real time', 
        col = tissueColor[c('BRAIN', 'GUT', 'FB')], border = 'black', 
        horizontal = T, pch = 19, cex.lab = 2, cex.axis = 1.5, 
        main = paste('Difference between predicted and real time for LASSO 7',
                     'leave one out, no outliers', sep = '\n'))
dev.off()

# LASSO, STEP 2: PHYSIOLOGICAL TIME BASED ON LASSO  ---------------------------
whiteDgrpProfnoOut <- whiteDgrpProf[!rownames(whiteDgrpProf) %in% outliers, ]
# 1st STEP: TRAIN MODEL FOR DETERMINE, TO WHICH PART OF THE DAY SAMPLE
# BELONGS
# run training
trainTest <- whiteDgrpProfnoOut[, complete.cases(t(whiteDgrpProfnoOut))]
trainTime <- INFOTissue[rownames(whiteDgrpProfnoOut)]$Time
trainTimeFact <- as.factor(ifelse(trainTime %% 24 > 12, 1, 0))
lasso_DayTime <- lassoOnTrainSet(as.matrix(trainTest), trainTimeFact, F,
                                 family = 'binomial')
lasso_DayTimeCoeffs <- coef(lasso_DayTime)[, 1][coef(lasso_DayTime)[, 1] != 0]

# 2nd STEP: TRAN MODELS FOR BEFORE AND AFTER 12
# create day time specific train sets
trainSet_Q1 <- trainTest[trainTime  %% 24 <= 12, ]
trainTime_Q1 <- INFOTissue[rownames(trainSet_Q1)]$Time %% 24
trainSet_Q2 <- trainTest[trainTime  %% 24 >= 12, ]
trainTime_Q2 <- INFOTissue[rownames(trainSet_Q2)]$Time %% 24
# fit 
lassoFit_Q1 <- lassoOnTrainSet(as.matrix(trainSet_Q1), trainTime_Q1, F)
coeffs_Q1 <- coef(lassoFit_Q1)[, 1][coef(lassoFit_Q1)[, 1] != 0]
lassoFit_Q2 <- lassoOnTrainSet(as.matrix(trainSet_Q2), trainTime_Q2, F)
coeffs_Q2 <- coef(lassoFit_Q2)[, 1][coef(lassoFit_Q2)[, 1] != 0]
# 3rd STEP: predict
# run prediction on test
dgrpNorm <- dgrpNorm[, complete.cases(t(dgrpNorm))]
dgrpNormQuadr <- predict(lasso_DayTime, newx = as.matrix(dgrpNorm), 
                         s = "lambda.min", family = 'binomial')
dgrpNormQuadrPred <- as.factor(ifelse(dgrpNormQuadr[, 1] > 0, 1, 0))
# divide them into quadrants
dgrpNorm_Q1 <- as.matrix(dgrpNorm[dgrpNormQuadrPred == 0, ])
dgrpNorm_Q2 <- as.matrix(dgrpNorm[dgrpNormQuadrPred != 0, ])
# make prediction of physiological time
dgrpNorm_Q1_PhysioTime <- predict(lassoFit_Q1, newx = dgrpNorm_Q1, 
                                  s = "lambda.min", interval = "predict")
dgrpNorm_Q2_PhysioTime <- predict(lassoFit_Q2, newx = dgrpNorm_Q2, 
                                  s = "lambda.min")
dgrpNormPhysioTime <- c(dgrpNorm_Q1_PhysioTime[, 1],
                        dgrpNorm_Q2_PhysioTime[, 1])
names(dgrpNormPhysioTime) <- c(rownames(dgrpNorm_Q1), rownames(dgrpNorm_Q2))
dgrpNormHarvTime <- INFOTissue[names(dgrpNormPhysioTime)]$Time %% 24
# put them to df for convinience
physioPredTimeDF <- data.frame(Sample = names(dgrpNormPhysioTime),
                               timeObs = dgrpNormHarvTime, 
                               timePred = dgrpNormPhysioTime)
physioPredTimeDF$timeError <- apply(physioPredTimeDF, 1, 
                                    function(x) diffInTime(as.numeric(x[2]), 
                                                           as.numeric(x[3])))

saveRDS(lasso_DayTime, paste0(saveRDSdir, tissue, '_lasso_DayTime_LASSO.Rds'))
saveRDS(lassoFit_Q1, paste0(saveRDSdir, tissue, '_lassoFit_Q1_LASSO.Rds'))
saveRDS(lassoFit_Q2, paste0(saveRDSdir, tissue, '_lassoFit_Q2_LASSO.Rds'))
saveRDS(physioPredTimeDF, paste0(saveRDSdir, tissue, '_physioPredTimeDF_LASSO.Rds'))

# BUBBLE PLOT: GENE-PREDICTORS AND DGRP-OUTLIERS TIME -------------------------
allTissues <- c('BRAIN', 'GUT', 'FB')
# restrict cycling genes only to significant ones
cyclingGenes <- lapply(cyclingGenes, function(x) x[x$ADJ.P < 0.05, ])
# names of all DGRP-s outliers
dgrpTimeOutliersAllTiss <- data.table(tissue = character(), 
                                      DGRPid = character(),
                                      numbHours = numeric())
# coefficients for the first model, which predicts quadrant of the day
dayTimeCoeffs <- data.table(tissue = character(), Gene = character(), 
                            Coeff = character(), Cycles = logical())
# same for the model for the first 12h
coeffQ1AllTiss <- data.table(tissue = character(), Gene = character(), 
                             Coeff = character(), Cycles = logical())
# same for the model for the second 12h
coeffsQ2AllTiss <- data.table(tissue = character(), Gene = character(), 
                              Coeff = character(), Cycles = logical())
for (tis in allTissues) {
  physioPredTimeDF <- readRDS(paste0(saveRDSdir, tis, '_physioPredTimeDF.Rds'))
  realPredTimeDiff <- readRDS(paste0(saveRDSdir, tis, '_realPredTimeDiff.Rds'))
  lasso_DayTime <- readRDS(paste0(saveRDSdir, tis, '_lasso_DayTime.Rds'))
  lassoFit_Q1 <- readRDS(paste0(saveRDSdir, tis, '_lassoFit_Q1.Rds'))
  lassoFit_Q2 <- readRDS(paste0(saveRDSdir, tis, '_lassoFit_Q2.Rds'))
  
  # find and plot outliers
  physioHarvTimeDiff <- apply(physioPredTimeDF, 1, 
                              function(x) diffInTime(x[1], x[2]))
  timeDiffCutOff <- max(realPredTimeDiff) + 1
  dgrpTimeOutliers <- sapply(rownames(physioPredTimeDF[physioHarvTimeDiff > 
                                                         timeDiffCutOff, ]), 
                             function(x) strsplit(x, "_")[[1]][4])
  dgrpTimeOutliersAllTiss <- rbind(dgrpTimeOutliersAllTiss,
                                   data.frame(tissue = tis, 
                                              DGRPid = dgrpTimeOutliers,
                                              numbHours = physioHarvTimeDiff[physioHarvTimeDiff > timeDiffCutOff]))
  
  svg(paste0(savePlotsDir, tis, '_physioHarvTimeOutliers.svg'),
      width = 8, height = 8)
  plot(physioPredTimeDF, pch = 20, bty = 'n',
       col = ifelse(physioHarvTimeDiff > timeDiffCutOff, tissueColor[tis], 
                    'black'), cex = 2, 
       xlim = c(0, 25), ylim = c(0, 25), cex.lab = 2, cex.axis = 2,
       xlab = 'Time of harversting, h', ylab = 'Physiological time, h',
       main = paste0('Physiological and harversting time, ', tis))
  abline(1, 1, lty = 2)
  abline(max(realPredTimeDiff) + 1, 1, lty = 2)
  abline(-max(realPredTimeDiff) - 1, 1, lty = 2)
  abline(2 * max(realPredTimeDiff) + 2, 1, lty = 2)
  abline(-2 * max(realPredTimeDiff) - 2, 1, lty = 2)
  text(physioPredTimeDF[physioHarvTimeDiff > timeDiffCutOff, ], 
       labels = dgrpTimeOutliers, pos = 3)
  grid()
  dev.off()
  
  # get genes - predictors
  dayTimeDF <- coef(lasso_DayTime)[, 1][coef(lasso_DayTime)[, 1] != 0]
  dayTimeDF <- data.table(tissue = tis, Gene = names(dayTimeDF), 
                          Coeff = dayTimeDF, 
                          Cycles = names(dayTimeDF) %in% 
                            rownames(cyclingGenes[[tis]])) 
  dayTimeCoeffs <- rbind(dayTimeCoeffs, dayTimeDF)
  
  Q1_df <- coef(lassoFit_Q1)[, 1][coef(lassoFit_Q1)[, 1] != 0]
  Q1_df <- data.table(tissue = tis, Gene = names(Q1_df), Coeff = Q1_df, 
                      Cycles = names(Q1_df) %in% 
                        rownames(cyclingGenes[[tis]]))
  coeffQ1AllTiss <- rbind(coeffQ1AllTiss, Q1_df)
  
  Q2_df <- coef(lassoFit_Q2)[, 1][coef(lassoFit_Q2)[, 1] != 0]
  Q2_df <- data.table(tissue = tis, Gene = names(Q2_df), Coeff = Q2_df,
                      Cycles = names(Q2_df) %in% 
                        rownames(cyclingGenes[[tis]])) 
  coeffsQ2AllTiss <- rbind(coeffsQ2AllTiss, Q2_df)
}

png(paste0(savePlotsDir, 'coeffsModel1.png'), width = 1150, height = 500)
print(plotModelCoeffs(dayTimeCoeffs, 'Model 1'))
dev.off()
png(paste0(savePlotsDir, 'coeffsModel2A.png'), width = 1150, height = 500)
print(plotModelCoeffs(coeffQ1AllTiss, 'Model 2A: 1 - 12h'))
dev.off()
png(paste0(savePlotsDir, 'coeffsModel2B.png'), width = 1150, height = 500)
print(plotModelCoeffs(coeffsQ2AllTiss, 'Model 2B: 12 - 24h'))
dev.off()

png(paste0(savePlotsDir, 'dgrpOutliers.png'), 1500, 500)
ggplot(dgrpTimeOutliersAllTiss, aes(DGRPid, tissue)) +
  geom_point(aes(size = abs(dgrpTimeOutliersAllTiss$numbHours))) + 
  mashaGgplot2Theme + scale_size_continuous(range = c(1, 15)) + 
  labs(size = "Physiological to\nharversting time\ndifference") +
  xlab("DGRP line") + ylab("Tissue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('Difference between physiological and harversting time')
dev.off()

# PLOT GENES-PREDICTORS -------------------------------------------------------
whiteDGRP <- readRDS(paste0(saveRDSdir, 'whiteDGRP.Rds'))
allTissues <- c('BRAIN', 'GUT', 'FB')
allTissCounts <- lapply(whiteDGRP[allTissues], 
                        function(x) x[, !grepl('_n_', colnames(x))])
names(allTissCounts) <- allTissues

# predictors for the first model
genesPred1stModel <- c('tim', 'CG13631', 'Cp190')
listToPlot <- list()
for (tiss in allTissues) {
  for (goi in genesPred1stModel) {
    INFOTiss <- INFO[Tissue == tiss]
    setkey(INFOTiss, rightGT_name)
    pl <- plotCountsTimecourse(goi, allTissCounts[[tiss]],
                               INFOTiss[colnames(allTissCounts[[tiss]])]$Time %% 24,
                               tissueColor[tiss], 
                               ggtitle(paste(goi, 'expression in', tiss)))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
}
png(paste0(savePlotsDir, 'genesPredModel1.png'), height = 1200,
    width = 1200)
multiplot(plotlist = listToPlot, cols = 3)
dev.off()

genesPred2AModel <- c('CG18493', 'CG11407', 'SP1173')
listToPlot <- list()
for (tiss in allTissues) {
  for (goi in genesPred2AModel) {
    INFOTiss <- INFO[Tissue == tiss]
    setkey(INFOTiss, rightGT_name)
    pl <- plotCountsTimecourse(goi, allTissCounts[[tiss]],
                               INFOTiss[colnames(allTissCounts[[tiss]])]$Time %% 24,
                               tissueColor[tiss],
                               ggtitle(paste(goi, 'expression in', tiss)))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
}
png(paste0(savePlotsDir, 'genesPredModel2A.png'), height = 1200,
    width = 1200)
multiplot(plotlist = listToPlot, cols = 3)
dev.off()

genesPred2BModel <- c('CG33127', 'CG11407', 'geko')
listToPlot <- list()
for (goi in genesPred2BModel) {
  for (tiss in allTissues) {
    INFOTiss <- INFO[Tissue == tiss]
    setkey(INFOTiss, rightGT_name)
    pl <- plotCountsTimecourse(goi, allTissCounts[[tiss]],
                               INFOTiss[colnames(allTissCounts[[tiss]])]$Time %% 24,
                               tissueColor[tiss],
                               ggtitle(paste(goi, 'expression in', tiss)))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
}
png(paste0(savePlotsDir, 'genesPredModel2B.png'), height = 1000,
    width = 1000)
multiplot(plotlist = listToPlot, cols = 3)
dev.off()

# OVERLAP BETWEEN GENES-PREDICTORS --------------------------------------------
# overlap between genes-predictors for the 1st model (daytime) in different
# tissue
forVenn <- list(BRAIN = dayTimeCoeffs[tissue == 'BRAIN']$Gene,
                GUT = dayTimeCoeffs[tissue == 'GUT']$Gene,
                FB = dayTimeCoeffs[tissue == 'FB']$Gene)
forVenn <- lapply(forVenn, function(x) x[x != '(Intercept)'])
venn(forVenn, ilab = T, zcolor = "style", cexil = 1.3)
# overlap between genes-predictors for the 2nd model (0 - 12h) in different
# tissue
forVenn <- list(BRAIN = coeffQ1AllTiss[tissue == 'BRAIN']$Gene,
                GUT = coeffQ1AllTiss[tissue == 'GUT']$Gene,
                FB = coeffQ1AllTiss[tissue == 'FB']$Gene)
forVenn <- lapply(forVenn, function(x) x[x != '(Intercept)'])
venn(forVenn, ilab = T, zcolor = "style", cexil = 1.3)
# overlap between genes-predictors for the 2nd model (12 - 24h) in different
# tissue
forVenn <- list(BRAIN = coeffsQ2AllTiss[tissue == 'BRAIN']$Gene,
                GUT = coeffsQ2AllTiss[tissue == 'GUT']$Gene,
                FB = coeffsQ2AllTiss[tissue == 'FB']$Gene)
forVenn <- lapply(forVenn, function(x) x[x != '(Intercept)'])
venn(forVenn, ilab = T, zcolor = "style", cexil = 1.3)
# overlap between genes-predictors for the same tissue, different models
for (tiss in allTissues) {
  forVenn <- list(Model1 = dayTimeCoeffs[tissue == tis]$Gene,
                  Model2A = coeffQ1AllTiss[tissue == tis]$Gene,
                  Model2B = coeffsQ2AllTiss[tissue == tis]$Gene)
  forVenn <- lapply(forVenn, function(x) x[x != '(Intercept)'])
  venn(forVenn, ilab = T, zcolor = "style", cexil = 1.3)
}

# MODELS, WHICH DIDN'T WORK ---------------------------------------------------
#' LASSO 1: train = white, test = DGRP  ---------------------------------------
pdf('~/Desktop/allLasso.pdf')
# 76% with 1st manual selection
CVerr <- c()
for (i in 1:100) {
  manualSel <- rownames(white)
  # divide into train and test dataset
  trainTestSets <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue,
                                          mixed = F, genesOI = manualSel)
  refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time %% 24
  # run training
  lassoFit <- lassoOnTrainSet(trainTestSets$train, refTrainTime)
  # run prediction on test
  testPredict <- predict(lassoFit, newx = trainTestSets$test, s = "lambda.min")
  refTestTime <- INFOTissue[rownames(trainTestSets$test)]$Time %% 24
  # get mean square error
  CVerr <- c(CVerr, MSEforClock(testPredict[, 1], refTestTime))
  # plot(INFOTissue[refTestTime, testPredict[, 1])
}
boxplot(CVerr, xlab = 'MSE CV', col = 'darkgrey', border = 'black', 
        names = 'CV error (MSE)', horizontal = T, 
        main = 'MSE for LASSO 1\n train = white, test = DGRP, all 24h')


#' LASSO 2: train = white&DGRPmix, test = white&DGRPmix -----------------------
CVerr <- c()
for (i in 1:100) {
  manualSel <- rownames(white)
  # divide into train and test dataset
  trainTestSets <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue,
                                          mixed = T, genesOI = manualSel)
  refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time %% 24
  # run training
  lassoFit <- lassoOnTrainSet(trainTestSets$train, refTrainTime)
  # run prediction on test
  testPredict <- predict(lassoFit, newx = trainTestSets$test, s = "lambda.min")
  refTestTime <- INFOTissue[rownames(trainTestSets$test)]$Time %% 24
  # get mean square error
  CVerr <- c(CVerr, MSEforClock(testPredict[, 1], refTestTime))
  # plot(INFOTissue[refTestTime, testPredict[, 1])
}
boxplot(CVerr, xlab = 'MSE CV', col = 'darkgrey', border = 'black', 
        names = 'CV error (MSE)', horizontal = T, 
        main = 'MSE for LASSO 2\n train = white&DGRPmix, test = white&DGRPmix, all 24h')


#' LASSO 3:train=white&DGRPmix, test=white&DGRPmix, time between 0 - 12--------
CVerr <- c()
for (i in 1:100) {
  manualSel <- rownames(white)
  # divide into train and test dataset
  trainTestSets <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue,
                                          mixed = T, genesOI = manualSel, 
                                          timePoint1 = 0, timePoint2 = 12,
                                          timeIntervalType = 'inner')
  refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time %% 24
  # run training
  lassoFit <- lassoOnTrainSet(trainTestSets$train, refTrainTime)
  # run prediction on test
  testPredict <- predict(lassoFit, newx = trainTestSets$test, s = "lambda.min")
  refTestTime <- INFOTissue[rownames(trainTestSets$test)]$Time %% 24
  # get mean square error
  CVerr <- c(CVerr, MSEforClock(testPredict[, 1], refTestTime))
  # plot(INFOTissue[refTestTime, testPredict[, 1])
}
boxplot(CVerr, xlab = 'MSE CV', col = 'darkgrey', border = 'black', 
        names = 'CV error (MSE)', horizontal = T, 
        main = 'MSE for LASSO 3\n train = white&DGRPmix, test = white&DGRPmix, 0 - 12h')


#' LASSO 4:train=white&DGRPmix, test=white&DGRPmix, time between 12 - 24--------
CVerr <- c()
for (i in 1:100) {
  manualSel <- rownames(white)
  # divide into train and test dataset
  trainTestSets <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue,
                                          mixed = T, genesOI = manualSel, 
                                          timePoint1 = 12, timePoint2 = 24,
                                          timeIntervalType = 'inner')
  refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time %% 24
  # run training
  lassoFit <- lassoOnTrainSet(trainTestSets$train, refTrainTime)
  # run prediction on test
  testPredict <- predict(lassoFit, newx = trainTestSets$test, s = "lambda.min")
  refTestTime <- INFOTissue[rownames(trainTestSets$test)]$Time %% 24
  # get mean square error
  CVerr <- c(CVerr, MSEforClock(testPredict[, 1], refTestTime))
  # plot(INFOTissue[refTestTime, testPredict[, 1])
}
boxplot(CVerr, xlab = 'MSE CV', col = 'darkgrey', border = 'black', 
        names = 'CV error (MSE)', horizontal = T, 
        main = 'MSE for LASSO 4\n train = white&DGRPmix, test = white&DGRPmix, 12-24h')


#' LASSO 5:train=white&DGRPmix, test=white&DGRPmix, get quadrant --------------
FPR <- c()
TPR <- c()
ACC <- c()
SPC <- c()
for (i in 1:200) {
  manualSel <- rownames(white)
  # divide into train and test dataset
  trainTestSets <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue,
                                          mixed = T, genesOI = manualSel)
  refTrainTime <- as.factor(ifelse(INFOTissue[rownames(trainTestSets$train)]$Time %% 24 > 12,
                                   1, 0))
  # run training
  lassoFit <- lassoOnTrainSet(trainTestSets$train, refTrainTime, F, 
                              family = 'binomial')
  # run prediction on test
  testPredict <- predict(lassoFit, newx = trainTestSets$test, s = "lambda.min",
                         family = 'binomial')
  refTestTime <-  ifelse(INFOTissue[rownames(trainTestSets$test)]$Time %% 24 > 12,
                         1, 0)
  # evaluate predictions
  predTestTime <- as.factor(ifelse(testPredict[, 1] > 0, 1, 0))
  testTime <- data.frame(predTestTime, refTestTime)
  TP <- sum(apply(testTime, 1, 
                  function(x) ifelse(x[1] == x[2] & x[2] == 1, 1, 0)))
  TN <- sum(apply(testTime, 1, 
                  function(x) ifelse(x[1] == x[2] & x[2] == 0, 1, 0)))
  FP <- sum(apply(testTime, 1, 
                  function(x) ifelse(x[1] != x[2] & x[2] == 0, 1, 0)))
  FN <- sum(apply(testTime, 1, 
                  function(x) ifelse(x[1] != x[2] & x[2] == 1, 1, 0)))
  FPR <- c(FPR, FP / (FP + TN))
  TPR <- c(TPR, TP / (TP + FN))
  ACC <- c(ACC, (TP + TN) / (TP + TN + FP + FN))
  SPC <- c(SPC, TN / (TN + FP))
}

#' LASSO 6:train=white&DGRPmix, test=white&DGRPmix, 2step model ---------------
CVerr <- c()
for (i in 1:100) {
  manualSel <- rownames(white)
  # divide into train and test dataset
  trainTestSets <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue,
                                          mixed = T, genesOI = manualSel)
  refTrainTime <- INFOTissue[rownames(trainTestSets$train)]$Time
  refTrainTimeFact <- as.factor(ifelse(refTrainTime %% 24 > 12, 1, 0))
  
  # 1st STEP: DETERMINE, TO WHICH PART OF THE DAY SAMPLE BELONGS
  # run training
  lassoFit_DayTime <- lassoOnTrainSet(trainTestSets$train, refTrainTimeFact, 
                                      F, family = 'binomial')
  # run prediction on test
  testPredict <- predict(lassoFit_DayTime, newx = trainTestSets$test, 
                         s = "lambda.min", family = 'binomial')
  predTestTime <- as.factor(ifelse(testPredict[, 1] > 0, 1, 0))
  
  # 2nd STEP: RUN MODEL ACCORDING TO THE TIME OF THE DAY
  # create day time specific train sets
  trainSet_1_12 <- trainTestSets$train[refTrainTime  %% 24 <= 12, ]
  refTrainTime_1_12 <- INFOTissue[rownames(trainSet_1_12)]$Time %% 24
  trainSet_12_24 <- trainTestSets$train[refTrainTime  %% 24 >= 12, ]
  refTrainTime_12_24 <- INFOTissue[rownames(trainSet_12_24)]$Time %% 24
  # create day time specific test sets
  testSet_1_12 <- trainTestSets$test[predTestTime == 0, ]
  testSet_12_24 <- trainTestSets$test[predTestTime == 1, ]
  # fit 
  lassoFit_1_12 <- lassoOnTrainSet(trainSet_1_12, refTrainTime_1_12, F)
  lassoFit_12_24 <- lassoOnTrainSet(trainSet_12_24, refTrainTime_12_24, F)
  # predict
  testPredict_1_12 <- predict(lassoFit_1_12, newx = testSet_1_12, 
                              s = "lambda.min")
  testPredict_12_24 <- predict(lassoFit_12_24, newx = testSet_12_24, 
                               s = "lambda.min")
  
  refTestTime_1_12 <- INFOTissue[rownames(testSet_1_12)]$Time %% 24
  refTestTime_12_24 <- INFOTissue[rownames(testSet_12_24)]$Time %% 24
  # get mean square error
  CVerr_1_12 <- MSEforClock(testPredict_1_12[, 1], refTestTime_1_12)
  CVerr_12_24 <- MSEforClock(testPredict_12_24[, 1], refTestTime_12_24)
  CVerr <- c(CVerr, CVerr_1_12, CVerr_12_24)
}
boxplot(CVerr, xlab = 'MSE CV', col = 'darkgrey', border = 'black', 
        names = 'CV error (MSE)', horizontal = T, 
        main = 'MSE for LASSO 6')
dev.off()

#' LASSO 8:leave one out, 2step model + 4 time intervals, 0 - 18 ---------------
realAndPredict <- data.frame()
for (i in 1:nrow(whiteDgrpProf)) {
  print(i)
  trainSet_0_18 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                          leaveOneOut = i, genesOI = manualSel,
                                          timePoint1 = 0, timePoint2 = 18, 
                                          timeIntervalType = 'inner')$train
  trainTime_0_18 <- INFOTissue[rownames(trainSet_0_18)]$Time %% 24
  lassoFit_0_18 <- lassoOnTrainSet(trainSet_0_18, trainTime_0_18, F)
  
  trainSet_6_24 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                          leaveOneOut = i, genesOI = manualSel, 
                                          timePoint1 = 6, timePoint2 = 24, 
                                          timeIntervalType = 'inner')$train
  trainTime_6_24 <- INFOTissue[rownames(trainSet_6_24)]$Time %% 24
  lassoFit_6_24 <- lassoOnTrainSet(trainSet_6_24, trainTime_6_24, F)
  
  trainSet_12_30 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                           leaveOneOut = i, 
                                           genesOI = manualSel, 
                                           timePoint1 = 6, timePoint2 = 12, 
                                           timeIntervalType = 'outer')$train
  trainTime_12_30 <- INFOTissue[rownames(trainSet_12_30)]$Time
  trainTime_12_30[trainTime_12_30 <= 6] <- trainTime_12_30[trainTime_12_30 <= 6] + 24
  lassoFit_12_30 <- lassoOnTrainSet(trainSet_12_30, trainTime_12_30, F)
  
  trainSet_18_36 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                           leaveOneOut = i, 
                                           genesOI = manualSel, 
                                           timePoint1 = 12, timePoint2 = 18, 
                                           timeIntervalType = 'outer')$train
  trainTime_18_36 <- INFOTissue[rownames(trainSet_18_36)]$Time
  trainTime_18_36[trainTime_18_36 <= 12] <- trainTime_18_36[trainTime_18_36 <= 12] + 24
  lassoFit_18_36 <- lassoOnTrainSet(trainSet_18_36, trainTime_18_36, F)
  
  testPredict_0_18 <- predict(lassoFit_0_18, 
                              newx = as.matrix(whiteDgrpProf[i, ]),
                              s = "lambda.min") %% 24
  testPredict_6_24 <- predict(lassoFit_6_24, 
                              newx = as.matrix(whiteDgrpProf[i, ]),
                              s = "lambda.min") %% 24
  testPredict_12_30 <- predict(lassoFit_12_30, 
                               newx = as.matrix(whiteDgrpProf[i, ]),
                               s = "lambda.min") %% 24
  testPredict_18_36 <- predict(lassoFit_18_36, 
                               newx = as.matrix(whiteDgrpProf[i, ]),
                               s = "lambda.min") %% 24
  realAndPredict <- rbind(realAndPredict, 
                          t(c(INFOTissue[rownames(whiteDgrpProf)[i]]$Time %% 24, 
                              testPredict_0_18, testPredict_6_24,
                              testPredict_12_30, testPredict_18_36)))
}

#' LASSO 9:leave one out, 2step model + 4 time intervals, 0 - 12 --------------
realAndPredict <- data.frame()
for (i in 1:nrow(whiteDgrpProf)) {
  print(i)
  trainSet_0_12 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                          leaveOneOut = i, genesOI = manualSel,
                                          timePoint1 = 0, timePoint2 = 18, 
                                          timeIntervalType = 'inner')$train
  trainTime_0_12 <- INFOTissue[rownames(trainSet_0_12)]$Time %% 24
  lassoFit_0_12 <- lassoOnTrainSet(trainSet_0_12, trainTime_0_12, F)
  
  trainSet_6_18 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                          leaveOneOut = i, genesOI = manualSel, 
                                          timePoint1 = 6, timePoint2 = 24, 
                                          timeIntervalType = 'inner')$train
  trainTime_6_18 <- INFOTissue[rownames(trainSet_6_18)]$Time %% 24
  lassoFit_6_18 <- lassoOnTrainSet(trainSet_6_18, trainTime_6_18, F)
  
  trainSet_12_24 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                           leaveOneOut = i, 
                                           genesOI = manualSel, 
                                           timePoint1 = 6, timePoint2 = 12, 
                                           timeIntervalType = 'outer')$train
  trainTime_12_24 <- INFOTissue[rownames(trainSet_12_24)]$Time
  trainTime_12_24[trainTime_12_24 <= 6] <- trainTime_12_24[trainTime_12_24 <= 6] + 24
  lassoFit_12_24 <- lassoOnTrainSet(trainSet_12_24, trainTime_12_24, F)
  
  trainSet_18_6 <- divideIntoTrainAndTest(whiteDgrpProf, INFOTissue, 
                                          leaveOneOut = i, 
                                          genesOI = manualSel, 
                                          timePoint1 = 12, timePoint2 = 18, 
                                          timeIntervalType = 'outer')$train
  trainTime_18_6 <- INFOTissue[rownames(trainSet_18_6)]$Time
  trainTime_18_6[trainTime_18_6 <= 12] <- trainTime_18_6[trainTime_18_6 <= 12] + 24
  lassoFit_18_6 <- lassoOnTrainSet(trainSet_18_6, trainTime_18_6, F)
  
  testPredict_0_12 <- predict(lassoFit_0_12, 
                              newx = as.matrix(whiteDgrpProf[i, ]),
                              s = "lambda.min") %% 24
  testPredict_6_18 <- predict(lassoFit_6_18, 
                              newx = as.matrix(whiteDgrpProf[i, ]),
                              s = "lambda.min") %% 24
  testPredict_12_24 <- predict(lassoFit_12_24, 
                               newx = as.matrix(whiteDgrpProf[i, ]),
                               s = "lambda.min") %% 24
  testPredict_18_6 <- predict(lassoFit_18_6, 
                              newx = as.matrix(whiteDgrpProf[i, ]),
                              s = "lambda.min") %% 24
  realAndPredict <- rbind(realAndPredict, 
                          t(c(INFOTissue[rownames(whiteDgrpProf)[i]]$Time %% 24, 
                              testPredict_0_12, testPredict_6_18, 
                              testPredict_12_24, testPredict_18_6)))
}