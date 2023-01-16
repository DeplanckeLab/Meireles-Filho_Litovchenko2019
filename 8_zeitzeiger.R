# FILE DESCRIPTION: 24_zeitzeiger ---------------------------------------------
#
# DESCRIPTION : Performs modelling of the physiological time based on gene 
# expression with use of zeitzeiger
#
# USAGE: 
#
# OPTIONS:  none
# REQUIREMENTS:  
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  06.04.2018
# REVISION: 06.04.2018

setwd('~/Desktop/BitBucket/AroundTheClock_Aug_Oct2017/')
setwd('~/Desktop/AroundTheClock_Aug_Oct2017/')
source('0_common_functions.R')

# INPUTS ----------------------------------------------------------------------
# info about samples
INFO <- readRDS('Rds/sampleInfo.Rds')
INFO <- as.data.table(INFO)
setkey(INFO, rightGT_name)
# list all tissues and corresponding colors
allTissues <- unique(INFO$Tissue)
tissueColor <- c('darkorchid3', 'brown3', 'darkgoldenrod1', 'darkgreen')
names(tissueColor) <- allTissues

# output dirs
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'
saveBedDir <- 'beds/'

tissueOutliers <- list(FB = c('FB43_Oct_2', 'FB30_Oct_2', 'BRB2_FB_ZT11_29655'),
                       GUT = c('BRB5_Gut_ZT11_29655', 'BRB5_Gut_ZT23_29655', 
                               'BRB5_Gut_ZT11_25205', 'BRB4_Gut_ZT9_28211', 
                               'GUT31_Jan_6', 'BRB4_Gut_ZT21_29655'),
                       BRAIN = c('BRAIN1_Oct_3', 'BRAIN13_Oct_3',
                                 'BRAIN28_Oct_1', 'BRAIN37_Oct_3',
                                 'BRB7_Brain_ZT21_28211', 
                                 'BRB8_Brain_ZT1_29655'))

# WHITE & DGRP: READ-IN COUNTS, NORMALIZE, BATCH CORRECT ----------------------
# in this case, I put dgrps and w- into one table for normalization and 
# batch correction, so there wouldn't be batch effect bacause of it
# prepare info table
tissue <- 'GUT'
whiteDGRP <- readRDS(paste0(saveRDSdir, 'whiteDGRP.Rds'))[[tissue]]
whiteDGRP <- whiteDGRP[, !colnames(whiteDGRP) %in% tissueOutliers[[tissue]]]

# separate into the training set + test set : white and profiled DGRPs
# and set of interest - DGRPs profiled once
whiteDgrpProf <- whiteDGRP[, !grepl('_n_', colnames(whiteDGRP))]
newRownames <- colnames(whiteDgrpProf)
whiteDgrpProf <- as.data.frame(t(whiteDgrpProf))
rownames(whiteDgrpProf) <- newRownames

dgrpNorm <- whiteDGRP[, grepl('_n_', colnames(whiteDGRP))]
newRownames <- colnames(dgrpNorm)
dgrpNorm <- as.data.frame(t(dgrpNorm))
rownames(dgrpNorm) <- newRownames

# time for both sets
whiteDgrpProfTime <- (INFO[rownames(whiteDgrpProf)]$Time %% 24 ) / 24
dgrpNormTime <- (INFO[rownames(dgrpNorm)]$Time %% 24 ) / 24

# STEP 1: best parameters -----------------------------------------------------
# To determine the best parameters for training a ZeitZeiger predictor, we can
# run cross-validation
registerDoParallel(cores = 4)
sumabsv <- c(1, 1.5, 3) # amount of regularization
nSpc <- 1:5 # how many SPCs are used for prediction
nFolds <- 10 # can only range from 1 to 10
foldid <- sample(rep(1:nFolds, length.out = nrow(whiteDgrpProf)))

# calculate cross-validation
# fit a periodic smoothing spline to the behavior of each feature as a function
# of time
fitCrossVal <- zeitzeigerFitCv(whiteDgrpProf, whiteDgrpProfTime, foldid)
spcCrossVal <- list()
predCrossVal <-  list()
for (oneSumabsv in 1:length(sumabsv)) {
  # fits to calculate sparse principal components (SPCs) for how the features
  # change over time
  spcCrossVal[[oneSumabsv]] <- zeitzeigerSpcCv(fitCrossVal, 
                                               sumabsv = sumabsv[oneSumabsv])
  # uses the training data and the SPCs to predict the corresponding time for
  # each test observation
  predCrossVal[[oneSumabsv]] <- zeitzeigerPredictCv(as.matrix(whiteDgrpProf), 
                                                    whiteDgrpProfTime, 
                                                    foldid, 
                                                    spcCrossVal[[oneSumabsv]],
                                                    nSpc = nSpc)
}

# reorganize the output, making a data.frame with the information for each 
# prediction.
timePredCrossVal <- lapply(predCrossVal, function(x) x$timePred)
crossValRes <- data.frame(do.call(rbind, timePredCrossVal),
                          timeObs = rep(whiteDgrpProfTime, length(sumabsv)),
                          sumabsv = rep(sumabsv, 
                                        each = length(whiteDgrpProfTime)),
                          obs = rep(rownames(whiteDgrpProf), length(sumabsv)),
                          stringsAsFactors = F)
crossValResGath <- gather(crossValRes, key = nSpc, value = timePred, -obs, 
                          -timeObs, -sumabsv)
crossValResGath$nSpc <- as.integer(sapply(as.character(crossValResGath$nSpc),
                                          function(x) substr(x, 2, nchar(x))))
crossValResGath$sumabsv <- factor(crossValResGath$sumabsv)
crossValResGath$timeError <- calcTimeDiff(crossValResGath$timeObs, 
                                          crossValResGath$timePred)
# median absolute error for each set of parameter values.
crossValResGathGroup <- crossValResGath %>% group_by(sumabsv, nSpc) %>%
                        summarize(medae = median(abs(timeError)))
# Plot the error for each set of parameter values
ggplot(crossValResGathGroup) +
  geom_point(aes(x = nSpc, y = medae, shape = sumabsv, color = sumabsv),
             size=2) +
  labs(x = 'Number of SPCs', y = 'Median absolute error') +
  theme_bw() + theme(legend.position = c(0.7, 0.7))
# !!! the best accuracy seems to be at sumabsv = 1.5 and nSpc = 2
optSumabsv <- 1.5
optSPC <- 2
timeDiffCutOff <- 5

# STEP 2: Train a model on white- and profiled DGRPs --------------------------
# fit a periodic smoothing spline to the behavior of each feature as a function
# of time
fitWhiteDgrpProf <- zeitzeigerFit(whiteDgrpProf, whiteDgrpProfTime)
# fits to calculate sparse principal components (SPCs) for how the features
# change over time
spcWhiteDgrpProf <- zeitzeigerSpc(fitWhiteDgrpProf$xFitMean, 
                                  fitWhiteDgrpProf$xFitResid, 
                                  sumabsv = optSumabsv)
# uses the training data and the SPCs to predict the corresponding time for
# each test observation
predfitWhiteDgrpProf <- zeitzeigerPredict(as.matrix(whiteDgrpProf), 
                                          whiteDgrpProfTime,
                                          as.matrix(whiteDgrpProf),
                                          spcWhiteDgrpProf, nSpc = optSPC)
# compare predicted and harversted time
harvPredTime <- data.frame(timeObs = whiteDgrpProfTime,
                           timePred = predfitWhiteDgrpProf$timePred,
                           timeError = calcTimeDiff(whiteDgrpProfTime, 
                                                    predfitWhiteDgrpProf$timePred))
plot(24 * harvPredTime[, 1:2], pch = 20, bty = 'n',
     col = ifelse(abs(harvPredTime$timeError) > timeDiffCutOff / 24, 
                  tissueColor[tissue], 'black'), cex = 2, 
     xlim = c(0, 25), ylim = c(0, 25), cex.lab = 2, cex.axis = 2,
     xlab = 'Time of harversting, h', ylab = 'Physiological time, h',
     main = paste0('Physiological and harversting time, ', tissue))
abline(1, 1, lty = 2)
text(24 * harvPredTime[abs(harvPredTime$timeError) > timeDiffCutOff / 24, 1:2],
     labels = rownames(whiteDgrpProf)[abs(harvPredTime$timeError) > 
                                        timeDiffCutOff / 24],
     pos = 3)
grid()

ggplot(harvPredTime, aes(x = 24 * timeObs, y = 24 * timeError)) +
  geom_point(size = 2, shape = 1) +
  scale_x_continuous(limits = 24 * c(0, 1)) +
  labs(x = 'Observed time', y = 'Error') + theme_bw()

# plot proporion of explained variance by each of SPC
dfVar <- data.frame(spc = 1:length(spcWhiteDgrpProf$d),
                    propVar = spcWhiteDgrpProf$d^2 / sum(spcWhiteDgrpProf$d^2))
ggplot(dfVar) + geom_point(aes(x = spc, y = propVar), size = 2, shape = 1) +
  scale_x_continuous(breaks = seq(1, 10)) +
  labs(x = 'SPC', y = 'Proportion of\nvariance explained') + theme_bw()

# genes which are in SPCs
genesPredictors <- colnames(whiteDgrpProf)[which(rowSums(spcWhiteDgrpProf$v[, 1:optSPC]) != 0)]
message(paste("Genes-predictors:", paste(genesPredictors, collapse = ',')))
spcNot0Matr <- spcWhiteDgrpProf$v[rowSums(spcWhiteDgrpProf$v[, 1:optSPC]) != 0, ][, 1:optSPC]
sapply(1:optSPC, 
       function(x) message(paste('SPC', x, ':', 
                           paste(genesPredictors[which(spcNot0Matr[, x] != 0)],
                               collapse = ', '))))
# See, which genes drive the predictions: project the observations from 
# feature-space to SPC-space, to look at how the optimal SPCs behave over time
spcMatr <- as.matrix(whiteDgrpProf) %*% spcWhiteDgrpProf$v[, 1:optSPC]
colnames(spcMatr) <- paste('SPC', 1:optSPC)
zGath <- gather(data.frame(spcMatr, obs = rownames(whiteDgrpProf), 
                           Time = whiteDgrpProfTime, check.names = F),
                key = SPC, value = Abundance, -obs, -Time)
ggplot(zGath) + facet_grid(SPC ~ ., scales = 'free_y') + 
  geom_point(aes(x = Time, y = Abundance), size = 2, shape = 1) + 
  ylab('Gene(s) expression') + theme_bw()
# Plot the coefficients of the genes for the SPCs
spcNot0Matr <- as.data.frame(spcNot0Matr)
colnames(spcNot0Matr) <- paste('SPC', 1:optSPC)
spcNot0Matr$Gene <- genesPredictors
spcNot0MatrGath <- gather(spcNot0Matr, key = spc, 
                          value = Coefficient, -Gene) %>% 
                   mutate(feature = factor(Gene, levels = rev(spcNot0Matr$Gene)))
ggplot(spcNot0MatrGath) + facet_wrap(~ spc, nrow = 1) + 
  geom_bar(aes(x = Gene, y = Coefficient), stat = 'identity') +
  labs(x = 'Gene') + coord_flip() +
  theme_bw() + theme(panel.spacing = unit(1.2, 'lines'))

# STEP 3: Predict time for DGRPs profiled once --------------------------------
optSumabsv <- 1.5
optSPC <- 2
timeDiffCutOff <- 5

whiteDGRP <- readRDS(paste0(saveRDSdir, 'whiteDGRP.Rds'))

for (tissue in allTissues) {
  whiteDGRPtiss <- whiteDGRP[[tissue]]
  whiteDGRPtiss <- whiteDGRPtiss[, !colnames(whiteDGRPtiss) %in% "" ]# tissueOutliers[[tissue]]]
  
  # separate into the training set + test set : white and profiled DGRPs
  # and set of interest - DGRPs profiled once
  whiteDgrpProftiss <- whiteDGRPtiss[, !grepl('_n_', colnames(whiteDGRPtiss))]
  newRownames <- colnames(whiteDgrpProftiss)
  whiteDgrpProftiss <- as.data.frame(t(whiteDgrpProftiss))
  rownames(whiteDgrpProftiss) <- newRownames
  
  dgrpNormtiss <- whiteDGRPtiss[, grepl('_n_', colnames(whiteDGRPtiss))]
  newRownames <- colnames(dgrpNormtiss)
  dgrpNormtiss <- as.data.frame(t(dgrpNormtiss))
  rownames(dgrpNormtiss) <- newRownames
  
  # time for both sets
  whiteDgrpProfTissTime <- (INFO[rownames(whiteDgrpProftiss)]$Time %% 24 ) / 24
  dgrpNormTissTime <- (INFO[rownames(dgrpNormtiss)]$Time %% 24 ) / 24
  
  # fit a periodic smoothing spline to the behavior of each feature as a function
  # of time
  fitWhiteDgrpProfTiss <- zeitzeigerFit(whiteDgrpProftiss, 
                                        whiteDgrpProfTissTime)
  # fits to calculate sparse principal components (SPCs) for how the features
  # change over time
  spcWhiteDgrpProfTiss <- zeitzeigerSpc(fitWhiteDgrpProfTiss$xFitMean, 
                                        fitWhiteDgrpProfTiss$xFitResid, 
                                        sumabsv = optSumabsv)
  # uses the training data and the SPCs to predict the corresponding time for
  # each test observation
  predfitDGRPsTiss <- zeitzeigerPredict(as.matrix(whiteDgrpProftiss), 
                                        whiteDgrpProfTissTime,
                                        as.matrix(dgrpNormtiss),
                                        spcWhiteDgrpProfTiss, nSpc = optSPC)
  # compare predicted and harversted time
  harvPredTimeDGRPsTiss <- data.frame(Sample = rownames(dgrpNormtiss), 
                                      timeObs = 24 * dgrpNormTissTime,
                                      timePred = 24 * predfitDGRPsTiss$timePred)
  harvPredTimeDGRPsTiss$timeError <- apply(harvPredTimeDGRPsTiss, 1, 
                                           function(x) diffInTime(as.numeric(x[2]), 
                                                                  as.numeric(x[3])))
  plot(harvPredTimeDGRPsTiss[, 2:3], pch = 20, bty = 'n',
       col = ifelse(abs(harvPredTimeDGRPsTiss$timeError) > timeDiffCutOff, 
                    tissueColor[tissue], 'black'), cex = 2, 
       xlim = c(0, 25), ylim = c(0, 25), cex.lab = 2, cex.axis = 2,
       xlab = 'Time of harversting, h', ylab = 'Physiological time, h',
       main = paste0('Physiological and harversting time, ', tissue))
  abline(1, 1, lty = 2)
  text(harvPredTimeDGRPsTiss[abs(harvPredTimeDGRPsTiss$timeError) > 
                                    timeDiffCutOff, 2:3],
       labels = rownames(dgrpNormtiss)[abs(harvPredTimeDGRPsTiss$timeError) > 
                                       timeDiffCutOff],
       pos = 3)
  grid()
  
  saveRDS(harvPredTimeDGRPsTiss, paste0(saveRDSdir, tissue, 
                                        '_physioPredTimeDF_ZZ.Rds'))
}