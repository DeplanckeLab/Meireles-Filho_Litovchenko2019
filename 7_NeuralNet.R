# FILE DESCRIPTION: 8_NeuralNet.R ---------------------------------------------
#
# DESCRIPTION : tries to predict physiological time with use of neural net
# https://www.r-bloggers.com/fitting-a-neural-network-in-r-neuralnet-package/
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

setwd("~/Desktop/AroundTheClock_Aug_Oct2017/")
setwd("~/Desktop/BitBucket/AroundTheClock_Aug_Oct2017/")
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
# many methods require matrix with samples in rows and genes in columns, here 
# you go
whiteDgrpProf <- lapply(whiteDgrpProf, function(x) as.data.frame(t(x)))

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
                                 'BRAIN28_Oct_1', 'BRB7_Brain_ZT21_28211'),
                       MT = c())

# -----------------------------------------------------------------------------
aaa <- c()
tissue <- 'BRAIN'
#for (tissue in allTissues) {
  outliers <- "" #tissueOutliers[[tissue]]
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
#  for (i in 1:nrow(whiteDgrpProfnoOut)) {
  for (i in 1:10) {
    message(paste0(tissue, ': ', i))
    # divide into train and test dataset
    trainTestSets <- divideIntoTrainAndTest(whiteDgrpProfnoOut, INFOTissue,
                                            leaveOneOut = i, genesOI = manualSel)
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
                                       #hidden = c(7, 5, 3),
                                       #hidden = c(5, 3),
                                       #hidden = c(12, 6),
                                       #hidden = c(12, 6, 3),
                                       hidden = 12,
                                       linear.output = T)
                            preditNN <- compute(neuralNet, testSet)}, 
                            warning = function(w) {#print('WARNING! trying again')
                                                   neuralNet = NULL},
                            error = function(e) {#print('ERROR! trying again')
                                                 neuralNet = NULL})
      #neuralNet <- neuralnet(modelFormula, data = trainSet, hidden = c(4, 2), 
      #                       linear.output = T)
    }
    aaa <- rbind(aaa, c(INFOTissue[rownames(trainTestSets$test)]$Time %% 24,
                        preditNN$net.result[, 1]))
    neuralNet <- NULL
    preditNN <- NULL
    #plot(neuralNet)
  }
#}

manualSel <- c('CG17562', 'vri', 'tim', 'CG18609', 'CG9459', 'CG17560', 'geko', 
               'CG17323', 'Clk', 'CG14893', 'CG32369')
#manualSel <- c('CG17562', 'vri', 'tim', 'Clk')
#manualSel <- rownames(white)