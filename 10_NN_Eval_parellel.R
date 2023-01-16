#!/usr/bin/env Rscript
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

#setwd('~/Desktop/AroundTheClock_Aug_Oct2017/')
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

cyclingGenes <- lapply(cyclingGenes, function(x) x[x$ADJ.P < 0.05, ])
cyclingGenes <- lapply(cyclingGenes, 
                       function(x) x[order(-x$AMP, x$ADJ.P), ])

args <- commandArgs(trailingOnly = T)
tissue <- args[1]
outliersName <- args[2]
neuronName <- args[3]
message(paste('Parameters:', outliersName, tissue, neuronName))
message(getwd())

# NEURAL NET ------------------------------------------------------------------
neuralNetPredTimeDF <- c()
#for (outliersName in names(tissueOutliers)) {
#  for (tissue in allTissues) {
    outliers <- tissueOutliers[[outliersName]][[tissue]]
    colorForTissue <- tissueColor[tissue]
    INFOTissue <- INFO[Tissue == tissue]
    setkey(INFOTissue, rightGT_name)
    
    cyclGenesTiss <- cyclingGenes[[tissue]]
    
    message(tissue)
    message(names(whiteDgrpProf))
    lapply(whiteDgrpProf, function(x) x[1:6, 1:6])
    whiteDgrpProfnoOut <- whiteDgrpProf[[as.character(tissue)]]
    
    head(whiteDgrpProfnoOut)
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[, colnames(whiteDgrpProfnoOut) %in% 
                                               rownames(cyclGenesTiss)]
    head(whiteDgrpProfnoOut)
    whiteDgrpProfnoOut <- whiteDgrpProfnoOut[!rownames(whiteDgrpProfnoOut) %in%
                                               outliers, ]
    head(whiteDgrpProfnoOut)
    neurons <- list(n12_6_3 = c(12, 6, 3), n12_6 = c(12, 6), n12 = 12)
    #for (neuronName in names(neurons)) {
      hiddenNeurons <- neurons[[neuronName]]
      for (i in 1:nrow(whiteDgrpProfnoOut)) {
        message(paste0('NEURAL NET ', tissue, ' ', neuronName, ': ', i))
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
#    }
#  }
#}

saveRDS(neuralNetPredTimeDF, 
        paste0(saveRDSdir, 'NeuralNet_realPredTimeDF_no25205', 
               '_', outliersName, '_', tissue, '_', neuronName, '.Rds'))
