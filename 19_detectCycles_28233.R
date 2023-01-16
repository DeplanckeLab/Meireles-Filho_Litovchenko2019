# FILE DESCRIPTION ------------------------------------------------------------
#
# FILE: 13_detectCycles.R
#
# DESCRIPTION : 28233 specific
# Detects cycling genes in every tissue and saves the results in Rds object
# Identifyes tissue-specific cycling genes, plots venn diagram for that
# Looks at expression of tissue specific cycling genes in other tissues
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
#  CREATED:  20.01.2018
#  REVISION: 20.01.2018

# INPUTS ----------------------------------------------------------------------
setwd('~/Desktop/AroundTheClock_Aug_Oct2017/')
source('0_common_functions.R')
# info about samples
INFO <- readRDS('Rds/sampleInfo_28233.Rds')
INFO <- as.data.table(INFO)

# list all tissues and corresponding colors
allTissues <- unique(INFO$Tissue)
tissueColor <- c('darkorchid3', 'darkgreen')
names(tissueColor) <- allTissues

# output dirs
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'
saveBedDir <- 'beds/'

# DETECT CYCLING GENES BY JTK in ALL TISSUES ----------------------------------
# not standartisezed
noStand28233 <- readRDS(paste0(saveRDSdir, 'noStand28233.Rds'))
# standartisezed
stand28233 <- readRDS(paste0(saveRDSdir, 'stand28233.Rds'))

cyclGenesList <- list() 
for (tissue in allTissues) {
  countsOneTiss <- stand28233[[tissue]]
  # prepare info table
  INFOTissue <- INFO[right_GT == '28233' & Tissue == tissue]
  setkey(INFOTissue, rightGT_name)
  
  INFOTissue <- INFOTissue[rightGT_name %in% colnames(countsOneTiss)]
  INFOTissue <- INFOTissue[order(Time)]
  countsOneTiss <- countsOneTiss[, INFOTissue$rightGT_name]
  
  # get cycling genes
  cyclGenes <- RunJTKCycle(countsOneTiss, table(INFOTissue$Time), 2)
  
  # save all results
  cyclGenesList[[length(cyclGenesList) + 1]] <- cyclGenes
}
names(cyclGenesList) <- allTissues
saveRDS(cyclGenesList, paste0(saveRDSdir, 'cyclGenesList_28233.Rds'))
# restrict to only genes passing p-value cutoff
cyclGenesList <- readRDS(paste0(saveRDSdir, 'cyclGenesList_28233.Rds'))
cyclGenesList <- lapply(cyclGenesList, function(x) x[x$ADJ.P < 0.05, ])

# DISTRIBUTION OF AMPLITUDE ---------------------------------------------------
# plot distribution of amplitude
for (tissue in allTissues) {
  png(paste0(savePlotsDir, tissue, '_cyclAmpDistr_28233.png'), width = 800,
      height = 800)
  hist(cyclGenesList[[tissue]]$AMP, xlab = 'Amplitude', ylab = 'Frequency', 
       main = paste('Amplitudes for the cycling genes in', tissue),
       bty = 'n', breaks = 25, col = tissueColor[tissue])
  grid()
  dev.off()
}
# plot typical gene with amplitude falling into quantile
for (tissue in allTissues) {
  # prepare info tab and counts
  countsOneTiss <- stand28233[[tissue]]
  # prepare info table
  INFOTissue <- INFO[right_GT == '28233' & Tissue == tissue]
  setkey(INFOTissue, rightGT_name)
  
  INFOTissue <- INFOTissue[rightGT_name %in% colnames(countsOneTiss)]
  INFOTissue <- INFOTissue[order(Time)]
  setkey(INFOTissue, rightGT_name)
  countsOneTiss <- countsOneTiss[, INFOTissue$rightGT_name]
  
  tisAmps <- cyclGenesList[[tissue]]$AMP
  ampQuant <- quantile(tisAmps, seq(0, 1, length.out = 10))
  listToPlot <- list()
  for (i in 2:length(ampQuant)) {
    quantMean <- mean(c(ampQuant[i], ampQuant[i - 1]))
    goiInd <- which(abs(tisAmps - quantMean) == min(abs(tisAmps - quantMean)))
    goi <- rownames(cyclGenesList[[tissue]])[goiInd]
    pl <- plotCountsTimecourse(goi, countsOneTiss,
                               INFOTissue[colnames(countsOneTiss)]$Time,
                               tissueColor[tissue],
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                              tissue)))
    pl <- pl + annotate("text", x = 10, y = 0.02, size = 8, 
                        label =   paste('AMP in range ', round(quantMean, 2)))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
  png(paste0(savePlotsDir, tissue, '_ampQuant_28233.png'), width = 1500, 
      height = 1500)
  multiplot(plotlist = listToPlot, cols = 3)
  dev.off()
}

# TISSUE SPECIFIC GENE CYCLING ------------------------------------------------
# plot Venn digram
png(paste0(savePlotsDir, 'tissueSpecCycl_28233.png'), width = 800, 
    height = 800)
forVenn <- lapply(cyclGenesList, function(x) unique(rownames(x)))
names(forVenn) <- names(cyclGenesList)
venn(forVenn, ilab = T, zcolor = "style", cexil = 1.3, cex.lab = 2)
dev.off()

# get overlaps between tissues
cyclOvrl <- Venn(forVenn)@IntersectionSets
# remove redundunt element
cyclOvrl$`0000` <- NULL
tissueIndex <- sapply(names(cyclOvrl), 
                      function(x) which(strsplit(x, '')[[1]] != "0"))
# rename, so names contain tissue
names(cyclOvrl) <- sapply(1:length(cyclOvrl), 
                          function(x) paste0(names(forVenn)[tissueIndex[[x]]], 
                                             collapse = '-'))

# PLOT GENES CYCLING IN ALL TISSUES -------------------------------------------
# plot normalized and standartizized counts
cycleAll4Tiss <- cyclOvrl$`BRAIN-GUT`
for (goi in cycleAll4Tiss) {
  listToPlot <- list()
  for (tiss in allTissues) {
    INFOTiss <- INFO[Tissue == tiss]
    setkey(INFOTiss, rightGT_name)
    pl <- plotCountsTimecourse(goi, stand28233[[tiss]], 
                               INFOTiss[colnames(stand28233[[tiss]])]$Time,
                               tissueColor[tiss],
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                      tiss)))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
  png(paste0(savePlotsDir, goi, '_CrossTissCycl_28233.png'), height = 425,
      width = 850)
  multiplot(plotlist = listToPlot, cols = 2)
  dev.off()
  
  listToPlot <- list()
  for (tiss in allTissues) {
    INFOTiss <- INFO[Tissue == tiss]
    setkey(INFOTiss, rightGT_name)
    pl <- plotCountsTimecourse(goi, noStand28233[[tiss]], 
                               INFOTiss[colnames(noStand28233[[tiss]])]$Time,
                               tissueColor[tiss],
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss)))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
  png(paste0(savePlotsDir, goi, '_CrossTissCyclNoStand_28233.png'), 
      height = 425, width = 850)
  multiplot(plotlist = listToPlot, cols = 2)
  dev.off()
}