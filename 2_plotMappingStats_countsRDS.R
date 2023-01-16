# FILE DESCRIPTION ------------------------------------------------------------
#
# FILE: 2_plotMappingStats_countsRDS.R
#
# DESCRIPTION : 
# QC check for RBR-seq: 
# 1) Number of RAW reads in each BRB libraries,
# 2) Percentage of demultiplexed reads in each BRB libraries
# 3) Distribution of reads between libraries
# 4) Percentage of uniqly mapped reads
# 
# Creating Rds with lists of counts for the following scripts:
# 1) list of tables for white- counts only, not standartized
# 2) list of tables for white- counts only, standartized - use it for the
#    detection of cycling genes
# 3) list of tables for white-&dgrp counts, not standartized
# 4) list of tables for white-&dgrp counts, standartized - use it for the
#    estimation of physiological time for dgrps
#
# USAGE: 
#
#  OPTIONS:  none
#  REQUIREMENTS:  ggplot2, data.table
#  BUGS: --
#  NOTES:  ---
#  AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
#  COMPANY:  EPFL, Lausanne, Switzerland
#  VERSION:  1
#  CREATED:  18.10.2017
#  REVISION: 18.10.2017

# Inputs ----------------------------------------------------------------------
source('0_common_functions.R')
savePlotsTo <- 'plots/'
saveRdsTo <- 'Rds/'
INFO <- readRDS('Rds/sampleInfo.Rds')

# PLOT STATS: MAPPED/UNMAPPED/DUPLICATED/DEDUPLICATED READS -------------------
# create df for plotting stacked bar plot
for (mon in unique(INFO$Month)) {
  statsMon <- INFO[Month == mon]
  statsMon[, Name := initial_name]
  for (tis in unique(statsMon$Tissue)) {
    png(paste0(savePlotsTo, '0_A_mapped_numb_', mon, '_', tis, '.png'),
        width = 1500, height = 750)
    print(plotMappedUnmapped(meltMappedUnmapped(statsMon[Tissue == tis])))
    dev.off()
    
    png(paste0(savePlotsTo, '0_B_mapped_perc_', mon, "_", tis, '.png'),
        width = 1500, height = 750)
    print(plotMappedUnmapped(meltMappedUnmapped(statsMon[Tissue == tis]), T))
    dev.off()
    
    png(paste0(savePlotsTo, '0_C_dedupl_perc_', mon, '_', tis, '.png'),
        width = 1500, height = 750)
    print(plotDupl(meltDupl(statsMon[Tissue == tis]), percentage = T))
    dev.off()
  }
}

# Create Rds with per tissue count tabs: white solo ---------------------------
allTissues <- unique(INFO$Tissue)
countsDirs <- c('~/Documents/AroundTheClock_Oct2017_white/2_counts/',
                '~/Documents/AroundTheClock_Aug2017/9_counts_GT_CORRECT/')
names(countsDirs) <- c('WHITE', 'DGRPs')

# not standartisezed
whiteNoStand <- lapply(allTissues, readCountsNormBatchCorr,countsDirs['WHITE'],
                       INFO, normWithVoom = T, combatBatch = T, 
                       standartize = F)
names(whiteNoStand) <- allTissues
saveRDS(whiteNoStand, paste0(saveRdsTo, 'whiteNoStand.Rds'))
# standartisezed
white <- lapply(allTissues, readCountsNormBatchCorr,countsDirs['WHITE'],
                INFO, normWithVoom = T, combatBatch = T, standartize = T)
names(white) <- allTissues
saveRDS(white, paste0(saveRdsTo, 'white.Rds'))

# Create Rds with per tissue count tabs: white solo with low expr genes -------
# Used for plotting purposes only
allTissues <- unique(INFO$Tissue)
countsDirs <- c('~/Documents/AroundTheClock_Oct2017_white/2_counts/',
                '~/Documents/AroundTheClock_Aug2017/9_counts_GT_CORRECT/')
names(countsDirs) <- c('WHITE', 'DGRPs')

# not standartisezed
whiteNoStand <- lapply(allTissues, readCountsNormBatchCorr,countsDirs['WHITE'],
                       INFO, normWithVoom = T, combatBatch = T, 
                       standartize = F, geneCut = 0)
names(whiteNoStand) <- allTissues
saveRDS(whiteNoStand, paste0(saveRdsTo, 'whiteNoStand_withLowGenes.Rds'))
# standartisezed
white <- lapply(allTissues, readCountsNormBatchCorr,countsDirs['WHITE'],
                INFO, normWithVoom = T, combatBatch = T, standartize = T,
                geneCut = 0)
names(white) <- allTissues
saveRDS(white, paste0(saveRdsTo, 'white_withLowGenes.Rds'))

# Create Rds with per tissue count tabs: white and DGRPs ----------------------
allTissues <- unique(INFO$Tissue)
countsDirs <- c('~/Documents/AroundTheClock_Oct2017_white/2_counts/',
                '~/Documents/AroundTheClock_Aug2017/9_counts_GT_CORRECT/')
names(countsDirs) <- c('WHITE', 'DGRPs')
# not standartisezed
whiteDGRPnoStand <- lapply(allTissues, readCountsNormBatchCorr, countsDirs,  
                           INFO, normWithVoom = T, combatBatch = T,
                           standartize = F)
names(whiteDGRPnoStand) <- allTissues
saveRDS(whiteDGRPnoStand, paste0(saveRdsTo, 'whiteDGRPnoStand.Rds'))
# standartisezed
whiteDGRP <- lapply(allTissues, readCountsNormBatchCorr, countsDirs, INFO, 
                    normWithVoom = T, combatBatch = T, standartize = T) 
names(whiteDGRP) <- allTissues
saveRDS(whiteDGRP, paste0(saveRdsTo, 'whiteDGRP.Rds'))

# White- for DE ---------------------------------------------------------------
allTissues <- unique(INFO$Tissue)
countsDir <- '~/Documents/AroundTheClock_Oct2017_white/2_counts/'

# restrict only to samples (time points) which are presented in all 4 tissues
INFO_w <- INFO[right_GT == 'w-']
setkey(INFO_w, Time_code, Month)
numbOfSampPerTimePoint <- INFO_w[,.N, by = .(Tissue, Time_code)]
numbOfSampPerTimePoint <- numbOfSampPerTimePoint[,.N, by = .(Time_code)]
numbOfSampPerTimePoint <- numbOfSampPerTimePoint[N == 4]
INFO_w <- INFO_w[Time_code %in% numbOfSampPerTimePoint$Time_code]

# now INFO_w contains the samples we need. However, some samples could 
# still be excluded, because they will not pass the cutoff of 300 000 reads 
# per sample.
allTissues <- unique(INFO_w$Tissue)
whiteDE <- lapply(allTissues, function(x) readCounts(countsDir, x,
                                                     removeNotMapped = T,
                                                     removeWhiteFromDGRPs = F))
names(whiteDE) <- allTissues
# restrict to the samples which are present in all 4 tissues and have coverage
# of > 300 000
for (tissue in allTissues) {
  whiteDE[[tissue]] <- whiteDE[[tissue]][, colnames(whiteDE[[tissue]]) %in%
                                           INFO_w$rightGT_name]
  whiteDE[[tissue]] <- cutoffCoverage(whiteDE[[tissue]], 300000)
}
# let's see, which time points still have 4 tissues presented
INFO_w <- INFO_w[rightGT_name %in% unlist(sapply(whiteDE, colnames)), ]
numbOfSampPerTimePoint <- INFO_w[,.N, by = .(Tissue, Time_code)]
numbOfSampPerTimePoint <- numbOfSampPerTimePoint[,.N, by = .(Time_code)]
numbOfSampPerTimePoint <- numbOfSampPerTimePoint[N == 4]
INFO_w <- INFO_w[Time_code %in% numbOfSampPerTimePoint$Time_code]
whiteDE <- lapply(whiteDE, 
                    function(x) x[, colnames(x) %in% INFO_w$rightGT_name])
# number of samples isn't the same for every tissue, but it's okay, the most 
# important is that all 4 tissues are present per time point. Some of the 
# tissues will have more or less replicates

# remove genes which are not expressed in any tissue
notExprGenes <- lapply(whiteDE, 
                       function(x) rownames(x)[apply(x, 1, 
                                                     function(y) sum(y == 0)) == ncol(x)])
notExprGenes <- unlist(notExprGenes)
# it should be low in all tissues to be excluded
notExprGenes <- table(notExprGenes)
notExprGenes <- names(notExprGenes[notExprGenes == 4])
whiteDE <- lapply(whiteDE, function(x) x[!rownames(x) %in% notExprGenes, ])

# save not normalized counts because DESeq2 loves so
saveRDS(whiteDE, paste0(saveRdsTo, 'whiteDE.Rds'))