# FILE DESCRIPTION ------------------------------------------------------------
#
# FILE: 14_compare_w_28233.R
#
# DESCRIPTION : compares cycling genes between w- and 28233
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
# INFO about 28233 samples
info28233 <- readRDS('Rds/sampleInfo_28233.Rds')
info28233 <- as.data.table(info28233)

# INFO about W samples
infoW <- readRDS('Rds/sampleInfo.Rds')
infoW <- as.data.table(infoW)
infoW <- infoW[Tissue %in% info28233$Tissue]
infoW <- infoW[initial_GT == 'w-']

INFO <- rbind(info28233, infoW)

# list all tissues and corresponding colors
allTissues <- unique(INFO$Tissue)
tissueColor <- c('darkorchid3', 'darkgreen')
names(tissueColor) <- allTissues

# output dirs
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'
saveBedDir <- 'beds/'

# cycling genes detection results from JTK
cyclGenes28233 <- readRDS(paste0(saveRDSdir, 'cyclGenesList_28233.Rds'))
# selecting significant ones 
cyclGenesList28233 <- lapply(cyclGenes28233, function(x) x[x$ADJ.P < 0.05, ])
# cycling genes detection results from JTK
cyclGenesW <- readRDS(paste0(saveRDSdir, 'cyclGenesList.Rds'))
# selecting significant ones 
cyclGenesListW <- lapply(cyclGenesW[names(cyclGenes28233)], 
                         function(x) x[x$ADJ.P < 0.05, ])

# LOAD COUNTS --------------------------------------------------------------
# The tables will contain both 28233 and w-
# not standartisezed
white28233noStand <- readRDS(paste0(saveRDSdir, 'white28233noStand.Rds'))
# standartisezed
white28233 <- readRDS(paste0(saveRDSdir, 'white28233.Rds'))

# remove crosses from count tables
white28233noStand <- lapply(white28233noStand, 
                            function(x) x[, colnames(x) %in%
                                            INFO[right_GT == '28233' |
                                                 right_GT == 'w-', ]$rightGT_name])
white28233 <- lapply(white28233,
                     function(x) x[, colnames(x) %in%
                                     INFO[right_GT == '28233' |
                                          right_GT == 'w-', ]$rightGT_name])

# COMPARE NUMBER OF DETECTED CYCLING GENES ------------------------------------
detectCyclGenesNumb <- data.frame(Genotype = c(rep('white-', 
                                                   length(allTissues)), 
                                               rep('28233', 
                                                   length(allTissues))),
                                  Tissue = rep(allTissues, 2))
detectCyclGenesNumb$detectCyclGenes <- c(sapply(allTissues,
                                               function(x) nrow(cyclGenesListW[[x]])),
                                        sapply(allTissues,
                                               function(x) nrow(cyclGenesList28233[[x]])))
png(paste0(savePlotsDir, 'cyclGenes_white_vs_28233.png'), height = 800,
    width = 900)
print(ggplot(detectCyclGenesNumb, aes(x = Tissue, y = detectCyclGenes)) +
      geom_bar(aes(fill = Genotype), position = "dodge", stat = "identity") +
      ylab("Cycling genes") + ggtitle("Number of detected cycling genes")+
      scale_fill_manual(values = c('hotpink3', 'lightblue3')) +
      geom_text(aes(label = detectCyclGenes), 
                position = position_dodge(width = .9)) +
      mashaGgplot2Theme)
dev.off()

# COMPARE CORE CIRCADIAN GENES ------------------------------------------------
coreCircGenes <- c('tim', 'per', 'Clk', 'vri', 'Pdp1', 'cry')
for (goi in coreCircGenes) {
  listToPlot <- list()
  for (tiss in allTissues) {
    infoTabTissW <- INFO[Tissue == tiss & right_GT == 'w-', ]
    infoTabTissW <- infoTabTissW[infoTabTissW$rightGT_name %in% 
                                 colnames(white28233[[tiss]]), ]
    setkey(infoTabTissW, rightGT_name)
    oneTissCountsW <- white28233[[tiss]][, infoTabTissW$rightGT_name]
    pl <- plotCountsTimecourse(goi, oneTissCountsW,
                               infoTabTissW[colnames(oneTissCountsW)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for W-')))
    listToPlot[[length(listToPlot) + 1]] <- pl + 
                                            scale_y_continuous(limits = c(0, 1))
  }
  for (tiss in allTissues) {
    infoTabTiss28233 <- INFO[Tissue == tiss  & right_GT == '28233', ]
    infoTabTiss28233 <- infoTabTiss28233[infoTabTiss28233$rightGT_name %in% 
                                           colnames( white28233[[tiss]]), ]
    setkey(infoTabTiss28233, rightGT_name)
    oneTissCounts28233 <- white28233[[tiss]][, infoTabTiss28233$rightGT_name]
    pl <- plotCountsTimecourse(goi, oneTissCounts28233,
                               infoTabTiss28233[colnames(oneTissCounts28233)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for 28233')))
    listToPlot[[length(listToPlot) + 1]] <- pl +
      scale_y_continuous(limits = c(0, 1))
  }
  png(paste0(savePlotsDir, goi, '_CrossTissCycl_28233.png'), height = 1200,
      width = 1200)
  multiplot(plotlist = listToPlot, cols = 2)
  dev.off()
  
  listToPlot <- list()
  for (tiss in allTissues) {
    infoTabTissW <- INFO[Tissue == tiss & right_GT == 'w-', ]
    infoTabTissW <- infoTabTissW[infoTabTissW$rightGT_name %in% 
                                   colnames(white28233noStand[[tiss]]), ]
    setkey(infoTabTissW, rightGT_name)
    oneTissCountsW <- white28233noStand[[tiss]][, infoTabTissW$rightGT_name]
    
    pl <- plotCountsTimecourse(goi, oneTissCountsW,
                               infoTabTissW[colnames(oneTissCountsW)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for W-')))
    listToPlot[[length(listToPlot) + 1]] <- pl
    
  }
  for (tiss in allTissues) {
    infoTabTiss28233 <- INFO[Tissue == tiss  & right_GT == '28233', ]
    infoTabTiss28233 <- infoTabTiss28233[infoTabTiss28233$rightGT_name %in% 
                                           colnames( white28233noStand[[tiss]]), ]
    setkey(infoTabTiss28233, rightGT_name)
    oneTissCounts28233 <- white28233noStand[[tiss]][, infoTabTiss28233$rightGT_name]
    pl <- plotCountsTimecourse(goi, oneTissCounts28233,
                               infoTabTiss28233[colnames(oneTissCounts28233)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for 28233')))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
  png(paste0(savePlotsDir, goi, '_CrossTissCyclNoStand_28233.png'), 
      height = 1200, width = 1200)
  multiplot(plotlist = listToPlot, cols = 2)
  dev.off()
}

# OVERLAP BETWEEN CYCLING IN W AND 28233  -------------------------------------
# plot Venn diagram of overlap
forVenn28233 <- lapply(cyclGenesList28233, function(x) unique(rownames(x)))
names(forVenn28233) <- paste0(names(forVenn28233), '_28233')
forVennW <- lapply(cyclGenesListW, function(x) unique(rownames(x)))
names(forVennW) <- paste0(names(forVennW), '_w-')
forVenn <- c(forVenn28233, forVennW)
png(paste0(savePlotsDir, 'cycl_w_vs_28233.png'), width = 800, height = 800)
venn(forVenn, ilab = T, zcolor = "style", cexil = 1.3, cex.lab = 2)
dev.off()

# get genes in overlaps between tissues
cyclOvrl <- Venn(forVenn)@IntersectionSets
# remove redundunt element
cyclOvrl$`0000` <- NULL
tissueIndex <- sapply(names(cyclOvrl), 
                      function(x) which(strsplit(x, '')[[1]] != "0"))
# rename, so names contain tissue
names(cyclOvrl) <- sapply(1:length(cyclOvrl), 
                          function(x) paste0(names(forVenn)[tissueIndex[[x]]], 
                                             collapse = '-'))
# translate them to flyBase IDs and do GO enrichment
cyclOvrlGO <- list()
# BG for GO
backGr <- lapply(white28233, function(x) getFlyBaseID(rownames(x))[[1]])
for (i in 1:length(cyclOvrl)) {
  if (length(cyclOvrl[[i]]) >= 10) {
    fbId <- getFlyBaseID(cyclOvrl[[i]])[[1]]
    # select proper tissue BG
    fbIdBG <- sapply(names(backGr), function(x) grepl(x, names(cyclOvrl)[1]))
    if (sum(fbIdBG) == 2) {
      fbIdBG <- intersect(backGr[[1]]$flybase_gene_id, 
                          backGr[[2]]$flybase_gene_id)
    } else {
      fbIdBG <- backGr[fbIdBG][[1]]$flybase_gene_id
    }
    cyclOvrlGO[[length(cyclOvrlGO) + 1]] <- GOenricment(fbId$flybase_gene_id,
                                                        fbIdBG, topNodes = 20)
  } else {
    message("Too little genes, no enrichment performed")
  }
}
names(cyclOvrlGO) <- names(cyclOvrl)[lapply(cyclOvrl, length) >= 10]
lapply(cyclOvrlGO, function(x) paste(x$Term, collapse = " "))

# AMPLITUDE BULK COMPARISON FOR 28233 and W -----------------------------------
ampDf <- data.frame(Tissue = character(), Genotype = character(), 
                    AMP = numeric())
for (tissue in allTissues) {
  ampDf <- rbind(ampDf, data.frame(Tissue = tissue, Genotype = 'w-',
                                   AMP = cyclGenesListW[[tissue]]$AMP[1:50]))
  ampDf <- rbind(ampDf, data.frame(Tissue = tissue, Genotype = '28233',
                                   AMP = cyclGenesList28233[[tissue]]$AMP[1:50]))
}

# plot 
png(paste0(savePlotsDir, 'AMP_w_vs_28233.png'), width = 800, height = 800)
ggplot(ampDf, aes(Tissue, AMP), fill = interaction(Genotype, Tissue)) +   
  geom_violin(aes(fill = Genotype), alpha = .5) +
  geom_boxplot(aes(fill = Genotype), alpha = .01, width = .2, 
               position = position_dodge(width = .9)) +
  scale_fill_manual(values = c("#E8F9A2", "darkgreen" )) + 
  ggtitle("Distribution of amplitude of the cycling genes") +
  mashaGgplot2Theme
dev.off()

for (tissue in unique(ampDf$Tissue)) {
  ampDfTissue <- ampDf[ampDf$Tissue == tissue, ]
  print(tissue)
  print(t.test(ampDfTissue[ampDfTissue$Genotype == 'w-', ]$AMP, 
               ampDfTissue[ampDfTissue$Genotype == '28233', ]$AMP)$p.value)
  
}

# PLOT TOP10 CYCLING IN 28233 GENES -------------------------------------------
cyclGenesList28233ord <- lapply(cyclGenesList28233, 
                                function(x) rownames(x[order(-x$AMP), ]))
for (goi in cyclGenesList28233ord$GUT[1:10]) {
  listToPlot <- list()
  for (tiss in allTissues) {
    infoTabTissW <- INFO[Tissue == tiss & right_GT == 'w-', ]
    infoTabTissW <- infoTabTissW[infoTabTissW$rightGT_name %in% 
                                   colnames(white28233[[tiss]]), ]
    setkey(infoTabTissW, rightGT_name)
    oneTissCountsW <- white28233[[tiss]][, infoTabTissW$rightGT_name]
    pl <- plotCountsTimecourse(goi, oneTissCountsW,
                               infoTabTissW[colnames(oneTissCountsW)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for W-')))
    listToPlot[[length(listToPlot) + 1]] <- pl + 
      scale_y_continuous(limits = c(0, 1))
  }
  for (tiss in allTissues) {
    infoTabTiss28233 <- INFO[Tissue == tiss  & right_GT == '28233', ]
    infoTabTiss28233 <- infoTabTiss28233[infoTabTiss28233$rightGT_name %in% 
                                           colnames( white28233[[tiss]]), ]
    setkey(infoTabTiss28233, rightGT_name)
    oneTissCounts28233 <- white28233[[tiss]][, infoTabTiss28233$rightGT_name]
    pl <- plotCountsTimecourse(goi, oneTissCounts28233,
                               infoTabTiss28233[colnames(oneTissCounts28233)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for 28233')))
    listToPlot[[length(listToPlot) + 1]] <- pl +
      scale_y_continuous(limits = c(0, 1))
  }
  #png(paste0(savePlotsDir, goi, '_CrossTissCycl_28233.png'), height = 1200,
  #    width = 1200)
  print(multiplot(plotlist = listToPlot, cols = 2))
  #dev.off()
  
  listToPlot <- list()
  for (tiss in allTissues) {
    infoTabTissW <- INFO[Tissue == tiss & right_GT == 'w-', ]
    infoTabTissW <- infoTabTissW[infoTabTissW$rightGT_name %in% 
                                   colnames(white28233noStand[[tiss]]), ]
    setkey(infoTabTissW, rightGT_name)
    oneTissCountsW <- white28233noStand[[tiss]][, infoTabTissW$rightGT_name]
    
    pl <- plotCountsTimecourse(goi, oneTissCountsW,
                               infoTabTissW[colnames(oneTissCountsW)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for W-')))
    listToPlot[[length(listToPlot) + 1]] <- pl
    
  }
  for (tiss in allTissues) {
    infoTabTiss28233 <- INFO[Tissue == tiss  & right_GT == '28233', ]
    infoTabTiss28233 <- infoTabTiss28233[infoTabTiss28233$rightGT_name %in% 
                                           colnames( white28233noStand[[tiss]]), ]
    setkey(infoTabTiss28233, rightGT_name)
    oneTissCounts28233 <- white28233noStand[[tiss]][, infoTabTiss28233$rightGT_name]
    pl <- plotCountsTimecourse(goi, oneTissCounts28233,
                               infoTabTiss28233[colnames(oneTissCounts28233)]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 'expression in', 
                                             tiss, 'for 28233')))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
  #png(paste0(savePlotsDir, goi, '_CrossTissCyclNoStand_28233.png'), 
  #    height = 1200, width = 1200)
  print(multiplot(plotlist = listToPlot, cols = 2))
  #dev.off()
}
