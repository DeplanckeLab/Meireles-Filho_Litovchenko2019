# FILE DESCRIPTION ------------------------------------------------------------
#
# FILE: 3_detectCycles.R
#
# DESCRIPTION : 
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
#  CREATED:  18.10.2017
#  REVISION: 18.10.2017
setwd('~/Desktop/AroundTheClock_Aug_Oct2017/')
source('0_common_functions.R')

# INPUTS ----------------------------------------------------------------------
# info about samples
INFO <- readRDS('Rds/sampleInfo.Rds')
INFO <- as.data.table(INFO)
# list all tissues and corresponding colors
allTissues <- unique(INFO$Tissue)
tissueColor <- c('#688B8A', '#A0B084', '#FAEFD4', '#A57C65')
names(tissueColor) <- allTissues

# peaks from Antonio's Cell paper
peaksDir <- 'inputs_for_Rscripts/'
headSpec <- readAndAnnotPeaks(paste0(peaksDir, "head_specific_peaks_dm3.bed"))
bodySpec <- readAndAnnotPeaks(paste0(peaksDir, "body_specific_peaks_dm3.bed"))
sharedPeaks <- readAndAnnotPeaks(paste0(peaksDir, "shared_peaks_dm3.bed"))
greyZone <- readAndAnnotPeaks(paste0(peaksDir, "grey_zone_peaks_dm3.bed"))

# output dirs
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'
saveBedDir <- 'beds/'

# counts - not standartisezed
whiteNoStand <- readRDS(paste0(saveRDSdir, 'whiteNoStand.Rds'))
# counts - standartisezed
white <- readRDS(paste0(saveRDSdir, 'white.Rds'))

# list of genes cycling in different tissues for baboon, table S6 from paper
# Diurnal transcriptome atlas of a primate across major neural and peripheral 
# tissues
baboonCyclPath <- paste0(saveRDSdir, 'aao0318_Mure_SM_Tables-S6.csv') 
baboonCycl <- fread(baboonCyclPath, skip = 1, header = T)
# to see, which genes cycle in all tissues:
# head(sort(table(unlist(baboonCycl)), decreasing = T), 10) / ncol(baboonCycl)

# genes cycling in mouse, derived from http://circadb.hogeneschlab.org/query
mouseCyclPath <- paste0(saveRDSdir, 'mouse_circaDB.csv') 
mouseCycl <- fread(mouseCyclPath, header = T)
mouseCycl$Symbol <- gsub('LOC100046232 /// Nfil3', "Nfil3", mouseCycl$Symbol)
# to see, which genes cycle in all tissues:
mouseCycl[,.N/ length(unique(mouseCycl$Tissue)), by = Symbol][order(V1)]

# DETECT CYCLING GENES BY JTK in ALL TISSUES ----------------------------------
cyclGenesList <- list() 
for (tissue in names(white)) {
  countsOneTiss <- white[[tissue]]
  # prepare info table
  INFOtiss <- INFO[right_GT == 'w-' & Tissue == tissue]
  setkey(INFOtiss, rightGT_name)
  
  INFOtiss <- INFOtiss[rightGT_name %in% colnames(countsOneTiss)]
  INFOtiss <- INFOtiss[order(Time)]
  countsOneTiss <- countsOneTiss[, INFOtiss$rightGT_name]
  
  # get cycling genes
  cyclGenes <- RunJTKCycle(countsOneTiss, table(INFOtiss$Time), 2)
  
  # save all results
  cyclGenesList[[length(cyclGenesList) + 1]] <- cyclGenes
}
names(cyclGenesList) <- names(white)
saveRDS(cyclGenesList, paste0(saveRDSdir, 'cyclGenesList.Rds'))

# restrict to only genes passing p-value cutoff
cyclGenesList <- lapply(cyclGenesList, function(x) x[x$ADJ.P < 0.05, ])

# PLOT CLASSICAL HEATMAP OF CYCLING GENES -------------------------------------
for (tissue in names(white)) {
  svg(paste0('plots/circadianHeatmap_', tissue, '.svg'), width = 10, 
      height = 20)
  plotCircadianHeatmap(cyclGenesList[[tissue]], white[[tissue]],
                       INFO[Tissue == tissue], tissueColor[tissue],
                       main = paste('All detected circadian genes in', tissue))
  dev.off()
}

# DISTRIBUTION OF AMPLITUDE ---------------------------------------------------
# plot distribution of amplitude
for (tissue in names(white)) {
  png(paste0(savePlotsDir, tissue, '_cyclAmpDistr.png'), width = 800,
      height = 800)
  hist(cyclGenesList[[tissue]]$AMP, xlab = 'Amplitude', ylab = 'Frequency', 
       main = paste('Amplitudes for the cycling genes in', tissue),
       bty = 'n', breaks = 25, col = tissueColor[tissue])
  grid()
  dev.off()
}
# plot typical gene with amplitude falling into quantile
for (tissue in names(white)) {
  # prepare info tab and counts
  countsOneTiss <- white[[tissue]]
  # prepare info table
  INFOtiss <- INFO[right_GT == 'w-' & Tissue == tissue]
  setkey(INFOtiss, rightGT_name)
  
  INFOtiss <- INFOtiss[rightGT_name %in% colnames(countsOneTiss)]
  INFOtiss <- INFOtiss[order(Time)]
  setkey(INFOtiss, rightGT_name)
  countsOneTiss <- countsOneTiss[, INFOtiss$rightGT_name]
  
  tisAmps <- cyclGenesList[[tissue]]$AMP
  ampQuant <- quantile(tisAmps, seq(0, 1, length.out = 10))
  listToPlot <- list()
  for (i in 2:length(ampQuant)) {
    quantMean <- mean(c(ampQuant[i], ampQuant[i - 1]))
    goiInd <- which(abs(tisAmps - quantMean) == min(abs(tisAmps - quantMean)))
    goi <- rownames(cyclGenesList[[tissue]])[goiInd]
    pl <- plotCountsTimecourse(goi, countsOneTiss,
                               INFOtiss[colnames(countsOneTiss)]$Time, 
                               tissueColor[tissue], 
                               ggtitle(paste('Timeline of', goi,
                                             'expression in',  tissue)))
    pl <- pl + annotate("text", x = 10, y = 0.02, size = 8,
                        label = paste('AMP in range ', round(quantMean, 2)))
    listToPlot[[length(listToPlot) + 1]] <- pl
  }
  png(paste0(savePlotsDir, tissue, '_ampQuant.png'), width = 1500, 
      height = 1500)
  print(multiplot(plotlist = listToPlot, cols = 3))
  dev.off()
}

# TISSUE SPECIFIC GENE CYCLING ------------------------------------------------
# plot Venn digram
png(paste0(savePlotsDir, 'tissueSpecCycl.png'), width = 800, height = 800)
forVenn <- lapply(cyclGenesList, function(x) unique(rownames(x)))
names(forVenn) <- names(cyclGenesList)
forVenn <- forVenn[sort(names(forVenn))]
venn(forVenn, ilab = T, zcolor = "style", cexil = 1.3, cex.lab = 2)
dev.off()

# get overlaps between tissues
cyclOvrl <- Vennerable::Venn(forVenn)@IntersectionSets
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
cycleAll4Tiss <- cyclOvrl$`BRAIN-FB-GUT-MT`
for (goi in cycleAll4Tiss) {
  listToPlot <- list()
  listToPlotNoStand <- list()
  for (tiss in allTissues) {
    INFOTiss <- INFO[Tissue == tiss]
    setkey(INFOTiss, rightGT_name)
    pl <- plotCountsTimecourse(goi, white[[tiss]], 
                               INFOTiss[colnames(white[[tiss]])]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 
                                             'expression in', tiss)))
    listToPlot[[length(listToPlot) + 1]] <- pl
    pl <- plotCountsTimecourse(goi, whiteNoStand[[tiss]], 
                               INFOTiss[colnames(white[[tiss]])]$Time,
                               tissueColor[tiss], 
                               ggtitle(paste('Timeline of', goi, 
                                             'expression in', tiss)))
    listToPlotNoStand[[length(listToPlotNoStand) + 1]] <- pl
  }
  png(paste0(savePlotsDir, goi, '_CrossTissCycl.png'), height = 1000,
      width = 1000)
  print(multiplot(plotlist = listToPlot, cols = 2))
  dev.off()
  png(paste0(savePlotsDir, goi, '_CrossTissCyclNoStand.png'), height = 1000,
      width = 1000)
  print(multiplot(plotlist = listToPlotNoStand, cols = 2))
  dev.off()
}

# CHECK ON HOMOLOGUES IN OTHER SPECIES: MOUSE, BABOON  ------------------------
# get flybase IDs
ensembl <- useMart("ensembl")
baboon <- useDataset("panubis_gene_ensembl", mart = ensembl)
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
mouse <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
fly <- useDataset("dmelanogaster_gene_ensembl", mart = ensembl)

toHuman <- getLDS(attributes = c('external_gene_name', 'flybase_gene_id'), 
                  filters = "external_gene_name", values = cycleAll4Tiss, 
                  mart = fly, martL = human,
                  attributesL = c("ensembl_gene_id", "hgnc_symbol"))
colnames(toHuman) <- c('flyName', 'FlyBaseID', 'HumanID', 'HumanName')
toHuman <- as.data.table(toHuman)

toMouse <- getLDS(attributes = c('external_gene_name', 'flybase_gene_id'),
                  filters = "external_gene_name", values = cycleAll4Tiss,
                  mart = fly, martL = mouse,
                  attributesL = c("ensembl_gene_id", "external_gene_name"))
colnames(toMouse) <- c('flyName', 'FlyBaseID', 'MouseID', 'MouseName')
toMouse <- as.data.table(toMouse)

toBaboon <- getLDS(attributes = c('external_gene_name', 'flybase_gene_id'),
                   filters = "external_gene_name", values = cycleAll4Tiss,
                   mart = fly, martL = baboon,
                   attributesL = c("ensembl_gene_id", "external_gene_name"))
colnames(toBaboon) <- c('flyName', 'FlyBaseID', "BaboonID", 'BaboonName')

toHuman <- as.data.table(toHuman)
toMouse <- as.data.table(toMouse)
toBaboon <- as.data.table(toBaboon)
setkey(toHuman, flyName)
setkey(toMouse, flyName)
setkey(toBaboon, flyName)
allSpecies <- merge(toHuman, toMouse, allow.cartesian = TRUE)
allSpecies <- merge(allSpecies, toBaboon, allow.cartesian = TRUE)
allSpecies[, FlyBaseID.y := NULL]
allSpecies[, FlyBaseID := NULL]
setnames(allSpecies, 'FlyBaseID.x', 'FlyBaseID')
noMatch <- setdiff(cycleAll4Tiss, allSpecies$flyName)
allSpecies <- rbind(allSpecies, 
                    data.frame(flyName = noMatch,
                               FlyBaseID = rep(NA, length(noMatch)),
                               HumanID = rep(NA, length(noMatch)),
                               HumanName = rep(NA, length(noMatch)),
                               MouseID = rep(NA, length(noMatch)),
                               MouseName = rep(NA, length(noMatch)),
                               BaboonID = rep(NA, length(noMatch)),
                               BaboonName = rep(NA, length(noMatch))))
saveRDS(allSpecies, paste0(saveRDSdir, 'cycl4tiss_Orthologs.Rds'))

# check overlap with genes cycling in baboon
cyclInBaboon <- apply(baboonCycl, 2, 
                      function(x) unique(allSpecies[which(allSpecies$BaboonName %in% x), ]$flyName))
cyclInBaboon <- table(unlist(cyclInBaboon))
cyclInBaboon <- data.table(cyclInBaboon)
names(cyclInBaboon) <- c('FlyName', 'NumbTissuesBaboon')
# check overlap in Mouse
mouseCycl <- mouseCycl[,.(length(unique(Tissue))), by = Symbol]
cyclInMouse <- data.table(Symbol = allSpecies$MouseName, 
                          FlyName = allSpecies$flyName)
setkey(mouseCycl, 'Symbol')
setkey(cyclInMouse, 'Symbol')
cyclInMouse <- merge(cyclInMouse, mouseCycl)
cyclInMouse <- cyclInMouse[!duplicated(cyclInMouse), ]
cyclInMouse <- cyclInMouse[,(max(V1)), by = FlyName]
# final table
setkey(cyclInBaboon, FlyName)
setkey(cyclInMouse, FlyName)
merge(cyclInMouse, cyclInBaboon, all = T)

# WHY TISSUE SPECIFIC CYCLING? Check on expression ----------------------------
# create a data table containing gene, tissue where it cycles, quantile of
# cycling amplitude, if it's expressed, quantile of expression
# It should contain ALL expressed genes detected in ALL tissues
tisSpecCycl <- data.table(Gene = unique(unlist(lapply(white, rownames))))
# info about where gene cycles
tisSpecCycl[, CyclesIn := sapply(tisSpecCycl$Gene, 
                                 function(x) names(cyclOvrl)[sapply(cyclOvrl,
                                                                    function(y) x %in% y)])]
tisSpecCycl[, CyclesIn := sapply(CyclesIn, function(x) ifelse(length(x) == 0,
                                                              '', x))]
# add information about the quantile of amplitude
cyclGenesList <- lapply(cyclGenesList, 
                        function(x) as.data.table(cbind(GENE = rownames(x), x)))
cyclGenesList <- lapply(cyclGenesList, 
                        function(x) x[, ampQuant := as.numeric(quantcut(-AMP,
                                                                        q = 5))])
cyclGenesList <- lapply(cyclGenesList, function(x) setkey(x, GENE))
tisSpecCycl[, Amp_Q_BRAIN := cyclGenesList$BRAIN[Gene]$ampQuant]
tisSpecCycl[, Amp_Q_FB := cyclGenesList$FB[Gene]$ampQuant]
tisSpecCycl[, Amp_Q_MT := cyclGenesList$MT[Gene]$ampQuant]
tisSpecCycl[, Amp_Q_GUT := cyclGenesList$GUT[Gene]$ampQuant]
# add information about the quantile of expression
white <- lapply(white, function(x) as.data.table(cbind(GENE = rownames(x), x)))
white <- lapply(white, function(x) x[, meanExpr := apply(x[, -1, with = F],
                                                         1, function(x) mean(as.numeric(x)))])
white <- lapply(white,
                function(x) x[, exprQuant := as.numeric(quantcut(-meanExpr,
                                                                q = 5))])
white <- lapply(white, function(x) setkey(x, GENE))
tisSpecCycl[, Expr_Q_BRAIN := white$BRAIN[Gene]$exprQuant]
tisSpecCycl[, Expr_Q_FB := white$FB[Gene]$exprQuant]
tisSpecCycl[, Expr_Q_MT := white$MT[Gene]$exprQuant]
tisSpecCycl[, Expr_Q_GUT := white$GUT[Gene]$exprQuant]

# let's see, how many of tissue-specific cycling genes also expressed in 1/2/3/
# 4 another tissues
for (tissue in allTissues) {
  tisSpecCycl_oneT <- tisSpecCycl[tisSpecCycl$CyclesIn == tissue, ]
  tisSpecCycl_oneT <- tisSpecCycl_oneT[, grepl('Expr', 
                                               colnames(tisSpecCycl_oneT)),
                                       with = F]
  exprCategory <- apply(tisSpecCycl_oneT, 1, 
                        function(x) c('BRAIN', 'FB', 'MT', 'GUT')[!is.na(x)])
  exprCategory <- sapply(exprCategory, paste, collapse = '&')
  
  #png(paste0(savePlotsDir, 'tissueSpecCyclExpr_', tissue, '.png'), width = 800,
  #    height = 800)
  par("mar" = c(5, 21, 2, 2))
  barplot(sort(table(exprCategory)), horiz = T, ylab = '', las = 2, 
          cex.axis = 2, cex.names = 2,
          main = paste0('Number of ', tissue, 
                        '-specific cycling genes expressed in other tissues'),
          xlab = '', col = tissueColor[tissue])
  #dev.off()
}

# chi square test on dependence of tissue specific cycling from tissue 
# specific expression
cyclOneTis <- tisSpecCycl[grepl(tissue, CyclesIn), ]
cyclOneTis <- cyclOneTis[, c('CyclesIn', 'Expr_Q_BRAIN', 'Expr_Q_FB',
                                       'Expr_Q_MT', 'Expr_Q_GUT')]
cyclOneTis[, Expr_Q_BRAIN := ifelse(is.na(Expr_Q_BRAIN), '', 'BRAIN')]
cyclOneTis[, Expr_Q_FB := ifelse(is.na(Expr_Q_FB), '', 'FB')]
cyclOneTis[, Expr_Q_GUT := ifelse(is.na(Expr_Q_GUT), '', 'GUT')]
cyclOneTis[, Expr_Q_MT := ifelse(is.na(Expr_Q_MT), '', 'MT')]
cyclOneTis[, abu := apply(cyclOneTis[, 2:5], 1, paste, collapse = '-')]
cyclOneTis[, abu := gsub('^-{1, }', '', abu)]
cyclOneTis[, abu := gsub('-{1, }$', '', abu)]
cyclOneTis[, abu := gsub('-{2,}', '-', abu)]
chisq.test(table(cyclOneTis$CyclesIn, cyclOneTis$abu))

# OVERLAP CYCLING GENES WITH ChIP-seq PEAKS -----------------------------------
# transform to data table for convinience
headSpec$SYMBOL <- unlist(headSpec$SYMBOL)
headSpec <- as.data.table(headSpec)
setkey(headSpec, SYMBOL)

bodySpec$SYMBOL <- unlist(bodySpec$SYMBOL)
bodySpec <- as.data.table(bodySpec)
setkey(bodySpec, SYMBOL)

greyZone$SYMBOL <- unlist(greyZone$SYMBOL)
greyZone <- as.data.table(greyZone)
setkey(greyZone, SYMBOL)

sharedPeaks$SYMBOL <- unlist(sharedPeaks$SYMBOL)
sharedPeaks <- as.data.table(sharedPeaks)
setkey(sharedPeaks, SYMBOL)

# Add info about peaks
tisSpecCycl[, headPeak := ifelse(Gene %in% headSpec$SYMBOL, T, NA)]
tisSpecCycl[, headPeakLoc := sapply(Gene, 
                                    function(x) ifelse(x %in% headSpec$SYMBOL,
                                                       paste(headSpec[x]$annotation,
                                                             collapse = ', '),
                                                       ''))]
tisSpecCycl[, bodyPeak := ifelse(Gene %in% bodySpec$SYMBOL, T, NA)]
tisSpecCycl[, bodyPeakLoc := sapply(Gene, 
                                    function(x) ifelse(x %in% bodySpec$SYMBOL,
                                                       paste(bodySpec[x]$annotation,
                                                             collapse = ', '),
                                                       ''))]
tisSpecCycl[, greyPeak := ifelse(Gene %in% greyZone$SYMBOL, T, NA)]
tisSpecCycl[, greyPeakLoc := sapply(Gene, 
                                    function(x) ifelse(x %in% greyZone$SYMBOL,
                                                       paste(greyZone[x]$annotation,
                                                             collapse = ', '),
                                                       ''))]
tisSpecCycl[, sharedPeak := ifelse(Gene %in% sharedPeaks$SYMBOL, T, NA)]
tisSpecCycl[, sharedPeakLoc := sapply(Gene, 
                                      function(x) ifelse(x %in% sharedPeaks$SYMBOL,
                                                         paste(sharedPeaks[x]$annotation,
                                                               collapse = ', '),
                                                         ''))]
saveRDS(tisSpecCycl, paste0(saveRDSdir, 'tisSpecCycl.Rds'))

# table for the presentation
cyclWithPeakToPrint <- c()
for (tiss in allTissues) {
  cyclWithPeak <- tisSpecCycl[CyclesIn == tiss & 
                              (headPeak != '' | bodyPeak != '' | 
                               greyPeak != '' | sharedPeak != ''), ]
  cyclWithPeak <- data.frame(Tissue = tiss,
                             genesWithPeaks = nrow(cyclWithPeak),
                             stringsAsFactors = F)
  cyclWithPeak$genesProm <- nrow(tisSpecCycl[headPeakLoc == 'Promoter' |
                                             bodyPeakLoc == 'Promoter' |
                                             greyPeakLoc == 'Promoter' |
                                             sharedPeakLoc == 'Promoter', ])
  genesPromSymb <- paste(tisSpecCycl[headPeakLoc == 'Promoter' |
                                     bodyPeakLoc == 'Promoter' |
                                     greyPeakLoc == 'Promoter' |
                                     sharedPeakLoc == 'Promoter', ]$gene,
                         collapse = ", ")
  cyclWithPeak$genesPromSymb <- genesPromSymb
  cyclWithPeakToPrint <- rbind(cyclWithPeakToPrint, cyclWithPeak)
}
cyclWithPeakToPrint

# PLOT TISSUE SPECIFIC CYCLING GENES ------------------------------------------
for (ttp in c('BRAIN', 'FB', 'GUT', 'MT')) {
  oneTisCycl <- tisSpecCycl[CyclesIn == ttp]
  oneTisCycl <- oneTisCycl[, c(1:10)]
  oneTisCycl <- oneTisCycl[, apply(oneTisCycl, 2, function(x) !all(is.na(x))),
                             with = F]
  oneTisCycl <- oneTisCycl[complete.cases(oneTisCycl)]
  colsToOrder <- which(grepl(ttp, colnames(oneTisCycl)) == T)
  
  oneTisCycl <- oneTisCycl[order(oneTisCycl[, colsToOrder[1], with = F],
                                 oneTisCycl[, colsToOrder[2], with = F])]
  print(oneTisCycl[1:10, ])
}

GTPs <- c(BRAIN = 'Dip-B', FB = 'CG17323', GUT = 'Cyp6a2', MT = 'Fas2')

for (gtp in GTPs) {
  listToPlot <- list()
  for (ttp in c('BRAIN', 'FB', 'GUT', 'MT')) {
    INFOTiss <- INFO[Tissue == ttp]
    setkey(INFOTiss, rightGT_name)
    plotToAdd <- plotCountsTimecourse(gtp, whiteNoStand[[ttp]],
                                      INFOTiss[colnames(white[[ttp]])]$Time,
                                      tissueColor[ttp], 
                                      ggtitle(paste('Timeline of', gtp, 
                                                    'expression in', ttp)))
    listToPlot[[length(listToPlot) + 1]] <- plotToAdd
  }
  print(multiplot(plotlist = listToPlot, cols = 2))
}

# GO-enrichment of tissue-specific cycling genes ------------------------------
ensembl91 <- useEnsembl("ensembl", version = 91, 
                        dataset = "dmelanogaster_gene_ensembl")
# select genes, which cycle in at least one tissue - this is our backkground
geneCycl <- tisSpecCycl[CyclesIn != '']
bg <- getFlyBaseID(geneCycl$Gene, ensembl91)[[1]]
bg <- as.data.table(bg)

goCycl1Tiss <- data.table()
for (tissue in allTissues) {
  # genes cycling in that tissue 
  cycl1T <- geneCycl[CyclesIn == tissue]
  # select genes which are also expressed in all tissues
  exprAllTiss <- complete.cases(cycl1T[, grepl('Expr', colnames(geneCycl)), 
                                       with = F])
  cycl1T <- cycl1T[exprAllTiss, ]
  
  # get flaybase IDs
  cycl1T <- getFlyBaseID(cycl1T$Gene, ensembl91)[[1]]
  cycl1T <- as.data.table(cycl1T)
  
  # perform GO enrichment
  goBP <- GOenricment(cycl1T$flybase_gene_id, bg$flybase_gene_id, ont = 'BP')
  goBP <- as.data.table(goBP)
  goBP[, Tissue := tissue]
  goBP[, ont := 'BP']
  goMF <- GOenricment(cycl1T$flybase_gene_id, bg$flybase_gene_id, ont = 'MF')
  goMF <- as.data.table(goMF)
  goMF[, Tissue := tissue]
  goMF[, ont := 'MF']
  goCC <- GOenricment(cycl1T$flybase_gene_id, bg$flybase_gene_id, ont = 'CC')
  goCC <- as.data.table(goCC)
  goCC[, Tissue := tissue]
  goCC[, ont := 'CC']
  
  # add to final table
  goCycl1Tiss <- rbind(goCycl1Tiss, goBP)
  goCycl1Tiss <- rbind(goCycl1Tiss, goMF)
  goCycl1Tiss <- rbind(goCycl1Tiss, goCC)
}

# GO-enrichment of genes cycling in several tissues ---------------------------
ensembl91 <- useEnsembl("ensembl", version = 91, 
                        dataset = "dmelanogaster_gene_ensembl")
# select genes, which cycle in at least one tissue - this is our backkground
geneCycl <- tisSpecCycl[CyclesIn != '']
bg <- getFlyBaseID(geneCycl$Gene, ensembl91)[[1]]
bg <- as.data.table(bg)

# patterns of gene cycling in different tissues
cyclTissPttrn <- tisSpecCycl[grep('-', CyclesIn)]
# select only ones which have > 10 genes
cyclTissPttrn <- cyclTissPttrn[,.N, by = CyclesIn][N >= 10]$CyclesIn

goCyclSevTiss <- data.table()
for (tissueCombo in cyclTissPttrn) {
  # genes cycling in that tissue 
  cyclSevTiss <- geneCycl[CyclesIn == tissueCombo]
  # select genes which are also expressed in all tissues
  exprAllTiss <- complete.cases(cyclSevTiss[, grepl('Expr', colnames(geneCycl)), 
                                            with = F])
  cyclSevTiss <- cyclSevTiss[exprAllTiss, ]
  
  # get flaybase IDs
  cyclSevTiss <- getFlyBaseID(cyclSevTiss$Gene, ensembl91)[[1]]
  cyclSevTiss <- as.data.table(cyclSevTiss)
  
  # perform GO enrichment
  goBP <- GOenricment(cyclSevTiss$flybase_gene_id, bg$flybase_gene_id, ont = 'BP')
  goBP <- as.data.table(goBP)
  goBP[, TissueCombo := tissueCombo]
  goBP[, ont := 'BP']
  goMF <- GOenricment(cyclSevTiss$flybase_gene_id, bg$flybase_gene_id, ont = 'MF')
  goMF <- as.data.table(goMF)
  goMF[, TissueCombo := tissueCombo]
  goMF[, ont := 'MF']
  goCC <- GOenricment(cyclSevTiss$flybase_gene_id, bg$flybase_gene_id, ont = 'CC')
  goCC <- as.data.table(goCC)
  goCC[, TissueCombo := tissueCombo]
  goCC[, ont := 'CC']
  
  # add to final table
  goCyclSevTiss <- rbind(goCyclSevTiss, goBP)
  goCyclSevTiss <- rbind(goCyclSevTiss, goMF)
  goCyclSevTiss <- rbind(goCyclSevTiss, goCC)
}