# FILE: 13_motifEnrichPostProc.R ----------------------------------------------
#
# DESCRIPTION : Post-processing of the motif enrichment done with Homer
#
# USAGE: 
#
# OPTIONS:  none
# REQUIREMENTS:  data.table, lubridate, stringr
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  02.08.2017
# REVISION: 02.08.2017

setwd('~/Desktop/AroundTheClock_Aug_Oct2017/')
setwd('~/Desktop/BitBucket/AroundTheClock_Aug_Oct2017/')
source('0_common_functions.R')

# FUNCTIONS -------------------------------------------------------------------
#' annotateMotifs
#' Annotates motifs with corresponding full gene name and GO. If it's mouse 
#' motif, gives fly ortholog
#' @param geneNames names of the genes/motifs to annotate
#' @param TFtoGene result of initNameToFBidTab
#' @param flyEnsembl fly ensembl biomart object
#' @param mouseEnsembl mouse ensembl biomart object
#' @return data table with FlybaseID, DrosoTFShort, DrosoTF, SpecieOfOrigin,
#' Full gene name, GO
annotateMotifs <- function(geneNames, TFtoGene, flyEnsembl, mouseEnsembl) {
  info <- data.table(DrosoTF = geneNames, 
                     SpecieOfOrigin = ifelse(grepl('not in Dmel', geneNames),
                                             'MMus', 'DMel'))
  info[, DrosoTFShort := gsub('\\(not in Dmel DB\\)', '', DrosoTF)]
  info[SpecieOfOrigin == 'MMus']$DrosoTFShort <- gsub('[[:punct:]].*', '', 
                                                      info[SpecieOfOrigin == 
                                                          'MMus']$DrosoTFShort)
  setkey(info, DrosoTFShort)
  
  # first of all, if it's a mouse gene, get corresponding ortholog flybase id
  mouseOrt <- getBM(attributes = c("external_gene_name", 
                                   'dmelanogaster_homolog_ensembl_gene'),
                    filters = "external_gene_name", 
                    values = info$DrosoTFShort, mouseEnsembl)
  colnames(mouseOrt) <- c('DrosoTFShort', 'FlybaseID')
  mouseOrt <- as.data.table(mouseOrt)
  
  # get flybase id for all fly genes
  flyGenes <- data.table(DrosoTFShort = info$DrosoTFShort,
                         FlybaseID = sapply(info$DrosoTFShort, 
                                            getGeneNamesIDsfromTab,
                                            TFtoGene, 'primary_FBgn'))
  flyGenes <- flyGenes[!grepl('not found in Dmel DB', FlybaseID), ]
  
  # add fly base IDs to info
  fbID <- rbind(mouseOrt, flyGenes)
  setkey(fbID, DrosoTFShort)
  info <- merge(info, fbID, all.x = T)
  
  # now get all the extra info: full gene name and GO terms
  flyGenes <- getBM(attributes = c('flybase_gene_id', 'description',
                                   'name_1006'),
                    filters = 'flybase_gene_id', 
                    values = info$FlybaseID, flyEnsembl)
  flyGenes <- as.data.table(flyGenes)
  flyGenes[, description := gsub('\\[Source.*', '', description)]
  flyGenes <- flyGenes[!name_1006 %in% c("DNA binding", "nucleus", "cytoplasm",
                                         "cytosol", "nucleoplasm")]
  flyGenes <- flyGenes[,. (unique(description), paste(unique(name_1006),
                                                      collapse = ', ')),
                       by = flybase_gene_id]
  colnames(flyGenes) <- c('FlybaseID', 'Full gene name', 'GO')
  
  # final merge
  setkey(info, FlybaseID)
  setkey(flyGenes, FlybaseID)
  info <- merge(info, flyGenes, all.x = T)
  info[, DrosoTFShort := NULL]
  info
}

#' compareDeNovoMotifs
#' @param deNovoRes1 results of de-novo to be subseted from
#' @param deNovoRes2 results of de-novo to be compared to
#' @param finalFileName file where to put results of tomtom
#' @param deNovoRes1 without motifs which are in deNovoRes2
compareDeNovoMotifs <- function(deNovoRes1, deNovoRes2, finalFileName) {
  if(nrow(deNovoRes2) != 0) {
    # write files for TomTom:
    deNovoRes1Name <- paste0('MotifEnrichment/deNovo1.meme')
    deNovoRes2Name <- paste0('MotifEnrichment/deNovo2.meme')
    printInMemeFormat(deNovoRes1$File, deNovoRes1Name)
    printInMemeFormat(deNovoRes2$File, deNovoRes2Name)
    
    # run TomTom
    system(paste0('tomtom ', deNovoRes1Name, ' ', deNovoRes2Name, ' -oc ',
                  'MotifEnrichment/deNovo1_vs_deNovo2/'))
    system(paste('mv MotifEnrichment/deNovo1_vs_deNovo2/tomtom.txt ', 
                 finalFileName))
    system(paste0('rm -r MotifEnrichment/deNovo1_vs_deNovo2/'))
    system(paste0('rm MotifEnrichment/deNovo1.meme'))
    system(paste0('rm MotifEnrichment/deNovo2.meme'))
    
    # read-in TomTom
    deNovo_res <- readInTomTom(finalFileName, 0.05, tfToGene)
    deNovo_res <- deNovoRes1[!Consensus %in% deNovo_res$Consensus]
    deNovo_res <- deNovo_res[,.(Consensus, group, `P-value`, PercTarg, PercBg, 
                                FC)]
  } else {
    deNovo_res <- deNovoRes1[,.(Consensus, group, `P-value`, PercTarg, PercBg,
                                FC)]
  }
  deNovo_res
}

#' plotHeatmapPerBg
#' Plots all heatmaps of motif enrichment per background
#' @param knownMotEnr enrichment of known motifs
#' @param deNovoMotEnr enrichment of de novo motifs
#' @return heatmaps in PDF
plotHeatmapPerBg <- function(knownMotEnr, deNovoMotEnr) {
  for (bgType in unique(knownMotEnr$BGname)) {
    known_TSC_bg <- knownMotEnr[BGname == bgType]
    setkey(known_TSC_bg, TargName, TargAmpQ, BGname, BGExprQ, upstrLen, tss50)
    deNovo_TSC_bg <- deNovoMotEnr[BGname == bgType]
    setkey(deNovo_TSC_bg, TargName, TargAmpQ, BGname, BGExprQ, upstrLen, tss50)
    
    numbOfMotifs <- numberOfMotifs(known_TSC_bg, deNovo_TSC_bg)
    numbOfMotifs <- numbOfMotifs[,.(sum(N)), by = .(TargName, TargAmpQ, BGname,
                                                    BGExprQ, upstrLen, tss50)]
    numbOfMotifs <- numbOfMotifs[order(-V1)]
    
    pdfName <- paste0('plots/', unique(known_TSC_bg$TargName), '_',
                      bgType, '.pdf')
    pdf(pdfName, width = 5, height = 10)
    for (i in 1:nrow(numbOfMotifs)) {
      forHeatmap <- numbOfMotifs[i, -ncol(numbOfMotifs), with = F]
      knownEnric <- known_TSC_bg[forHeatmap]
      deNovoEnric <- deNovo_TSC_bg[forHeatmap]
      # plot title
      plotTitle <- paste0('Motif enrichment of ', knownEnric$TargFullName[1],
                          '(Q-amp=', knownEnric$TargAmpQ[1], ')\n over ',
                          knownEnric$BgFullName[1], 
                          '(Q-expr=', knownEnric$BGExprQ[1], ')\n on ', 
                          knownEnric$upstrLen[1], 'bp around TSS with ',
                          knownEnric$tss50[1], 'bp excluded')
      plotTitle <- tolower(plotTitle)
      
      plotMotifHeatmap(knownEnric, deNovoEnric, main = plotTitle)
    }
    dev.off()
  }
}

#' readInTomTom
#' Reads in and parces tomtom output
#' @param filePath path to tomtom output file
#' @param eValCut cut off on e value
#' @param TFtoGeneDB result of initNameToFBidTab
#' @return data table with MotifName, Consensus, DrosoTF, E-value, 
#' DrosoTF_consensus DrosoTF_orientation
readInTomTom <- function(filePath, eValCut, TFtoGeneDB) {
  tomtom <- fread(filePath, header = T)
  tomtom <- tomtom[`E-value` < eValCut]
  colnames(tomtom)[1] <- 'Consensus'
  tomtom <- tomtom[order(Consensus, `E-value`)]
  tomtom <- tomtom[,.SD[1], Consensus]
  tomtom[, `Target ID` := gsub('_.*', '', `Target ID`)]
  
  tomtom[, DrosoTF := sapply(tomtom$`Target ID`, getGeneNamesIDsfromTab,
                             TFtoGeneDB)]
  # do it twice because sometimes FBpp show up
  if (any(sapply(tomtom$DrosoTF, function(x) grepl('FBpp', x)))) {
    FBpp <- sapply(tomtom$DrosoTF, function(x) grepl('FBpp', x))
    tomtom$DrosoTF[FBpp] <- sapply(tomtom$DrosoTF[FBpp], 
                                   getGeneNamesIDsfromTab,
                                   TFtoGeneDB)
  }
  tomtom <- tomtom[,.(Consensus, DrosoTF, `E-value`, `Target consensus`, 
                      Orientation)]
  tomtom$DrosoTF <- gsub('E-box(not in Dmel DB)', 'E-box', 
                         tomtom$DrosoTF)
  tomtom$DrosoTF <- gsub('TATA-box(not in Dmel DB)', 'TATA-box',
                         tomtom$DrosoTF)
  colnames(tomtom)[4:5] <- c('DrosoTF_consensus', 'DrosoTF_orientation')
  tomtom
}

# INPUTS ----------------------------------------------------------------------
# info about samples
infoTab <- readRDS('Rds/sampleInfo.Rds')
infoTab <- as.data.table(infoTab)

# list all tissues and corresponding colors
allTissues <- unique(infoTab$Tissue)
tissueColor <- c('#688B8A', '#A0B084', '#FAEFD4', '#A57C65')
names(tissueColor) <- allTissues

# all tested combinations of targets and backgrounds
TSE_targBgCombs <- readRDS('Rds/targBgCombs_TSE.Rds')
TSC_targBgCombs <- readRDS('Rds/targBgCombs_TSC.Rds')
targBgCombs <- rbind(TSE_targBgCombs, TSC_targBgCombs)

# results for the length 8
# ME_len8 <- 'MotifEnrichment/len_8/'
# results for the length 10
ME_len10 <- 'MotifEnrichment/len_10/'

# motifs of interest and cutoff on q value
MOI <- c('E-box', 'GATA', 'Unknown2', 'BMAL1', 'NPAS2', 'bHLH', 'CLOCK', 'Myc',
         'IBL1', 'PIF', 'Max', 'SPCH', 'USF1', 'HIF', 'Cbf')
qValueCutoff <- 0.1

# database of TF to gene
tfToGene <- initNameToFBidTab()

# output dirs
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'

# HOMER MOTIF DISCOVERY read-in -----------------------------------------------
#known_10 <- readHomerKnownAllIn_v2(ME_len10, qValueCutoff)
#known_10[, motifLen := 10]
#setkey(known_10, TargName, TargAmpQ, BGname, BGExprQ, upstrLen, tss50)
# save RDS, because it takes long time to read-in
#saveRDS(known_10, 'Rds/known_10.Rds')
known_10 <- readRDS('Rds/known_10.Rds')

#deNovo_10 <- readHomerDeNovoAllIn_v2(ME_len10) 
#deNovo_10[, motifLen := 10]
#setkey(deNovo_10, TargName, TargAmpQ, BGname, BGExprQ, upstrLen, tss50)
# save RDS, because it takes long time to read-in
#saveRDS(deNovo_10, 'Rds/deNovo_10.Rds')
deNovo_10 <- readRDS('Rds/deNovo_10.Rds')

# ANNOTATE TO DROSO MOTIFS ----------------------------------------------------
# known - PWM are the same for the same MotifName
#printInMemeFormat(known_10[!duplicated(MotifName)]$File, 
#                  "MotifEnrichment/known_len10.meme")
#system(paste('tomtom MotifEnrichment/known_len10.meme',
#             'MotifEnrichment/tomtom_db_fly/*.meme -oc', 
#             'MotifEnrichment/known_len10_tomtom'))
#system(paste('mv MotifEnrichment/known_len10_tomtom/tomtom.txt', 
#             'MotifEnrichment/known_len10_tomtom.txt'))
#system('rm -r MotifEnrichment/known_len10_tomtom')
known_10_tomtom <- readInTomTom('MotifEnrichment/known_len10_tomtom.txt', 0.05,
                                tfToGene) 
setkey(known_10_tomtom, Consensus)
setkey(known_10, Consensus)
known_10 <- merge(known_10, known_10_tomtom, all = T)
known_10[, DrosoTF := sapply(1:nrow(known_10),
                             function(x) if(is.na(known_10$DrosoTF[x])) {
                               paste0(known_10$MotifName[x], '(not in Dmel DB)')
                               } else {known_10$DrosoTF[x]})]
saveRDS(known_10, 'Rds/known_10_withDrosoTF.Rds')
  
# de novo : have to run on everything for the unbiased matching
#printInMemeFormat(deNovo_10$File, "MotifEnrichment/deNovo_len10.meme")
#system(paste('tomtom MotifEnrichment/deNovo_len10.meme',
#             'MotifEnrichment/tomtom_db_fly/*.meme -oc', 
#             'MotifEnrichment/deNovo_len10_tomtom'))
#system(paste('mv MotifEnrichment/deNovo_len10_tomtom/tomtom.txt', 
#             'MotifEnrichment/deNovo_len10_tomtom.txt'))
#system('rm -r MotifEnrichment/deNovo_len10_tomtom')
deNovo_10_tomtom <- readInTomTom('MotifEnrichment/deNovo_len10_tomtom.txt', 0.05,
                                 tfToGene) 
setkey(deNovo_10_tomtom, Consensus)
setkey(deNovo_10, Consensus)
deNovo_10 <- merge(deNovo_10, deNovo_10_tomtom, all = T)
deNovo_10[, DrosoTF := sapply(1:nrow(deNovo_10),
                             function(x) if(is.na(deNovo_10$DrosoTF[x])) {
                               deNovo_10$Consensus[x]}
                             else {deNovo_10$DrosoTF[x]})]
saveRDS(deNovo_10, 'Rds/deNovo_10_withDrosoTF.Rds')

# HIST OF NUMBER OF DETECTED MOTIFS (LEN 10) ----------------------------------
numbMotifs_10 <- numberOfMotifs(known_10, deNovo_10)
pdf('plots/MotifEnrich_TSC_Numb.pdf', height = 12, width = 25)
for (targQ in unique(numbMotifs_10$TargAmpQ)) {
  print(plotNumberMotifs(numbMotifs_10[TargName == 'TSC' & TargAmpQ == targQ],
                         paste('Overview of the motif enrichment for', 
                               'tissue-specific cycling genes, amp. quantile =',
                               targQ)))
}
dev.off()
pdf('plots/MotifEnrich_TSC_noDeNovoFP_Numb.pdf', height = 12, width = 25)
for (targQ in unique(numbMotifs_10$TargAmpQ)) {
  print(plotNumberMotifs(numbMotifs_10[TargName == 'TSC' & TargAmpQ == targQ &
                                       Type != 'de-novo_FP'],
                         paste('Overview of the motif enrichment for', 
                               'tissue-specific cycling genes, amp. quantile =',
                               targQ)))
}
dev.off()

pdf('plots/MotifEnrich_TSE_Numb.pdf', height = 12, width = 16)
for (targQ in unique(numbMotifs_10$TargAmpQ)) {
  print(plotNumberMotifs(numbMotifs_10[TargName == 'TSE' & TargAmpQ == targQ],
                         paste('Overview of the motif enrichment for', 
                               'tissue-specific expressed genes, expr. quantile =',
                               targQ)))
}
dev.off()
pdf('plots/MotifEnrich_TSE_noDevoFP_Numb.pdf', height = 12, width = 16)
for (targQ in unique(numbMotifs_10$TargAmpQ)) {
  print(plotNumberMotifs(numbMotifs_10[TargName == 'TSE' & TargAmpQ == targQ &
                                       Type != 'de-novo_FP'],
                         paste('Overview of the motif enrichment for', 
                               'tissue-specific expressed genes, expr. quantile =',
                               targQ)))
}
dev.off()

# PRINT HEATMAPS (LEN 10) ------------------------------------------------------
plotHeatmapPerBg(known_10[TargName == 'TSC'], deNovo_10[TargName == 'TSC'])
plotHeatmapPerBg(known_10[TargName == 'TSE'], deNovo_10[TargName == 'TSE'])

# PRODUCE FIGURES -------------------------------------------------------------
known_10 <- readRDS('Rds/known_10_withDrosoTF.Rds')
deNovo_10 <- readRDS('Rds/deNovo_10_withDrosoTF.Rds')

# for unknown 2
pdf('plots/TSC_over_RANDOM.pdf', width = 5, height = 10)
known_unknown2 <- known_10[TargName == 'TSC' & TargAmpQ == 5 & 
                           BGname == 'RANDOM' & BGExprQ == 5 & 
                           upstrLen == 500 & tss50 == 0]
deNovo_unknown2 <- deNovo_10[TargName == 'TSC' & TargAmpQ == 5 & 
                             BGname == 'RANDOM' & BGExprQ == 5 & 
                             upstrLen == 500 & tss50 == 0 & FP == F]
unknown2HeatMatrix <- createHeatMatrix(known_unknown2, deNovo_unknown2, 'FC',
                                       'P-value')
par(cex.main = 0.75)
colorScheme <- c('white', "#8E9B97", "#5B7065", "#537072", "#304040",
                 "#2C4A52")
heatColors <- colorRampPalette(colorScheme)(n = 100)
colBreaks <- seq(0, max(unknown2HeatMatrix$heatMatrix), 
                 length.out = length(heatColors) + 1)
heatmap.2(unknown2HeatMatrix$heatMatrix, density.info = "none", trace = "none",
          margins = c(6, 10), dendrogram = "none",  
          colRow = unknown2HeatMatrix$rowColors, #Colv = "NA", Rowv = "NA",
          key = T, notecol = 'black', col = heatColors, 
          breaks = colBreaks, lmat = rbind(c(3, 0), c(1, 2), c(4, 0)),
          lwid = c(5, 0.1), lhei = c(0.5, 5, 0.85),
          main = paste0('Motif enrichment of tissue-specific cycling genes',
                        '\n over random background on 500bp around TSS'))
dev.off()

# for TSC vs ECOTG
known_tsc <- known_10[TargName == 'TSC' & TargAmpQ == 5 & BGname == 'NCOTG' & 
                      BGExprQ == 5 & upstrLen == 500 & tss50 == 0]
deNovo_tsc <- deNovo_10[TargName == 'TSC' & TargAmpQ == 5 & BGname == 'NCOTG' & 
                        BGExprQ == 5 & upstrLen == 500 & tss50 == 0 & FP == F]
# for TSE vs TSEOT
known_tse <- known_10[TargName == 'TSE' & TargAmpQ == 5 & BGname == 'TSEOT' &
                      BGExprQ == 5 & upstrLen == 500 & tss50 == 0]
deNovo_tse <- deNovo_10[TargName == 'TSE' & TargAmpQ == 5 & BGname == 'TSEOT' &
                        BGExprQ == 5 & upstrLen == 500 & tss50 == 0 & FP == F]

flyEnsembl <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
mouseEnsembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

known_tsc_tse <- c()
deNovo_tsc_tse <- c()
for (tissue in allTissues) {
  known_tscTis <- known_tsc[group == tissue]
  known_tseTis <- known_tse[group == tissue]
  
  # get known motifs and annotate them
  known_res <- setdiff(known_tscTis$DrosoTF, known_tseTis$DrosoTF)
  known_res <- known_tscTis[DrosoTF %in% known_res]
  message(paste0('Length of known results in ', tissue, ': ',
                 length(unique(known_res$DrosoTF))))
  if (nrow(known_res) != 0) {
    knownInfo <- annotateMotifs(unique(known_res$DrosoTF), tfToGene, 
                                flyEnsembl, mouseEnsembl)
    known_tscTisReduced <- known_res[, .(DrosoTF, group, `q-value`, PercTarg,
                                         PercBg, FC)]
    setkey(knownInfo, DrosoTF)
    setkey(known_tscTisReduced, DrosoTF)
    knownInfo <- merge(known_tscTisReduced, knownInfo, all = T)
  }

  # get de-novo motifs and compare them with tomtom
  deNovo_tscTis <- deNovo_tsc[group == tissue]
  message(paste0('Length of de novo results in ', tissue, ': ',
                 nrow(deNovo_tscTis)))
  deNovo_tseTis <- deNovo_tse[group == tissue]
  deNovoInfo <- compareDeNovoMotifs(deNovo_tscTis, deNovo_tseTis,
                                    paste0('MotifEnrichment/', tissue,
                                           '_TSC_vs_TSE_deNovo_tomtom.txt'))
  
  # put in result
  known_tsc_tse <- rbind(known_tsc_tse, knownInfo)
  deNovo_tsc_tse <- rbind(deNovo_tsc_tse, deNovoInfo)
}
write.table(known_tsc_tse, 'known_TSC_over_NCOTG_vs_TSE_over_TSEOT.csv', 
            quote = F, sep = '\t', row.names = F, col.names = T)
write.table(deNovo_tsc_tse, 'deNovo_TSC_over_NCOTG_vs_TSE_over_TSEOT.csv', 
            quote = F, sep = '\t', row.names = F, col.names = T)

# plot heatmap
tsc_tseHeatMatrix <- createHeatMatrix(known_tsc, deNovo_tsc, 'FC', 'P-value')
rowColors <- sapply(rownames(tsc_tseHeatMatrix$heatMatrix),
                    function(x) {
                      x <- gsub('de-novo:', '', x)
                      result <- '#F9D423'
                      if (x %in% known_tsc_tse$DrosoTF) {result <- '#107FC9'}
                      if (x %in% deNovo_tsc_tse$Consensus) {result <- '#A0C55F'}
                      result
                    })
pdf('plots/TSC_over_NCOTG_vs_TSE_over_TSEOT.pdf', width = 5, height = 10)
par(cex.main = 0.75)
colorScheme <- c('white', "#8E9B97", "#5B7065", "#537072", "#304040",
                 "#2C4A52")
heatColors <- colorRampPalette(colorScheme)(n = 100)
colBreaks <- seq(0, max(tsc_tseHeatMatrix$heatMatrix), 
                 length.out = length(heatColors) + 1)
heatmap.2(tsc_tseHeatMatrix$heatMatrix, density.info = "none", trace = "none",
          margins = c(6, 10), dendrogram = "none",  
          colRow = rowColors, #Colv = "NA", Rowv = "NA",
          key = T, notecol = 'black', col = heatColors, 
          breaks = colBreaks, lmat = rbind(c(3, 0), c(1, 2), c(4, 0)),
          lwid = c(5, 0.1), lhei = c(0.5, 5, 0.85),
          main = paste0('Motif enrichment of tissue-specific cycling genes',
                        '\n overnot cycling in any tissue genes \n on 500bp around TSS'))
dev.off()