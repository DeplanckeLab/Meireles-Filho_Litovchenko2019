# FILE: 0_common_functions.R --------------------------------------------------
#
# DESCRIPTION : 
# Contains functions used in other scripts
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
# CREATED:  15.11.2017
# REVISION: 15.11.2017

source('JTK_CYCLEv3.1.R')
library(zeitzeiger)
library(biomaRt)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(ChIPseeker)
library(circular)
library(circacompare)
library(colorspace)
library(data.table)
library(DESeq2)
library(Directional)
#library(doParallel)
library(GENIE3)
library(GenomicRanges)
library(igraph)
library(ggplot2)
#library(ggplus)
library(glmnet)
library(gtools)
library(gplots)
library(lattice)
library(limma)
library(MotIV)
#library(NPCirc)
library(rGADEM)
library(refGenome)
library(stringr)
#library(sva)
library(topGO)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(VariantAnnotation)
library(venn)
library(XML)

# LOAD GENOME ANNOTATION ------------------------------------------------------
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene

currentDir <- getwd()
setwd("~/Documents/RefGen/dm3/")
# create USCS object to store the genome
uscsRefGen <- ensemblGenome()
# read GTF file into uscsRefGen object
read.gtf(uscsRefGen, "dm3refGene.srt.gtf")
uscsRefGenCorrds <- as.data.table(getGenePositions(uscsRefGen))
setkey(uscsRefGenCorrds, gene_name)
setwd(currentDir)

# READS INPUT, NORMALIZATIONS, BATCH CORRECTION -------------------------------
# readCounts
# cutoffCoverage
# removeLowExpressed
# normVoom
# removeBatch

#' cutoffCoverage
#' Removes samples with the coverage less than cutoff 
#' @param inputTab input count table, columns are samples and rows are genes
#' gene names are in the row names
#' @param cutOff cutoff on coverage
#' @return count table only containing samples which passed coverage cutoff
cutoffCoverage <- function(inputTab, cutOff) {
  samplesToKeep <- which(colSums(inputTab) > cutOff)
  message(paste('Removed', ncol(inputTab) -  length(samplesToKeep), 'samples'))
  message(cat(colnames(inputTab)[colSums(inputTab) <  cutOff], sep = '\n'))
  outputTab <- inputTab[, samplesToKeep]
  outputTab
}

#' normVoom
#' normalizes counts with Voom
#' @param inputTab input count table, columns are samples and rows are genes
#' gene names are in the row names
#' @return normalized count table
#' @keywords normalize, counts
normVoom <- function(inputTab) {
  # categories for the future fake model matrix
  if (ncol(inputTab) %% 2 == 0) {
     categories <- c(ncol(inputTab) / 2, ncol(inputTab) / 2)
  } else {
     categories <- c(ncol(inputTab) %/% 2, ncol(inputTab) %/% 2 + 1)
  }
  # create fake design matrix
  fakeDesign <- model.matrix(~ c(rep('1', categories[1]), 
                                 rep('2', categories[2])))
  # normalize by voom
  outputTab <- voom(inputTab, design = fakeDesign, 
                    normalize.method = "quantile")
  outputTab <- outputTab$E
  outputTab
}

#' readAndAnnotPeaks
#' Reads-in bed file with peaks and annotates it with Chipseeker 
#' @param pathToBed path to bed file
#' @details IMPORTANT NOTE: when annotatePeak annotates peaks, there's a bug!
#' It's a bug only for introns/exons. In annotation field it gives REAL 
#' Flybase ID of the gene where peak falls, but other fiels (i.e. geneID, 
#' distance to TSS, etc) correspond to the closest genes by TSS! Since Droso
#' genome is very packed, it might not be the same gene in which intron the
#' peak falls. This function takes cake of that and corrects it.
#' @return Granges with annotated peaks
readAndAnnotPeaks <- function(pathToBed) {
  peaks <- fread(pathToBed)
  colnames(peaks) <- c('chr', 'start', 'end')
  peaks <- makeGRangesFromDataFrame(peaks)
  annotPeaks <- annotatePeak(peaks, tssRegion = c(-2000, 500), TxDb = txdb,
                             level = 'gene', annoDb = "org.Dm.eg.db")@anno
  
  # take cake of what is written in the details
  allGenes <- genes(txdb)
  allGenesTss <- ifelse(strand(allGenes) == '+', start(allGenes), 
                        end(allGenes))
  names(allGenesTss) <- names(allGenes)
  ENSEMBL2EG <- as.list(org.Dm.egENSEMBL2EG)
  EG2SYMBOL <- as.list(org.Dm.egSYMBOL)
  EG2GENENAME <- as.list(org.Dm.egGENENAME)
  
  # for INTRONS
  annotPeaksIntrons <- annotPeaks[grepl('Intron', annotPeaks$annotation), ]
  annotPeaksIntrons$rightGene <- gsub(",.*", "", 
                                      gsub(".*FBgn", "FBgn",
                                           annotPeaksIntrons$annotation))
  message(paste(sum(annotPeaksIntrons$rightGene != annotPeaksIntrons$geneId),
          'of peaks were annotated to introns and have the bug.',
          'Correcting it.'))
  # correction of coordinates
  annotPeaksIntrons$geneStart <- start(allGenes[annotPeaksIntrons$rightGene])
  annotPeaksIntrons$geneEnd <- end(allGenes[annotPeaksIntrons$rightGene])
  annotPeaksIntrons$geneLength <- abs(annotPeaksIntrons$geneEnd - 
                                      annotPeaksIntrons$geneStart) + 1
  # of strand
  annotPeaksIntrons$geneStrand <- strand(allGenes[annotPeaksIntrons$rightGene])
  annotPeaksIntrons$geneStrand <- ifelse(annotPeaksIntrons$geneStrand == '-', 2, 1)
  annotPeaksIntrons$geneId <- annotPeaksIntrons$rightGene
  annotPeaksIntrons$distanceToTSS <- ifelse(strand(annotPeaksIntrons) == 1,
                                            start(annotPeaksIntrons) -
                                            allGenesTss[annotPeaksIntrons$rightGene],
                                            allGenesTss[annotPeaksIntrons$rightGene] -
                                            start(annotPeaksIntrons))
  # of gene info
  annotPeaksIntrons$ENTREZID <- sapply(annotPeaksIntrons$geneId,
                                       function(x) ENSEMBL2EG[[x]])
  annotPeaksIntrons$SYMBOL <- sapply(annotPeaksIntrons$ENTREZID,
                                     function(x) ifelse(is.null(x), NA,
                                              EG2SYMBOL[as.character(x)]))
  annotPeaksIntrons$GENENAME <- sapply(annotPeaksIntrons$ENTREZID,
                                       function(x) ifelse(is.null(x), NA,
                                                          EG2GENENAME[as.character(x)]))
  annotPeaksIntrons$annotation <- gsub(' \\(.*', '', annotPeaksIntrons$annotation)
  annotPeaksIntrons$rightGene <- NULL
  
  annotPeaksExons <- annotPeaks[grepl('Exon', annotPeaks$annotation), ]
  annotPeaksExons$rightGene <- gsub(",.*", "", 
                                    gsub(".*FBgn", "FBgn",
                                         annotPeaksExons$annotation))
  message(paste(sum(annotPeaksExons$rightGene != annotPeaksExons$geneId),
                'of peaks were annotated to Exons and have the bug.',
                'Correcting it.'))
  # correction of coordinates
  annotPeaksExons$geneStart <- start(allGenes[annotPeaksExons$rightGene])
  annotPeaksExons$geneEnd <- end(allGenes[annotPeaksExons$rightGene])
  annotPeaksExons$geneLength <- abs(annotPeaksExons$geneEnd - 
                                    annotPeaksExons$geneStart) + 1
  # of strand
  annotPeaksExons$geneStrand <- strand(allGenes[annotPeaksExons$rightGene])
  annotPeaksExons$geneStrand <- ifelse(annotPeaksExons$geneStrand == '-', 2, 1)
  annotPeaksExons$geneId <- annotPeaksExons$rightGene
  annotPeaksExons$distanceToTSS <- ifelse(strand(annotPeaksExons) == 1,
                                          start(annotPeaksExons) -
                                            allGenesTss[annotPeaksExons$rightGene],
                                          allGenesTss[annotPeaksExons$rightGene] -
                                            start(annotPeaksExons))
  # of gene info
  annotPeaksExons$ENTREZID <- sapply(annotPeaksExons$geneId,
                                     function(x) ENSEMBL2EG[[x]])
  annotPeaksExons$SYMBOL <- sapply(annotPeaksExons$ENTREZID,
                                   function(x) ifelse(is.null(x), NA,
                                                      EG2SYMBOL[as.character(x)]))
  annotPeaksExons$GENENAME <- sapply(annotPeaksExons$ENTREZID,
                                     function(x) ifelse(is.null(x), NA,
                                                        EG2GENENAME[as.character(x)]))
  annotPeaksExons$annotation <- gsub(' \\(.*', '', annotPeaksExons$annotation)
  annotPeaksExons$rightGene <- NULL
  
  annotPeaks <- annotPeaks[!grepl('Exon|Intron', annotPeaks$annotation), ]
  annotPeaks <- c(annotPeaks, annotPeaksIntrons, annotPeaksExons)
  annotPeaks
}

#' readCounts
#' Reads counts (result of HTSeq) for the selected tissue and/or line
#' Names should contain tissue name and line name
#' @param pathToFolder - full or relative path to folder with counts
#' @param tissueName - name of the tissue as its encoded in the files
#' @param lineName - name of the line as its encoded in the files
#' @param removeNotMapped - whatever or not last five lines will be removed
#' @param removeWhiteFromDGRPs - whatever or not white- samples should be 
#'        removed in case we read folder with DGRPs
#' @return data frame with rownames = genes and columns = samples
#' @keywords read, counts
readCounts <- function(pathToFolder, tissueName = NA, lineName = NA, 
                       removeNotMapped = T, removeWhiteFromDGRPs = T) {
  # list all files
  countFiles <- list.files(pathToFolder, full.names = T)
  # select tissue
  if (!is.na(tissueName)) {
     countFiles <- countFiles[grepl(toupper(tissueName), toupper(countFiles))]
  }
  # select line
  if (!is.na(lineName)) {
     countFiles <- countFiles[grepl(toupper(tissueName), toupper(countFiles))]
  }

  # combine all files for different samples together
  countTable <- read.table(countFiles[1], header = F, sep = '\t',
                           stringsAsFactors = F)
  if (nrow(countTable) != 17000) {
    stop('Row number != 17000')
  }
  for (i in 2:length(countFiles)) {
    countFile <- read.table(countFiles[i], header = F, sep = '\t',
                            stringsAsFactors = F)
    if (nrow(countFile) != 17000) {
      stop('Row number != 17000')
    }
    countTable <- cbind(countTable, countFile[, 2])
  }

  # add rowname and colnames
  rownames(countTable) <- countTable[, 1]
  countTable <- countTable[, -1]
  colnames(countTable) <- tools::file_path_sans_ext(basename(countFiles))
  
  # remove last 5 lines: '__no_feature', '__ambiguous', '__too_low_aQual', 
  # '__not_aligned', '__alignment_not_unique'
  if (removeNotMapped) {
    noFeature <- countTable['__no_feature', ]
    ambiguous <- countTable['__ambiguous', ]
    numbOfCounts <- colSums(countTable)
    message(paste('no feature counts:',
                  paste(round(100 * range(noFeature/numbOfCounts)), 
                              collapse = '-'), '%', collapse = ' '))
    message(paste('ambiguous counts:',
                  paste(round(100 * range(ambiguous/numbOfCounts)), 
                              collapse = '-'), '%', collapse = ' '))
    countTable <- countTable[!grepl('^__', rownames(countTable)), ]
  }
  
  # remove white- samples from the table in case we read DGRP folder
  if (removeWhiteFromDGRPs) {
    if (sum(grepl('_n_', colnames(countTable))) > 0) {
      countTable <- countTable[, !grepl('w-', colnames(countTable))]
    }
  }
  countTable
}

#' readCountsNormBatchCorr
#' Reads-in counts, normalizes them, batch correct and standartizize expression
#' to 1. All is based on the standart functions from 0_common_functions.R
#' @param tissueName name of the tissue
#' @param countsDirPaths vector of paths to 
#' @param infoTable info table for the samples
#' @param coverageCut cut off on coverage of the sample to be included
#' @param geneCut cut off on percentage of samples in which gene needs to be
#'                expressed to be incuded into ananlysis
#' @param normWithVoom - whatever normalization with VOOM should be performed
#' @param combatBatch - whatever batch correction with Combat should be done
#' @param standartize - whatever standartization to 0 - 1 range should be done
#' @param selectSamples if vector, then should contain names of the samples to
#'                      select (they will be in the table, all the rest - no)
#' @return data frame 
readCountsNormBatchCorr <- function(tissueName, countsDirPaths, 
                                    infoTable, coverageCut = 300000,
                                    geneCut = 0.8, normWithVoom = T,
                                    combatBatch = T, standartize = T,
                                    selectSamples = NULL) {
  # read counts from every folder
  allDirCounts <- lapply(countsDirPaths, readCounts, tissueName,
                         removeNotMapped = T, removeWhiteFromDGRPs = T)
  allCounts <- do.call(cbind, allDirCounts)
  colnames(allCounts) <- gsub('.*\\.', '', colnames(allCounts))
  # restrict to the selected samples, if desired
  if (is.vector(selectSamples) & !is.null(selectSamples)) {
    allCounts <- allCounts[, colnames(allCounts) %in% selectSamples]
  }
  
  # cutoff on coverage
  allCounts <- cutoffCoverage(allCounts, coverageCut)
  numbOfreads <- colSums(allCounts)
  # remove Low Expressed genes
  allCounts <- removeLowExpressed(allCounts, geneCut)
  # normalize with voom
  if (normWithVoom) {
    allCounts <- normVoom(allCounts)
  }
  # allCounts <- allCounts[!grepl('^__', rownames(allCounts)), ]
  # remove batch with ComBat
  if (combatBatch) {
    infoTableTissue <- infoTable[Tissue == tissueName, ]
    setkey(infoTableTissue, rightGT_name)
    allLibs <- infoTableTissue[colnames(allCounts)]$Lib
    allMonths <- infoTableTissue[colnames(allCounts)]$Month
    allCounts <- removeBatch(allCounts, as.factor(paste0(allLibs, allMonths)))
  }
  # standartize expression to 1
  if (standartize) {
    allCounts <- t(apply(allCounts, 1, standartizeTo1))
  }
  # remove genes with NaN or NA
  message(paste('Removed', sum(rowSums(is.na(allCounts)) != 0), 'genes',
                'because of the NAs'))
  allCounts <- allCounts[rowSums(is.na(allCounts)) == 0, ]
  allCounts
}

#' removeBatch
#' Remove batch effect from the data with use of combat
#' @param inputTab input count table, columns are samples and rows are genes
#' gene names are in the row names
#' @param batchVector for every sample in the input tab, which batch it 
#' correponds to
#' @return data frame with corrected for the batch effect values
#' @keywords batch effect
removeBatch <- function(inputTab, batchVector) {
  normCountTab <- ComBat(dat = as.matrix(inputTab), 
                         batch = batchVector, 
                         mod = model.matrix(~1, data = batchVector), 
                         par.prior = T, prior.plots = F)
  normCountTab
}

#' removeLowExpressed
#' Removes low expressed genes
#' @param inputTab input count table, columns are samples and rows are genes
#' gene names are in the row names
#' @param removeGenes cutoff on how many genes sample should be expressed in 
#' order to not to be excluded ( 0 < x < 1)
#' @return count table without low expressed genes
#' @keywords low expressed genes
removeLowExpressed <- function(inputTab, removeGenes) {
  # lowely expressed genes
  #lowExpr <- inputTab[rowSums(inputTab) < removeGenes * ncol(inputTab), ]
  lowExpr <- inputTab[apply(inputTab, 1, 
                            function(x) sum(x > 1)) < removeGenes * 
                                                      ncol(inputTab), ]
  # remove low expressed genes
  message(paste('Removed', nrow(lowExpr), 'low expressed genes(', 
                round(100 * nrow(lowExpr) / nrow(inputTab)) , '%)'))
  outputTab <- inputTab[!rownames(inputTab) %in% rownames(lowExpr), ]
  outputTab
}

#' standartizeTo1
#' Standartizes numerical vector x to the range of 0 to 1
#' @param x numerical vector
#' @return numerical vector with range of 0 to 1
standartizeTo1 <- function(x) {
  result <- (x - min(x))/(max(x) - min(x))
  result
}

# READ VCFs -------------------------------------------------------------------
#' readInVCFgeno
#' @param pathToVCF path to VCF file
#' @param dfOfCovars data frame of covariats, to change of ids from rals to
#' bloomington
#' @return matrix with genotypes
readInVCFgeno <- function(pathToVCF, dfOfCovars) {
  vcf <- readVcf(pathToVCF, 'dm3')
  vcfGeno <- geno(vcf)$GT
  colnames(vcfGeno) <- dfOfCovars[colnames(vcfGeno), ]$ID
  vcfGeno
}

# CYCLING GENES DETECTION -----------------------------------------------------
#' jtkLittle
#' Basic JTK function
#' @param z row from normalized count table
#' @return JTK result
jtkLittle <- function(z) {
  jtkx(z)
  c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
}

#' RunJTKCycle
#' Runs JTKCycle on normalized count table
#' @param normCountTab normalized count table
#' @param replPerTimePoint number of replicates per time point
#' @param numbHoursBetween number of hours between time points
#' @param periodToFind period of the oscillations to find
#' @return JTKCycle results
#' @author Maria Litovchenko
#' @example  
RunJTKCycle <- function(normTabOfCounts, replPerTimePoint, numbHoursBetween,
                        periodToFind = 24) {
  # first agr - total time points, 
  # second - replicates per time point
  jtkdist(length(replPerTimePoint), replPerTimePoint)

  # These number are NOT measured in hours, but instead based on the number of
  # time points per cycle.  For example, if your time points are spaced 2 hours
  # apart and you’re looking for circadian genes (20-28h) Lines 14 and 15 will
  # read: “periods <- 10:14” and “jtk.init(periods,2)”.
  periods <- c(((periodToFind/numbHoursBetween) - 2) : 
               ((periodToFind/numbHoursBetween) + 2))
  jtk.init(periods, numbHoursBetween)
  
  message(paste("JTK analysis started on", date()))
  
  result <- apply(normTabOfCounts, 1, jtkLittle)
  result <- as.data.frame(t(result))
  # adjust p values with benjamini hochberg
  bhq <- p.adjust(unlist(result[, 1]), "BH")
  result <- cbind(bhq, result)
  colnames(result) <- c("BH.Q", "ADJ.P", "PER", "LAG", "AMP")

  result <- result[order(result$ADJ.P, -result$AMP),]
  message(paste('Number of cycling genes (padj < 0.05)', 
                sum(result$ADJ.P < 0.05)))
  result
}

# HOMER preparation and postprocessing ----------------------------------------
#' readHomerKnownIn
#' Reads homer results in
#' @param filePath path to the file
#' @param qvalueCut cut off on q value
#' @param toRemove what to remove from file name before parcing
#' @return data frame with results from homer, q-value < 0.05, 
#' motif name parsed
readHomerKnownIn <- function(filePath, qvalueCut = 0.05) {
  homer <- fread(filePath)
  toParse <- strsplit(filePath, '/')[[1]]
  
  # group = tissue
  homer[, group := strsplit(grep('BRAIN|FB|GUT|MT', toParse, value = T), 
                            '_')[[1]][1]]
  
  # motif len, doesn't work for v2
  homer[, motifLen := grep('^\\d', toParse, value = T)]
  
  # select needed columns
  homer <- homer[, c(10:11, 1, 2, 3, 5, 7, 9), with = F]
  colnames(homer) <- c('group', 'motifLen', 'MotifName', 'Consensus', 
                       'P-value', 'q-value', 'PercTarg', 'PercBg')
  homer$PercTarg <- as.numeric(gsub('%', '', homer$PercTarg))
  homer$PercBg <- as.numeric(gsub('%', '', homer$PercBg))
  homer$MotifName <- gsub('/.*', '', homer$MotifName)
  homer$MotifName <- gsub('\\(.*', '', homer$MotifName)
  # add info about motif file location
  homer$File <- paste0(gsub('.txt', '/known', filePath), 1:nrow(homer), 
                       '.motif')
  homer$File <- sapply(homer$File, function(x) if(file.exists(x)) {x} 
                                   else {NA})
  
  # select motifs passing q-value cut-off
  homer <- homer[`q-value` < qvalueCut & `P-value` < qvalueCut]
  if (nrow(homer) != 0) {
    homer <- homer[order(`q-value`), ]
  }
  
  # sometimes, then the same motifs come from different sourses, we see 
  # duplication of the entries
  homer <- homer[!duplicated(homer[, -9, with = F])]
  # remove all seq biases
  homer <- homer[!grepl('SeqBias', MotifName)]
  # add rank
  homer[, Rank := 1:nrow(homer)]
  # add fold change
  homer[, FC := PercTarg / PercBg]
  # rearrange
  homer[, .(group, motifLen, MotifName, Consensus, `P-value`, `q-value`, 
            PercTarg, PercBg, Rank, FC, File)]
}

#' readHomerKnownAllIn
#' Reads in all knownMotifs.txt from all folders and subfolders of homerResDir
#' @param homerResDir directory full of homer results
#' @param enrichCode code for the enrichment performed with homer
#' @param qValCut cutoff on q-value
#' @return data table with homer results
readHomerKnownAllIn <- function(homerResDir, enrichCode, qValCut) {
  # read all the results files from all the folders
  resultDirs <- list.dirs(homerResDir, recursive = F)
  knownResultsFiles <- sapply(resultDirs, 
                              function(x) list.files(x, 'knownResults.txt',
                                                     recursive = T, 
                                                     full.names = T))
  knownResultsFiles <- unlist(knownResultsFiles)
  knownResultsFiles <- knownResultsFiles[lapply(knownResultsFiles, 
                                                length) != 0]
  knownResults <- lapply(knownResultsFiles, readHomerKnownIn, qValCut)
  names(knownResults) <- knownResultsFiles
  
  # result = summary all knownMotifs.txt
  allEnrichMotifs <- c()
  # add info about upstream, amplitude quantile and BG quantile
  for (resultDir in names(knownResults)) {
    knownMotifsOneDir <- knownResults[[resultDir]]
    resultDir <- strsplit(resultDir, '//')[[1]][2]
    if (grepl('/', resultDir)) {
      resultDir <- strsplit(resultDir, '/')[[1]][1]
    }
    resultDir <- strsplit(resultDir, '_')[[1]]
    knownMotifsOneDir[, upstr := grep('^\\d', resultDir, value = T)]
    if (enrichCode == '_specCyclGenes_') {
      knownMotifsOneDir[, ampQT := resultDir[3]]
      knownMotifsOneDir[, BG := resultDir[4]]
    } else {
      knownMotifsOneDir[, ampQT := "all"]
      knownMotifsOneDir[, BG := resultDir[3]]
    }
    allEnrichMotifs <- rbind(allEnrichMotifs, knownMotifsOneDir)
  }
  # sometimes there are several motif lenght in one folder, so upstr get confused
  if (max(as.integer(allEnrichMotifs$upstr) %% 100) != 0) {
    allEnrichMotifs[, upstr := ifelse(as.integer(upstr) %% 100 != 0,
                                      as.integer(gsub('.$', '', upstr)),
                                      as.integer(upstr))]
  }
  allEnrichMotifs
}

#' deCodeFileName
#' Decoded information from the folder name, such as argFullName TargName 
#' TargAmpQ BgFullName BGname BGExprQ upstrLen tss50
#' @param filePath full path
#' @return vector with TargFullName TargName TargAmpQ BgFullName BGname BGExprQ
#'  upstrLen tss50
deCodeFileName <- function(filePath) {
  # parse the code
  code <- strsplit(filePath, '/')[[1]]
  code <- code[sapply(code, function(x) grepl('BRAIN|GUT|FB|MT', x))]
  code <- code[sapply(code, function(x) grepl('_', x))]
  code <- strsplit(code, '_')[[1]]
  
  # create a result vector 
  result <- data.frame(TargFullName = '', TargName = code[2],
                       TargAmpQ = as.integer(code[3]), BgFullName = '', 
                       BGname = code[4], BGExprQ = as.integer(code[5]),
                       upstrLen = as.integer(code[6]), 
                       tss50 = as.integer(code[7]))
  if (result$TargName == 'TSC') {
    result$TargFullName <- 'TISSUE-SPECIFIC CYCLING GENES'
  }
  if (result$TargName == 'TSE') {
    result$TargFullName <- 'TISSUE SPECIFIC EXPRESSED'
  }
  if (result$BGname == 'RANDOM') {
    result$BgFullName <- 'RANDOM'
  }
  if (result$BGname == 'AEG') {
    result$BgFullName <- 'ALL EXPRESSED GENES'
  }
  if (result$BGname == 'TSEAT') {
    result$BGname <- 'TSEOT'
    result$BgFullName <- 'TISSUE SPECIFIC EXPRESSED OTHER TISSUES'
  }
  if (result$BGname == 'ENCG') {
    result$BgFullName <- 'EXPRESSED IN THE TISSUE, BUT NOT CYCLING'
  }
  if (result$BGname == 'ENCATG') {
    result$BgFullName <- 'EXPRESSED IN THE TISSUE, BUT NOT CYCLING IN ANY TISSUE'
  }
  if (result$BGname == 'NCATG') {
    result$BgFullName <- 'NOT CYCLING IN ANY TISSUE'
  }
  if (result$BGname == 'ECOTG') {
    result$BgFullName <- 'CYCLING IN OTHER TISSUES AND EXPRESSED IN THIS'
  }
  if (result$BGname == 'COTG') {
    result$BgFullName <- 'CYCLING IN OTHER TISSUES'
  }
  result
}

#' readHomerKnownAllIn_v2
#' Reads in all knownMotifs.txt from all folders and subfolders of homerResDir
#' Version 2, not compatible with version 1
#' @param homerResDir directory full of homer results
#' @param qValCut cutoff on q-value
#' @return data table with homer results
readHomerKnownAllIn_v2 <- function(homerResDir, qValCut) {
  # read all the results files from all the folders
  resultDirs <- list.dirs(homerResDir, recursive = F)
  knownResultsFiles <- sapply(resultDirs, 
                              function(x) list.files(x, 'knownResults.txt',
                                                     recursive = T, 
                                                     full.names = T))
  knownResultsFiles <- unlist(knownResultsFiles)
  knownResultsFiles <- knownResultsFiles[lapply(knownResultsFiles, 
                                                length) != 0]
  knownResults <- lapply(knownResultsFiles, readHomerKnownIn, qValCut)
  
  # decode the file names to get info about upstream, amplitude quantile and BG
  # quantile
  knownResultsInfo <- lapply(knownResultsFiles, deCodeFileName)
  
  # merge
  knownResults <- lapply(1:length(knownResults), 
                         function(x) if (nrow(knownResults[[x]]) != 0) {
                                     cbind(knownResults[[x]], 
                                           knownResultsInfo[[x]])
                         } else {''})
  knownResults <- knownResults[sapply(knownResults, 
                                      function(x) is.data.frame(x))]
  knownResults <- do.call(rbind, knownResults)
}

#' parseHomerDeNovoHTML
#' Parses HTML output of de-novo motifs from homer
#' @param filePath path to html 
#' @return data frame of if the de novo motif is false positive
parseHomerDeNovoHTML <- function(filePath) {
  html <- xmlParse(filePath)
  html <- xmlChildren(html)$HTML[[2]][['TABLE']]
  isFP <- c()
  i = 2
  while(!is.null(html[[i]])) {
    if (is.null(html[[i]][1]$TD[['FONT']])) {
      isFP <- c(isFP, F)
    } else {
      if (xmlAttrs(html[[i]][1]$TD[['FONT']])[['color']] == 'red') {
        isFP <- c(isFP, T)
      } else {
        isFP <- c(isFP, NA)
      }
    }
    i <- i + 1
  }
  # add identifiers for subsequent merge: group, motif len and upstr len
  toParse <- strsplit(filePath, '/|_')[[1]]
  group <- strsplit(grep('BRAIN|FB|GUT|MT', toParse, value = T), '_')[[1]][1]
  motifLen <- min(as.integer(grep('^\\d', toParse, value = T)))
  upstr <- max(as.integer(grep('^\\d', toParse, value = T)))
  result <- data.table(group = group, motifLen = motifLen, upstr = upstr,
                       Rank = 1:length(isFP), FP = isFP)
  result
}

#' parseHomerDeNovoHTML_v2
#' Parses HTML output of de-novo motifs from homer. Not compatible with v1.
#' @param filePath path to html 
#' @return data frame of if the de novo motif is false positive
parseHomerDeNovoHTML_v2 <- function(filePath) {
  # parse html 
  html <- xmlParse(filePath)
  # get to the table
  html <- xmlChildren(html)$HTML[[2]][['TABLE']]
  isFP <- c() # vector which will contain if it's false positive
  motifPath <- c() # path to motif file, need for merging with rest of data
  i = 2
  while(!is.null(html[[i]])) {
    # get motif paths
    oneMotifPath <- xmlAttrs(html[[i]][2]$TD[['IMG']])
    oneMotifPath <- gsub('homerResults.html', oneMotifPath, filePath)
    oneMotifPath <- gsub('logo.png', 'motif', oneMotifPath)
    motifPath <- c(motifPath, oneMotifPath)
    
    # getting to FP
    if (is.null(html[[i]][1]$TD[['FONT']])) {
      isFP <- c(isFP, F)
    } else {
      if (xmlAttrs(html[[i]][1]$TD[['FONT']])[['color']] == 'red') {
        isFP <- c(isFP, T)
      } else {
        isFP <- c(isFP, NA)
      }
    }
    i <- i + 1
  }
  result <- data.table(filePath = motifPath, FP = isFP)
  result
}

#' readHomerDeNovoIn
#' Reads-in one .motif file (result of de-novo motif discovery from Homer)
#' @param filePath path to the file
#' @return data table with columns motifLen, MotifName, Consensus, P-value,
#' q-value, PercTarg, PercBg, Rank
readHomerDeNovoIn <- function(filePath) {
  # read-in and leave only 1st line starting with ">" because it contains
  # info about motif
  motifInfo <- read.table(filePath, stringsAsFactors = F, nrows = 1)
  motifInfo <- motifInfo[c(1, 2, length(motifInfo))]
  motifInfo <- c(motifInfo, gsub('', '', gsub('%.*', '', motifInfo[3])))
  Rank <- gsub('.motif', '', gsub('.*/motif', '', filePath))
  Rank <- as.integer(Rank)
  group <- strsplit(filePath, '_')[[1]]
  group <- group[sapply(group, function(x) grepl('BRAIN|FB|GUT|MT', x))]
  group <- gsub('.*/', '' , group)
  motifInfoDT <- data.table(group = group,
                            motifLen = nchar(gsub('>', '', motifInfo[1])),
                            MotifName = gsub('>', '', motifInfo[1]), 
                            Consensus = gsub('>', '', motifInfo[1]),
                            Pvalue = as.numeric(gsub('.*:', '', 
                                                     motifInfo[3])),
                            qvalue = NA,
                            PercTarg = as.numeric(gsub('%.*', '', 
                                                       gsub('.*\\(', '', 
                                                            motifInfo[4]))),
                            PercBg = as.numeric(gsub('%.*', '', 
                                                     gsub('.*\\(', '', 
                                                          motifInfo[3]))),
                            Rank = Rank,
                            BestMatch = gsub('\\(.*', '', gsub('.*:', '',
                                                               motifInfo[2])))
  motifInfoDT[, MotifName := BestMatch]
  motifInfoDT[, MotifName := gsub('/.*', '', MotifName)]
  motifInfoDT[, FC := PercTarg / PercBg]
  colnames(motifInfoDT)[5:6] <- c('P-value', 'q-value')
  motifInfoDT$filePath <- filePath
  
  motifInfoDT
}

#' readHomerDeNovoAllIn
#' Reads in all .motif from all folders and subfolders of homerResDir
#' @param homerResDir directory full of homer results
#' @param enrichCode code for the enrichment performed with homer
#' @return data table with homer results
readHomerDeNovoAllIn <- function(homerResDir, enrichCode) {
  # read all the results files from all the folders
  resultDirs <- list.dirs(homerResDir, recursive = T, full.names = T)
  resultDirs <- resultDirs[grepl('homerResults', resultDirs)]
  deNovoResultsFiles <- sapply(resultDirs, 
                               function(x) list.files(x, '\\d.motif$',
                                                      recursive = T, 
                                                      full.names = T))
  deNovoResultsFiles <- unlist(deNovoResultsFiles) 
  deNovoResultsFiles <- deNovoResultsFiles[lapply(deNovoResultsFiles, 
                                                  length) != 0]
  deNovoResults <- lapply(deNovoResultsFiles, readHomerDeNovoIn)
  names(deNovoResults) <- deNovoResultsFiles
  deNovoResults <- deNovoResults[sapply(deNovoResults, 
                                        function(x) !is.na(x$Rank))]
 
  # result = summary all knownMotifs.txt
  allDeNovoMotifs <- c()
  # add info about upstream, amplitude quantile and BG quantile
  for (resultDir in names(deNovoResults)) {
    deNovoMotifsOne <- deNovoResults[[resultDir]]
    resultDir <- strsplit(resultDir, '//')[[1]][2]
    if (grepl('/', resultDir)) {
      resultDir <- strsplit(resultDir, '/')[[1]]
    }
    resultDir <- strsplit(resultDir, '_')[[1]]
    deNovoMotifsOne[, upstr := grep('^\\d', resultDir, value = T)]
    if (enrichCode == '_specCyclGenes_') {
      deNovoMotifsOne[, ampQT := resultDir[3]]
      deNovoMotifsOne[, BG := resultDir[4]]
    } else {
      deNovoMotifsOne[, ampQT := "all"]
      deNovoMotifsOne[, BG := resultDir[3]]
    }
    allDeNovoMotifs <- rbind(allDeNovoMotifs, deNovoMotifsOne)
  }
  # sometimes there are several motif lenght in one folder, so upstr get confused
  if (max(as.integer(allDeNovoMotifs$upstr) %% 100) != 0) {
    allDeNovoMotifs[, upstr := ifelse(as.integer(upstr) %% 100 != 0,
                                      as.integer(gsub('.$', '', upstr)),
                                      as.integer(upstr))]
  }
  
  # add info about if motif is false-positive
  allHtml <- sapply(resultDirs, 
                    function(x) list.files(x, 'homerResults.html',
                                           recursive = T, 
                                           full.names = T))
  deNovoHtml <- unlist(allHtml)
  names(deNovoHtml) <- deNovoHtml
  deNovoFP <- lapply(deNovoHtml, parseHomerDeNovoHTML)
  deNovoFP <- do.call(rbind, deNovoFP)
  deNovoFP <- deNovoFP[order(group, motifLen, upstr, Rank), ]
  allDeNovoMotifs <- allDeNovoMotifs[order(group, motifLen, upstr, Rank), ]
  allDeNovoMotifs <- cbind(allDeNovoMotifs, FP = deNovoFP$FP)
  
  allDeNovoMotifs 
}

#' readHomerDeNovoAllIn_v2
#' Reads in all .motif from all folders and subfolders of homerResDir. Not 
#' compatible with v1
#' @param homerResDir directory full of homer results
#' @param enrichCode code for the enrichment performed with homer
#' @return data table with homer results
readHomerDeNovoAllIn_v2 <- function(homerResDir) {
  # read all the results files from all the folders
  resultDirs <- list.dirs(homerResDir, recursive = T, full.names = T)
  resultDirs <- resultDirs[grepl('homerResults', resultDirs)]
  deNovoResultsFiles <- sapply(resultDirs, 
                               function(x) list.files(x, 'f\\d+.motif$',
                                                      recursive = T, 
                                                      full.names = T))
  deNovoResultsFiles <- unlist(deNovoResultsFiles) 
  deNovoResultsFiles <- deNovoResultsFiles[lapply(deNovoResultsFiles, 
                                                  length) != 0]
  deNovoResults <- lapply(deNovoResultsFiles, readHomerDeNovoIn)
  
  # decode the file names to get info about upstream, amplitude quantile and BG
  # quantile
  deNovoResultsInfo <- lapply(deNovoResultsFiles, deCodeFileName)
  
  # merge
  deNovoResults <- lapply(1:length(deNovoResults), 
                         function(x) if (nrow(deNovoResults[[x]]) != 0) {
                           cbind(deNovoResults[[x]], 
                                 deNovoResultsInfo[[x]])
                         } else {''})
  deNovoResults <- deNovoResults[sapply(deNovoResults, 
                                      function(x) is.data.frame(x))]
  deNovoResults <- do.call(rbind, deNovoResults)
  
  # if Best match is empty, then it's seq bias
  deNovoResults <- deNovoResults[BestMatch != '']
  
  # add info about if motif is false-positive
  allHtml <- sapply(homerResDir, 
                    function(x) list.files(x, 'homerResults.html',
                                           recursive = T, 
                                           full.names = T))
  deNovoHtml <- unlist(allHtml)
  names(deNovoHtml) <- deNovoHtml
  deNovoFP <- lapply(deNovoHtml, parseHomerDeNovoHTML_v2)
  deNovoFP <- do.call(rbind, deNovoFP)
  
  # merge
  setkey(deNovoFP, filePath)
  setkey(deNovoResults, filePath)
  deNovoResults <- merge(deNovoResults, deNovoFP)
  setnames(deNovoResults, "filePath", 'File')
  
  deNovoResults 
}

#' readHomerMotif
#' Reads-in result of homer motif map 
#' @param filePath path to the file
#' @return data table with columns FASTA, ID, Offset, Sequence, Motif, Name, 
#' Strand, MotifScore, chr, regStart, regEnd, start, end
#' where start and end is an absolute start and end of motif
readHomerMotif <- function(filePath) {
  motifDF <- fread(filePath, sep = '\t', header = T)
  motifDF[, chr := sapply(motifDF$`FASTA ID`, 
                          function(x) strsplit(x, ':')[[1]][1])]
  motifDF[, regStart := sapply(`FASTA ID`,
                               function(x) as.integer(strsplit(gsub('.*:', '',
                                                           x), '-')[[1]][1]))]
  motifDF[, regEnd := sapply(`FASTA ID`,
                              function(x) as.integer(strsplit(gsub('.*:', '',
                                                           x), '-')[[1]][2]))]
  motifDF[, start := ifelse(Strand == '+', 
                            regStart + ((regEnd - regStart) / 2) + Offset + 1,
                            regStart + ((regEnd - regStart) / 2) + Offset + 2 - 
                            nchar(Sequence))]
  motifDF[, end := start + nchar(Sequence) - 1]
  motifDF
}

#' writeBedForHomer
#' @param GOI vector with the names of genes of interest
#' @param lenUpstr upstream of the gene length
#' @param lenDownstr downstream of the gene length
#' @param excludeAroundTssUp # of bp to to exlude upsteam TSS
#' @param excludeAroundTssDown # of bp to to exlude downstrem TSS
#' @param outName output bed file name
writeBedForHomer <- function(GOI, lenUpstr, lenDownstr, 
                             excludeAroundTssUp = 50,
                             excludeAroundTssDown = 50, outName) {
  COI <- getGeneCoords(GOI) # Coordinates Of Interest
  message(paste('Lost ', sum(is.na(COI$chr))/nrow(COI), 'of genes'))
  COI <- COI[complete.cases(COI), ]
  # get region of interest excluding 50bp around TSS
  TSS_OI <- getTSScoords(COI, lenUpstr, lenDownstr, excludeAroundTssUp, 
                         excludeAroundTssDown)
  TSS_OI$chr <- gsub('chr', '', TSS_OI$chr)
  TSS_OI <- TSS_OI[order(TSS_OI$gene_name), ]
  TSS_OI[, 1] <- paste0('chr', TSS_OI[, 1])
  # target regions are ready
  write.table(TSS_OI[, 1:3], col.names = F, row.names = F, sep = '\t',
              file = outName, quote = F)
}

#' numberOfMotifs
#' Collects statistics about number of detected motifs
#' @param knownDF data table with known detected motifs
#' @param deNovoDF data table with de novo detected motifs
#' @return data table with number and type of the dected motifs
numberOfMotifs <- function(knownDF, deNovoDF) {
  numbKnown <- knownDF[, .N, by = .(group, TargName, TargAmpQ, BGname, BGExprQ,
                                    upstrLen, tss50)]
  numbKnown <- numbKnown[, Type := 'known']
  numbDeNovo <- deNovoDF[, .N, by = .(group, TargName, TargAmpQ, BGname,
                                      BGExprQ, upstrLen, tss50, FP)]
  numbDeNovo <- numbDeNovo[, FP := sapply(FP, 
                                          function(x) ifelse(x, '_FP', '_TP'))]
  numbDeNovo <- numbDeNovo[, Type := paste0('de-novo', FP)]
  numbDeNovo <- numbDeNovo[, FP := NULL]
  
  numbMotifs <- rbind(numbKnown, numbDeNovo)
  numbMotifs
}

#' motifToMeme
#' Converts a Homer .motif PWM/PFM into MEME format
#' @param inFile path to input file
#' @return list ready to be printed in MEME format
motifToMeme <- function(filePath) {
  # Reading the input file
  motif <- scan(file = filePath, character(0), sep = "\n", quote = NULL)
  motifName <- motif[1]
  
  # parcing of the matrix
  motifMatr <- motif[2:length(motif)] 
  motifMatr <- strsplit(motifMatr, split = "\t")
  motifMatr <- do.call(rbind, lapply(motifMatr, as.numeric))
  
  # geting name and p-value
  pval <- gsub('.*:', '', motifName)
  concensus <- gsub('>', '', strsplit(motifName, split = "\t")[[1]][1])
  motifName <- strsplit(motifName, split = "\t")[[1]][2]
  if (!grepl('BestGuess', motifName)) {
    motifName <- gsub('/.*', '', motifName)
    motifName <- gsub('\\(.*', '', motifName)
  } else {
    motifName <- gsub('\\(.*', '', gsub('.*:', '', motifName))
    motifName <- gsub('/.*', '', motifName)
  }
  
  result <- list(header = paste0('MOTIF ', concensus, ' ', motifName, '\n',
                                 'letter-probability matrix: alength= 4 w=',
                                 nrow(motifMatr), ' nsites= 20 E= ', pval),
                 matrix = motifMatr)
  result
}

#' printInMemeFormat
#' Prints all motifs in MEME format to be used with TomTom
#' @param motifFiles path to motif files
#' @param outFile path to the output file
printInMemeFormat <- function(motifFiles, outFile) {
  write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n", outFile)
  for(motifPath in motifFiles) {
    memeFormat <- motifToMeme(motifPath)
    write(memeFormat$header, outFile, append = T)
    write.table(memeFormat$matrix, outFile, append = T, sep = '\t', quote = F,
                col.names = F, row.names = F)
    write("\n", outFile, append = T)
  }
}

# PLOT MAPPING STATS ----------------------------------------------------------
shadesOfGreen <- "#0C8954"
shadesOfRed <- c("#FFD9D9", '#FFBABA', '#C48484', '#991D1D', '#730202', 
                 'black')

#' meltMappedUnmapped
#' Melts statsDF data frame to suitable for ggplot2 format to plot mapped/
#' unmapped stats
#' @param data.frame statsDF with statistics, columns: Name, NR_total, 
#' NR_mapped, NR_mapped_multLoci, NR_mapped_tooManyLoci, 
#' Perc_unmapped_mismatch, Perc_unmapped_short, Perc_unmapped_other, 
#' NR_dedupl
#' @return melted data frame
meltMappedUnmapped <- function(statsDF) {
  result <- data.frame(Name = rep(statsDF$Name, 6),
                       numbReads = c(unlist(statsDF[, 10:15, with = F])),
                       type = c(rep(c('uniquely mapped', 
                                      'mapped to mult. loci', 
                                      'mapped to too many loci',
                                      'NOTmapped - mismatches',
                                      'NOTmapped - too short',
                                      'NOTmapped - other'),
                                    each = nrow(statsDF))),
                       NR_total = rep(statsDF$NR_total, 6))
  result
}

#' plotMappedUnmapped
#' plots statistics about mapped and unmapped reads
#' @param statsMeltedDF result of function meltMappedUnmapped
#' @param percentage boolean, indicating, whatever percentage should be 
#' calculated
#' @return ggplot2
plotMappedUnmapped <- function(statsMeltedDF, percentage = F) {
  if (percentage) {
    yToPlot <- 100 * statsMeltedDF$numbReads/statsMeltedDF$NR_total
    axisBreaks <-  seq(0, 100, 10)
  } else {
    yToPlot <- statsMeltedDF$numbReads/1000000
    axisBreaks <-  seq(0, 1.1 * max(yToPlot), 1)
  }
  
  result <- ggplot(statsMeltedDF, aes(x = Name, y = yToPlot, fill = type)) + 
    geom_bar(stat = "identity") + 
    ggtitle('Number of uniquely mapped/unmapped reads') +
    xlab("Sample") + ylab("Number of reads, mlns") + mashaGgplot2Theme +
    scale_fill_manual("legend", 
                      values = c("NOT mapped" = shadesOfRed[1], 
                                 "uniquely mapped" = shadesOfGreen[1], 
                                 "mapped to mult. loci" = shadesOfRed[2],
                                 'mapped to too many loci' = shadesOfRed[6],
                                 "NOTmapped - mismatches" = shadesOfRed[3],
                                 'NOTmapped - too short' = shadesOfRed[4],
                                 'NOTmapped - other' = shadesOfRed[5])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = axisBreaks)
  result
}

#' meltDupl
#'  Melts statsDF data frame to suitable for ggplot2 format to plot 
#'  duplicated/not duplicated stats
#'  @param data.frame statsDF with statistics, columns: Name, NR_total, 
#' NR_mapped, NR_mapped_multLoci, NR_mapped_tooManyLoci, 
#' Perc_unmapped_mismatch, Perc_unmapped_short, Perc_unmapped_other, 
#' NR_dedupl
#' @return melted data frame
meltDupl <- function(statsDF) {
  statsDF$NR_dedupl <- as.numeric(as.character(statsDF$NR_dedupl))
  result <-  data.frame(Name = rep(statsDF$Name, 2),
                        numbReads = c(statsDF$NR_dedupl, 
                                      statsDF$NR_mapped - statsDF$NR_dedupl),
                        type = c(rep('not duplicated', nrow(statsDF)),
                                 rep('duplicated', nrow(statsDF))),
                        NR_mapped = rep(statsDF$NR_mapped, 2))
  result
}

#' plotDupl
#' plots statistics about duplicated and not duplicate reads
#' @param statsMeltedDF result of function meltDupl
#' @param percentage boolean, indicating, whatever percentage should be 
#' calculated
#' @return ggplot2
plotDupl <- function(statsMeltedDF, percentage = F) {
  if (percentage) {
    yToPlot <- 100 * statsMeltedDF$numbReads / statsMeltedDF$NR_mapped
    axisBreaks <-  seq(0, 100, 10)
  } else {
    yToPlot <- statsMeltedDF$numbReads / 1000000
    axisBreaks <-  seq(0, 1.1 * max(yToPlot), 1)
  }
  
  result <- ggplot(statsMeltedDF, aes(x = Name, y = yToPlot, fill = type)) +
    geom_bar(stat = "identity") + 
    ggtitle('Number of deduplicated reads') +
    xlab("Sample") + ylab("Number of reads, mlns") + mashaGgplot2Theme +
    scale_fill_manual("legend", 
                      values = c("duplicated" = shadesOfRed[5],
                                 "not duplicated" = shadesOfGreen[1]))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = axisBreaks)
  result
}

# PLOT COUNTS -----------------------------------------------------------------
mashaGgplot2Theme <- list(
  theme_classic(base_size = 20) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size = 0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "grey", size = 0.5,
                                          linetype = 2))
)

# color-blind friendly pallete
colfunc <- colorRampPalette(c("#D7191C", "#FDAE61", "#FFFFBF", "#ABD9E9", 
                              "#2C7BB6"))
colfunc <- heat.colors

#' LightenDarkenColor
#' @param col hex color value
#' @param value value to perform multiplicative decrease of lightness
#' @return darker/lighter color in hex
LightenDarkenColor <- function(col, value = 0.75) {
  cols1 <- readhex(file = textConnection(col), class = "RGB")
  # transform to hue/lightness/saturation colorspace
  cols1 <- as(cols1, "HLS")
  cols2 <- cols1
  # multiplicative decrease of lightness
  cols2@coords[, "L"] <- cols2@coords[, "L"] * value
  cols2 <- as(cols2, "RGB")
  cols2 <- hex(cols2)
}

#' Multiple plot function
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot 
#' objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' plotCountsTimecourse
#' Plots timeline of expression for a gene
#' @param geneName name of the gene to plot
#' @param inputTab normalized and batch corrected input count table, columns 
#' are samples and rows are genes gene names are in the row names
#' @param timePoints vector of the time then samples indicated in the column
#' names were taken
#' @param pointColor color of the points
#' @return plot
plotCountsTimecourse <- function(geneName, inputTab, timePoints, pointColor, 
                                 ...) {
  if (geneName %in% rownames(inputTab)) {
    dfToPlot <- data.frame(Time = timePoints, GeneExpr = inputTab[geneName, ],
                           stringsAsFactors = F)
  } else {
    dfToPlot <- data.frame(Time = timePoints, GeneExpr = 0,
                           stringsAsFactors = F)
    pointColor <- 'white'
  }
  
  thePlot <- ggplot(dfToPlot, aes(x = Time, y = GeneExpr)) + 
                    geom_point(shape = 16, colour = pointColor, 
                               fill = pointColor) + 
                    geom_smooth(span = 0.3, color = pointColor,
                                fill = pointColor) +
                    xlab('Time, h') + ylab('Normalized expression') +
                    # geom_text(label = colnames(inputTab)) +
                    mashaGgplot2Theme + ...
  thePlot
}

#' plotCCA
#' Plots projection of samples into first 2 canonical variables
#' @param ccaObj result of cc function
#' @param whatToPlot either white- or DGRP
#' @param timeOfSamples time point when sample were taken
#' @return plot
plotCCA <- function(ccaObj, whatToPlot = 'white-', timeOfSamples, ...) {
  dfToPlot <- data.frame(CanVar1 = ccaObj$scores$corr.X.xscores[, 1], 
                         CanVar2 = ccaObj$scores$corr.X.xscores[, 2])
  if (whatToPlot == 'DGRP') {
    dfToPlot <- data.frame(CanVar1 = ccaObj$scores$corr.Y.xscores[, 1], 
                           CanVar2 = ccaObj$scores$corr.Y.xscores[, 2])
  }
  # assign color for the samples
  samplColor <- colfunc(24)[round(timeOfSamples %% 24)]
  
  plot(dfToPlot, xlab = 'Canonical variable 1', ylab = 'Canonical variable 2',
       pch = 20, bty = 'n', col = samplColor, cex = 3, ...)
  text(dfToPlot$CanVar1, dfToPlot$CanVar2, labels = timeOfSamples %% 24, 
       pos = 3)
}


plotGeneExprAcrossTissues <- function(geneName, tissCountsList, 
                                      tissueColorVect) {
  exprDF <- data.frame(Tissue = character(), Expression = numeric())
  for (tissue in names(tissCountsList)) {
    tissueExpr <- tissCountsList[tissue]
    if (geneName %in% rownames(tissueExpr)) {
      toAdd <- data.frame(Tissue = tissue, Expression = tissueExpr[geneName, ])
    } else {
      toAdd <- data.frame(Tissue = tissue, Expression = NA)
    }
    exprDF <- rbind(exprDF, toAdd)
  }
  boxplot(Expression ~ Tissue, data = exprDF, lwd = 2, 
          ylab = 'Normalized expression', xlab = '',  bty = 'n', 
          cex.axis = 1.5, cex.main = 1.5, 
          main = paste('Expression of', geneName, 'across tissues'),
          col = tissueColorVect[names(tissCountsList)], border = 'black')
  stripchart(Expression ~ Tissue, vertical = TRUE, data = exprDF, 
             method = "jitter", add = T, pch = 20, col = 'black')
}

#' plotEQTL
#' @param geneName name of the gene CG42489
#' @param varCode id of the variant, i.e. 2L_5372_SNP
#' @param genotypeTab genotype table
#' @param geMatr gene expression matrix, genes in rows, samples in columns
#' @param infoTab info table
#' @return ggplot2
plotEQTL <- function(geneName, varCode, genotypeTab, geMatr, infoTab, 
                     tissColor) {
  bloomCodeGE <- sapply(colnames(geMatr), function(x) strsplit(x, '_')[[1]][4])
  commonSamples <- intersect(colnames(genotypeTab), bloomCodeGE)
  genotypeTab <- genotypeTab[, commonSamples]
  geMatr <- geMatr[, bloomCodeGE %in% commonSamples]
  
  exprDF <- data.frame(Expression = as.numeric(as.character(geMatr[geneName, ])),
                       Tissue = infoTab[colnames(geMatr), ]$Tissue)
  exprDF$Allele <- genotypeTab[varCode, 
                               sapply(colnames(geMatr), 
                                      function(x) strsplit(x, '_')[[1]][4])]
  exprDF <- exprDF[exprDF$Allele != '.', ]
  
  tissColor <- tissColor[as.character(unique(exprDF$Tissue))]
  colorToPlot <- as.vector(c(LightenDarkenColor(tissColor, 1.75), tissColor))
  result <- ggplot(exprDF, aes(x = Tissue, y = Expression, 
                               fill = interaction(Tissue, Allele), 
                               dodge = Allele)) + 
    geom_boxplot() + geom_jitter(position = position_dodge(0.8)) +
    scale_fill_manual(values = colorToPlot) + 
    ylab("Normalized expression") +
    ggtitle(paste('eQTL:', geneName, '&', varCode)) +
    mashaGgplot2Theme + theme(legend.position = "none")
  result
}

#' meanAcrossTimeReps
#' Calculates mean of gene expression across replicates of one time point
#' Helper function for plotCircadianHeatmap
#' @param exprOneTimePoint matrix/vector of gene expression across replicates
#'                         of one time point. Genes in rows, replicates in 
#'                         columns
#' @return vector of means
meanAcrossTimeReps <- function(exprOneTimePoint) {
  if (is.matrix(exprOneTimePoint)) {
    result <- apply(exprOneTimePoint, 1, mean)
  } else {
    result <- exprOneTimePoint
  }
  result
}

#' plotCircadianHeatmap
#' Plots classical circadian heatmap for 1 tissue
#' @param JTKres results of JTK 
#' @param exprMatr ecpression matrix, use standartized one
#' @param INFOTissMatr info table for 1 tissue
#' @param lightColor start color of gradient (other color is black)
#' @return heatmap
plotCircadianHeatmap <- function(JTKres, exprMatr, INFOTissMatr, lightColor, 
                                 ...) {
  # to get the pattern - order by phase (LAG)
  JTKres <- JTKres[order(JTKres$LAG), ]
  # get expression ordered by time
  exprInTime <- exprMatr[rownames(JTKres), ]
  setkey(INFOTissMatr, rightGT_name)
  INFOTissMatr <- INFOTissMatr[colnames(exprInTime)]
  INFOTissMatr <- INFOTissMatr[order(Time)]
  # take mean of different replicates across time points
  exprInTimeMean <- lapply(unique(INFOTissMatr$Time), 
                           function(x) meanAcrossTimeReps(exprInTime[, INFOTissMatr[Time == x]$rightGT_name]))
  exprInTimeMean <- do.call(rbind, exprInTimeMean)
  
  colorPallet <- colorRampPalette(c(lightColor, "black"))(n = 12)
  colorBreaks <- seq(0, 1, length.out = length(colorPallet) + 1)
  heatmap.2(t(exprInTimeMean), density.info = "none", trace = "none", 
            margins = c(4, 1), Colv = "NA", Rowv = "NA",
            col = colorPallet, breaks = colorBreaks, dendrogram = "none", 
            labRow = F, labCol = unique(INFOTissMatr$Time_code),
            lhei = c(1, 15), ...)
}

# PLOT MOTIFS -----------------------------------------------------------------
#' plotNumberMotifs
#' Plots number of motifs detected under each background
#' @param numbMotifsDT result of function numberOfMotifs
#' @return ggplot
plotNumberMotifs <- function(numbMotifsDT, plotTitle) {
  if (length(unique(numbMotifsDT$Type)) == 2) {
    colorScheme <- c("#B3DBC0", "#67BACA")
  } else {
    colorScheme <- c("#FD3C3C", "#B3DBC0", "#67BACA")
  }
  ggplot(numbMotifsDT, aes(x = as.factor(upstrLen), y = N, fill = Type)) +
    geom_bar(stat = 'identity', position = 'stack') + 
    xlab('Upstream TSS, bp') +  ylab('Number of motifs') +
    facet_grid(group + tss50 ~ BGname + BGExprQ) + mashaGgplot2Theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colorScheme) + ggtitle(plotTitle)
}

#' motifEnrichToMatr
#' Creates matrix to be ploted in heatmap
#' @param enrichMotif data table of enriched motifs
#' @param heatColName - name of the column, which will be in color, 
#' i.e. FC or q-value
#' @param noteColName -  name of the column, which will be in color, 
#' i.e. FC or q-value/perctargets
#' @return list with heatMatrix, cellNoteMatrix, labColours
motifEnrichToMatr <- function(enrichMotif, heatColName, noteColName) {
  # reduce the table to the columns which we will need
  enrReduc <- enrichMotif[, c('group', 'DrosoTF', heatColName, noteColName),
                          with = F]
  colnames(enrReduc)[3] <- 'heat'
  colnames(enrReduc)[4] <- 'note'
  if ('FP' %in% colnames(enrichMotif)) {
    enrReduc <- cbind(enrReduc, enrichMotif$FP)
    colnames(enrReduc)[5] <- 'FP'
  }
  enrReduc <- enrReduc[!is.na(heat)]
  setkey(enrReduc, group, DrosoTF)
  
  # number of rows/columns in the future 
  heatRow <- length(unique(enrReduc$DrosoTF))
  heatCol <- length(unique(enrReduc$group))
  
  # for heatmap matrix - fold change
  heatMatrix <- data.frame(matrix(NA, heatRow, heatCol))
  colnames(heatMatrix) <- unique(enrReduc$group)
  rownames(heatMatrix) <- unique(enrReduc$DrosoTF)
  
  # for cell not matrix - percentage of target
  cellNoteMatrix <- data.frame(matrix(NA, heatRow, heatCol))
  colnames(cellNoteMatrix) <- unique(enrReduc$group)
  rownames(cellNoteMatrix) <- unique(enrReduc$DrosoTF)
  
  # colors to label false positives
  labColours <- c()
  
  for (i in 1:heatRow) {
    for (j in 1:heatCol) {
      # sometimes in the same tissue same motif is identified twice
      value <- enrReduc[.(colnames(heatMatrix)[j], rownames(heatMatrix)[i]),
                           heat]
      # in case motif isn't found in the tissue
      if (!all(is.na(value))) { 
        valueMax <- which(value == max(value))[1]
        targetsPerc <- enrReduc[.(colnames(cellNoteMatrix)[j],
                                  rownames(cellNoteMatrix)[i]),
                                   note]
        heatMatrix[i, j] <- value[valueMax]
        cellNoteMatrix[i, j] <- round(targetsPerc[valueMax], 2)
      } else {
        heatMatrix[i, j] <- NA
        cellNoteMatrix[i, j] <- NA
      }
    }
    if ('FP' %in% colnames(enrReduc)) {
      fpStatus <- enrReduc[.(colnames(heatMatrix), rownames(heatMatrix)[i]),
                           FP]
      fpStatus <- fpStatus[!is.na(fpStatus)]
      if (all(fpStatus == F)) {
        labColours <- c(labColours, "#6FB98F")
      } else {
        labColours <- c(labColours, '#E4535E')
      }
    } else {
      if (rownames(heatMatrix)[i] == 'E-box') {
        labColours <- c(labColours, "#1995AD")
      } else {
        labColours <- c(labColours, 'black')
      }
    }
  }
  heatMatrix[is.na(heatMatrix)] <- 0
  heatMatrix <- as.matrix(heatMatrix)
  heatMatrix[is.infinite(heatMatrix)] <- max(na.omit(heatMatrix[is.finite(heatMatrix)]))
  result <- list(heatMatrix = as.data.frame(heatMatrix), 
                 cellNoteMatrix = cellNoteMatrix, labColours = labColours)
  result
}

#' createHeatMatrix
#' Creates matrix ready to be plotted as a heatmap from known and de-novo 
#' enriched motifs
#' @param knownMotifs known motif data table
#' @param deNovoMotifs de novo motif data table
#' @param heatColName - name of the column, which will be in color, 
#' i.e. FC or q-value
#' @param noteColName -  name of the column, which will be in color, 
#' i.e. FC or q-value/perctargets
createHeatMatrix <- function(knownMotifs, deNovoMotifs, heatColName, 
                             noteColName) {
  # in case known enriched motifs were found
  knownHeatMatr <- NULL
  if (!all(is.na(knownMotifs$group))) {
    knownHeatMatr <- motifEnrichToMatr(knownMotifs, heatColName, noteColName)
    message('Heatmap matrix for the known motifs was created')
  } 
  # in case de-novo enriched motifs were found
  deNovoHeatMatr <- NULL
  if (!all(is.na(deNovoMotifs$group))) { 
    deNovoHeatMatr <- motifEnrichToMatr(deNovoMotifs, heatColName, noteColName)
    # to distinguish from known
    rownames(deNovoHeatMatr$heatMatrix) <- paste0('de-novo:',
                                                  rownames(deNovoHeatMatr$heatMatrix))
    message('Heatmap matrix for the de novo motifs was created')
  }
  
  if (!is.null(knownHeatMatr) & !is.null(deNovoHeatMatr)) {
    heatMatrix <- plyr::rbind.fill(knownHeatMatr$heatMatrix,
                                   deNovoHeatMatr$heatMatrix)
    rownames(heatMatrix) <- c(rownames(knownHeatMatr$heatMatrix),
                              rownames(deNovoHeatMatr$heatMatrix))
    cellNoteMatrix <- plyr::rbind.fill(knownHeatMatr$cellNoteMatrix,
                                       deNovoHeatMatr$cellNoteMatrix)
    rowColors <- c(knownHeatMatr$labColours, deNovoHeatMatr$labColours)
  }
  
  if (is.null(knownHeatMatr)) {
    heatMatrix <- deNovoHeatMatr$heatMatrix
    cellNoteMatrix <-  deNovoHeatMatr$cellNoteMatrix
    rowColors <- deNovoHeatMatr$labColours
  }
  if (is.null(deNovoMotifs$group)) {
    heatMatrix <- knownHeatMatr$heatMatrix
    cellNoteMatrix <- knownHeatMatr$cellNoteMatrix
    rowColors <- knownHeatMatr$labColours
  }
  
  if (exists("heatMatrix")) {
    # one needs at least 2 rows and 2 columns to do a heatmap
    if (ncol(heatMatrix) == 1) {
      heatMatrix <- cbind(0, heatMatrix)
      cellNoteMatrix <- cbind(0, cellNoteMatrix)
    }
    if (nrow(heatMatrix) == 1) {
      heatMatrix <- rbind(0, heatMatrix)
      cellNoteMatrix <- rbind(0, cellNoteMatrix)
    }
    
    heatMatrix[is.na(heatMatrix)] <- 0
    heatMatrix <- as.matrix(heatMatrix)
    
    # sometimes inf may arise, replace it with the max not inf value
    maxVal <- as.vector(heatMatrix)
    maxVal <- maxVal[is.finite(maxVal)]
    maxVal <- max(maxVal)
    heatMatrix[is.infinite(heatMatrix)] <- maxVal
    
    result <- list(heatMatrix = heatMatrix, 
                   cellNoteMatrix = as.matrix(cellNoteMatrix),
                   rowColors = rowColors)
  } else {
    result <- NULL
  }
  result
}

#' plotMotifHeatmap
#' Plots heatmap of the known and de novo motif enrichment
#' @param knownMotifs known motif data table
#' @param deNovoMotifs de novo motif data table
#' @param heatColName name of the column, which will be in color, i.e. FC or 
#' q-value
#' @param noteColName name of the column, which will be in color, i.e. FC or 
#' q-value/perctargets
#' @param withCellNote if cell notes should be displayd
#' @return heatmap
plotMotifHeatmap <- function(knownMotifs, deNovoMotifs, heatColName = 'FC',
                             noteColName = 'PercTarg', withCellNote = T, 
                             colorScheme = c('white', "#8E9B97", "#5B7065",
                                             "#537072", "#304040", 
                                             "#2C4A52"), ...) {
  heatMatrix <- createHeatMatrix(knownMotifs, deNovoMotifs, heatColName, 
                                 noteColName)
  if (!is.null(heatMatrix)) {
    if (all(dim(heatMatrix$heatMatrix) > c(1, 1))) {
      heatColors <- colorRampPalette(colorScheme)(n = 100)
      colBreaks <- seq(0, max(heatMatrix$heatMatrix), 
                       length.out = length(heatColors) + 1)
      
      if (!withCellNote) {
        heatMatrix$cellNoteMatrix <- matrix(NA, 
                                            nrow(heatMatrix$cellNoteMatrix),
                                            ncol(heatMatrix$cellNoteMatrix))
      }
      
      par(cex.main = 0.75)
      heatmap.2(heatMatrix$heatMatrix, density.info = "none", trace = "none",
                margins = c(6, 10), dendrogram = "none", 
                cellnote = heatMatrix$cellNoteMatrix, 
                colRow = heatMatrix$rowColors, #Colv = "NA", Rowv = "NA",
                key = T, notecol = 'black', col = heatColors, 
                breaks = colBreaks, lmat = rbind(c(3, 0), c(1, 2), c(4, 0)),
                lwid = c(5, 0.1), lhei = c(0.5, 5, 0.85), ...)
    }
  }
  
  #heatColors <- c('white', '#ffffcc', '#ffeda0', '#fed976', '#feb24c',
  #                '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026')
}

#' plotMotifsGviz
#' Plots motif enrichment meta figure with use of Gviz
#' @param enrichMotsFltr filtered enriched motifs (do selection of motifs to 
#' plot for amp quantile, bg and len yourself)
#' @param numbTopMotifs number of top motifs to restrict the motifs to
#' @return Gviz plot
plotMotifsGviz <- function(enrichMotsFltr, numbTopMotifs) {
  enrichMotsFltr <- enrichMotsFltr[!duplicated(enrichMotsFltr[, -9,
                                                              with = F]), ]
  enrichMotsFltr[, MotifName := ifelse(MotifName == 'CTCF-SatelliteElement',
                                       'CTCF-SE', MotifName)]
  enrichMotsFltr[, MotifName := ifelse(MotifName == 'NF1-halfsite', 
                                       'NF1-HS', MotifName)]
  
  enrichMotsFltr[, FE := round(as.numeric(gsub('%', '', PercTarg)) / 
                                 as.numeric(gsub('%', '', PercBg)), 2)]
  # GATA and gata-like has a lot of motif hits, I want to replace them with one
  # best on q-value one
  enrichMotsFltrGATA <- enrichMotsFltr[grepl('GATA|Gata|Unknown5|PQM',
                                             MotifName), ]
  enrichMotsFltrGATA <- enrichMotsFltrGATA[order(group, `q-value`)][, .SD[1],
                                                                    group]
  enrichMotsFltrGATA$MotifName <- 'GATA-family'
  enrichMotsFltr <- rbind(enrichMotsFltrGATA,
                          enrichMotsFltr[!grepl('GATA|Gata|Unknown5|PQM',
                                                MotifName)])
  
  # Number of motifs occurence 
  enrichMotsFltrNumb <- enrichMotsFltr[,.N, by = MotifName]
  enrichMotsFltrNumb <- enrichMotsFltrNumb[order(-N), ]
  # make sure to include HLH-1
  enrichMotsFltrNumb <- rbind(enrichMotsFltrNumb[1:(numbTopMotifs - 1), ],
                              enrichMotsFltrNumb[MotifName == 'HLH-1', ])
  enrichMotsFltrNumb <- enrichMotsFltrNumb[complete.cases(enrichMotsFltrNumb)]
  # creating tracks for GVIZ
  trackList <- list()
  for (tissue in unique(enrichMotsFltr$group)) {
    enrichTiss <- enrichMotsFltr[group == tissue, ]
    setkey(enrichTiss, MotifName)
    enrichTiss <- enrichTiss[enrichMotsFltrNumb$MotifName]
    
    # if color is avaible - use it, if not, put grey
    if (is.na(tissueColor[tissue])) {
      trackColor <- 'darkgrey'
    } else {
      trackColor <- tissueColor[tissue]
    }
    aTrack.groups <- AnnotationTrack(start = c(20 * which(!is.na(enrichTiss$group)),
                                               20 * (topMotifs + 2)),
                                     width = c(rep(10, sum(!is.na(enrichTiss$group))), 
                                               100),
                                     chromosome = 'chr2L',
                                     group = c(enrichTiss[which(!is.na(enrichTiss$group))]$FE,
                                               paste0(tissue, '-specific cycling genes')),
                                     feature = c(enrichTiss[which(!is.na(enrichTiss$group))]$FE,
                                                 paste0(tissue, '-specific cycling genes')),
                                     genome = 'dm3', 
                                     strand = c(rep("*", sum(!is.na(enrichTiss$group))),
                                                "+"),
                                     name = tissue, fontcolor.item = 'black',
                                     background.title = trackColor,
                                     fill = trackColor)
    
    trackList[[length(trackList) + 1]] <- aTrack.groups
  }
  motifNamesTrack <- AnnotationTrack(start = 20*(1:nrow(enrichMotsFltrNumb)),
                                     width = 10, chromosome = 'chr2L',
                                     group = enrichMotsFltrNumb$MotifName,
                                     feature = enrichMotsFltrNumb$MotifName,
                                     genome = 'dm3', strand = '*', name = 'Motifs',
                                     fontcolor.item = 'black')
  trackList[[length(trackList) + 1]] <- motifNamesTrack
  plotTracks(trackList, groupAnnotation = "none", shape = 'arrow',
             featureAnnotation = "feature",
             fontsize.group = 65)
}

# GENE ANNOTATION -------------------------------------------------------------
#' initNameToFBidTab
#' initialize FBgn ID to the gene name database
#' files are downloaded from http://flybase.org/cgi-bin/get_static_page.pl?file=bulkdata7.html&title=Current%20Release
#' @return data table
initNameToFBidTab <- function() {
  # flybase FBgn to gene names
  flybaseFBgn <- fread('inputs_for_Rscripts/fbgn_annotation_ID_fb_2018_03.tsv',
                       header = T)
  colnames(flybaseFBgn) <- gsub('#', '', colnames(flybaseFBgn))
  colnames(flybaseFBgn) <- gsub('\\(.*', '', colnames(flybaseFBgn))
  setkey(flybaseFBgn, primary_FBgn)
  
  # flybase FBpp to FBgn
  flybaseFBpp <- fread('inputs_for_Rscripts/fbgn_fbtr_fbpp_fb_2018_03.tsv',
                       header = T, fill = T)
  colnames(flybaseFBpp)[1] <- 'primary_FBgn'
  setkey(flybaseFBpp, primary_FBgn)
  
  # merge them
  flybase <- merge(flybaseFBgn, flybaseFBpp)
  
  onTheFly <- fread('inputs_for_Rscripts/onTheFlyCodes.csv', header = T)
  colnames(onTheFly) <- cbind('primary_FBgn', 'gene_symbol')
  TFtoGene <- rbind(flybase, onTheFly, fill = T)
  TFtoGene
}

#' getGeneNamesIDsfromTab
#' Gets gene name based on Fbgn, FBpp, etc
#' @param gene something like Fbgn, FBpp, gene name
#' @param NameToFBidDB result of initNameToFBidTab
#' @param colToGet column(s) to return: primary_FBgn, gene_symbol, 
#' organism_abbreviation, secondary_FBgn, annotation_ID, 
#' secondary_annotation_ID, FlyBase_FBtr, FlyBase_FBpp
#' @return gene name if found, tfID + (not in Dmel DB) else
getGeneNamesIDsfromTab <- function(gene, NameToFBidDB, 
                                   colToGet = 'gene_symbol') {
  # scan though all columns to see which ID was submitted and if it's in the DB
  idType <- apply(NameToFBidDB, 2, function(x) sum(gene %in% x))
  if (sum(idType) == 0) {
    result <- paste0(gene, '(not found in Dmel DB)')
  } else {
    idType <- names(idType[which(idType != 0)[1]])
    result <- NameToFBidDB[which(NameToFBidDB[, idType, with = F] == gene)]
    if (length(colToGet) == 1) {
      result <- unique(as.character(result[, colToGet, with = F]))
      if (length(result) != 1) {
        message(paste('For', gene, 'multiple entries were found in the table,',
                      'taking first one'))
        result <- result[1]
      }
    } else {
      result <- result[, colToGet, with = F]
      result <- result[!duplicated(result)]
      if (nrow(result) != 1) {
        message(paste('For', gene, 'multiple entries were found in the table'))
      }
    }
  }
  result
}

#' getFlyBaseID 
#' @param externName vector of external gene names
#' @param ensembl91 ensembl object, the higher the version, the more genes will
#' be found
#' @return list, first element is data frame with genes for which flybase ID 
#' found and second one - vector for which it wasnt'
getFlyBaseID <- function(externName, ensembl91) {
  # found 
  egs <- getBM(attributes = c('flybase_gene_id', 'external_gene_name'),
               filters = 'external_gene_name', values = externName, 
               mart = ensembl91)
  # not found
  notFound <- externName[!externName %in% egs$external_gene_name]
  message(paste("Not translated to FlyBase ID "), 
          100 * round(length(notFound) / nrow(egs), 4), "% genes: ",
          paste(notFound, collapse = ", "))
  result <- list(egs, notFound)
  result
}

#' getGeneCoords
#' Returns coordinates of the gene and nothing else!
#' @param geneNames vector of gene names
#' @return data.table whith seqid, start, end, strand, gene_name
getGeneCoords <- function(geneNames) {
  result <- uscsRefGenCorrds[geneNames][, c('seqid', 'start', 'end', 'strand',
                                            'gene_name')]
  colnames(result)[1] <- c('chr')
  setkey(result, 'gene_name')
  result
}

#' getGenesStructure
#' Returns gene structure (promoter, 3'UTR, 5'UTR, exons, introns) for the 
#' input genes
#' @param genes genes to get structure for
#' @param txdbObj txdb made from GTF data was aligned to
#' @param promotUpstLen
#' @return GRanges object with gene structure
getGenesStructure <- function(genes, txdbObj, promotUpstLen = 2000) {
  # get corresponding transcript IDs
  trIDs <- select(txdbObj, keys = genes, columns = "TXNAME", 
                  keytype = "GENEID")
  trIDs <- as.data.table(trIDs)
  setkey(trIDs, TXNAME)
  # get all exons
  exonsForGene <- exonsBy(txdbObj, by = "gene")[genes]
  exonsForGene <- unlist(exonsForGene)
  mcols(exonsForGene)$tr_id <- gsub('\\..*', '', mcols(exonsForGene)$exon_name)
  message('Got exons')
  
  # get introns
  intronsForGene <- intronsByTranscript(txdbObj, use.names = T)[trIDs$TXNAME]
  tr_id <- rep(names(intronsForGene), sapply(intronsForGene, length))
  intronsForGene <- unlist(intronsForGene)
  mcols(intronsForGene)$tr_id <- tr_id
  names(intronsForGene) <- trIDs[names(intronsForGene)]$GENEID
  message('Got introns')
  
  # get 5'UTRs
  fiveUTRsForGene <- fiveUTRsByTranscript(txdbObj, use.names = T)
  fiveUTRsForGene <- unlist(fiveUTRsForGene)
  names(fiveUTRsForGene) <- fiveUTRsForGene$exon_name
  fiveUTRsForGene <- fiveUTRsForGene[intersect(exonsForGene$exon_name,
                                               names(fiveUTRsForGene))]
  names(fiveUTRsForGene) <- names(exonsForGene)[match(names(fiveUTRsForGene),
                                                      exonsForGene$exon_name)]
  mcols(fiveUTRsForGene)$tr_id <- gsub('\\..*', '',
                                       mcols(fiveUTRsForGene)$exon_name)
  message("Got 5' UTR")
  
  # get 3'UTRs
  threeUTRsForGene <- threeUTRsByTranscript(txdbObj, use.names = T)
  threeUTRsForGene <- unlist(threeUTRsForGene)
  names(threeUTRsForGene) <- threeUTRsForGene$exon_name
  threeUTRsForGene <- threeUTRsForGene[intersect(exonsForGene$exon_name,
                                                 names(threeUTRsForGene))]
  names(threeUTRsForGene) <- names(exonsForGene)[match(names(threeUTRsForGene),
                                                       exonsForGene$exon_name)]
  mcols(threeUTRsForGene)$tr_id <- gsub('\\..*', '',
                                        mcols(threeUTRsForGene)$exon_name)
  message("Got 3' UTR")
  
  # proximal & distal promoters
  proxPromForGenes <- promoters(txdbObj, upstream = promotUpstLen / 2, 
                                downstream = 0)
  proxPromForGenes <- proxPromForGenes[proxPromForGenes$tx_name %in%
                                         trIDs$TXNAME]
  names(proxPromForGenes) <- trIDs[proxPromForGenes$tx_name]$GENEID
  distPromForGenes <- promoters(txdbObj, upstream = promotUpstLen,
                                downstream = 0)
  distPromForGenes <- distPromForGenes[distPromForGenes$tx_name %in% trIDs$TXNAME]
  distPromForGenes <- sapply(1:length(proxPromForGenes), 
                             function(x) if (as.character(strand(proxPromForGenes[x])) == '+') {
                               disjoin(c(proxPromForGenes[x], 
                                         distPromForGenes[x]))[1]
                             } else {
                               disjoin(c(proxPromForGenes[x], 
                                         distPromForGenes[x]))[2]})
  distPromForGenes <- do.call(c, distPromForGenes)
  names(distPromForGenes) <- names(proxPromForGenes)
  mcols(distPromForGenes)$tx_name <- proxPromForGenes$tx_name
  message('Got promoters')
  
  mcols(exonsForGene) <- data.frame(annotation = 'Exon', 
                                    gene = names(exonsForGene), 
                                    tr_id = exonsForGene$tr_id,
                                    stringsAsFactors = F)
  mcols(intronsForGene) <- data.frame(annotation = 'Intron', 
                                      gene = names(intronsForGene), 
                                      tr_id = intronsForGene$tr_id,
                                      stringsAsFactors = F)
  mcols(fiveUTRsForGene) <- data.frame(annotation = "5' UTR", 
                                       gene = names(fiveUTRsForGene), 
                                       tr_id = fiveUTRsForGene$tr_id, 
                                       stringsAsFactors = F)
  mcols(threeUTRsForGene) <- data.frame(annotation = "3' UTR", 
                                        gene = names(threeUTRsForGene), 
                                        tr_id = threeUTRsForGene$tr_id, 
                                        stringsAsFactors = F)
  mcols(proxPromForGenes) <- data.frame(annotation = "Proximal promoter", 
                                        gene = names(proxPromForGenes), 
                                        tr_id = proxPromForGenes$tx_name, 
                                        stringsAsFactors = F)
  mcols(distPromForGenes) <- data.frame(annotation = "Distal promoter", 
                                        gene = names(distPromForGenes), 
                                        tr_id = distPromForGenes$tx_name, 
                                        stringsAsFactors = F)
  
  result <- c(exonsForGene, intronsForGene, fiveUTRsForGene,
              threeUTRsForGene, proxPromForGenes, distPromForGenes)
  names(result) <- NULL
  result <- result[!duplicated(as.data.frame(result))]
  names(result) <- result$gene
  result
}

#' getFeatureBoundary
#' Returns start / end of the granges object(s) 
#' @param gr Granges
#' @param bound what to return, "start" or "stop"
#' @return Granges 
getFeatureBoundary <- function(gr, bound = 'start') {
  negStr <- factor('-', levels = c('+', '-', '*')) # to simplify comparison
  # strands for all ranges
  allStrands <- as.factor(strand(gr))
  # get tss coordinates
  if (bound == 'start') {
    boundCoords <- sapply(1:length(gr),
                          function(x) if(allStrands[x] == negStr) {
                            end(gr[x])} else {start(gr[x])})
  } else {
    boundCoords <- sapply(1:length(gr),
                          function(x) if(allStrands[x] == negStr) {
                            start(gr[x])} else {end(gr[x])})
  }
  result <- data.frame(chr = as.character(seqnames(gr)),
                       start = boundCoords, end = boundCoords,
                       strand = allStrands)
  result <- makeGRangesFromDataFrame(result)
  result
}

#' getLocationInGene
#' @param geneStr structure for the gene of interest 
#' @param vars variants - eQTLs for this gene
#' @return Granges object with variants as ranges, genes as names and mcols
#' annotation, distToAnnotation, distToTSS
#' @note Transcript ID is very important! Without it you can't calculate proper
#' distance to TSS
getLocationInGene <- function(geneStr, vars) {
  # annotate variants which fall into gene features
  ovrl <- findOverlaps(vars, geneStr)
  if (length(ovrl) != 0) {
    inFeatures <- vars[ovrl@from]
    mcols(inFeatures) <- data.frame(annotation = geneStr[ovrl@to]$annotation, 
                                    tr_id = geneStr[ovrl@to]$tr_id,
                                    distToAnnotation = 0,
                                    stringsAsFactors = F)
  } else {
    inFeatures <- GRanges()
  }
  
  # ... and the ones which are not
  # note: do not use setdiff(vars, inFeatures), it mergers adjusent variants in
  # one!
  outFeatures <- vars[setdiff(1:length(vars), ovrl@from)]
  if (length(outFeatures) != 0) {
    outFeaturesDist <- distanceToNearest(outFeatures, geneStr)
    if (length(outFeaturesDist) != length(outFeatures)) {
      message(paste('Possibly variant from another chr!', names(geneStr)[1]))
    }
    mcols(outFeatures) <- data.frame(annotation = 'Downstream', 
                                     tr_id = geneStr[outFeaturesDist@to]$tr_id,
                                     stringsAsFactors = F)
    mcols(outFeatures)$annotation[grepl("promoter|5'UTR", 
                                        geneStr[outFeaturesDist@to]$annotation)] <- 'Upstream'
    mcols(outFeatures)$distToAnnotation <- outFeaturesDist@elementMetadata$distance
  } else {
    outFeatures <- GRanges()
  }
  
  # put them together
  result <- c(inFeatures, outFeatures)
  names(result) <- rep(names(geneStr)[1], length(result))
  strand(result) <- strand(geneStr)[1]
  
  # add distance to TSS
  tssCoord <- getFeatureBoundary(geneStr[geneStr$annotation == 
                                           'Proximal promoter'], bound = 'end')
  names(tssCoord) <- geneStr[geneStr$annotation == 'Proximal promoter']$tr_id
  mcols(result)$distToTSS <- sapply(result, 
                                    function(x) distance(x, tssCoord[x$tr_id]))
  
  result
}

#' getKnownCircGenes
#' Returns names of all known circadian genes
getKnownCircGenes <- function() {
  ensembl = useEnsembl("ensembl", version = 78, 
                       dataset = "dmelanogaster_gene_ensembl")
  egs = getBM(attributes = c('flybase_gene_id', 'external_gene_name', 
                             'name_1006'), values = '*', 
              mart = ensembl)
  CircGenes <- egs[grepl('circad|rhyt', egs$name_1006), ]
  CircGenes <- CircGenes[!duplicated(CircGenes[, c('flybase_gene_id',
                                                  'external_gene_name')]), ]
  CircGenes
}

#' getTSScoords
#' Returns coordinates of region surrounding TSS
#' @param coordsDT data table with coordinates, columns, chr, start, end, 
#' gene_name
#' @param before number of bp before TSS (upstream)
#' @param after  number of bp after TSS (downstream)
#' @param excludeBefore number of bp to exlude before TSS
#' @param excludeAfter number of bp to exlude after TSS
#' @return data.frame with coordiantes
getTSScoords <- function(coordsDT, before, after, excludeBefore = 0, 
                         excludeAfter = 0) {
  # get before bp tss excludung excludeBefore bp around TSS
  upstream <- data.frame(chr = coordsDT$chr,
                         start = ifelse(coordsDT$strand == "+", 
                                        coordsDT$start - before,
                                        coordsDT$end + excludeBefore),
                         end = ifelse(coordsDT$strand == "+", 
                                      coordsDT$start - excludeBefore,
                                      coordsDT$end + before))
  upstream <- cbind(upstream, gene_name = coordsDT$gene_name, 
                    strand = coordsDT$strand, type = 'upstream')
  # get bp after tss excludung bp around TSS
  downstream <- data.frame(chr = coordsDT$chr,
                           start = ifelse(coordsDT$strand == "+",
                                          coordsDT$start + excludeAfter,
                                          coordsDT$end - after),
                           end = ifelse(coordsDT$strand == "+",
                                        coordsDT$start + after, 
                                        coordsDT$end - excludeAfter))
  downstream <- cbind(downstream, gene_name = coordsDT$gene_name, 
                      strand = coordsDT$strand, type = 'downstream')
  result <- rbind(upstream, downstream)
  result
}

#' parseAnno
#' Simple parser for the annotation
#' @param anno string, SiteClass info for a variant
#' @param varID variant ID
parseAnno <- function(anno, varID) {
  anno <- gsub('SiteClass\\[', '', anno)
  anno <- gsub(']$', '', anno)
  oneVarAnno <- strsplit(anno, ';')[[1]]
  oneVarAnno <- lapply(oneVarAnno, function(x) strsplit(x, '\\|')[[1]])
  oneVarAnno <- do.call(rbind, oneVarAnno)
  oneVarAnno <- cbind(varID, oneVarAnno)
  if (ncol(oneVarAnno) != 5) {
    oneVarAnno <- cbind(oneVarAnno, "")
  }
  oneVarAnno <- as.data.table(oneVarAnno)
  colnames(oneVarAnno) <- c('SNP', 'flybaseID', 'gene', 'location',
                            'distance')
  oneVarAnno
}

#' getSiteClass
#' Parses annotation from fb57 for dgrp2 to extract SiteClass annotation
#' @param VOI vector with variants of interest
#' @param annoFilePath path to annotation file from dgrp2 
#' @return data table with columns VarID, flybaseID, geneName, 
#'         location(Intron/Exon/etc), distance
getSiteClass <- function(VOI, 
                         annoFilePath = 'inputs_for_Rscripts/dgrp.fb557.annot.txt') {
  # read in annotation file
  dgrpFB57anno <- fread(annoFilePath, sep = '\t', header = F)
  # restrict to the variants of the interest only
  dgrpFB57anno[, V1 := gsub('_SNP|_DEL|_INS|_MNP', '', V1)]
  dgrpFB57anno[, V1 := paste0('chr', V1)]
  dgrpFB57anno <- dgrpFB57anno[V1 %in% VOI]
  # remove TranscriptAnnot
  SiteClassAnno <- dgrpFB57anno$V3
  SiteClassAnno <- lapply(SiteClassAnno, 
                          function(x) gsub(',TranscriptAnnot.*', '', x))
  names(SiteClassAnno) <- dgrpFB57anno$V1
  
  # parse
  SiteClassAnnoTab <- lapply(1:length(SiteClassAnno), 
                             function(m) parseAnno(SiteClassAnno[[m]],
                                                   names(SiteClassAnno)[m]))
  SiteClassAnnoTab <- do.call(rbind, SiteClassAnnoTab)
  SiteClassAnnoTab <- SiteClassAnnoTab[!duplicated(SiteClassAnnoTab), ]
  SiteClassAnnoTab
}

#' getTranscriptAnnot
#' Parses annotation from fb57 for dgrp2 to extract TranscriptAnnot annotation
#' @param oneLineVars a data frame, directly read from the annotated with fb57
#'                    dgrp2 vcf (download from dgrp2 website and merged by sh 
#'                    script)
#' @param ralID RAL id of the line
#' @param bloomNumb bloomington number
#' @return data frame with columns RAL, bloomID, VarID, effect, geneName
getTranscriptAnnot <- function(oneLineVars, ralID, bloomNumb) {
  # extract column with annotation
  if ('INFO' %in% colnames(oneLineVars)) {
    TranscriptAnnotAnno <- oneLineVars$INFO
    # extract column with variant ID
    names(TranscriptAnnotAnno) <- oneLineVars$ID
  } else {
    TranscriptAnnotAnno <- oneLineVars$V8
    # extract column with variant ID
    names(TranscriptAnnotAnno) <- oneLineVars$V3
  }
  
  # remove SiteClass
  TranscriptAnnotAnno <- lapply(TranscriptAnnotAnno, 
                                function(x) gsub('.*TranscriptAnnot', '', x))
  TranscriptAnnotAnno <- lapply(TranscriptAnnotAnno, 
                                function(x) gsub('].*', '', x))
  TranscriptAnnotAnno <- lapply(TranscriptAnnotAnno,
                                function(x) gsub('\\[', '', x))
  
  TranscriptAnnotAnnoTab <- unlist(lapply(1:length(TranscriptAnnotAnno), 
                                          function(m) parseAnno(TranscriptAnnotAnno[[m]],
                                                                names(TranscriptAnnotAnno)[m],
                                                                c(1, 6))))
  TranscriptAnnotAnnoTab <- matrix(TranscriptAnnotAnnoTab,  ncol = 3, 
                                   byrow = T)
  TranscriptAnnotAnnoTab <- as.data.frame(TranscriptAnnotAnnoTab, 
                                          stringsAsFactors = F)
  
  colnames(TranscriptAnnotAnnoTab) <- c('VarID', 'effect', 'geneName')
  TranscriptAnnotAnnoTab <- TranscriptAnnotAnnoTab[TranscriptAnnotAnnoTab$geneName != '', ]
  TranscriptAnnotAnnoTab <- TranscriptAnnotAnnoTab[!duplicated(TranscriptAnnotAnnoTab), ]
  TranscriptAnnotAnnoTab$RAL <- ralID
  TranscriptAnnotAnnoTab$bloomID <- bloomNumb
  TranscriptAnnotAnnoTab
}

#' annoToTable
#' Converts annotation to table like format
#' @param vcfTab vcf table
#' @param loiID line of interest ID, in dgrp or bloomington form
#' @param RALtoBloomConvert ral to bloomington convertion table
annoToTable <- function(vcfTab, loiID, RALtoBloomConvert) {
  # if something, except numbers, is present - it's not Bloomington
  if (grepl('\\D', loiID) ) {
    loiRAL <- loiID
    loiBloom <- RALtoBloomConvert[line == loiRAL]$stock
  } else {
    loiBloom <- loiID
    loiRAL <- RALtoBloomConvert[stock == loiBloom]$line
  }
  
  # extract SiteClass info
  vcfSiteClass <- getSiteClass(vcfTab, loiRAL, loiBloom)
  vcfSiteClass <- as.data.table(vcfSiteClass)
  setkey(vcfSiteClass, RAL, bloomID, VarID, geneName)
  
  # extract TranscriptAnnot info
  vcfTranscriptAnnot <-  getTranscriptAnnot(vcfTab, loiRAL, loiBloom)
  vcfTranscriptAnnot <- as.data.table(vcfTranscriptAnnot)
  setkey(vcfTranscriptAnnot, RAL, bloomID, VarID, geneName)
  # merge TranscriptAnnot and SiteClass info
  vcfAnnot <- merge(vcfSiteClass, vcfTranscriptAnnot, all = T)
  vcfAnnot[, location := NULL]
  setkey(vcfTab, ID)
  vcfAnnot[, REF := vcfTab[vcfAnnot$VarID]$REF]
  vcfAnnot[, ALT := vcfTab[vcfAnnot$VarID]$ALT]
  vcfAnnot
}

#' GOenrichment
#' Performs GO enrichment of the selected gene set with TopGO
#' @param selectedGenes FlybaseIDs for the genes of interest
#' @param allGenesList background list of genes, selectedGenes should be
#' part of this list
#' @param ont ontology, BP, MF, CC
GOenricment <- function (selectedGenes, allGenesList, ont = 'BP',
                         topNodes = 10) {
  allGenesList_bg <- rep(1, length(allGenesList))
  names(allGenesList_bg) <- allGenesList
  allGenesList_bg[selectedGenes] <- 0.01
  
  tg.1 <- new("topGOdata", description = 'GO analysis',
              ontology =  ont, allGenes = allGenesList_bg,
              geneSel = function (x) {return (x < 0.05)},
              annot = annFUN.org ,
              nodeSize = 5 , # minimum number of genes in a GO categorie
              ID = "ENSEMBL", mapping = "org.Dm.eg.db") 
  GO.res <- runTest(tg.1, algorithm = "elim", statistic = "fisher")
  
  result <- GenTable(tg.1, Fisher = GO.res, orderBy = "Fisher", 
                     ranksOf = "Fisher", topNodes = topNodes)
  result
}

# LASSO -----------------------------------------------------------------------
#' diffInTime
#' Computes difference in time
#' @param timePoint1 integer, 0 <= x <= 24
#' @param timePoint2 integer, 0 <= x <= 24 
#' @return integer, difference in time
diffInTime <- function(timePoint1, timePoint2) {
  biggerTime <- max(timePoint1, timePoint2)
  smallerTime <- min(timePoint1, timePoint2)
  
  clockWise <- biggerTime - smallerTime
  counterClockWise <- 24 - biggerTime + smallerTime
  result <- min(clockWise, counterClockWise)
  result
}

#' divideIntoTrainAndTest
#' Divides data frame into test and train data sets
#' @param whiteDgrpProfDF initial data frame with white and profiled DGRPs
#' @param infoDF data frame containing time of sampling
#' @param leaveOneOut index of the sample to test on in case of leave one out
#' strategy. Not compatible with mixed!
#' @param mixed boolean, indicating, whatever white and dgrps should be mixed
#' @param genesOI genes of interest, if any
#' @param timePoint1 rigth time point of the interval of interes, if any
#' @param timePoint2 left time point of the interval of interest, if any
#' @param timeIntervalType string, "inner" or "outer". Between or outside time
#' points
#' @return list with 1) train data set 2) test data set
divideIntoTrainAndTest <- function(whiteDgrpProfDF, infoDF, leaveOneOut = 0, 
                                   mixed = F, genesOI = c(), timePoint1 = NA,
                                   timePoint2 = NA, 
                                   timeIntervalType = 'inner') {
  whiteDgrpProfDF <- whiteDgrpProfDF[, complete.cases(t(whiteDgrpProfDF))]
  # leave one ou dominates everything else
  if (leaveOneOut != 0) {
    trainingSet <- whiteDgrpProfDF[-leaveOneOut, ]
    testSet <- whiteDgrpProfDF[leaveOneOut, ]
    # now, select time interval, if needed
    if (!is.na(timePoint1) & !is.na(timePoint2)) {
      if (timeIntervalType == 'inner') {
        selectSamp <- (infoDF[rownames(trainingSet)]$Time %% 24 <= timePoint2) &
          (infoDF[rownames(trainingSet)]$Time %% 24 >= timePoint1)
      }
      if (timeIntervalType == 'outer') {
        selectSamp <- (infoDF[rownames(trainingSet)]$Time %% 24 >= timePoint2) |
          (infoDF[rownames(trainingSet)]$Time %% 24 <= timePoint1)
      }
      trainingSet <- trainingSet[selectSamp, ]
    }
  } else {
    # select time interval
    if (!is.na(timePoint1) & !is.na(timePoint2)) {
      if (timeIntervalType == 'inner') {
        selectSamp <- (infoDF[rownames(whiteDgrpProfDF)]$Time %% 24 <= timePoint2) &
          (infoDF[rownames(whiteDgrpProfDF)]$Time %% 24 >= timePoint1)
      }
      if (timeIntervalType == 'outer') {
        selectSamp <- (infoDF[rownames(whiteDgrpProfDF)]$Time %% 24 >= timePoint2) |
          (infoDF[rownames(whiteDgrpProfDF)]$Time %% 24 <= timePoint1)
      }
      whiteDgrpProfDF <- whiteDgrpProfDF[selectSamp, ]
    }
    if (mixed) {
      index <- sample(1:nrow(whiteDgrpProfDF),
                      round(0.75 * nrow(whiteDgrpProfDF)))
    } else {
      index <- grepl('\\d{5}', rownames(whiteDgrpProfDF))
    }
    trainingSet <- whiteDgrpProfDF[index, ]
    testSet <- whiteDgrpProfDF[-index, ]
  }
  result <- list(as.matrix(trainingSet), as.matrix(testSet))
  names(result) <- c('train', 'test')
  result
}

#' lassoOnTrainSet
#' Performs LASSO regression on train set
#' @param trainDF train data frame, samples in rows, genes in columns
#' @param trainTimeLine time when samples were taken
#' @param doPlots whatever plots shouls be done
#' https://jasdumas.github.io/tech-short-papers/glmnet_lasso_tutorial.html
lassoOnTrainSet <- function(trainDF, trainTimeLine, doPlots = F, ...) {
  # different values of alpha return different estimators
  # alpha = 1 is the lasso
  fit <- glmnet(x = trainDF, y = trainTimeLine, alpha = 1, ...)
  # crossvalidation
  crossval <- cv.glmnet(x = trainDF, y = trainTimeLine, ...)
  penalty <- crossval$lambda.min #optimal lambda
  #estimate the model with optimal lambda
  result <- glmnet(x = trainDF, y = trainTimeLine, alpha = 1, lambda = penalty,
                   ...)
  if (doPlots) {
    plot(fit, xvar = "lambda", main = 'Visualizing of the lasso coefficients')
    plot(crossval, main = 'Crossvalidation')
  }
  result
}

#' MSEforClock
#' Calculates mean square error for the clock setting
#' @param estimatedTime time given by the model
#' @param realTime real time of harversting
#' @return MSE vector
MSEforClock <- function(estimatedTime, realTime) {
  df <- data.frame(estimatedTime, realTime)
  timeDiff <- apply(df, 1, function(x) diffInTime(x[1], x[2]))
  result <- sum(timeDiff ^ 2)/length(timeDiff)
  result
}

#' plotModelCoeffs
#' Plots coefficients in front of the genes for models
#' @param allModelCoeffs data frame containing columns tissue, Gene, Coeff, 
#' Cycles
#' @param modelName model name
#' @return ggplot2
plotModelCoeffs <- function(allModelCoeffs, modelName) {
  # plot genes with the biggest impact on model
  allModelCoeffs$Coeff <- as.numeric(allModelCoeffs$Coeff)
  modelGenRank <- allModelCoeffs[,.(sum(abs(Coeff))), by = Gene]
  modelGenRank <- modelGenRank[order(-abs(V1))][Gene != '(Intercept)']
  modelGenRank <- modelGenRank[1:30]
  allModelCoeffs <- allModelCoeffs[Gene %in% modelGenRank$Gene]
  allModelCoeffs$Gene <- factor(allModelCoeffs$Gene, 
                                levels = modelGenRank$Gene)
  
  result <- ggplot(allModelCoeffs, aes(Gene, tissue, colour = Cycles)) +
    geom_point(aes(size = abs(allModelCoeffs$Coeff))) + 
    mashaGgplot2Theme + scale_size_continuous(range = c(1, 15)) +
    scale_color_manual(values = c("#C02942", "#3B8686")) +
    labs(size = "Coefficient") + xlab("Gene") + ylab("Tissue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste('Coefficients for the', modelName))
  result
}

# Molecular Time Table ---------------------------------------------------------
#' cosineBundle
#' Returns bundle of cosine waves with designated period on timePoints with
#' incriment as matrix
#' @param timePoints time points on which the wave should be generated, no need
#' to order, duplicated values are fine
#' @param period period of cosine wave
#' @param incriment incriment, default = 10 min
#' @return matrix
cosineBundle <- function(timePoints, period = 24, incriment = 1/6) {
  perMult <- 2 * pi / period
  result <- lapply(0:((period / incriment) - 1), 
                   function(x) cos(perMult * (timePoints - x * incriment)))
  result <- do.call(rbind, result)
}

#' plotCosineBundle
#' Plots cosine bundle
#' @param timePoints time points on which the wave should be generated, no need
#' to order, duplicated values are fine
#' @param period period of cosine wave
#' @param incriment incriment, default = 10 min
#' @return lots of plots
plotCosineBundle <- function(timePoints, period = 24, incriment = 1/6) {
  cosBundlMatr <- cosineBundle(timePoints, period, incriment)
  cosBundlPlot <- data.frame(Time = timePoints, t(cosBundlMatr))
  phases <- round(incriment*(1:nrow(cosBundlMatr)), 2)
  colnames(cosBundlPlot) <- c('Time', paste0('P', phases))
  cosBundlPlot <- data.frame(Time = rep(cosBundlPlot$Time, 
                                        ncol(cosBundlPlot) - 1),
                             melt(cosBundlPlot[, -1]))
  result <- ggplot(cosBundlPlot, aes(x = Time, y = value, color = variable)) + 
    geom_point(shape = 20, cex = 2) + mashaGgplot2Theme +
    theme(legend.position = "none")
  print(facet_multiple(plot = result, facets = 'variable', ncol = 12, 
                       nrow = 12))
}

#' corCosBundleMolPeakTime
#' Computes best corelation with cosine bundle and molecular peak time
#' @param timePoints time points on which the wave should be generated, no need
#' to order, duplicated values are fine
#' @param geneExpr matrix of gene expression
#' @param period period of cosine wave
#' @param incriment incriment, default = 10 min
#' @return vector with maximum correlation and molecular peak time (MPT)
corCosBundleMolPeakTime <- function(timePoints, geneExprMatr, period = 24,
                                    incriment = 1/6) {
  # create cos bundle matrix
  cosBundlMatr <- cosineBundle(timePoints, period, incriment)
  # correlation
  cosCor <- apply(geneExprMatr, 2, function(x) apply(cosBundlMatr,
                                                     1, cor, x))
  # maximum correlation for every gene
  cosCorMax <- apply(cosCor, 2, max)
  # calculate Molecular Peak Time
  mptPhase <- apply(cosCor, 2, function(x) which(x == max(x))[1])
  mptPhase <- incriment * mptPhase# because incriment = shift to the right
  
  result <- data.frame(Gene = colnames(geneExprMatr), maxCor = cosCorMax,
                       MPT = mptPhase)
  result
}

#' calcBodyTime
#' Calculates body time 
#' @param MPTs Molecular Peak Times
#' @param TIGexprSamp time indicating genes expression in the sample
#' @param period period of cosine wave
#' @param incriment incriment, default = 10 min
#' @return body time of the sample, according to time indicating genes
calcBodyTime <- function(MPTs, TIGexprSamp, period = 24, incriment = 1/6) {
  # create cos bundle matrix
  cosBundlMatr <- cosineBundle(MPTs, period, incriment)
  # correlation
  cosCor <- apply(cosBundlMatr, 1, cor, TIGexprSamp)
  # maximum correlation over all cos curve
  mptPhase <- which(cosCor == max(cosCor))[1]
  # calculate body time
  bodyTime <- incriment * mptPhase# because incriment = shift to the right
  bodyTime
}

# Dirertional statistics ------------------------------------------------------
#' hoursToRadians
#' Converts time in hours to time in radians
#' @param x vector or numeric, time in hours, can be > 24h
#' @return vector of radians
hoursToRadians <- function(x) {
  # hours in 24h system
  y <- x %% 24
  # hours in degrees
  y <- y * (360 / 24)
  # degrees in radians
  y <- (pi / 180) * y
  y
}

#' Converts time in hours to time in radians
#' @param x vector or numeric, time in hours, can be > 24h
#' @return vector of radians
radiansToHours <- function(x) {
  # radians to degrees
  y <- x * (180 / pi)
  # degrees to hours
  y <- y * (24 / 360)
  y
}

fitDirectionalModel <- function(trainGeneMatr, trainPointsOfTime, testGeneVect, 
                                testPointsOfTime) {
  # time points to radians
  pointsOfTimeRads <- hoursToRadians(trainPointsOfTime)
  # Forward Backward Early Dropping selection for circular data using the SPML 
  # regression
  selectGenes <- spml.fbed(pointsOfTimeRads, trainGeneMatr)
  selectGenes <- selectGenes$res[, 1]
  # Fit model
  directModel <- spml.reg(pointsOfTimeRads, trainGeneMatr[, selectGenes],
                          rads = T, xnew = t(testGeneVect[selectGenes]))
  # predicted time is now in radians, convert it to hours 
  estTime <- radiansToHours(directModel$est)
  # result 
  result <- c(testPointsOfTime, estTime)
}


# RcisTarget ------------------------------------------------------------------
#' applyRcisTarget
#' Applies RcisTarget to tha named list inList
#' @param inList input NAMED list
#' @param returnFull if T, returns list with motifsAUC, motifEnrichTab
applyRcisTarget <- function(inList, returnFull = F) {
  # 1. Calculate AUC
  motifsAUC <- calcAUC(inList, motifRankings)
  
  # 2. Select significant motifs, add TF annotation & format as table
  motifEnrichTab <- addMotifAnnotation(motifsAUC, 
                                       motifAnnot = motifAnnotations_dmel_v8)
  
  # 3. Identify significant genes for each motif
  # (i.e. genes from the gene set in the top of the ranking)
  # Note: Method 'iCisTarget' instead of 'aprox' is more accurate, but slower
  motifEnrichTab <- addSignificantGenes(motifEnrichTab, geneSets = inList,
                                        rankings = motifRankings, nCores = 1,
                                        method = "iCisTarget")
  motifEnrichTab <- addLogo(motifEnrichTab)
  
  if (returnFull) {
    result <- list(motifsAUC, motifEnrichTab)
  } else {
    result <- motifEnrichTab
  }
  result
}

# GENIE3 ----------------------------------------------------------------------
#' formatExprTSForDynGENIE3
#' Prepares list of gene expression and time points to be submitted to 
#' dynGENIE3. DynSCENIC can not handle replicates then they are in the same 
#' matrix, we need to put them in the different matrices and create lists out
#' of matrices
#' @param inCounts list of expression matrices, one matrix per tissue
#' @param infoDT table with support information
#' @param repNumb number of replicates
#' @return list with list of expression matrices, and time points
formatExprTSForDynGENIE3 <- function(inCounts, infoDT, repNumb) {
  allTiss <- sort(unique(infoDT$Tissue))
  setkey(infoDT, rightGT_name)
  
  timePoints <- lapply(inCounts, function(x) infoDT[colnames(x)]$Time)
  timePoints <- lapply(timePoints, function(x) x %% 24)
  # dynSCENIC can not handle replicates then they are in the same matrix, we need
  # to put them in the different matrices and create lists out of matrices
  timeSeriesMatr <- list()
  timeSeriesTPs <- list()
  for (tis in allTiss) {
    # get time points for the tissue with replicates
    tisTimePoints <- timePoints[[tis]]
    # get time points for the tissue without replicates
    timeTicks <- sort(unique(tisTimePoints))
    # get samples which correspond to the replicates of time points
    timeTicksReps <- lapply(timeTicks, function(x) which(tisTimePoints == x))
    # restrict to the number of replicates per time point set before
    timeTicksReps <- lapply(timeTicksReps, 
                            function(x) if(length(x) > repNumb) {x[1:repNumb]}
                            else {x})
    names(timeTicksReps) <- timeTicksReps
    # result list for the tissue time series matrices
    tisTimeSeriesMatr <- list()
    for (i in 1:repNumb) {
      timeSeries <- matrix(nrow = nrow(inCounts[[tis]])) # expression
      for (j in 1:length(timeTicksReps)) {
        if (length(timeTicksReps[[j]]) >= i) {
          timeSeries <- cbind(timeSeries, 
                              inCounts[[tis]][, timeTicksReps[[j]][i]])
          sampleName <- colnames(inCounts[[tis]])[timeTicksReps[[j]][i]]
          colnames(timeSeries)[ncol(timeSeries)] <- sampleName
          # I don't  fill in time points here, because it's going to be one more
          # sanity check
        }
      }
      timeSeries <- timeSeries[, -1] # remove NAs
      tisTimeSeriesMatr[[length(tisTimeSeriesMatr) + 1]] <- timeSeries
    }
    # time points 
    tisTimeSeriesTPs <- lapply(tisTimeSeriesMatr,
                               function(x) infoDT[colnames(x)]$Time %% 24)
    
    # add to the big list of lists of matrices
    timeSeriesMatr[[length(timeSeriesMatr) + 1]] <- tisTimeSeriesMatr
    timeSeriesTPs[[length(timeSeriesTPs) + 1]] <- tisTimeSeriesTPs
  }
  names(timeSeriesMatr) <- allTiss
  names(timeSeriesTPs) <- allTiss
  
  result <- list(tsExpr = timeSeriesMatr, ts = timeSeriesTPs)
}

#' plotCircNet
#' Plots network of tissue-specific circadian expression, inferred with GENIE3
#' @param tissue tissue
#' @param weightMatr weighted matrix
#' @param numbOfLinks number of top links to include in the network
#' @param cyclExprInf data table with info about tissue-specific cycling and
#'                    expression
#' @param netColors colors to give to the net
#' @return network 
plotCircNet <- function(tissue, weightMatr, numbOfLinks, cyclExprInf, 
                        tisBaseColors, genesTFs) {
  # set key 
  setkey(cyclExprInf, Gene)
  
  # get colors for all possible combinations
  netColors <- c(tisBaseColors, rep("#F2DA8A", 6), rep("#E6AC27", 4),
                 'snow4')
  names(netColors)[-4:-1] <- unlist(sapply(2:4, 
                                           function(x) apply(combn(names(tisBaseColors), x),
                                                             2, paste, collapse = '-')))
  
  # get link lists from GENIE3 result
  linkList <- getLinkList(weightMatr, reportMax = numbOfLinks)
  linkList$regulatoryGene <- as.character(linkList$regulatoryGene)
  linkList$targetGene <- as.character(linkList$targetGene)
  
  # create color and shape info about the network
  netInfo <- data.table(name = unique(c(linkList$regulatoryGene,
                                        linkList$targetGene)))
  netInfo[, color := sapply(name, function(x) cyclExprInf[x]$CyclesIn)]
  netInfo[, color := netColors[color]]
  netInfo[, shape := sapply(name, 
                            function(x) if(cyclExprInf[x]$ExprIn == tissue) {
                              'triangle'
                            } else {
                              if(x %in% genesTFs$external_gene_name) 
                              {'square'} else {'circle'}
                            })]
  netInfo[, nameModified := sapply(name, 
                                   function(x) if(x %in% genesTFs$external_gene_name) {x}
                                   else {NA})]
  
  # create net itself
  net <- graph_from_data_frame(linkList, directed = T, vertices = netInfo)
  layout <- layout_with_fr(net)
  
  # plot
  add_shape("triangle", clip = shapes("circle")$clip, plot = mytriangle)
  plot(net, layout = layout, vertex.size = 3, vertex.label.degree = 0,
       vertex.label.color = "black", vertex.label.cex = 1.2,
       vertex.color = V(net)$color,  edge.arrow.size = 0.1,
       edge.width = 10 * E(net)$weight, vertex.shape = V(net)$shape, 
       edge.curved = .4, vertex.label = V(net)$nameModified,
       edge.color = 'slategrey', main = tissue)
}

#' triangle vertex shape for igraph
#' @param coords
mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}

# MISC-------------------------------------------------------------------------
#' calcMaf
#' Calculates minor allele frequnecy (MAF) for one variant
#' @param oneVarFromGenoTab one line of GT from vcf
#' @return named vector with mac and maf 
calcMaf <- function(oneVarFromGenoTab) {
  # alternative allele count
  altCount <- length(oneVarFromGenoTab[oneVarFromGenoTab != '.' &
                                         oneVarFromGenoTab != '0/0'])
  # count also reference allele, to circumvent troubles with NAs
  refCount <- length(oneVarFromGenoTab[oneVarFromGenoTab == '0/0'])
  # minor allele count
  mac <- min(altCount, refCount)
  result <- mac/(altCount + refCount)
  result
}

#' displayColorPallet
#' @param colVect vector with HEX colors
#' @return color pallet dysplay
displayColorPallet <- function(colVect) {
  image(1:length(colVect), 1, col = colVect, axes = F,
        as.matrix(1:length(colVect)), xlab = "",
        ylab = "")
}
