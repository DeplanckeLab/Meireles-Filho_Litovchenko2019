# FILE: 12_inputsForHomer.R ---------------------------------------------------
#
# DESCRIPTION : 
# Prints beds for promoters of tissue specific cycling genes  and tissue 
# specific expressed genes for motif enrichment with Homer
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

# INPUTS ----------------------------------------------------------------------
# list all tissues and corresponding colors
# info about samples
infoTab <- readRDS('Rds/sampleInfo.Rds')
infoTab <- as.data.table(infoTab)
allTissues <- unique(infoTab$Tissue)

# output dirs
saveRDSdir <- 'Rds/'
savePlotsDir <- 'plots/'
saveTSEbedDir <- 'MotifEnrichment/input_TSE_beds/'
dir.create(saveTSEbedDir)
saveTSCbedDir <- 'MotifEnrichment/input_TSC_beds/'
dir.create(saveTSCbedDir)
saveCiTbedDir <- 'MotifEnrichment/input_CiT_beds/'
dir.create(saveCiTbedDir)

# data frame with info about where gene cycles and if it has peak
tisSpecCycl <- readRDS(paste0(saveRDSdir, 'tisSpecCycl.Rds'))

# Tissue specific expresion
TSE <- readRDS('Rds/TSE.Rds')

# all lengths to scan
upstrVect <- c(300, 500, 1000)
downstrVect <- c(150, 500, 500)
skipTSSLen <- c(0, 50)
ampNumbQuant <- 5
exprNumbQuant <- 5

# BEDs FOR HOMER MOTIF ENRICHMENT: list of all combinations -------------------
# HOMER - 1:
#   Target: TISSUE-SPECIFIC EXPRESSED GENES
#   BGs: RANDOM
#   EXPECTATIONS: TISSUE-SPECIFIC MOTIFS
#   QUANTILES OF AMPLITUDE (TARGETS): Q1-Q5 (this case it's quantiles of 
#                                     expression)
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT

# HOMER - 2:
#   Target: TISSUE-SPECIFIC EXPRESSED GENES
#   BGs: ALL EXPRESSED GENES
#   EXPECTATIONS: TISSUE-SPECIFIC MOTIFS
#   QUANTILES OF AMPLITUDE (TARGETS): Q1-Q5 (this case it's quantiles of 
#                                     expression)
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT

# HOMER - 3:
#   Target: TISSUE-SPECIFIC EXPRESSED GENES
#   BGs: GENES TISSUE-SPECIFICALLY EXPRESSED FOR OTHER TISSUES
#   EXPECTATIONS: TISSUE-SPECIFIC MOTIFS
#   QUANTILES OF AMPLITUDE (TARGETS): Q1-Q5 (this case it's quantiles of 
#                                     expression)
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT

# HOMER - 4:
#   Target: TISSUE-SPECIFIC CYCLING GENES
#   BGs: RANDOM
#   EXPECTATIONS: COMMON CIRCADIAN MOTIFS, tissue-specific circadian motifs 
#                 could be found if they are ubiquitously present
#   QUANTILES OF AMPLITUDE (TARGETS): Q1 - Q5
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT

# HOMER - 5:
#   Target: TISSUE-SPECIFIC CYCLING GENES
#   BGs: ALL EXPRESSED GENES
#   EXPECTATIONS: COMMON CIRCADIAN MOTIFS, tissue-specific circadian motifs 
#                 could be found because there will be genes cycling in other 
#                 tissues in the background too
#   QUANTILES OF AMPLITUDE (TARGETS): Q1 - Q5
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT
#   Note: circadian motifs might not be found because there will be genes in 
#         the bg, which cycling in other tissues

# HOMER - 6:
#   Target: TISSUE-SPECIFIC CYCLING GENES
#   BGs: GENES EXPRESSED IN THE TISSUE, BUT NOT CYCLING
#   EXPECTATIONS: COMMON CIRCADIAN MOTIFS, tissue-specific circadian motifs 
#                 could be found if they are ubiquitously present
#   QUANTILES OF AMPLITUDE (TARGETS): Q1 - Q5
#   QUANTILES OF EXPRESSION (BGs): Q1 - Q5
#   50bp TSS: WITH AND WITHOUT
#   Note: circadian motifs might not be found because there will be genes in 
#         the bg, which cycling in other tissues

# HOMER - 7:
#   Target: TISSUE-SPECIFIC CYCLING GENES
#   BGs: GENES EXPRESSED IN THE TISSUE, BUT NOT CYCLING IN ANY TISSUE
#   EXPECTATIONS: COMMON CIRCADIAN MOTIFS, tissue-specific circadian motifs 
#                 could be found if they are ubiquitously present
#   QUANTILES OF AMPLITUDE (TARGETS): Q1 - Q5
#   QUANTILES OF EXPRESSION (BGs): Q1 - Q5
#   50bp TSS: WITH AND WITHOUT
#   Note: that should be a good one

# HOMER - 8:
#   Target: TISSUE-SPECIFIC CYCLING GENES
#   BGs: GENES NOT CYCLING IN ANY TISSUE (they might not be expressed in this
#        tissue)
#   EXPECTATIONS: COMMON CIRCADIAN MOTIFS, tissue-specific circadian motifs 
#                 could be found if they are ubiquitously present
#   QUANTILES OF AMPLITUDE (TARGETS): Q1 - Q5
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT
#   Note: that should be a good one

# HOMER - 9:
#   Target: TISSUE-SPECIFIC CYCLING GENES
#   BGs: GENES CYCLING IN OTHER TISSUES AND EXPRESSED IN THIS ONE
#   EXPECTATIONS: COMMON CIRCADIAN MOTIFS, tissue-specific circadian motifs 
#   QUANTILES OF AMPLITUDE (TARGETS): Q1 - Q5
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT
#   Note: that should be a good one

# HOMER - 10:
#   Target: TISSUE-SPECIFIC CYCLING GENES
#   BGs: GENES CYCLING IN OTHER TISSUES (they might not be expressed in this
#        tissue)
#   EXPECTATIONS: COMMON CIRCADIAN MOTIFS, tissue-specific circadian motifs 
#   QUANTILES OF AMPLITUDE (TARGETS): Q1 - Q5
#   QUANTILES OF EXPRESSION (BGs): NO
#   50bp TSS: WITH AND WITHOUT
#   Note: that should be a good one

# OUTPUTS----------------------------------------------------------------------
# Except beds which would be printed, I would like also to have a table with
# all target - bgs combinations which I would like to run
targBgCombs <- data.table(Tissue = character(), TargFullName = character(),
                          TargName = character(), TargAmpQ = integer(), 
                          TargNumb = integer(), BgFullName = character(), 
                          BGname = character(), BGExprQ = integer(), 
                          BGnumb = integer(), upstrLen = integer(), 
                          downstrLen = integer(), tss50 = factor(), 
                          TargetBed = character(),  BGbed = character())

# TARGETs: TISSUE SPECIFIC EXPRESSION (TSE) -----------------------------------
targetsList <- list()
for (tissue in allTissues) {
  # produce target beds: get genes which expression only in this tissue
  targetGenes <- TSE[[tissue]]
  targetGenes <- tisSpecCycl[Gene %in% targetGenes]
  # add information about amplitude quantile
  ampColName <- grep(paste0('Expr.*', tissue), colnames(targetGenes), value = T)
  targetGenes <- targetGenes[, c('Gene', ampColName), with = F]
  # rename column for simplicity
  colnames(targetGenes)[2] <- 'Expr_Q'
  
  # select the amplitude and write bed
  for (ampQuant in 1:ampNumbQuant) {
    targetGenesAmp <- targetGenes[Expr_Q <= ampQuant]$Gene
    for (i in 1:length(upstrVect)) {
      writeBedForHomer(targetGenesAmp, upstrVect[i], downstrVect[i], 50, 50,
                       paste0(saveTSEbedDir, tissue, '_TSE_A', ampQuant, '_', 
                              upstrVect[i], '_50.bed'))
      writeBedForHomer(targetGenesAmp, upstrVect[i], downstrVect[i], 0, 0,
                       paste0(saveTSEbedDir, tissue, '_TSE_A', ampQuant, '_',
                              upstrVect[i], '_0.bed'))
    }
    targetsList[[length(targetsList) + 1]] <- targetGenesAmp
    names(targetsList)[length(targetsList)] <- paste(tissue, ampQuant)
  }
}

# BGs (1): RANDOM -------------------------------------------------------------
for (tissue in allTissues) {
  for (ampQuant in 1:ampNumbQuant) {
    for (i in 1:length(upstrVect)) {
      for (avoidLen in skipTSSLen) {
        targetBedPath <- paste0(saveTSEbedDir, tissue, '_TSE_A', ampQuant, '_', 
                                upstrVect[i], '_', avoidLen, '.bed')
        addToTargBgCombs <- data.table(Tissue = tissue, 
                                       TargFullName = 'TISSUE SPECIFIC EXPRESSED',
                                       TargName = 'TSE', TargAmpQ = ampQuant,
                                       TargNumb = length(targetsList[[paste(tissue, ampQuant)]]),
                                       BgFullName = 'RANDOM', BGname = 'RANDOM', 
                                       BGExprQ = exprNumbQuant, BGnumb = NA,
                                       upstrLen = upstrVect[i], 
                                       downstrLen = downstrVect[i],
                                       tss50 = avoidLen, 
                                       TargetBed = targetBedPath, BGbed = NA)
        targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
      }
    }
  }
}

# BGs (2): ALL EXPRESSED GENES ------------------------------------------------
for (tissue in allTissues) {
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    
    # get background genes
    BGgenesAmp <- setdiff(tisSpecCycl$Gene, 
                          targetsList[[paste(tissue, ampQuant)]])
    # no need to check that gene is expressed in at least one tissue - it is 
    # already the case in the tisSpecCycl
    for (i in 1:length(upstrVect)) {
      for (avoidLen in skipTSSLen) {
        targetBedPath <- paste0(saveTSEbedDir, tissue, '_TSE_A', ampQuant, '_',
                                upstrVect[i], '_', avoidLen, '.bed')
        bgBedPath <- paste0(saveTSEbedDir, tissue, '_AEG_A', ampQuant, '_', 
                            upstrVect[i], '_', avoidLen, '.bed')
        writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], avoidLen,
                         avoidLen, bgBedPath)
        
        # add to the data table with all possile combinations
        addToTargBgCombs <- data.table(Tissue = tissue, 
                                       TargFullName = 'TISSUE SPECIFIC EXPRESSED',
                                       TargName = 'TSE', TargAmpQ = ampQuant,
                                       TargNumb = length(targetGenesAmp), 
                                       BgFullName = 'ALL EXPRESSED GENES',
                                       BGname = 'AEG', BGExprQ = exprNumbQuant,
                                       BGnumb = length(BGgenesAmp), 
                                       upstrLen = upstrVect[i], 
                                       downstrLen = downstrVect[i],
                                       tss50 = avoidLen, 
                                       TargetBed = targetBedPath, 
                                       BGbed = bgBedPath)
        targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
      }
    }
  }
}

# BGs (3): GENES TISSUE-SPECIFICALLY EXPRESSED FOR OTHER TISSUES --------------
for (tissue in allTissues) {
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    
    # get background genes
    BGgenesAmp <- unlist(TSE[names(TSE) != tissue])
    for (i in 1:length(upstrVect)) {
      for (avoidLen in skipTSSLen) {
        targetBedPath <- paste0(saveTSEbedDir, tissue, '_TSE_A', ampQuant, '_',
                                upstrVect[i], '_', avoidLen, '.bed')
        bgBedPath <- paste0(saveTSEbedDir, tissue, '_TSEAT_A', ampQuant, '_', 
                            upstrVect[i], '_', avoidLen, '.bed')
        writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], avoidLen,
                         avoidLen, bgBedPath)
        
        # add to the data table with all possile combinations
        addToTargBgCombs <- data.table(Tissue = tissue, 
                                       TargFullName = 'TISSUE SPECIFIC EXPRESSED',
                                       TargName = 'TSE', TargAmpQ = ampQuant,
                                       TargNumb = length(targetGenesAmp), 
                                       BgFullName = 'TISSUE SPECIFIC EXPRESSED OTHER TISSUES',
                                       BGname = 'TSEAT', BGExprQ = exprNumbQuant,
                                       BGnumb = length(BGgenesAmp), 
                                       upstrLen = upstrVect[i], 
                                       downstrLen = downstrVect[i],
                                       tss50 = avoidLen, 
                                       TargetBed = targetBedPath, 
                                       BGbed = bgBedPath)
        targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
      }
    }
  }
}

# WRITE HOMER SCRIPT ----------------------------------------------------------
saveRDS(targBgCombs, 'Rds/targBgCombs_TSE.Rds')
# path to the reference genome
refGenPath <- '/scratch/el/monthly/mlitovch/RefGen/dm3/dm3.Wolb.fa'

# add name of the output folder
targBgCombs[, outputCode := apply(targBgCombs, 1, 
                                  function(x) paste(x['Tissue'], x['TargName'],
                                                    x['TargAmpQ'], x['BGname'],
                                                    x['BGExprQ'], 
                                                    x['upstrLen'],
                                                    x['tss50']))]
targBgCombs[, outputCode := gsub(' ', '_', outputCode)]
targBgCombs[, outputCode := gsub('__', '_', outputCode)]

# assemble commands
homerCommands <- paste('"findMotifsGenome.pl', targBgCombs$TargetBed,
                       refGenPath, targBgCombs$refGenPath, 
                       targBgCombs$outputCode, '-size given -len 10 -bg',
                       targBgCombs$BGbed, '"')
homerCommands <- gsub(' -bg NA', '', homerCommands)
write(homerCommands, 'homer_TSE.sh', ncolumns = 1)

# now just past what is below in the file

#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e motifs_TSE_len10.%I.err
#BSUB -o motifs_TSE_len10.%I.out
#BSUB -J motifs_TSE_len10[1-360]
#BSUB -M 24000000
#BSUB -R rusage[mem=24000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#export PATH=/software/bin:$PATH;
#module use /software/module/;
#module add UHTS/Analysis/homer/4.2;

# commands=( zero

# HERE ARE ALL THE COMMANDS

# toDo=${commands[${LSB_JOBINDEX}]};
# $toDo
# exit 0;

# OUTPUTS----------------------------------------------------------------------
# Except beds which would be printed, I would like also to have a table with
# all target - bgs combinations which I would like to run
targBgCombs <- data.table(Tissue = character(), TargFullName = character(),
                          TargName = character(), TargAmpQ = integer(), 
                          TargNumb = integer(), BgFullName = character(), 
                          BGname = character(), BGExprQ = integer(), 
                          BGnumb = integer(), upstrLen = integer(), 
                          downstrLen = integer(), tss50 = factor(), 
                          TargetBed = character(),  BGbed = character())

# TARGETs: TISSUE SPECIFIC CYCLING (TSC) --------------------------------------
targetsList <- list()
for (tissue in allTissues) {
  # produce target beds: get genes which cycling only in this tissue
  targetGenes <- tisSpecCycl[CyclesIn == tissue, ]
  # add information about amplitude quantile
  ampColName <- grep(paste0('Amp.*', tissue), colnames(targetGenes), value = T)
  targetGenes <- targetGenes[, c('Gene', ampColName), with = F]
  # rename column for simplicity
  colnames(targetGenes)[2] <- 'Amp_Q'
  
  # select the amplitude and write bed
  for (ampQuant in 1:ampNumbQuant) {
    targetGenesAmp <- targetGenes[Amp_Q <= ampQuant]$Gene
    for (i in 1:length(upstrVect)) {
      writeBedForHomer(targetGenesAmp, upstrVect[i], downstrVect[i], 50, 50,
                       paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                              upstrVect[i], '_50.bed'))
      writeBedForHomer(targetGenesAmp, upstrVect[i], downstrVect[i], 0, 0,
                       paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                              upstrVect[i], '_0.bed'))
    }
    targetsList[[length(targetsList) + 1]] <- targetGenesAmp
    names(targetsList)[length(targetsList)] <- paste(tissue, ampQuant)
  }
}

# BGs (3): RANDOM -------------------------------------------------------------
for (tissue in allTissues) {
  for (ampQuant in 1:ampNumbQuant) {
    for (i in 1:length(upstrVect)) {
      for (avoidLen in skipTSSLen) {
        targetBedPath <- paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                                upstrVect[i], '_', avoidLen, '.bed')
        addToTargBgCombs <- data.table(Tissue = tissue, 
                                       TargFullName = 'TISSUE SPECIFIC CYCLING',
                                       TargName = 'TSC', TargAmpQ = ampQuant,
                                       TargNumb = length(targetsList[[paste(tissue, ampQuant)]]),
                                       BgFullName = 'RANDOM', BGname = 'RANDOM', 
                                       BGExprQ = exprNumbQuant, BGnumb = NA,
                                       upstrLen = upstrVect[i], 
                                       downstrLen = downstrVect[i],
                                       tss50 = avoidLen, 
                                       TargetBed = targetBedPath, BGbed = NA)
        targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
      }
    }
  }
}

# BGs (4): ALL EXPRESSED GENES ------------------------------------------------
# BG = all expressed genes in all tissues - the target ones
# EXPECTATIONS: circadian motifs (i.e. E-box) maybe tissue-specific cycling
# motifs
for (tissue in allTissues) {
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    
    # get background genes
    BGgenesAmp <- setdiff(tisSpecCycl$Gene, 
                          targetsList[[paste(tissue, ampQuant)]])
    # no need to check that gene is expressed in at least one tissue - it is 
    # already the case in the tisSpecCycl
    for (i in 1:length(upstrVect)) {
      for (avoidLen in skipTSSLen) {
        targetBedPath <- paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                                upstrVect[i], '_', avoidLen, '.bed')
        bgBedPath <- paste0(saveTSCbedDir, tissue, '_AEG_A', ampQuant, '_', 
                            upstrVect[i], '_', avoidLen, '.bed')
        writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], avoidLen,
                         avoidLen, bgBedPath)
        
        # add to the data table with all possile combinations
        addToTargBgCombs <- data.table(Tissue = tissue, 
                                       TargFullName = 'TISSUE SPECIFIC CYCLING',
                                       TargName = 'TSC', TargAmpQ = ampQuant,
                                       TargNumb = length(targetGenesAmp), 
                                       BgFullName = 'ALL EXPRESSED GENES',
                                       BGname = 'AEG', BGExprQ = exprNumbQuant,
                                       BGnumb = length(BGgenesAmp), 
                                       upstrLen = upstrVect[i], 
                                       downstrLen = downstrVect[i],
                                       tss50 = avoidLen, 
                                       TargetBed = targetBedPath, 
                                       BGbed = bgBedPath)
        targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
      }
    }
  }
}

# BGs (5) : GENES EXPRESSED IN THE TISSUE, BUT NOT CYCLING---------------------
# BG = all genes expressed in this tissue - the target ones
# EXPECTATIONS: circadian motifs (i.e. E-box) maybe tissue-specific cycling
# motifs. However, there will be genes which are expressed in this tissue, but
# cycle only in the other tissues, which might remove circadian motifs. Here it
# actually doesn't make sence to do splits on quantiles of amplitude in targets
# because then rest will be in BG and dilute circadian motifs
for (tissue in allTissues) {
  # select genes expressed in the tissue of choise
  exprColName <- grep(paste0('Expr.*', tissue), colnames(tisSpecCycl), 
                      value = T)
  exprGenes <- tisSpecCycl[, c('Gene', exprColName), with = F]
  colnames(exprGenes)[2] <- 'Expr_Q'
  exprGenes <- exprGenes[complete.cases(exprGenes)]
  
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    
    for (exprQuant in 1:exprNumbQuant) {
      # get BG genes
      BGgenesAmp <- exprGenes[Expr_Q == exprQuant]$Gene
      BGgenesAmp <- setdiff(BGgenesAmp, targetGenesAmp)
      
      for (i in 1:length(upstrVect)) {
        for (avoidLen in skipTSSLen) {
          targetBedPath <- paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                                  upstrVect[i], '_', avoidLen, '.bed')
          bgBedPath <- paste0(saveTSCbedDir, tissue, '_ENCG_A', 
                              ampQuant, '_E', exprQuant, '_', upstrVect[i], 
                              '_', avoidLen, '.bed')
          writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], avoidLen,
                           avoidLen, bgBedPath)
          
          # add to the data table with all possile combinations
          addToTargBgCombs <- data.table(Tissue = tissue, 
                                         TargFullName = 'TISSUE SPECIFIC CYCLING',
                                         TargName = 'TSC', TargAmpQ = ampQuant,
                                         TargNumb = length(targetGenesAmp),
                                         BgFullName = 'EXPRESSED IN THE TISSUE, BUT NOT CYCLING',
                                         BGname = 'ENCG', BGExprQ = exprQuant,
                                         BGnumb = length(BGgenesAmp),
                                         upstrLen = upstrVect[i], 
                                         downstrLen = downstrVect[i],
                                         tss50 = avoidLen, 
                                         TargetBed = targetBedPath, 
                                         BGbed = bgBedPath)
          targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
        }
      }
    }
  }
}

# BGs (6) : GENES EXPRESSED IN THE TISSUE, BUT NOT CYCLING IN ANY TISSUE ------
# BG = genes expressed in this tissue - the target ones - cycling in other 
# tissues

# select genes which are not cycling in any tissue
notCycling <- tisSpecCycl[CyclesIn == '']
for (tissue in allTissues) {
  # select genes expressed in the tissue of choise
  exprColName <- grep(paste0('Expr.*', tissue), colnames(notCycling), 
                      value = T)
  exprGenes <- notCycling[, c('Gene', exprColName), with = F]
  colnames(exprGenes)[2] <- 'Expr_Q'
  exprGenes <- exprGenes[complete.cases(exprGenes)]
  
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    
    for (exprQuant in 1:exprNumbQuant) {
      # get BG genes
      BGgenesAmp <- exprGenes[Expr_Q == exprQuant]$Gene
      BGgenesAmp <- setdiff(BGgenesAmp, targetGenesAmp)
      
      for (i in 1:length(upstrVect)) {
        for (avoidLen in skipTSSLen) {
          targetBedPath <- paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                                  upstrVect[i], '_', avoidLen, '.bed')
          bgBedPath <- paste0(saveTSCbedDir, tissue, '_ENCATG_A', 
                              ampQuant, '_E', exprQuant, '_', upstrVect[i],
                              '_', avoidLen, '.bed')
          writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], 50, 50,
                           bgBedPath)
          
          # add to the data table with all possile combinations
          addToTargBgCombs <- data.table(Tissue = tissue, 
                                         TargFullName = 'TISSUE SPECIFIC CYCLING',
                                         TargName = 'TSC', TargAmpQ = ampQuant,
                                         TargNumb = length(targetGenesAmp),
                                         BgFullName = 'EXPRESSED IN THE TISSUE, BUT NOT CYCLING IN ANY TISSUE',
                                         BGname = 'ENCATG', BGExprQ = exprQuant,
                                         BGnumb = length(BGgenesAmp),
                                         upstrLen = upstrVect[i], 
                                         downstrLen = downstrVect[i],
                                         tss50 = avoidLen, 
                                         TargetBed = targetBedPath, 
                                         BGbed = bgBedPath)
          targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
        }
      }
    }
  }
}

# BGs (7) : GENES NOT CYCLING IN ANY TISSUE (might not be expressed) ----------
# select genes which are not cycling in any tissue
notCycling <- tisSpecCycl[CyclesIn == '']$Gene
for (tissue in allTissues) {
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    # get bg genes
    BGgenesAmp <- setdiff(notCycling, targetGenesAmp)
    
    for (i in 1:length(upstrVect)) {
      for (avoidLen in skipTSSLen) {
        targetBedPath <- paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                                upstrVect[i], '_', avoidLen, '.bed')
        bgBedPath <- paste0(saveTSCbedDir, tissue, '_NCATG_A',  ampQuant, '_E',
                            exprNumbQuant, '_', upstrVect[i], '_', avoidLen, 
                            '.bed')
        writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], 50, 50,
                         bgBedPath)
        
        # add to the data table with all possile combinations
        addToTargBgCombs <- data.table(Tissue = tissue, 
                                       TargFullName = 'TISSUE SPECIFIC CYCLING',
                                       TargName = 'TSC', TargAmpQ = ampQuant,
                                       TargNumb = length(targetGenesAmp),
                                       BgFullName = 'NOT CYCLING IN ANY TISSUE',
                                       BGname = 'NCATG', 
                                       BGExprQ = exprNumbQuant,
                                       BGnumb = length(BGgenesAmp),
                                       upstrLen = upstrVect[i], 
                                       downstrLen = downstrVect[i],
                                       tss50 = avoidLen, 
                                       TargetBed = targetBedPath, 
                                       BGbed = bgBedPath)
        targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
      }
    }
  }
}

# BGs(8) : GENES CYCLING IN OTHER TISSUES AND EXPRESSED IN THIS ONE -----------
for (tissue in allTissues) {
  # genes which cycle in other tissues
  cyclOtherTiss <- tisSpecCycl[!grepl(tissue, CyclesIn)]
  cyclOtherTiss <- cyclOtherTiss[CyclesIn != '']
  # select genes expressed in the tissue of choise
  exprColName <- grep(paste0('Expr.*', tissue), colnames(cyclOtherTiss), 
                      value = T)
  cyclOtherTissExpr <- cyclOtherTiss[, c('Gene', exprColName), with = F]
  colnames(cyclOtherTissExpr)[2] <- 'Expr_Q'
  cyclOtherTissExpr <- cyclOtherTissExpr[complete.cases(cyclOtherTissExpr)]
  
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    for (exprQuant in 1:exprNumbQuant) {
      # get BG genes
      BGgenesAmp <- cyclOtherTissExpr[Expr_Q == exprQuant]$Gene
      BGgenesAmp <- setdiff(BGgenesAmp, targetGenesAmp)
      
      for (i in 1:length(upstrVect)) {
        for (avoidLen in skipTSSLen) {
          targetBedPath <- paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                                  upstrVect[i], '_', avoidLen, '.bed')
          bgBedPath <- paste0(saveTSCbedDir, tissue, '_ECOTG_A',  ampQuant, '_E',
                              exprQuant, '_', upstrVect[i], '_', avoidLen, '.bed')
          writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], avoidLen,
                           avoidLen, bgBedPath)
          
          # add to the data table with all possile combinations
          addToTargBgCombs <- data.table(Tissue = tissue, 
                                         TargFullName = 'TISSUE SPECIFIC CYCLING',
                                         TargName = 'TSC', TargAmpQ = ampQuant,
                                         TargNumb = length(targetGenesAmp),
                                         BgFullName = 'CYCLING IN OTHER TISSUES AND EXPRESSED IN THIS',
                                         BGname = 'ECOTG', BGExprQ = exprQuant,
                                         BGnumb = length(BGgenesAmp),
                                         upstrLen = upstrVect[i], 
                                         downstrLen = downstrVect[i],
                                         tss50 = avoidLen, 
                                         TargetBed = targetBedPath, 
                                         BGbed = bgBedPath)
          targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
        }
      }
    }
  }
}

# BGs(9) : GENES CYCLING IN OTHER TISSUES (might not be expressed in this) ----
for (tissue in allTissues) {
  # genes which cycle in other tissues
  cyclOtherTiss <- tisSpecCycl[!grepl(tissue, CyclesIn)]
  cyclOtherTiss <- cyclOtherTiss[CyclesIn != '']$Gene
  
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    # get BG genes
    BGgenesAmp <- setdiff(cyclOtherTiss, targetGenesAmp)
    
    for (i in 1:length(upstrVect)) {
      for (avoidLen in skipTSSLen) {
        targetBedPath <- paste0(saveTSCbedDir, tissue, '_TSC_A', ampQuant, '_', 
                                upstrVect[i], '_', avoidLen, '.bed')
        bgBedPath <- paste0(saveTSCbedDir, tissue, '_COTG_A', ampQuant, '_E', 
                            exprNumbQuant, '_', upstrVect[i], '_', avoidLen,
                            '.bed')
        writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], avoidLen,
                        avoidLen, bgBedPath)
          
        # add to the data table with all possile combinations
        addToTargBgCombs <- data.table(Tissue = tissue, 
                                       TargFullName = 'TISSUE SPECIFIC CYCLING',
                                       TargName = 'TSC', TargAmpQ = ampQuant,
                                       TargNumb = length(targetGenesAmp),
                                       BgFullName = 'CYCLING IN OTHER TISSUES',
                                       BGname = 'COTG', 
                                       BGExprQ = exprNumbQuant,
                                       BGnumb = length(BGgenesAmp),
                                       upstrLen = upstrVect[i], 
                                       downstrLen = downstrVect[i],
                                       tss50 = avoidLen, 
                                       TargetBed = targetBedPath, 
                                       BGbed = bgBedPath)
        targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
      }
    }
  }
}

# WRITE HOMER SCRIPT ----------------------------------------------------------
saveRDS(targBgCombs, 'Rds/targBgCombs_TSC.Rds')
# path to the reference genome
refGenPath <- '/scratch/el/monthly/mlitovch/RefGen/dm3/dm3.Wolb.fa'

# add name of the output folder
targBgCombs[, outputCode := apply(targBgCombs, 1, 
                                  function(x) paste(x['Tissue'], x['TargName'],
                                                    x['TargAmpQ'], x['BGname'],
                                                    x['BGExprQ'], 
                                                    x['upstrLen'],
                                                    x['tss50']))]
targBgCombs[, outputCode := gsub(' ', '_', outputCode)]
targBgCombs[, outputCode := gsub('__', '_', outputCode)]

# assemble commands
homerCommands <- paste('"findMotifsGenome.pl', targBgCombs$TargetBed,
                       refGenPath, targBgCombs$refGenPath, 
                       targBgCombs$outputCode, '-size given -len 10 -bg',
                       targBgCombs$BGbed, '"')
homerCommands <- gsub(' -bg NA', '', homerCommands)
write(homerCommands, 'homer_TSC.sh', ncolumns = 1)

# now just past what is below in the file

#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e motifs_len10.%I.err
#BSUB -o motifs_len10.%I.out
#BSUB -J motifs_len10[1-2280]
#BSUB -M 24000000
#BSUB -R rusage[mem=24000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#export PATH=/software/bin:$PATH;
#module use /software/module/;
#module add UHTS/Analysis/homer/4.2;

# commands=( zero

# HERE ARE ALL THE COMMANDS

# toDo=${commands[${LSB_JOBINDEX}]};
# $toDo
# exit 0;

# TARGETs: GENES CYCLING IN TISSUE --------------------------------------------
targetsList <- list()
for (tissue in allTissues) {
  # produce target beds: get genes which cycling only in this tissue
  targetGenes <- tisSpecCycl[grepl(tissue, CyclesIn) ]
  # add information about amplitude quantile
  ampColName <- grep(paste0('Amp.*', tissue), colnames(targetGenes), value = T)
  targetGenes <- targetGenes[, c('Gene', ampColName), with = F]
  # rename column for simplicity
  colnames(targetGenes)[2] <- 'Amp_Q'
  
  # select the amplitude and write bed
  for (ampQuant in 1:ampNumbQuant) {
    targetGenesAmp <- targetGenes[Amp_Q <= ampQuant]$Gene
    for (i in 1:length(upstrVect)) {
      writeBedForHomer(targetGenesAmp, upstrVect[i], downstrVect[i], 50, 50,
                       paste0(saveCiTbedDir, tissue, '_CiT_A', ampQuant, '_', 
                              upstrVect[i], '_50.bed'))
      writeBedForHomer(targetGenesAmp, upstrVect[i], downstrVect[i], 0, 0,
                       paste0(saveCiTbedDir, tissue, '_CiT_A', ampQuant, '_', 
                              upstrVect[i], '_0.bed'))
    }
    targetsList[[length(targetsList) + 1]] <- targetGenesAmp
    names(targetsList)[length(targetsList)] <- paste(tissue, ampQuant)
  }
}

# BGs(1) : GENES CYCLING IN OTHER TISSUES AND EXPRESSED IN THIS ONE -----------
for (tissue in allTissues) {
  # genes which cycle in other tissues
  cyclOtherTiss <- tisSpecCycl[!grepl(tissue, CyclesIn)]
  cyclOtherTiss <- cyclOtherTiss[CyclesIn != '']
  # select genes expressed in the tissue of choise
  exprColName <- grep(paste0('Expr.*', tissue), colnames(cyclOtherTiss), 
                      value = T)
  cyclOtherTissExpr <- cyclOtherTiss[, c('Gene', exprColName), with = F]
  colnames(cyclOtherTissExpr)[2] <- 'Expr_Q'
  cyclOtherTissExpr <- cyclOtherTissExpr[complete.cases(cyclOtherTissExpr)]
  
  for (ampQuant in 1:ampNumbQuant) {
    # get target genes 
    targetGenesAmp <- targetsList[[paste(tissue, ampQuant)]]
    for (exprQuant in 1:exprNumbQuant) {
      # get BG genes
      BGgenesAmp <- cyclOtherTissExpr[Expr_Q == exprQuant]$Gene
      BGgenesAmp <- setdiff(BGgenesAmp, targetGenesAmp)
      
      for (i in 1:length(upstrVect)) {
        for (avoidLen in skipTSSLen) {
          targetBedPath <- paste0(saveCiTbedDir, tissue, '_CiT_A', ampQuant, '_', 
                                  upstrVect[i], '_', avoidLen, '.bed')
          bgBedPath <- paste0(saveCiTbedDir, tissue, '_ECOTG_A',  ampQuant, '_E',
                              exprQuant, '_', upstrVect[i], '_', avoidLen, '.bed')
          writeBedForHomer(BGgenesAmp, upstrVect[i], downstrVect[i], avoidLen,
                           avoidLen, bgBedPath)
          
          # add to the data table with all possile combinations
          addToTargBgCombs <- data.table(Tissue = tissue, 
                                         TargFullName = 'CYCLING IN TISSUE',
                                         TargName = 'CiT', TargAmpQ = ampQuant,
                                         TargNumb = length(targetGenesAmp),
                                         BgFullName = 'CYCLING IN OTHER TISSUES AND EXPRESSED IN THIS',
                                         BGname = 'ECOTG', BGExprQ = exprQuant,
                                         BGnumb = length(BGgenesAmp),
                                         upstrLen = upstrVect[i], 
                                         downstrLen = downstrVect[i],
                                         tss50 = avoidLen, 
                                         TargetBed = targetBedPath, 
                                         BGbed = bgBedPath)
          targBgCombs <- rbind(targBgCombs, addToTargBgCombs)
        }
      }
    }
  }
}

# WRITE HOMER SCRIPT ----------------------------------------------------------
saveRDS(targBgCombs, 'Rds/targBgCombs_CiT.Rds')
# path to the reference genome
refGenPath <- '/scratch/el/monthly/mlitovch/RefGen/dm3/dm3.Wolb.fa'

# add name of the output folder
targBgCombs[, outputCode := apply(targBgCombs, 1, 
                                  function(x) paste(x['Tissue'], x['TargName'],
                                                    x['TargAmpQ'], x['BGname'],
                                                    x['BGExprQ'], 
                                                    x['upstrLen'],
                                                    x['tss50']))]
targBgCombs[, outputCode := gsub(' ', '_', outputCode)]
targBgCombs[, outputCode := gsub('__', '_', outputCode)]

# assemble commands
homerCommands <- paste('"findMotifsGenome.pl', targBgCombs$TargetBed,
                       refGenPath, targBgCombs$refGenPath, 
                       targBgCombs$outputCode, '-size given -len 10 -bg',
                       targBgCombs$BGbed, '"')
homerCommands <- gsub(' -bg NA', '', homerCommands)
write(homerCommands, 'homer_CiT.sh', ncolumns = 1)

# now just past what is below in the file

#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e motifs_len10.%I.err
#BSUB -o motifs_len10.%I.out
#BSUB -J motifs_len10[1-2280]
#BSUB -M 24000000
#BSUB -R rusage[mem=24000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#export PATH=/software/bin:$PATH;
#module use /software/module/;
#module add UHTS/Analysis/homer/4.2;

# commands=( zero

# HERE ARE ALL THE COMMANDS

# toDo=${commands[${LSB_JOBINDEX}]};
# $toDo
# exit 0;


