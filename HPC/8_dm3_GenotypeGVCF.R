#!/usr/bin/env Rscript

#  FILE: genotyping.R
#
#  USAGE: change paths to files, run the script
#
#  DESCRIPTION: genotypes in-house sequnced lines to reference DGRP lines
#
#  OPTIONS:  none
#  REQUIREMENTS:
#  BUGS: --
#  NOTES:  ---
#  AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
#  COMPANY:  EPFL, Lausanne, Switzerland
#  VERSION:  2
#  CREATED:  09.08.2016
#  REVISION: 09.08.2016

library(data.table)

# Inputs ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = T)

refSNPsPath <- args[1]
deplSNPsPath <- args[2]
indexInDepl <- as.integer(args[3])

refSNPs <- fread(refSNPsPath, header = T)
deplSNPs <- fread(deplSNPsPath, header = T)

message('Started processing', colnames(deplSNPs)[indexInDepl], 'at ',
        Sys.time())

# Functions -------------------------------------------------------------------
#' annotateDGRPmatch
#' Selects 3 best mathes for the each of the sequenced lines based on 
#' similarity matrix
#' @param similarities similarity matrix
#' @return string translated similarity to one of the terms: MATCH, WEAK MATCH
#' SWAP or I DONT KNOW
annotateDGRPmatch <- function(similarities) {
  deplLineName <- rownames(similarities)
  similarities <- as.vector(similarities[1, ])
  similaritiesSort <- sort(similarities, decreasing = T)
  
  result <- similaritiesSort[1:3]
  result <- c(names(result), result)
  bestDGRPmatch <- gsub('line', 'DGRP', result[1])
  qual <- ''
  if (deplLineName == bestDGRPmatch) {
    if (as.numeric(result[4]) >= 0.99) {
      qual <- 'MATCH'
    }
    if (as.numeric(result[4]) > 0.91 & as.numeric(result[4]) < 0.99) {
      qual <- 'WEAK MATCH'
    }
    if (as.numeric(result[4]) < 0.91) {
      qual <- 'I DONT KNOW'
    }
  } else {
    if (as.numeric(result[4]) > 0.99) {
      qual <- 'SWAP'
    }
    if (as.numeric(result[4]) > 0.91 & as.numeric(result[4]) < 0.99) {
      qual <- 'WEAK SWAP'
    }
    if (as.numeric(result[4]) < 0.91) {
      qual <- 'I DONT KNOW'
    }
  }
  result <- as.vector(unlist(c(deplLineName, result, qual)))
  names(result) <- c()
  result
}

#' calcSimilarityToDGRPs
#' Calculates similarity of one in-house sequenced line to all DGRP lines
#' @param deplanckeOneSeqLine genotype of one in-house sequenced line
#' @param dgrpOneSeqLine genotype of all DGRP lines
#' @return percentage of variants which overlap between two lines
calcSimilarityToDGRP <- function(deplanckeOneSeqLine, dgrpOneSeqLine) {
  names(dgrpOneSeqLine) <- c('chr', 'start', 'GT')
  dgrpOneSeqLine <- dgrpOneSeqLine[GT %in% c('0', '2'), ]
  
  setkey(deplanckeOneSeqLine, chr, start)
  setkey(dgrpOneSeqLine, chr, start)
  
  deplanckeOneSeqLine <- deplanckeOneSeqLine[dgrpOneSeqLine]
  dgrpOneSeqLine <- NULL
  deplanckeOneSeqLine <- deplanckeOneSeqLine[complete.cases(deplanckeOneSeqLine)]
  
  nrow(deplanckeOneSeqLine[GT == i.GT]) / nrow(deplanckeOneSeqLine)
}

# Perform processing ----------------------------------------------------------
deplOneLine <- deplSNPs[, c(1, 2, indexInDepl), with = F]
names(deplOneLine) <- c('chr', 'start', 'GT')
deplOneLine <- deplOneLine[GT %in% c('0', '2'), ]
  
similarityVector <- sapply(3:ncol(refSNPs),
                           function(x) calcSimilarityToDGRP(deplOneLine,
                                                            refSNPs[, c(1,2,x),
                                                                    with = F]))
  
print('Finished', colnames(deplSNPs)[indexInDepl], 'at ', Sys.time())
 
# Save result and clean up the memory -----------------------------------------
delpOneLine <- NULL
refSNPs <- NULL
deplSNPs <- NULL

saveRDS(similarityVector, 
        paste0('8_dm3_similarityVectors/similarityVector_', 
               indexInDepl, '.Rds'))