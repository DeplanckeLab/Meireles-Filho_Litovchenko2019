#!/usr/bin/env Rscript

#  FILE: genotyping_part2.R
#
#  USAGE: change paths to files, run the script
#
#  DESCRIPTION: genotypes in-house sequnced lines to reference DGRP lines,
#		part 2, producing best hits
#
#  OPTIONS:  none
#  REQUIREMENTS: data.table
#  BUGS: --
#  NOTES:  ---
#  AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
#  COMPANY:  EPFL, Lausanne, Switzerland
#  VERSION:  2
#  CREATED:  15.08.2016
#  REVISION: 15.08.2016

library(data.table)

# Functions -------------------------------------------------------------------
#' bloomToRalAndBack
#' translates bloomington ids to rals(dgrps) and back
#' @param input vector with ids to convert
#' @param transTabPath translate table path, if not standard
#' @return vector of translated ids
bloomToRalAndBack <- function(input, replaceDGRPbyLine = F,
                              transTabPath = 'RefGen/RAL_to_Bloomington.txt') {
  translateTable <- read.table(transTabPath, header = T,  stringsAsFactors = F)
  input <- sapply(input, function(x) gsub('line', 'DGRP', x))
  result <- c()
  for (id in input) {
    coordsInTransTab <- as.vector(which(translateTable == id, arr.ind = T))
    if (grepl('DGRP', id)) {
      result <- c(result, translateTable[coordsInTransTab[1], 1])
    } else {
      result <- c(result, translateTable[coordsInTransTab[1], 2])
    }
  }
  if (replaceDGRPbyLine) {
    result <- sapply(result, function(x) gsub('DGRP', 'line', x))
  }
  result
}

#' annotateDGRPmatch
#' Selects 3 best mathes for the each of the sequenced lines based on 
#' similarity matrix
#' @param similarities similarity matrix
#' @return vector containing IDs of 3 DGRP lines which are the closest
#' gentically
annotateDGRPmatch <- function(similarities) {
  initialName <-  rownames(similarities)
  deplLineName <- strsplit(rownames(similarities), '_')[[1]][4]
  
  similarities <- as.vector(similarities[1, ])
  if (!grepl('^[0-9]', deplLineName)) {
    result <- c(initialName, rep('NA', 3), rep(0, 3), 'NA')
  } else {
    similaritiesSort <- sort(similarities, decreasing = T)
    
    result <- similaritiesSort[1:3]
    result <- c(names(result), result)
    bestDGRPmatch <- bloomToRalAndBack(gsub('line', 'DGRP', result[1]))
    qual <- ''
    
    if (is.na(deplLineName) | is.na(bestDGRPmatch)) {
      qual <- 'NOT IN DB'
    } else {
      if (as.numeric(gsub("\\D", "", deplLineName)) == as.numeric(bestDGRPmatch)) {
        if (as.numeric(result[4]) - as.numeric(result[5]) > 0.1) {
          if (as.numeric(result[4]) >= 0.90) {
            qual <- 'MATCH'
          } 
          if (as.numeric(result[4]) >= 0.85 & as.numeric(result[4]) < 0.90) {
            qual <- 'WEAK MATCH'
          } 
          if (as.numeric(result[4]) < 0.85) {
            qual <- 'VERY WEAK MATCH' 
          } 
        } else {
          qual <- 'I DONT KNOW'
        } 
      } else {
        if (as.numeric(result[4]) - as.numeric(result[5]) > 0.1) {
          if (as.numeric(result[4]) >= 0.90) {
            qual <- 'SWAP'
          }
          if (as.numeric(result[4]) >= 0.85 & as.numeric(result[4]) < 0.90) {
            qual <- 'WEAK SWAP'
          }
          if (as.numeric(result[4]) < 0.85) {
            qual <- 'VERY WEAK SWAP'
          }
        } else {
          qual <- 'I DONT KNOW'
        }
      }
    }
    resultBloom <- bloomToRalAndBack(result)
    resultBloom[is.na(resultBloom)] <- result[is.na(resultBloom)]
    result <- as.vector(unlist(c(initialName, resultBloom, qual)))
  }
  names(result) <- c()
  result
}

# Inputs ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
refSNPsPath <- args[1]
deplSNPsPath <- args[2]

refSNPs <- unlist(fread(refSNPsPath, nrows = 1, header = F))
deplSNPs <- unlist(fread(deplSNPsPath, nrows = 1, header = F))

# Create similarity matrix out of individual similarity vectors ---------------
similarityMatrix <- c()
for (i in 3:ncol(deplSNPs)) {
   similarityVector <- readRDS(paste0('8_dm3_similarityVectors/similarityVector_', 
                                      i, '.Rds'))
   similarityMatrix <- rbind(similarityMatrix, similarityVector)
} 

rownames(similarityMatrix) <- c()
similarityMatrix <- as.data.frame(similarityMatrix)
colnames(similarityMatrix) <- refSNPs[-2:-1] 
rownames(similarityMatrix) <- deplSNPs[-2:-1]

# Compute best hits -----------------------------------------------------------
bestHits <- c()
for (i in 1:nrow(similarityMatrix)) {
  bestHits <- rbind(bestHits, annotateDGRPmatch(similarityMatrix[i, ]))
}
colnames(bestHits) <- c('DEPL_line', 'DGRP_1stMatch', 'DGRP_2ndMatch', 
                        'DGRP_3rdMatch', 'DGRP_1stMatch_perc', 
                        'DGRP_2ndMatch_perc', 'DGRP_3rdMatch_perc', 'Decision')

# Output ----------------------------------------------------------------------
write.table(similarityMatrix, 'AroundTheClock_similarityMatrix_DGRP2_dm3.csv',
            sep = '\t', quote = F)
write.table(bestHits, 'AroundTheClock_genotyping_DGRP2_dm3.csv',  sep = '\t',
            quote = F, row.names = F)
saveRDS(bestHits, 'AroundTheClock_genotyping_DGRP2_dm3.Rds')