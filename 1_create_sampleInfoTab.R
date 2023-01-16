# FILE DESCRIPTION ------------------------------------------------------------
# : 1_create_sampleInfoTab.R
#
# DESCRIPTION : 
# This code puts all the information about sequenced samples, such as
# initial_name, rightGT_name, initial_GT, Lib, Tissue, Time_code, right_GT,
# Time, NR_total, NR_mapped, Month into one table
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
#  REVISION: 22.01.2018

# GENERAL NOTE: samples were taken 3 times: in January, August and October.
# JANUARY:
#     1. samples from DGRPs (one sample every 9 minutes, 150 lines in total,
#                            3 tissues)
#     2. one replicate of w- profiled every 2 h for the 1st day and every 4h
#        for the second day. First 12h light, all the rest - dark
#     3. lines 25205, 25174, 28211, 29655 were profiled for 24h every 2h
# AUGUST: 
#     w- was sampled as in January, but every 2h in duplicate and sequenced
#     spoiler alert: only gut and MT worked
# OCTOBER:
#     w- Brain and FB were resequenced

# NOTE: 1. samples from January, were realigned in August, because I lost them
#       due to Vital-IT
#       2. samples from August(w-) and January, ONLY W-, were realigned again 
#          in October

# SO BASICALLY, ALL SAMPLES FOR W- ARE IN OCTOBER, AND ALL SAMPLES FOR DGRPS 
# ARE IN AUGUST

# LIBs and FUNCTIONS ----------------------------------------------------------
library(data.table)
library(lubridate)
library(stringr)

#' GetNumericTime
#' Converts ZT/CT times into numberic time
#' @param timeAndDGRP vector of ZT/CT values, 
#'                    i.e. c("ZT5_w-", CT17_w-", "CT9_w-", "BRB2_FB_ZT1_25205")
#' @return numeric time corresponding to time points c(5, 41, 33)
#' @details DO NOT HANDLE MESSED UP SAMPLES! NO BRB7_Brain_ZT5_w-A please
#' @keywords sort, ZT, CT, time
#' @author Maria Litovchenko
GetNumericTime <- function(timeAndDGRP) {
  timePoints <- str_extract(timeAndDGRP, "[CZ]T\\d*")
  CTs <- grepl('CT', timePoints)
  # this is time points, I can't sort time points with ZTs/CTs, because 
  # they will be sorted lexicographically
  timePointsNumeric <- as.numeric(str_extract(timePoints, "[0-9]+"))
  timePointsNumeric[CTs == T] <- 24 + timePointsNumeric[CTs == T]
  timePointsNumeric
}

# INPUTS ----------------------------------------------------------------------
RdsDir <- 'Rds/'

# AUGUST
# read in samples name to real genotype map
initGTtoRightGTpath <- 'inputs_for_Rscripts/initialGT_to_rightGT_map.txt'
timeOfSampInitGTpath <- 'inputs_for_Rscripts/timeOfSampling_InitGT.txt'
mappingStatsAUGPath <- 'inputs_for_Rscripts/align_to_dm3_stats_AUG.txt'
initGTtoRightGT <- fread(initGTtoRightGTpath, stringsAsFactors = F)
names(initGTtoRightGT)[1] <- 'initial_name'
timeOfSampInitGT <- fread(timeOfSampInitGTpath, stringsAsFactors = F)
names(timeOfSampInitGT)[1] <- 'initial_GT'
timeOfSampInitGT$initial_GT <- as.character(timeOfSampInitGT$initial_GT)
mappingStatsAUG <- fread(mappingStatsAUGPath, stringsAsFactors = F)
names(mappingStatsAUG)[1] <- 'initial_name'

# OCTOBER
# collection codes - time of samples and replicates
collectCodesPath <- 'inputs_for_Rscripts/white_collection_codes_OCT.txt'
mappingStatsOCTPath <- 'inputs_for_Rscripts/align_to_dm3_stats_OCT.txt'

collectCodes <- fread(collectCodesPath, stringsAsFactors = F)
mappingStatsOCT <- fread(mappingStatsOCTPath, stringsAsFactors = F)
names(mappingStatsOCT)[1] <- 'initial_name'

# initiate future big table with all the info about samples
INFO <- initGTtoRightGT
w_INFO <- data.table(initial_name = mappingStatsOCT$initial_name,
                     rightGT_name = mappingStatsOCT$initial_name,
                     initial_GT = 'w-', right_GT = 'w-')
  
# 1:AUGUST (DGRPs) LIB, TISSUE, TIME CODE AND RIGHT GENOTYPE ------------------
# GET INFO ABOUT LIBRARY, TISSUE, TIME CODE AND RIGHT GENOTYPE

# extract initial GT
INFO$initial_GT <- sapply(INFO$initial_name, 
                          function(x) strsplit(x, '_')[[1]][4])
# add right genotype
INFO$right_GT <- sapply(INFO$rightGT_name,  
                        function(x) strsplit(x, '_')[[1]][4])

# add brb-library, it's the same for initial and right genotype
INFO$Lib <- sapply(INFO$rightGT_name, 
                   function(x) strsplit(x, '_')[[1]][1])
# add tissue, it's the same for initial and right genotype
INFO$Tissue <- toupper(sapply(INFO$rightGT_name,
                              function(x) strsplit(x, '_')[[1]][2]))
# add time point code, it's the same for initial and right genotype
INFO$Time_code <- sapply(INFO$rightGT_name, 
                         function(x) strsplit(x, '_')[[1]][3])

# remove all w-, they will come from October
INFO <- INFO[!initial_GT == 'w-']

# 2:OCTOBER (w-) LIB, TISSUE, TIME CODE AND RIGHT GENOTYPE --------------------
# GET INFO ABOUT LIBRARY, TISSUE, TIME CODE AND RIGHT GENOTYPE
w_INFO$Lib <- paste0('BRB', sapply(w_INFO$initial_name,  
                                   function(x) strsplit(x, '_')[[1]][3]))
w_INFO$Tissue <- gsub('\\d', '', sapply(w_INFO$initial_name,
                                        function(x) strsplit(x, '_')[[1]][1]))
# get corresponding time code in usual (ZT/CT) format
octoberCodes <- as.integer(gsub('\\D', '', 
                                sapply(w_INFO$initial_name, 
                                       function(x) strsplit(x, '_')[[1]][1])))
w_INFO$Time_code<- collectCodes[octoberCodes, TimeCode]


# 3:AUGUST ADD TIME INFO TO DGRPs ---------------------------------------------
# NOTE: I think, it's unlikely that the swap occured at the level of flies in 
# the tube. It's rather likely to occur on the level of test tubes, so I'll 
# assign time according to right genotype

setkey(timeOfSampInitGT, initial_GT)
INFO$Time <- timeOfSampInitGT[INFO$right_GT, CollectionTime]

# right down time for time-profiled lines
profiledTime <- GetNumericTime(INFO[grepl("[ZC]T", initial_name)]$initial_name)
# for future convertion to time class
profiledTime <- paste0(profiledTime, ':00:00')

INFO[grepl("[ZC]T", initial_name)]$Time <- profiledTime
# convert to time class. CT will be NAs, I'll deal with them later
timeClass <- strptime(INFO$Time, "%H:%M:%S")
numTime <- hour(timeClass) + minute(timeClass)/60
# deal with CTs
numTime[is.na(numTime)] <- as.numeric(gsub(':00:00', '',
                                           INFO[grepl("CT",
                                                      initial_name)]$Time))
INFO$Time <- numTime

# 4:OCTOBER ADD TIME INFO TO WHITE --------------------------------------------
w_INFO$Time <- paste0(GetNumericTime(w_INFO$Time_code), ':00:00')
# convert to time class. CT will be NAs, I'll deal with them later
w_timeClass <- strptime(w_INFO$Time, "%H:%M:%S")
w_numTime <- hour(w_timeClass) + minute(w_timeClass)/60
# deal with CTs
w_numTime[is.na(w_numTime)] <- as.numeric(gsub(':00:00', '', 
                                               w_INFO[grepl("CT",
                                                            Time_code)]$Time))
w_INFO$Time <- w_numTime

# 5: AUGUST ADD MAPPING INFO --------------------------------------------------
# recalculate percentages to the real number of reads
mappingStatsAUG[, Perc_unmapped_mismatch := 0.01 * Perc_unmapped_mismatch * 
                  NR_total]
mappingStatsAUG[, Perc_unmapped_short := 0.01 * Perc_unmapped_short * NR_total]
mappingStatsAUG[, Perc_unmapped_other := 0.01 * Perc_unmapped_other * NR_total]
setnames(mappingStatsAUG, grep('Perc', names(mappingStatsAUG), value = T),
         gsub('Perc', 'NR', grep('Perc', names(mappingStatsAUG), value = T)))

# Mapping info will be added according to initial name
setkey(INFO, initial_name)
setkey(mappingStatsAUG, initial_name)
mappingStatsAUG$Run_Numb <- NULL
mappingStatsAUG$BRBLib <- NULL
INFO <- merge(INFO, mappingStatsAUG)
INFO$NR_dedupl <- NA

# 6: OCTOBER ADD MAPPING INFO --------------------------------------------------
# recalculate percentages to the real number of reads
mappingStatsOCT[, Perc_unmapped_mismatch := 0.01 * Perc_unmapped_mismatch * 
                  NR_total]
mappingStatsOCT[, Perc_unmapped_short := 0.01 * Perc_unmapped_short * NR_total]
mappingStatsOCT[, Perc_unmapped_other := 0.01 * Perc_unmapped_other * NR_total]
setnames(mappingStatsOCT, grep('Perc', names(mappingStatsOCT), value = T),
         gsub('Perc', 'NR', grep('Perc', names(mappingStatsOCT), value = T)))

# Mapping info will be added according to initial name
setkey(w_INFO, initial_name)
setkey(mappingStatsOCT, initial_name)
mappingStatsOCT$Run_Numb <- NULL
mappingStatsOCT$BRBlib <- NULL
w_INFO <- merge(w_INFO, mappingStatsOCT)

# 7: ADD MONTH ----------------------------------------------------------------
INFO$Month = 'Jan'
w_INFO$Month <- sapply(w_INFO$initial_name, 
                       function(x) strsplit(x, '_')[[1]][2])

# 8: SAVE ---------------------------------------------------------------------
INFO <- rbind(INFO, w_INFO)
saveRDS(INFO, paste0(RdsDir, 'sampleInfo.Rds'))