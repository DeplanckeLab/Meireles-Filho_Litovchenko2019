source('~/Desktop/AroundTheClock_Aug_Oct2017/0_common_functions.R')

# INPUTS ----------------------------------------------------------------------
dgrpPath <- '~/Desktop/eQTL_input/dgrp2_fb57_variable_R.vcf'
RALtoBloomPath <- '~/Desktop/AroundTheClock_HiSeq_dm3/data_files/RAL_to_Bloomington.txt'
allCircGenesPath <- '~/Desktop/AroundTheClock_HiSeq_dm3/all_circadian_genes.txt'
locomotDefPath <- '~/Desktop/AroundTheClock_HiSeq_dm3/locomotor_rhythm_defective_avaible_lines.txt'
photorecDefPath <- '~/Desktop/AroundTheClock_HiSeq_dm3/photoperiod_response_defective_avaible_lines.txt'

# read all variable sites in dgrps
dgrp <- fread(dgrpPath)
dgrpHigh <- dgrp[grepl('HIGH', dgrp$INFO)]

# read RAL to bloomington convertion table
RALtoBloom <- fread(RALtoBloomPath, header = T, stringsAsFactors = F, 
                    sep = '\t')
RALtoBloom$line <- gsub('DGRP', 'line', RALtoBloom$line)
setkey(RALtoBloom, line)

# read all known circadian genes
allCircGenes <- read.table(allCircGenesPath, header = T, stringsAsFactors = F, 
                           sep = '\t')
allCircGenesCollapse <- paste0(allCircGenes$flybase_gene_id, collapse = '|')

# locomoter rhythm and photoreceptor defective
locomotDef <- fread(locomotDefPath, header = T, sep = '\t')
photorecDef <- fread(photorecDefPath, header = T, sep = '\t')

# QUICK ANALYSIS OF THE HIGH-IMPACT STUFF -------------------------------------
# all disruptive variants in this line, disruptive variants unique to this line
# disruptive variants unique to this line in circadian genes
linesOfInterest <- paste0('line_', c(88, 91, 158, 176, 340, 796, 832))

loisHighVars <- list()
loisHighUniqVars <- list()
loisHighUniqCircVars <- list()

for (LOI in linesOfInterest) {
  # all disruptive variants in this line
  LOI_high <- dgrpHigh[unlist(dgrpHigh[, which(colnames(dgrpHigh) == LOI), 
                                       with = F]) == 1]
  LOI_high_Annot <- annoToTable(LOI_high, LOI, RALtoBloom)
  affectLoco <- sapply(LOI_high_Annot$geneName, 
                       function(x) grepl(x, locomotDef$downSymbolup))
  affectPhoto <- sapply(LOI_high_Annot$geneName, 
                        function(x) grepl(x, photorecDef$downSymbolup))
  if (sum(affectLoco) != 0) {
    message(LOI_high_Annot$geneName[affectLoco])
  }
  if (sum(affectLoco) != 0) {
    message(LOI_high_Annot$geneName[affectPhoto])
  }
  
  loisHighVars[[length(loisHighVars) + 1]] <- LOI_high_Annot
  
  
  #  disruptive variants unique to this line
  LOI_high_uniq <- LOI_high[rowSums(LOI_high[, 10:ncol(LOI_high), with = F], 
                                    na.rm = T) == 1, ]
  LOI_high_uniqAnnot <- annoToTable(LOI_high_uniq, LOI, RALtoBloom)
  loisHighUniqVars[[length(loisHighUniqVars) + 1]] <- LOI_high_uniqAnnot
  
  # disruptive variants unique to this line in circadian genes
  LOI_high_uniq_circad <- LOI_high_uniq[grepl(allCircGenesCollapse, 
                                              LOI_high_uniq$INFO)]
  if (nrow(LOI_high_uniq_circad) != 0 ) {
    LOI_high_uniq_circadAnnot <- annoToTable(LOI_high_uniq_circad, LOI,
                                             RALtoBloom)
    loisHighUniqCircVars[[length(loisHighUniqCircVars) + 1]] <- LOI_high_uniq_circadAnnot
  } else {
    loisHighUniqCircVars[[length(loisHighUniqCircVars) + 1]] <- NA
  }
  
  #write.table(LOI_high_Annot, paste0('~/Desktop/', LOI, '_high.csv'), sep = '\t',
  #            col.names = T, row.names = F, quote = F)
  #write.table(LOI_high_uniqAnnot, paste0('~/Desktop/', LOI, '_high_uniq.csv'), 
  #            sep = '\t', col.names = T, row.names = F, quote = F)
  #write.table(LOI_high_uniq_circadAnnot, 
  #            paste0('~/Desktop/', LOI, '_high_uniq_circad.csv'),
  #            sep = '\t', col.names = T, row.names = F, quote = F)
}
names(loisHighVars) <- linesOfInterest
names(loisHighUniqVars) <- linesOfInterest
names(loisHighUniqCircVars) <- linesOfInterest

# UNIQUE COMBINATION OF DISRUPTIVE alleles in CIRCADIAN genes -----------------
# it's useless to look not in circadian genes, because then LOI_high is the 
# unique combination, because there are variants specific to this LOI
LOI_high_circad <- LOI_high[grepl(paste0(allCircGenes$flybase_gene_id, 
                                         collapse = '|'), 
                                  LOI_high$INFO)]
LOI_high_circadAnnot <- annoToTable(LOI_high_circad, LOI, RALtoBloom)
LOI_high_circadAnnot <- LOI_high_circadAnnot[grepl('HIGH', effect)]
#write.table(LOI_high_circadAnnot, 
#            paste0('~/Desktop/', LOI, '_high_circad.csv'),
#            sep = '\t', col.names = T, row.names = F, quote = F)

# which lines have as well this combination
alsoAffected <- RALtoBloom[colnames(LOI_high_circad)[9 + 
                                                       which(colSums(LOI_high_circad[, -1:-9,
                                                                                     with = F],
                                                                     na.rm = T) == 
                                                               ncol(LOI_high_circad))]]
print(alsoAffected)
for (varid in LOI_high_circad$ID) {
  message(paste(varid, sum(dgrp[ID == varid][, -1:-9, with = F], na.rm = T)))
}

# variants UNIQUE to this line, in CIRCADIAN genes ----------------------------
LOI_uniq <- dgrp[unlist(dgrp[, which(colnames(dgrp) == LOI), with = F]) == 1]
LOI_uniq <- LOI_uniq[rowSums(LOI_uniq[, 10:ncol(LOI_uniq), with = F], 
                             na.rm = T) == 1, ]
LOI_uniq_circad <- LOI_uniq[grepl(paste0(allCircGenes$flybase_gene_id, 
                                         collapse = '|'), 
                                  LOI_uniq$INFO)]
LOI_uniq_circadAnnot <- annoToTable(LOI_uniq_circad, LOI, RALtoBloom)
LOI_uniq_circadAnnot <- LOI_uniq_circadAnnot[order(effect)]
write.table(LOI_uniq_circadAnnot, 
            paste0('~/Desktop/', LOI, '_uniq_circad.csv'),
            sep = '\t', col.names = T, row.names = F, quote = F)

for (circgene in unique(LOI_uniq_circadAnnot$geneName)) {
  affectLoco <- locomotDef[grepl(circgene, locomotDef$downSymbolup), ]
  if (nrow(affectLoco) != 0) {
    message(circgene)
    #print(affectLoco)
  }
}

for (circgene in unique(LOI_uniq_circadAnnot$geneName)) {
  affectPhoto <- photorecDef[grepl(circgene, photorecDef$downSymbolup), ]
  if (nrow(affectPhoto) != 0) {
    message(circgene)
    #print(affectPhoto)
  }
}

# ALL variants in this line, in CIRCADIAN genes -------------------------------
LOI_circad <- dgrp[unlist(dgrp[, which(colnames(dgrp) == LOI), with = F]) == 1]
LOI_circad <- LOI_circad[grepl(paste0(allCircGenes$flybase_gene_id, 
                                      collapse = '|'), 
                               LOI_circad$INFO)]
LOI_circadAnnot <- annoToTable(LOI_circad, LOI, RALtoBloom)
#write.table(LOI_circadAnnot, 
#            paste0('~/Desktop/', LOI, '_circad.csv'),
#            sep = '\t', col.names = T, row.names = F, quote = F)