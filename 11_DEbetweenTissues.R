# FILE DESCRIPTION 11_DEbetweenTissues ----------------------------------------
#
# DESCRIPTION : calculates differential expression between tissues
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
# CREATED:  17.03.2018
# REVISION: 17.03.2018

setwd('~/Desktop/BitBucket/AroundTheClock_Aug_Oct2017/')
source('0_common_functions.R')

# Inputs ----------------------------------------------------------------------
saveRdsTo <- 'Rds/'

# counts NOT normalized because DESeq2 wants integer counts 
# genes not expressed in any tissue are removed
whiteDE <- readRDS(paste0(saveRdsTo, 'whiteDE.Rds'))

# collapse list into one table
whiteDE_DF <- do.call(cbind, whiteDE)
colnames(whiteDE_DF) <- gsub('.*[.]', '', colnames(whiteDE_DF))

# info about samples
INFO <- as.data.table(readRDS(paste0(saveRdsTo, 'sampleInfo.Rds')))
# leave in only samples in out table
setkey(INFO, rightGT_name)
INFO <- INFO[colnames(whiteDE_DF)]

# list all tissues and corresponding colors
allTissues <- unique(INFO$Tissue)
tissueColor <- c('#688B8A', '#A0B084', '#FAEFD4', '#A57C65')
names(tissueColor) <- allTissues

# tissue-specific cycling genes
tisSpecCycl <- readRDS(paste0(saveRdsTo, 'tisSpecCycl.Rds'))
# only the ones cycling in one tissue
tisSpecCycl <- tisSpecCycl[CyclesIn %in% allTissues][, .(Gene, CyclesIn)]

# PCA to determine if there are any outliers ----------------------------------
# quantile normalization with voom
whiteDEvoom <- voom(counts = whiteDE_DF, normalize.method = "quantile")
whiteDEvoom <- whiteDEvoom$E

whiteDEvoomPCA <- prcomp(t(whiteDEvoom))
# calculate how much of variation (in percentage term) every PC explains
percPCexpl <- round(whiteDEvoomPCA$sdev / sum(whiteDEvoomPCA$sdev) * 100, 2)
names(percPCexpl) <- paste0('PC', 1:length(percPCexpl))

# Check with elbow plot to see how many PCs to include into consideration
plot(percPCexpl, xlab = 'PC index', ylab = '% of the variation explained', 
     main = 'Elbow plot for the PCA', pch = 20, type = 'o', bty = 'n',
     ylim = c(0, 1.1 * max(percPCexpl)))
# it's ugly, but it does the job. You need to look where the line has an 
# "elbow" - it's your altimate number of PCs.

# create matrix for plotting
whiteDEvoomPCA <- as.data.frame(whiteDEvoomPCA$x)
# add to this matrix some factors which may explain separation, such as tissue,
# depth of coverage, etc. INFO and whiteDEvoomPCA are sorted in the same way
identical(INFO$rightGT_name, rownames(whiteDEvoomPCA))
pcaToPlot <- cbind(whiteDEvoomPCA, INFO)

ggplot(pcaToPlot, aes(x = PC1, y = PC2, color = Tissue)) + 
  geom_point() + #geom_text() + 
  xlab(paste0('PC1 (', percPCexpl['PC1'], '%)')) + 
  ylab(paste0('PC2 (', percPCexpl['PC2'], '%)')) + 
  ggtitle("PCA based on voom normalized data") + theme_classic()

# one brain sample is an outlier and cluster with gut, remove it from the 
# consideration
ggplot(pcaToPlot[rownames(pcaToPlot) != 'BRAIN1_Oct_3', ],
       aes(x = PC1, y = PC2, color = Tissue)) + 
  geom_point() + #geom_text() + 
  xlab(paste0('PC1 (', percPCexpl['PC1'], '%)')) + 
  ylab(paste0('PC2 (', percPCexpl['PC2'], '%)')) + 
  ggtitle("PCA based on voom normalized data") + theme_classic()

whiteDE_DF <- whiteDE_DF[, colnames(whiteDE_DF) != 'BRAIN1_Oct_3']
setkey(INFO, rightGT_name)
INFO <- INFO[colnames(whiteDE_DF)]

# TISSUE-SPECIFIC EXPRESSED GENES as 80% --------------------------------------
# I count a gene to be tissue-specifically expressed if it's not expressed in
# 3 other tissues at all.
# Gene is expressed in a tissue if it has (> 1 count) in at least 80% of
# samples
# get expressed genes per tissue
exprGenes <- lapply(whiteDE, 
                    function(x) rownames(x)[apply(x, 1, 
                                                 function(y) sum(y > 1)) > 
                                                 0.8 * ncol(x)])
TSE <- Vennerable::Venn(exprGenes)@IntersectionSets
# remove redundunt element
TSE$`0000` <- NULL
tissueIndex <- sapply(names(TSE),
                      function(x) which(strsplit(x, '')[[1]] != "0"))
# rename, so names contain tissue
names(TSE) <- sapply(1:length(TSE), 
                     function(x) paste0(names(exprGenes)[tissueIndex[[x]]], 
                                              collapse = '-'))
TSE <- TSE[allTissues]
saveRDS(TSE, 'Rds/TSE.Rds')

# Differential expression between tissues with DESeq2 -------------------------
# info table for DE
setkey(INFO, rightGT_name)
de_INFO <- INFO[colnames(whiteDE_DF)]
de_INFO <- data.frame(de_INFO)
rownames(de_INFO) <- de_INFO$rightGT_name
de_INFO$Lib <- as.factor(de_INFO$Lib)
de_INFO$Tissue <- as.factor(de_INFO$Tissue)

# perform DE
dds <- DESeqDataSetFromMatrix(whiteDE_DF, colData = de_INFO,
                              design = ~ Lib + Tissue)
#dds <- DESeq(dds)
#saveRDS(dds, paste0(saveRdsTo, 'DE_between_tissues_DESeq2.Rds'))
#dds <- readRDS(paste0(saveRdsTo, 'DE_between_tissues_DESeq2.Rds'))

## Iterate over all possible pairs of tissues to get ALL genes UP regulated in
## a tissue of interest in comparison to at least one other tissue
#pairWiseDE <- list()
#for (toi in allTissues) {
#  deg_for_toi <- list()
#  for (tissueToCompare in allTissues[allTissues != toi]) {
#    onePair_de <- results(dds, contrast = c('Tissue', toi, tissueToCompare))
#    
#    # create data frame with Gene, TOI = tissue of interest, 
#    # Compared_with = tissue compared with and results of DE
#    onePair_de <- data.frame(Gene = rownames(onePair_de), TOI = toi,
#                             Compared_with = tissueToCompare, onePair_de)
#    onePair_de <- as.data.table(onePair_de)
#    
#    # then DESeq2 can't estimate FDR/pval/logFC, it put's NA, let's replace 
#    # them with 0 for log2FC, and 1 for pval and FDR respectively
#    onePair_de[, log2FoldChange := ifelse(is.na(log2FoldChange), 0, 
#                                                    log2FoldChange)]
#    onePair_de[, pvalue := ifelse(is.na(pvalue), 1, pvalue)]
#    onePair_de[, padj := ifelse(is.na(padj), 1, padj)]
#    
#    deg_for_toi[[length(deg_for_toi) + 1]] <- onePair_de
#  }
#  deg_for_toi <- do.call(rbind, deg_for_toi)
#  pairWiseDE[[length(pairWiseDE) + 1]] <- deg_for_toi
#}
#saveRDS(pairWiseDE, paste0(saveRdsTo, 'pairWiseDEbetweenTissues.Rds'))
pairWiseDE <- readRDS(paste0(saveRdsTo, 'pairWiseDEbetweenTissues.Rds'))
names(pairWiseDE) <- sapply(pairWiseDE, 
                            function(x) as.character(unique(x$TOI)))

# Append to the table S2 with info about cycling in tissues -------------------
toAppend <- data.table()

# tableS2 - just gene, tissue, index
tableS2 <- fread('../../ForOutput.csv', header = T, stringsAsFactors = F)
for (i in 1:nrow(tableS2)){
  complTissues <- sort(setdiff(names(pairWiseDE), tableS2[i]$Tissue))
  
  toAdd <- pairWiseDE[[tableS2[i]$Tissue]][Gene == tableS2[i]$Gene]
  toAdd <- toAdd[, .(Gene, TOI, baseMean, Compared_with, log2FoldChange, 
                     stat, padj)]
  
  if (length(unique(toAdd$Gene)) == 1 & length(unique(toAdd$TOI)) == 1 &
      length(unique(toAdd$baseMean)) == 1) {
    res <- data.table(Gene = unique(toAdd$Gene), TOI = unique(toAdd$TOI), 
                      baseMean = unique(toAdd$baseMean),
                      toAdd[Compared_with == complTissues[1]][, 4:7, with = F],
                      toAdd[Compared_with == complTissues[2]][, 4:7, with = F],
                      toAdd[Compared_with == complTissues[3]][, 4:7, with = F])
    toAppend <- rbind(toAppend, res)
    
  } else {
    if (nrow(toAdd) == 0) {
      res <- data.table(Gene = NA, TOI = NA, baseMean = NA, 
                        Compared_with = NA, log2FoldChange = NA, 
                        stat = NA, padj = NA,
                        Compared_with = NA, log2FoldChange = NA, 
                        stat = NA, padj = NA,
                        Compared_with = NA, log2FoldChange = NA, 
                        stat = NA, padj = NA)
      toAppend <- rbind(toAppend, res)
    } else {
      print(paste0('olala', i))
      stop()
    }
  }
}

# Assign DE status to tissue specific cycling genes ---------------------------
tisSpecCycl_DE <- lapply(names(pairWiseDE), 
                         function(x) pairWiseDE[[x]][Gene %in% 
                                                       tisSpecCycl[CyclesIn == x]$Gene])
names(tisSpecCycl_DE) <- names(pairWiseDE)


# Create bar plot showing how many tissue specific cycling genes are de -------
# TSC genes which are DE with at least one tissue
cyclAndDEpatterns <- lapply(tisSpecCycl_DE, 
                            function(x) x[abs(log2FoldChange) > 2 & 
                                          padj < 0.05])
cyclAndDEpatterns <- lapply(cyclAndDEpatterns, 
                            function(x) x[,.(Gene, Compared_with, 
                                             log2FoldChange)])
lapply(cyclAndDEpatterns, 
       function(x) x[, log2FoldChange := ifelse(log2FoldChange > 0, 
                                                'Up', 'Down')])

# patterns of gene de directions: if gene is de with 2+ tissues, it can be up 
# in one tissue and down in the other
directPttrn <- lapply(cyclAndDEpatterns, 
                      function(x) x[,.(paste0(unique(log2FoldChange),
                                              collapse = '')), 
                                    by = Gene])
lapply(directPttrn, function(x) x[, V1 := ifelse(V1 == 'Up' | V1 == 'Down',
                                                 V1, 'Multi')])
lapply(directPttrn, setnames, 'V1', 'deDirect')

# merge cyclAndDEpatterns and directPttrn
cyclAndDEpatterns <- lapply(cyclAndDEpatterns, 
                            function(x) x[,.(paste(sort(as.character(Compared_with)), 
                                                   collapse = '-')), 
                                          by = Gene])
lapply(cyclAndDEpatterns, setnames, 'V1', 'Compared_with')
lapply(directPttrn, setkey, Gene)
lapply(cyclAndDEpatterns, setkey, Gene)
cyclAndDEpatterns <- lapply(names(cyclAndDEpatterns), 
                            function(x) merge(cyclAndDEpatterns[[x]], 
                                              directPttrn[[x]], all = T))
names(cyclAndDEpatterns) <- names(directPttrn)

# summary data table for barplot
cyclAndDEpatterns_N <- lapply(cyclAndDEpatterns,
                              function(x) x[,.N, by = .(Compared_with, 
                                                       deDirect)])
cyclAndDEpatterns_N <- lapply(names(cyclAndDEpatterns_N), 
                              function(x) data.table(Tissue = x, 
                                                     cyclAndDEpatterns_N[[x]]))
cyclAndDEpatterns_N <- do.call(rbind, cyclAndDEpatterns_N)

# TSC genes which are NOT DE with none of tissues
notDE <- lapply(tisSpecCycl_DE, 
                function(x) x[!(abs(log2FoldChange) > 2 &  padj < 0.05)])
notDE <- lapply(notDE, function(x) x[,.N, by = Gene][N == 3])

notDE_N <- lapply(notDE, function(x) length(unique(x$Gene)))
notDE_N <- lapply(names(notDE_N), 
                  function(x) data.table(Tissue = x, Compared_with = "Not DE",
                                         deDirect = 'Up', N = notDE_N[[x]]))
notDE_N <- do.call(rbind, notDE_N)
cyclAndDEpatterns_N <- rbind(cyclAndDEpatterns_N, notDE_N)

# make the plot
cyclAndDEpatterns_N[, Compared_with := gsub('BRAIN', 'Brain', Compared_with)]
cyclAndDEpatterns_N[, Compared_with := gsub('GUT', 'Gut', Compared_with)]
cyclAndDEpatterns_N[, Tissue := gsub('BRAIN', 'Brain', Tissue)]
cyclAndDEpatterns_N[, Tissue := gsub('GUT', 'Gut', Tissue)]
cyclAndDEpatterns_N[, Tissue := gsub('FB', 'Fat body', Tissue)]
cyclAndDEpatterns_N[, Tissue := gsub('MT', 'Malpighian tubules', Tissue)]
cyclAndDEpatterns_N$Compared_with <- factor(cyclAndDEpatterns_N$Compared_with,
                                     c("Not DE", "Brain", "FB", "Gut", "MT",
                                       "Brain-FB", "Brain-Gut", "Brain-MT",
                                       "FB-Gut", "FB-MT", "Gut-MT", 
                                       "Brain-FB-Gut", "Brain-FB-MT", 
                                       "Brain-Gut-MT", "FB-Gut-MT"))
cyclAndDEpatterns_N$Tissue <- factor(cyclAndDEpatterns_N$Tissue,
                                     levels = sort(unique(cyclAndDEpatterns_N$Tissue)))

pdf(paste0('figures/edited/Fig1F_', 'v13', '.pdf'),
    height = 0.8 * 200 * 2 / 90, width = 250 * 2 / 90)
ggplot(cyclAndDEpatterns_N, aes(x = Compared_with, y = N, fill = deDirect)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label = N), vjust = -0.5,
            position = position_dodge(0.9), size = 2) +
  scale_fill_manual(values = c('#910E08', "black", "#2A987B")) +
  facet_wrap(~ Tissue, scales = 'free') +
  xlab('Tissue(s) compared with') + ylab("# of genes") +
  theme_classic(base_size = 10) +  
  theme(axis.line.x = element_line(colour = 'black', size = 0.5, 
                                   linetype = 'solid'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black',
                                   size = 6),
        axis.text.y = element_text(colour = 'black', size = 6),
        axis.line.y = element_line(colour = 'black', size = 0.5,
                                   linetype = 'solid'),
        strip.text.x = element_text(size = 8),
        panel.background = element_blank(), plot.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), legend.position = "none",
        strip.background = element_rect(colour = "white", fill = "white")) 
dev.off()

# Level of expression for DE and TSC genes ------------------------------------
# I want to see if genes which are up-regulated between tissue of interest and
# all of the rest tissue (i.e. Brain vs FB-Gut-MT) are expressed on the both 
# ends of comparison or for example, genes are not expressed in FB-Gut-MT

# normalized by voom count
normCounts <-  voom(counts = whiteDE_DF, normalize.method = "quantile")$E

# number of tissues in which de gene is expressed
numbOfTissGeneExpr <- data.table()

for (toi in names(cyclAndDEpatterns)) {
  # complement tissues
  complTiss <- paste(sort(setdiff(names(cyclAndDEpatterns), toi)), 
                     collapse = '-')
  
  # tsc genes which are de 
  deGenesToi <- cyclAndDEpatterns[[toi]]
  deGenesToi <- deGenesToi[Compared_with == complTiss & deDirect == 'Up']$Gene

  # all tsc
  #deGenesToi <- unique(tisSpecCycl_DE[[toi]]$Gene)
  
  deGenesToi <- as.character(deGenesToi)
  
  # norm. counts for those genes
  normCounts_deGenesToi <- as.data.table(melt(normCounts[deGenesToi, ]))
  setnames(normCounts_deGenesToi, colnames(normCounts_deGenesToi),
           c('Gene', 'rightGT_name', 'Expr'))
  normCounts_deGenesToi[, Tissue := gsub('[0-9].*', '', rightGT_name)]
  
  # mean expression in tissue
  meanNormCounts_deGenesToi <- normCounts_deGenesToi[,.(tissMean = mean(Expr)),
                                                     by = .(Tissue, Gene)]
  # binarize it
  meanNormCounts_deGenesToi[, tissMean := ifelse(tissMean > 0, 1, 0)]
  # see, if gene is expressed not only in tissue of interest
  meanNormCounts_deGenesToi <- meanNormCounts_deGenesToi[,.(sum(tissMean)), 
                                                         by = Gene]

  toAdd <- melt(meanNormCounts_deGenesToi[,.("More than 1" = sum(V1 > 1), 
                                             "TOI only" = sum(V1 == 1))])
  toAdd[, Tissue := toi]
  numbOfTissGeneExpr <- rbind(numbOfTissGeneExpr, toAdd)
}

numbOfTissGeneExpr[order(Tissue)]

# GO of NOT DE genes ----------------------------------------------------------
ensembl91 <- useEnsembl("ensembl", version = 91, 
                        dataset = "dmelanogaster_gene_ensembl")

goTSC_NotDE <- data.table()
goTSC_DE <- data.table()

for (toi in names(cyclAndDEpatterns)) {
  # not DE
  target_notDE <- getFlyBaseID(notDE[[toi]]$Gene, ensembl91)[[1]]
  bg <- getFlyBaseID(tisSpecCycl[CyclesIn == toi]$Gene, ensembl91)[[1]]
  
  toAdd <- lapply(c('BP', 'MF', 'CC'), 
                  function(x) GOenricment(target_notDE$flybase_gene_id,
                                          bg$flybase_gene_id, ont = x, 
                                          topNodes = 30))
  toAdd <- lapply(toAdd, as.data.table)
  names(toAdd) <- c('BP', 'MF', 'CC')
  lapply(names(toAdd), function(x) toAdd[[x]][, ont := x])
  toAdd <- do.call(rbind, toAdd)
  toAdd[, Tissue := toi]
  goTSC_NotDE <- rbind(goTSC_NotDE, toAdd)
  
  # DE
  target_DE <- getFlyBaseID(cyclAndDEpatterns[[toi]]$Gene, ensembl91)[[1]]
  bg <- getFlyBaseID(tisSpecCycl[CyclesIn == toi]$Gene, ensembl91)[[1]]
  
  toAdd <- lapply(c('BP', 'MF', 'CC'), 
                  function(x) GOenricment(target_DE$flybase_gene_id,
                                          bg$flybase_gene_id, ont = x, 
                                          topNodes = 30))
  toAdd <- lapply(toAdd, as.data.table)
  names(toAdd) <- c('BP', 'MF', 'CC')
  lapply(names(toAdd), function(x) toAdd[[x]][, ont := x])
  toAdd <- do.call(rbind, toAdd)
  toAdd[, Tissue := toi]
  goTSC_DE <- rbind(goTSC_DE, toAdd)
}
goTSC_NotDE[, padj := p.adjust(Fisher, method = 'fdr')]
View(goTSC_NotDE[ont == 'BP'][order(padj)][,.SD[1:8], Tissue])
goTSC_NotDE[ont == 'BP'][order(padj)][,.SD[1:8], Tissue][grepl('intra', Term)][order(Tissue)]

goTSC_DE[, padj := p.adjust(Fisher, method = 'fdr')]
View(goTSC_DE[ont == 'MF'][order(padj)][,.SD[1:8], Tissue])
goTSC_DE[ont == 'BP'][order(padj)][,.SD[1:8], Tissue]
goTSC_DE[ont == 'CC'][order(padj)][,.SD[1:8], Tissue]
