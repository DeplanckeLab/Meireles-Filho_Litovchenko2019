# Extensive tissue-specific expression variation and novel regulators underlying circadian behavior

Source code to support the paper: "Extensive mitochondrial population structure and haplotype-specific
variation in metabolic phenotypes in the Drosophila Genetic Reference Panel" published in 
[Science Advances](https://advances.sciencemag.org/content/7/5/eabc3781). 

The data are freely avaible at [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126018).


| File                             | Description            |
| :------------                    | :---                   |
| HPC                              |                        |
| 0_common_functions.R             | Functions used though out the study                                                          | 
| 1_create_sampleInfoTab.R         | Script which created information table holding all meta data about samples                   |
| 2_plotMappingStats_countsRDS.R   | Script which plots mapping statistics of samples and creates Rds with integer counts         |
| JTK_CYCLEv3.1.R                  | |
| 3_detectCycles.R                 | Script which performs detection of cycling genes in all tissues with use of JTK CYCLE        |
| 6_LASSO.R                        | Script to perform LASSO regression on data                                                   |
| 7_NeuralNet.R                    | Script to perform NeuralNet regression on data                                               |
| 8_zeitzeiger.R                   | Script to run ZeitZeiger on data                                                             |
| 9_MolecularTimeTable.R           | Script to run molecular time table on data                                                   |
| 10_NN_Eval_parellel.R            | Script to evaluate various configurations of neural net to detected physiological time       | 
| 10_NN_Eval_parellel.sh           | Script - wrapper to run 10_NN_Eval_parellel.R                                                |
| 10_PhysioTimePredModelsEval.R    | Script to perform comparison of various models to detected physiological time                |
| 10_PhysioTimePredMTTonDGRP.R     | Estimation of physiological time of the static DGRP samples with use of molecular time table |
| 11_DEbetweenTissues.R            | Script to perform detection of differential expression between different tissues             |
| 11_DGRPoutliersVars.R            | Script to 
| 12_inputsForHomer.R              | Script to create input files for run of Homer                                                |
| 13_motifEnrichPostProc.R         | Script to perform postprocessing of Homer run                                                |
| 17_create_sampleInfoTab_28233.R  | |
| 18_plotMappingStats_28233.R      | |
| 19_detectCycles_28233.R          | |
| 20_compare_w_28233.R             | |
| 33_dynGENIE3.R                   | |
| 35_compareNets_white_28233.R ||
| 36A_activity_fromRawData.R ||
| 37_variantsInDGRP-796.R ||
| 39_down_Cell_ChiP-seq.sh ||
| 40_trimMapDedupl_Cell_ChiP-seq.sh ||
| 41_mergeReps_Cell_ChIP-seq.sh || 
| 42_bigwigs_Cell_ChiP-seq.sh  ||
| 43_geneBodyBed.R  ||
| 44_plotPeaks_ChIP-seq.R  ||
| 45_plotBaboonMouseCyclOrth.R  ||
| 46_phaseShift.R  ||
| drawNets2.R  ||
| figures.R  ||
| model_the_expression.R  ||

 
