## Scipts for data pre-processing on HPC: QC/mapping/genotyping/read counting
Since the data uploaded to GEO were already demultiplexed, demultiplexing is not presented 
in the scripts below. Adding modules could be different on your cluster. fastqc was done.
genotypes were checked and re-labelled.

| File                             | Description            |
| :------------                    | :---                   |
| 4_STAR_2pass.sh  | |
| 5_map_dm3.sh  | |
| 7_GenotypeGVCF.sh  | |
| 8_dm3_GenotypeGVCF.R  | |
| 8_dm3_GenotypeWithR.sh  | |
| 9_dm3_CreateSimilarityMatrix.R  | |
