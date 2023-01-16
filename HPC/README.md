## Scipts for data pre-processing on HPC: QC/mapping/genotyping/read counting

### Obtaining raw sequencing data 
In order to download the raw SRA data, install [SRA toolkit](), and use `prefetch`
command, i.e.:
```
prefetch 
```
To convert data from SRA to fastq format, use `fastq-dump`, i.e.:
```
fastq-dump 
```
You can find a full list of SRA accession codes in []().

### Data processing
Raw sequencing data uploaded to SRA are already demultiplexed. Therefore, 
demultiplexing is not presented in the scripts below. Also, the genotype match
to the original genotypes determined by VCF of DGRP2 panel was ensured. Files
were relabelled in case genotype swap was detected and uploaded to SRA/GEO 
under their proper molecular genotype name

| File                             | Description                                                              |
| :------------                    | :---                                                                     |
| 0_prepare_reference_genome.sh    | Downloads and indexes reference genome (dm3) for future use in the study |
| 
| 5_map_dm3.sh  | |
| 7_GenotypeGVCF.sh  | |
| 8_dm3_GenotypeGVCF.R  | |
| 8_dm3_GenotypeWithR.sh  | |
| 9_dm3_CreateSimilarityMatrix.R  | |
