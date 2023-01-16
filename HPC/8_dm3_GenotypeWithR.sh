#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 8_GenoR.%I.err
#BSUB -o 8_GenoR.%I.out
#BSUB -J 8_GenoR[1-844]
#BSUB -M 100000000
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

export PATH=/software/bin:$PATH;
module use /software/module/;
module add R/3.3.2;

#------------------------------------------------------------------------------
# PERFORM GENOTYPING WITH USE OF CUSTOM R SCRIPT
#------------------------------------------------------------------------------
# number of columns (DGRP files)
# head -n 1 MitoSeq_AllRuns_dm3.gvcf.snps.fltr.GT.FORMAT | awk '{print NF}'
# substract 2 from this number

outputdir=/scratch/el/monthly/mlitovch/AroundTheClock_Aug2017;

samp=($(seq 3 846))
zeroArr=( zero )
samples=("${zeroArr[@]}" "${samp[@]}")

sample=${samples[${LSB_JOBINDEX}]};

echo $sample
# run R script
mkdir $outputdir/8_dm3_similarityVectors
Rscript --vanilla $outputdir/8_dm3_GenotypeGVCF.R \
                  /scratch/el/monthly/mlitovch/RefGen/dm3/dgrp2.SNPs.GT.FORMAT \
                  $outputdir/AroundTheClock_Aug2017_dm3.snps.flrt.GT.FORMAT \
                  $sample

exit 0;
