#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 7_GVCFs.err
#BSUB -o 7_GVCFs.out
#BSUB -J 7_GVCFs
#BSUB -M 100000000
#BSUB -R rusage[mem=100000]
#BSUB -n 4
#BSUB -u maria.litovchenko@epfl.ch

export PATH=/software/bin:$PATH;
module use /software/module/;
module add Development/java_jdk/1.8.0_66;
module add UHTS/Analysis/vcftools/0.1.12b;
module add R/3.3.2;

#------------------------------------------------------------------------------
# PERFORM MERGING OF RAW GVCF WITH GenotypeGVCFs
#------------------------------------------------------------------------------
outputdir='.';
GATKjar=tools/GenomeAnalysisTK.jar
ref=RefGen/dm3/dm3.Wolb.fa
rawGVCFdir=5-rawGVCF

# list samples to merge and remove bad samples
toMerge=`ls -d -1 $rawGVCFdir/* | grep noWolbNoM.SNPs.vcf$ | xargs | sed "s@/scratch@--variant /scratch@g"`

# try first use mergeGVCF to make it faster
java  -Xms64g -Xmx64g -jar $GATKjar -T CombineGVCFs \
      -R $ref $toMerge -o $outputdir/AroundTheClock_Aug2017_dm3_merged.g.vcf

# do GenotypeGVCFs, it improves genotype calls
java -Xms64g -Xmx64g -jar $GATKjar -T GenotypeGVCFs \
     -allSites -stand_call_conf 30 -stand_emit_conf 10 \
     -R $ref -nt 4 \
     --variant $outputdir/AroundTheClock_Aug2017_dm3_merged.g.vcf \
     -o $outputdir/AroundTheClock_Aug2017_dm3_geno.gvcf.vcf
echo "Finished GenotypeGVCFs"

# select only SNPs
java -Xms64g -Xmx64g -jar $GATKjar -T SelectVariants \
     -R $ref -V $outputdir/AroundTheClock_Aug2017_dm3_geno.gvcf.vcf \
     -selectType SNP --selectTypeToExclude MNP -select "DP > 5" \
     -o $outputdir/AroundTheClock_Aug2017_dm3.snps.gvcf.vcf
echo "Finished SelectVariants"

# filter out bad quality
java -Xms64g -Xmx64g -jar $GATKjar -T VariantFiltration \
     -R $ref -V $outputdir/AroundTheClock_Aug2017_dm3.snps.gvcf.vcf \
     -filterName FS -filter "FS > 30.0" -filterName QD \
     -filter "QD < 2.0" \
     -o $outputdir/AroundTheClock_Aug2017_dm3.snps.flrt.gvcf.vcf
echo "Finished VariantFiltration"

# extract only genotype info
vcftools --vcf $outputdir/AroundTheClock_Aug2017_dm3.snps.flrt.gvcf.vcf \
         --extract-FORMAT-info GT \
         --out $outputdir/AroundTheClock_Aug2017_dm3.snps.flrt
echo "Finished GT extraction"

# create table suitable for follow up R script
sed -i 's@\.\/\.@-@g' $outputdir/AroundTheClock_Aug2017_dm3.snps.flrt.GT.FORMAT
sed -i 's@0\/0@0@g' $outputdir/AroundTheClock_Aug2017_dm3.snps.flrt.GT.FORMAT
sed -i 's@1\/1@2@g' $outputdir/AroundTheClock_Aug2017_dm3.snps.flrt.GT.FORMAT
echo "Finished preparation of file for R"

exit 0;
