#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 4_map.%I.err
#BSUB -o 4_map.%I.out
#BSUB -J 4_map[1-778]
#BSUB -M 16000000
#BSUB -R rusage[mem=16000]
#BSUB -n 4
#BSUB -u maria.litovchenko@epfl.ch

# WARNING!
# Create directory called tools in your project directory and ensure presence
# of picard.jar and GenomeAnalysisTK.jar in it before running this script

export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Aligner/STAR/2.5.0b;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Analysis/HTSeq/0.6.1;

fastqDir=2_demultiplex/
mappedDir=3_mapped/
countsDir=4_counts/
cleanSamDir=5-a-cleanSam/
replaceRGdir=5-b-rgReplaced/
markedDuplDir=5-c-MerkedDupl/
SplitAtrimDir=5-d-SplitAtrim/
targetDir=5-e-target/
realignDir=5-f-realign/
vcfDir=5-rawGVCF/
mkdir $mappedDir
mkdir $countsDir
mkdir $cleanSamDir
mkdir $replaceRGdir
mkdir $markedDuplDir
mkdir $SplitAtrimDir
mkdir $targetDir
mkdir $realignDir
mkdir $vcfDir

genomeDir=RefGen/dm3_2pass/;
genomeFa=RefGen/dm3/dm3.Wolb.fa
inputGTF=RefGen/dm3/dm3refGene.srt.gtf;
GATKjar=tools/GenomeAnalysisTK.jar
picardJar=tools/picard.jar
end="_Aligned.sortedByCoord.out.bam"
endClean="_Aligned.sortedByCoord.out.Clean.bam"
endRGfixed="_Aligned.sortedByCoord.out.RGrepl.bam"
endMarkDupl="_Aligned.sortedByCoord.out.MD.bam"
endSplitAtrim="_Aligned.sortedByCoord.out.SaT.bam"
bamFilterJar=tools/BAMFiltering-1.0.jar

toMap=($(ls -d -1 $fastqDir/*/* | grep .fastq.gz$ ))
zeroArr=( zero )
firstInPair=("${zeroArr[@]}" "${toMap[@]}")

R1=${firstInPair[${LSB_JOBINDEX}]};
sample=`echo $R1 | sed 's@.*x/@@g' | sed 's@_/@_@g' | sed 's/.fastq.gz//g'`

#------------------------------------------------------------------------------
# 1) Map data to the reference genome with STAR
STAR --runMode alignReads --runThreadN 4 \
    --genomeDir $genomeDir \
    --outFilterScoreMinOverLread 0.20 \
    --outFilterMatchNminOverLread 0.20 \
    --outFilterMultimapNmax 1 \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $mappedDir$sample"_" \
    --readFilesIn $R1
samtools index $mappedDir$sample"_"Aligned.sortedByCoord.out.bam
samtools flagstat $mappedDir$sample"_"Aligned.sortedByCoord.out.bam > $mappedDir$sample"_"Aligned.sortedByCoord.out.flagstat

#------------------------------------------------------------------------------
# CLEAN BAM WITH bamFilterJar
java -jar $bamFilterJar -f $mappedDir$sample$end  -maxLength 1000 -filter \
          -o $mappedDir
end="_Aligned.sortedByCoord.out.filtered.bam"

#------------------------------------------------------------------------------
# COUNT READS OVERLAPPING WITH GENES
htseq-count -s no -m union \
            -f bam $mappedDir$sample$end \
             $inputGTF > $countsDir$sample".counts"

#------------------------------------------------------------------------------
# PERFORM A GENOTYPING CALLING ACCORDING TO GATK PIPELINE
# https://software.broadinstitute.org/gatk/guide/article?id=3891
# CLEAN SAM, JUST IN CASE...
java -jar $picardJar CleanSam I=$mappedDir$sample$end \
                              O=$cleanSamDir$sample$endClean

#------------------------------------------------------------------------------
# ADD READ GROUPS
java -jar $picardJar AddOrReplaceReadGroups \
          I=$cleanSamDir$sample$endClean \
          O=$replaceRGdir$sample$endRGfixed \
          SO=coordinate RGID=id RGLB=library \
          RGPL=platform RGPU=machine RGSM=sample

#------------------------------------------------------------------------------
# MARK DUPLICATES
java -jar $picardJar MarkDuplicates I=$replaceRGdir$sample$endRGfixed  \
                                    O=$markedDuplDir$sample$endMarkDupl \
          CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
          M=$mappedDir$sample".metrics"

#------------------------------------------------------------------------------
# Split'N'Trim and reassign mapping qualities
java -jar $GATKjar -T SplitNCigarReads \
     -R $genomeFa \
     -I $markedDuplDir$sample$endMarkDupl \
     -o $SplitAtrimDir$sample$endSplitAtrim \
     -rf ReassignOneMappingQuality -RMQF 255 \
     -RMQT 60 -U ALLOW_N_CIGAR_READS

#------------------------------------------------------------------------------
# Local realignment around indels
java -jar $GATKjar -T RealignerTargetCreator \
     -R $genomeFa -I $SplitAtrimDir$sample$endSplitAtrim \
     -o $targetDir$sample.intervals

#------------------------------------------------------------------------------
# Perform realignment
java -jar $GATKjar -T IndelRealigner -R $genomeFa \
     -I $SplitAtrimDir$sample$endSplitAtrim \
     -targetIntervals $targetDir$sample.intervals \
     -o $realignDir$sample.realign.bam

#------------------------------------------------------------------------------
# call variants
java -jar $GATKjar -T HaplotypeCaller -dontUseSoftClippedBases \
     -nct 8 -R $genomeFa -I $realignDir$sample.realign.bam \
     -stand_call_conf 30 -stand_emit_conf 10 \
     --emitRefConfidence GVCF -o $vcfDir$sample.raw.g.vcf
toReplace=sample
sed -i "s@$toReplace@$sample@g" $vcfDir$sample.raw.g.vcf

sed '/Wolbachia/d' $vcfDir$sample.raw.g.vcf > $vcfDir$sample.raw.g.noWolbNoM.vcf
sed -i '/dmel_mitochondrion_genome/d' $vcfDir$sample.raw.g.noWolbNoM.vcf

#------------------------------------------------------------------------------
#select no indels
java -jar $GATKjar -T SelectVariants -R $genomeFa \
     --variant $vcfDir$sample.raw.g.noWolbNoM.vcf \
     -o $vcfDir$sample.raw.g.noWolbNoM.SNPs.vcf \
     --selectTypeToExclude INDEL \
     --selectTypeToExclude MNP -select "DP > 5"

exit 0;
