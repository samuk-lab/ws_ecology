#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J align
#SBATCH -o tmp/bwa_align_out-%j.out

#Pipeline for using GATK and BWA to call snps
# KMS April 2018

# folders
proj_home="/work/kms173/wht_stbk"
ref="$proj_home/ref/revisedAssemblyMasked.fa"
demulti="$proj_home/data/trimmed"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf"

# plate, barcode and project ids
project="wht_stbk"

# tools
bin_folder="/dscrhome/kms173/bin"
java_1_8="$bin_folder/jdk1.8.0_101/bin/java"
tabix="$bin_folder/tabix-0.2.6"
vcftools="$bin_folder/vcftools_0.1.11/bin"
picardtools="$bin_folder/picard.jar"
gatk="$bin_folder/GATK"
bwa="$bin_folder/bwa"

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

project="iso_seq"


# sample list is an argument (see previous scripts)
# I split the sample list for deploying script on 2 nodes of the cluster
samplelist=$1


###Align using BWA. Turn from sam to bam. Sort by coordinate and add read group data.
while read prefix
do
  
  #convert ubam to fastq, align, merge bams
  
  ~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar SamToFastq \
  I=$bam/"$prefix"_markadapters.bam \
  FASTQ=/dev/stdout \
  CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
  TMP_DIR=./tmp | \
  $bwa mem -M -t 18 -p $ref /dev/stdin | \
  ~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar MergeBamAlignment \
  R=$ref \
  UNMAPPED_BAM=$bam/"$prefix"_unmapped.bam \
  ALIGNED_BAM=/dev/stdin \
  OUTPUT=$bam/"$prefix"_merged.bam \
  CREATE_INDEX=true ADD_MATE_CIGAR=true \
  CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
  INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
  PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
  TMP_DIR=./tmp

  #mark pcr duplicates
  #not done for GBS
  
# ~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar MarkDuplicates \
#  INPUT=$bam/"$prefix"_merged.bam \
#  OUTPUT=$bam/"$prefix"_dedup.bam \
#  METRICS_FILE=tmp/"$prefix"_dedup_metrics.txt \
#  CREATE_INDEX=true \
#  TMP_DIR=./tmp

#done < bam.list
done < $samplelist
