#!/bin/bash

#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH -p noor
#SBATCH -J hap-mini-%j
#SBATCH -o log/hap-mini/haplotype-mini-%j.out

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
project="iso_seq"

# tools
bin_folder="/dscrhome/kms173/bin"
java_1_8="$bin_folder/jdk1.8.0_101/bin/java"
tabix="$bin_folder/tabix-0.2.6"
vcftools="$bin_folder/vcftools_0.1.11/bin"
picardtools="$bin_folder/picard.jar"
gatk="$bin_folder/GenomeAnalysisTK.jar"
bwa="$bin_folder/bwa"

#export _JAVA_OPTIONS="-Xmx20g"
#export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

project="iso_seq"

# sample is an argument
prefix=$1

#Call GATK HaplotypeCaller

 #Call GATK HaplotypeCaller
gatk --java-options "-Xmx4g" HaplotypeCaller \
   -R $ref \
   -I $bam/"$prefix"_merged.bam \
   -O $gvcf/$prefix.g.vcf.gz \
   --native-pair-hmm-threads 1 \
   -ERC GVCF
 