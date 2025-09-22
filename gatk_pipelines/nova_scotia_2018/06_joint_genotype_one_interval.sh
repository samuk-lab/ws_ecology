#!/bin/bash

#SBATCH --mem=9000
#SBATCH --cpus-per-task=2
#SBATCH -p scavenger
#SBATCH -J joint_geno-%j
#SBATCH -o tmp/joint-geno/joint_geno-%j.out

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
gdb="$proj_home/data/gdb"
vcf="$proj_home/data/vcf"

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

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

project="iso_seq"

#Call GATK joint genotyper

 gatk --java-options "-Xmx8g" GenotypeGVCFs \
   -R $ref \
   --intervals $1 \
   -V gendb://$gdb/"$1" \
   -new-qual true \
   -O $vcf/"$1".vcf.gz
