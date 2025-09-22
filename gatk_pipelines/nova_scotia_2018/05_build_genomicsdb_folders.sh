#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH -p noor
#SBATCH -J genomics_db-%j
#SBATCH -o tmp/genomics_db-%j.out

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

# sample list for genomics db import
ls $gvcf/*.gz  > tmp/gvcf_list.list

# array of contig names
#declare -a arr=("chrI" "chrII" "chrIII" "chrIV" "chrV" "chrVI" "chrVII" "chrVIII" "chrIX" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI" "chrXVII" "chrXVIII" "chrXIX" "chrXX" "chrXXI" "chrUn" "chrM")

while read interval
do
sbatch 05_one_genomics_db.sh $interval
done < tmp/stickleback_intervals_string.txt




