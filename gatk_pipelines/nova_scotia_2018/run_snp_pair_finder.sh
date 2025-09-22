#!/bin/bash

#SBATCH --mem=40G
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J bam_views
#SBATCH -o tmp/bam_views-%j.out

#Pipeline for using GATK and BWA to call snps
# KMS April 2018

# folders
proj_home="/dscrhome/kms173/noor2/kms173/iso_seq"
ref="$proj_home/data/ref_genome/dpse-all-chromosome-r3.04.fasta"
raw="$proj_home/data/raw_unzip"
demulti="$proj_home/data/demulti"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf2"

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

Rscript r_scripts/04_create_bam_views.R