#!/bin/bash

#SBATCH --mem=40G
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J hard_filter-%j
#SBATCH -o tmp/hard_filter-%j.out

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

# sample list is an argument (see previous scripts)
# I split the sample list for deploying script on 2 nodes of the cluster
ls $gvcf/*.vcf > tmp/gvcf.list

#Call GATK filtration tools


# retain only snps
#$java_1_8 -Xmx32g -jar $gatk \
#    -T SelectVariants \
#    --disable_auto_index_creation_and_locking_when_reading_rods \
#    -R $ref \
#    -nt 16 \
#    -V $vcf/dpse_iso_seq_raw.vcf\
#    -selectType SNP \
#    -o $vcf/dpse_iso_seq_snps.vcf

$java_1_8 -Xmx32g -jar $gatk \
    -T VariantFiltration \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -R $ref \
    -V $vcf/dpse_iso_seq_snps.vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "best_prac_filters" \
    --setFilteredGtToNocall \
    -o $vcf/dpse_iso_seq_hard_filtered_snps.vcf
