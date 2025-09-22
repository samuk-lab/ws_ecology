#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J haplotype-%j
#SBATCH -o tmp/haplotype-caller-%j.out

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

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

project="iso_seq"

# sample list is an argument (see previous scripts)
# I split the sample list for deploying script on 2 nodes of the cluster
samplelist=$1

#Realign around local indels
while read prefix
do

	#Call GATK HaplotypeCaller
	$java_1_8 -Xmx32g -jar $gatk \
	-nct 20 \
	-R $ref \
	-log $log/"$prefix"_HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I $bam/"$prefix"_dedup.bam \
	--emitRefConfidence GVCF \
	--max_alternate_alleles 4 \
	-o $gvcf/$prefix.g.vcf
	
done < $samplelist