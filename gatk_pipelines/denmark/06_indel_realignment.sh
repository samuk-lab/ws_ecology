#!/bin/bash

# Demultiplexes paired-end GBS reads (fastq) using a perl script
# !!!!assumes PstI enzyme!!!!
# USAGE: ./SB_gatk_pipe.sh <plate identifier> <path to barcode file>

# January 2016 - Kieran Samuk
# Modified from stickleback meta analysis pipeline (github.com/ksamuk/gene_flow_linkage)

# folders
proj_home="/home/samuk/north_sea"
ref="/home/samuk/gbs2015/ref/revisedAssemblyMasked.fa"
raw="$proj_home/data/raw"
fastq="$proj_home/data/raw/fastq"
cleandata="$proj_home/data/combined"
trimmeddata="$proj_home/data/trimmed"
unpaired="$proj_home/data/unpaired"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
bam_realigned="$proj_home/data/bam_realigned"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf"

# plate, barcode and project ids
plate="north_sea"
project="north_sea"

# tools
demultiplex="$proj_home/bin/GBS_fastq_Demultiplexer_v8.GO.pl"
tabix="/home/samuk/bin/tabix-0.2.6"
vcftools="/home/samuk/bin/vcftools_0.1.11/bin"
picardtools="/home/samuk/bin/picard-tools-1.96"
gatk="/home/samuk/bin/GATK"
trim="/home/samuk/bin/Trimmomatic-0.32"
bwa="/home/samuk/bin/bwa/bwa.kit/bwa"

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

#Make bam.list for GATK
ls -d $bam/*.* | grep sortrg.bam  > $proj_home/bamlist.list

#identify local indels
java -Xmx10g -jar $gatk/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $ref \
   -I $proj_home/bamlist.list \
   -nt 3 \
   -log $log/RealignerTargetCreator.log \
   -allowPotentiallyMisencodedQuals\
   -o $proj_home/realign.intervals

#Realign around local indels
while read prefix
do
java -Xmx20g -jar $gatk/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $ref \
	-I $bam/$prefix.sortrg.bam \
	-targetIntervals $proj_home/realign.intervals \
	-o $bam_realigned/$prefix.realign.bam \
	-log $log/$prefix.IndelRealigner.log 

done < $proj_home/Samplelist.txt
