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
ls -d $bam_realigned/*.* | grep realign.bam  > $proj_home/bamlist.list

while read prefix
do
#Call GATK HaplotypeCaller
java -Xmx10g -jar $gatk/GenomeAnalysisTK.jar \
	-nct 12 \
	-l INFO \
	-R $ref \
	-log $log/$plate.$prefix.HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I  $bam_realigned/$prefix.realign.bam \
	--emitRefConfidence GVCF \
	--variant_index_type LINEAR \
    --variant_index_parameter 128000 \
	--max_alternate_alleles 2 \
	-o $gvcf/${prefix}.GATK.g.vcf
	
done < $proj_home/Samplelist.txt
