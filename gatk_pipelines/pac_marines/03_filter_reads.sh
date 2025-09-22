#!/bin/bash

# Demultiplexes paired-end GBS reads (fastq) using a perl script
# !!!!assumes PstI enzyme!!!!
# USAGE: ./SB_gatk_pipe.sh <plate identifier> <path to barcode file>

# January 2016 - Kieran Samuk
# Modified from stickleback meta analysis pipeline (github.com/ksamuk/gene_flow_linkage)

project="pac_marine"

# folders
proj_home="/home/samuk/pac_marines"
ref="/home/samuk/gbs2015/ref/revisedAssemblyMasked.fa"
raw="$proj_home/data/raw"
fastq="$proj_home/data/raw/fastq"
cleandata="$proj_home/data/raw/fastq"
trimmeddata="$proj_home/data/trimmed"
unpaired="$proj_home/data/unpaired"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf"

# tools
demultiplex='~/bin/GBS_fastq_Demultiplexer_v8.GO.pl'
barcodes=''
tabix="/home/samuk/bin/tabix-0.2.6"
vcftools="/home/samuk/bin/vcftools_0.1.11/bin"
picardtools="/home/samuk/bin/picard-tools-1.96"
gatk="/home/samuk/bin/GATK"
trim="/home/samuk/bin/Trimmomatic-0.32"

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

# make list of samples

ls $cleandata | grep -v "nobar" | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > $proj_home/Samplelist.txt

# #Trim the data using Trimmomatic. Removes bad reads and illumina adapter contamination.
 
while read prefix
do
echo "Trimming ${prefix}..."
#note: TruSeq3-PE.fa should be TruSeq2-PE.fa, actually. Used 3 for meta though.
java -Xmx20g -jar $trim/trimmomatic-0.32.jar PE -threads 8 -phred33 $cleandata/"$prefix"_R1.fastq $cleandata/"$prefix"_R2.fastq $trimmeddata/"$prefix"_R1.fastq $unpaired/"$prefix"_unR1.fastq $trimmeddata/"$prefix"_R2.fastq $unpaired/"$prefix"_unR2.fastq SLIDINGWINDOW:4:15 MINLEN:36
done < Samplelist.txt

