#!/bin/bash

# USAGE: ./SB_gatk_pipe.sh <plate identifier> 
# January 2016 - Kieran Samuk
# Modified from the stickleback meta-analysis pipeline (github.com/ksamuk/gene_flow_linkage)
# Pipeline for using GATK and BWA to call snps on GBS data

# plate and project ids
plate_id="1"
project="north_sea"

# folders
proj_home="/home/samuk/north_sea"
data="$proj_home/data"
ref="/home/samuk/gbs2015/ref/revisedAssemblyMasked.fa"
raw="$proj_home/data/raw"
fastq="$proj_home/data/raw/fastq"
cleandata="$proj_home/data/combined"
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

# make directories if they don't exist
if [ ! -d "$data" ]; then
	mkdir $data
fi
if [ ! -d "$raw" ]; then
	mkdir $raw
fi
if [ ! -d "$fastq" ]; then
	mkdir $fastq
fi
if [ ! -d "$sam" ]; then
	mkdir $sam
fi
if [ ! -d "$bam" ]; then
	mkdir $bam
fi
if [ ! -d "$log" ]; then
	mkdir $log
fi
if [ ! -d "$gvcf" ]; then
        mkdir $gvcf
fi
if [ ! -d "$trimmeddata" ]; then
        mkdir $trimmeddata
fi
if [ ! -d "$unpaired" ]; then
        mkdir $unpaired
		
fi

