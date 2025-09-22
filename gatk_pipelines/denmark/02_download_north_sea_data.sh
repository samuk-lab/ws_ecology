#!/bin/bash

# USAGE: ./SB_gatk_pipe.sh <plate identifier> 
# January 2016 - Kieran Samuk
# Modified from the stickleback meta-analysis pipeline (github.com/ksamuk/gene_flow_linkage)
# Pipeline for using GATK and BWA to call snps on GBS data

# plate and project ids
plate_id="1"
project="pac_marines"

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
sra="/SciBorg/array0/samuk/ncbi/public/sra"
tmp="$proj_home/data/tmp"

# tools
demultiplex='~/bin/GBS_fastq_Demultiplexer_v8.GO.pl'
barcodes=''
tabix="/home/samuk/bin/tabix-0.2.6"
vcftools="/home/samuk/bin/vcftools_0.1.11/bin"
picardtools="/home/samuk/bin/picard-tools-1.96"
gatk="/home/samuk/bin/GATK"
trim="/home/samuk/bin/Trimmomatic-0.32"
sratoolkit="/home/samuk/bin/sratoolkit/bin"

echo "Donwloading SRA files..."

#North Sea Accessions
srx_accessions="ERS123916 ERS123917 ERS123918 ERS123919 ERS123920 ERS123921"

for srx in $srx_accessions
	do
	wget -qO- "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$srx" | grep $srx | cut -f1 -d"," > srr_accessions.txt

	while read srr
		do
		echo "Downloading $srr..."
		$sratoolkit/prefetch $srr
		echo "fastq-dumping $srr..."
		$sratoolkit/fastq-dump --split-files /SciBorg/array0/samuk/ncbi/public/sra/$srr.sra --outdir $fastq/$srx
	done < srr_accessions.txt
	
	ls $fastq/$srx | grep "_1" > tmp/read_1.txt
	ls $fastq/$srx | grep "_2" > tmp/read_2.txt
	
	while read r1
		do 
		cat $fastq/$srx/$r1 >> $fastq/${srx}_R1.fastq
	done < tmp/read_1.txt
	
	while read r2
		do 
		cat $fastq/$srx/$r2 >> $fastq/${srx}_R2.fastq
	done < tmp/read_2.txt
done

