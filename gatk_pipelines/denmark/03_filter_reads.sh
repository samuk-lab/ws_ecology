#!/bin/bash

# Demultiplexes paired-end GBS reads (fastq) using a perl script
# !!!!assumes PstI enzyme!!!!
# USAGE: ./SB_gatk_pipe.sh <plate identifier> <path to barcode file>

# January 2016 - Kieran Samuk
# Modified from stickleback meta analysis pipeline (github.com/ksamuk/gene_flow_linkage)

plate_id="1"
project="north_sea"

# folders
proj_home="/home/samuk/north_sea"
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

### NOTE THIS HAS BEEN MODIFIED TO RERUN ON TRIMMED DATA
### SEE ORIIGNAL PIPEFILE FOR CORRECTION

# make list of samples

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

#ls $cleandata *.fastq | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > $proj_home/Samplelist.txt
#ls $trimmeddata | grep -v "nobar" | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > $proj_home/Samplelist.txt


# #Trim the data using Trimmomatic. Removes bad reads and illumina adapter contamination.
 
while read prefix
do
echo "Trimming ${prefix}..."
#note: genome quebec reads have adapters pre-removed, so we don't use the ILLUMINACLIP: argument as normal
java -Xmx20g -jar $trim/trimmomatic-0.32.jar PE -threads 8 $cleandata/"$prefix"_R1.fastq $cleandata/"$prefix"_R2.fastq $trimmeddata/"$prefix"_R1.fastq $unpaired/"$prefix"_unR1.fastq $trimmeddata/"$prefix"_R2.fastq $unpaired/"$prefix"_unR2.fastq SLIDINGWINDOW:4:20 MINLEN:50
done < $proj_home/Samplelist_trim.txt

