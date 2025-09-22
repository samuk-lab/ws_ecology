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

# assumes data demultiplexed, trimmed with trimmomatic, combined 
# build list of files to align
ls  $trimmeddata/BS*.fastq | xargs -n1 basename |sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | grep -v sai | uniq  > Samplelist.txt

align_samplelist=Samplelist.txt

# make a list of bam files
ls $bam | sed s/.sortrg.bai//g | sed s/.sortrg.bam//g | uniq  > bam_list.txt

# compare to the input file list
grep -v -f bam_list.txt Samplelist.txt > rerun_samples.txt

if [ -s $rerun_lines ]; then
	echo "The following files failed to align:"
	cat rerun_samples.txt
fi

#allows reruns on failed alignments, just keep running the script
#rerun_lines=$(wc -l rerun_samples.txt | awk '{print $1}')

if [ -s $rerun_lines  ]; then
	echo "Rerunning ${rerun_lines} samples..."
	align_samplelist="rerun_samples.txt"
fi

###Align using BWA. Turn from sam to bam. Sort by coordinate and add read group data.
echo "Aligning files in ${align_samplelist}..."
while read prefix
do
echo $prefix
	$bwa mem -M -t 8 $ref $trimmeddata/"$prefix"_R1.fastq $trimmeddata/"$prefix"_R2.fastq > $sam/$prefix.sam
	samtools view -Sb $sam/$prefix.sam > $bam/$prefix.bam	
	echo "Cleaning..."
	java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar $picardtools/CleanSam.jar INPUT=$bam/$prefix.bam OUTPUT=$bam/$prefix.clean.bam 2> $log/$prefix.cleansam.log
	echo "Sorting..."
	java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar $picardtools/SortSam.jar INPUT=$bam/$prefix.clean.bam OUTPUT=$bam/$prefix.sort.bam SORT_ORDER=coordinate 2> $log/$prefix.sortsam.log
	echo "Adding read groups..."
	java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar $picardtools/AddOrReplaceReadGroups.jar I=$bam/$prefix.sort.bam O= $bam/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True 2> $log/$prefix.addRG.log
	rm $sam/$prefix.sam
	rm $bam/$prefix.bam
	rm $bam/$prefix.clean.bam
	rm $bam/$prefix.sort.bam
done < $align_samplelist

##CHECK IF ALIGNMENT WAS SUCCESSFUL




