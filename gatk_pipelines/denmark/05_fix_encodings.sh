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
bam_fixed="$proj_home/data/bam_fixed"
bam_fixed_rg="$proj_home/data/bam_fixed_rg"
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

while read file
do
java -Xmx20g -jar $gatk/GenomeAnalysisTK.jar \
   -T PrintReads \
   -R $ref \
   -I $file \
   -allowPotentiallyMisencodedQuals \
   -o $file.fixed.bam \
   
mv $bam/*fixed* $bam_fixed
   
done < bamlist.list

#Make bam.list for GATK
ls  $trimmeddata/BR*.fastq | xargs -n1 basename |sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | grep -v sai | uniq  > Samplelist.txt

while read prefix
do
echo "fixing $prefix..."
echo "cleaning..."
	java -jar $picardtools/CleanSam.jar INPUT=$bam_fixed/$prefix.sortrg.bam.fixed.bam OUTPUT=$bam_fixed/$prefix.clean.bam 2> $log/$prefix.cleansam.log
	echo "sorting..."
	java -jar $picardtools/SortSam.jar INPUT=$bam_fixed/$prefix.clean.bam OUTPUT=$bam_fixed/$prefix.sort.bam SORT_ORDER=coordinate 2> $log/$prefix.sortsam.log
	echo "adding read groups..."
	java -jar $picardtools/AddOrReplaceReadGroups.jar I=$bam_fixed/$prefix.sort.bam O= $bam_fixed_rg/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True 2> $log/$prefix.addRG.log
	rm $bam_fixed/$prefix.clean.bam
	rm $bam_fixed/$prefix.sort.bam
   
done < Samplelist.txt








