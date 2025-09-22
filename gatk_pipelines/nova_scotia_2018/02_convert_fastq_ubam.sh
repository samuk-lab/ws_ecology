#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH -J convert
#SBATCH -p noor
#SBATCH -o tmp/convert.out

#Pipeline for using GATK and BWA to call snps
# KMS March 2017

# use samplelist as input
samplelist=$1

# folders
proj_home="/work/kms173/wht_stbk"
ref="$proj_home/data/ref/revisedAssemblyMasked.fa"
demulti="$proj_home/data/trimmed"
sam="$proj_home/data/sam"
bam="$proj_home/data/bam"
log="$proj_home/log"
gvcf="$proj_home/data/gvcf"

# plate, barcode and project ids
project="wht_stbk"


# tools
bin_folder="/dscrhome/kms173/bin"
java_1_8="$bin_folder/jdk1.8.0_101/bin/java"
tabix="$bin_folder/tabix-0.2.6"
vcftools="$bin_folder/vcftools_0.1.11/bin"
picardtools="$bin_folder/picard.jar"
gatk="$bin_folder/GATK"
bwa="$bin_folder/bwa"

export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"

project="iso_seq"

#Prerequisites:
#Index the reference genome, we are using the Broad 75 genome

#ensure the folders exist.

#rm $proj_home/tmp/convert.out

echo "Generating sample list..."
#ls $raw | sed 's/whtstbk_gbs_2012_brds_/2012./g' | sed 's/_.*//g' | uniq > $proj_home/tmp/Samplelist.txt

#ls $demutli | grep 2012 | sed 's/whtstbk_gbs_2012_brds_//g' | sed 's/_.*//g' | uniq > $proj_home/tmp/Samplelist_2012.txt
#ls $demulti | grep -v whtstbk_gbs_2012 | sed 's/_.*//g' | uniq > $proj_home/tmp/Samplelist_2014.txt

cat $proj_home/tmp/Samplelist_2012.txt


# clean + redo
#rm $sam/*
#rm $bam/*

# test mode

#/dscrhome/kms173/bin/bwa mem -t 1 -v 3 -a -M data/ref_genome/dpse-all-chromosome-r3.04.fasta \
#   <(gunzip -c data/demulti/S1_R1.fastq.gz) \
#   <(gunzip -c data/demulti/S1_R2.fastq.gz) > data/sam/S1.sam

#head -n 1 tmp/Samplelist.txt > tmp/Samplelist2.txt
#rm tmp/Samplelist.txt 
#mv tmp/Samplelist2.txt tmp/Samplelist.txt


### perform fastq to ubam conversion
#while read prefix
#do
#	echo "$prefix..."
#
#  ~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar FastqToSam \
#    FASTQ=$demulti/"$prefix"_R1_combined.fastq \
#    FASTQ2=$demulti/"$prefix"_R2_combined.fastq \
#    OUTPUT=$bam/"$prefix"_unmapped.bam \
#    READ_GROUP_NAME=$prefix \
#    SAMPLE_NAME=$prefix \
#    LIBRARY_NAME=$project \
#    PLATFORM=illumina \
#    SEQUENCING_CENTER=DUKEGCB \
#    RUN_DATE=2018-04-20T00:00:00-0400
#    
#  ~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar MarkIlluminaAdapters \
#	I=$bam/"$prefix"_unmapped.bam \
#	O=$bam/"$prefix"_markadapters.bam \
#	M=$proj_home/log/"$prefix"_markadapters.txt \
#	TMP_DIR=./tmp
#  
#done < tmp/Samplelist_2014.txt

while read prefix
do
	echo "$prefix..."

  ~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar FastqToSam \
    FASTQ=$demulti/whtstbk_gbs_2012_brds_"$prefix"_R1.fastq \
    FASTQ2=$demulti/whtstbk_gbs_2012_brds_"$prefix"_R2.fastq \
    OUTPUT=$bam/"$prefix".2012_unmapped.bam \
    READ_GROUP_NAME="$prefix".2012 \
    SAMPLE_NAME="$prefix".2012 \
    LIBRARY_NAME=$project \
    PLATFORM=illumina \
    SEQUENCING_CENTER=DUKEGCB \
    RUN_DATE=2018-04-20T00:00:00-0400
    
  ~/bin/jdk1.8.0_101/bin/java -Djava.io.tmpdir=./tmp -jar ~/bin/picard.jar MarkIlluminaAdapters \
	I=$bam/"$prefix".2012_unmapped.bam \
	O=$bam/"$prefix".2012_markadapters.bam \
	M=$proj_home/log/"$prefix".2012_markadapters.txt \
	TMP_DIR=./tmp
  
done < tmp/Samplelist_2012.txt
