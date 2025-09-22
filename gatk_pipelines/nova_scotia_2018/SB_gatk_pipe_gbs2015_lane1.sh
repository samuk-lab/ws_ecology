#!/bin/bash

#USAGE: ./SB_gatk_pipe.sh <plate identifier> 
#Sept 2013 GLO
#Pipeline for using GATK and BWA to call snps
plate='whtstbk_gbs_2015_L1'
Gpath='/home/samuk/bin'
Hpath='/home/samuk/gbs2015'
rawdata='/home/samuk/gbs2015/lane1/raw'
bwa='/home/samuk/bin/bwa/bwa-0.7.10'
stampy='/home/samuk/bin/stampy-1.0.23/stampy.py'
demultiplex='/home/samuk/bin/GBS_fastq_Demultiplexer_v8.GO.pl'
barcodes='/home/samuk/gbs2015/whtsbk_2015_barcodes_lib1.txt'
ref='/home/samuk/review/ref/GA.broad1_75_dna_rm.fa'
stampyref='/home/samuk/review/ref/GA.broad1_75_dna_rm.fa.stampyref'
cleandata="/home/samuk/gbs2015/lane1/combined"
cleandata_chopped="/home/samuk/gbs2015/lane1/combined_chopped"
trimmeddata="/home/samuk/gbs2015/lane1/$plate.trimmed_data_paired"
unpaired="/home/samuk/gbs2015/lane1/$plate.trimmed_data_unpaired"
sam="/home/samuk/gbs2015/lane1/$plate.sam"
bam="/home/samuk/gbs2015/lane1/$plate.bam"
log="/home/samuk/gbs2015/lane1/$plate.log"
gvcf="/home/samuk/gbs2015/lane1/gvcf"
tabix='/home/samuk/bin/tabix-0.2.6'
vcftools='/home/samuk/bin/vcftools_0.1.11/bin'
picardtools='/home/samuk/bin/picard-tools-1.96'
trim="/home/samuk/bin/Trimmomatic-0.32"
project='whtstbk_gbs_2015_L1'
#Prerequisites:
#Index the reference for GATK and BWA
# 
# #/home/samuk/bin/stampy-1.0.23/stampy.py -G GA.broad1_75_dna_rm.fa.stampyref ./GA.broad1_75_dna_rm.fa
# #/home/samuk/bin/stampy-1.0.23/stampy.py -g GA.broad1_75_dna_rm.fa.stampyref -H GA.broad1_75_dna_rm.fa.stampyref
# 
# #Change the variables to fit your file structure and ensure the folders exist.
# 
# #set java opts
export _JAVA_OPTIONS="-Xmx12g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=4G"
# 
# echo "Starting GBS processing pipeline on $plate."
# #Make directories if they don't exist
# if [ ! -d "$sam" ]; then
# 	mkdir $sam
# fi
# if [ ! -d "$bam" ]; then
# 	mkdir $bam
# fi
# if [ ! -d "$log" ]; then
# 	mkdir $log
# fi
# if [ ! -d "$gvcf" ]; then
#         mkdir $gvcf
# fi
# if [ ! -d "$cleandata" ]; then
#         mkdir $cleandata
# fi
# if [ ! -d "$trimmeddata" ]; then
#         mkdir $trimmeddata
# fi
# if [ ! -d "$unpaired" ]; then
#         mkdir $unpaired
# fi
# 
# #Demultiplex
# #perl $demultiplex $barcodes $rawdata/"$plate"_R1.fastq $rawdata/"$plate"_R2.fastq $cleandata/
# 
# 
# ls $cleandata | grep -v "nobar" | grep -v "POND" | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > $Hpath/Samplelist.$plate.txt
# 
# # chop reads with fastx_trimmer
# # only necessary because of messed up second read (first 3 bases on 2nd are NNN and ### quality, trimmomatic sliding window fails to run at all)
# while read prefix
# do
# echo "Chopping ${prefix}_R2..."
# fastx_trimmer -f 4 -Q33 -i $cleandata/"$prefix"_R2.fastq -o $cleandata_chopped/"$prefix"_R2.fastq
# done < Samplelist.$plate.txt
# 
# # #Trim the data using Trimmomatic. Removes bad reads and illumina adapter contamination.
#  
# while read prefix
# do
# echo "Trimming ${prefix}..."
# java -jar $trim/trimmomatic-0.32.jar PE -phred33 $cleandata/"$prefix"_R1.fastq $cleandata_chopped/"$prefix"_R2.fastq $trimmeddata/"$prefix"_R1.fastq $unpaired/"$prefix"_unR1.fastq $trimmeddata/"$prefix"_R2.fastq $unpaired/"$prefix"_unR2.fastq ILLUMINACLIP:$trim/adapters/TruSeq3-PE.fa:2:30:10:8:T SLIDINGWINDOW:4:15 MINLEN:36
# done < Samplelist.$plate.txt

# ls  $trimmeddata | grep -v "nobar" | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > Samplelist.$plate.txt

# align_samplelist=Samplelist.$plate.txt

# # allows reruns on failed alignments, just keep running the script
# rerun_lines=$(wc -l rerun_samples_lane1.txt | awk '{print $1}')

# if [ $rerun_lines -gt 0 ]; then
	# echo "Rerunning ${rerun_lines} samples..."
	# align_samplelist="rerun_samples_lane1.txt"
# fi

# ###Align using BWA. Turn from sam to bam. Sort by coordinate and add read group data.
# echo "Aligning files in ${align_samplelist}..."
# while read prefix
# do
	# $bwa/bwa aln -t 4 $ref $trimmeddata/"$prefix"_R1.fastq 1> $trimmeddata/"$prefix"_R1.sai
	# $bwa/bwa aln -t 4 $ref $trimmeddata/"$prefix"_R2.fastq 1> $trimmeddata/"$prefix"_R2.sai
	# $bwa/bwa sampe $ref $trimmeddata/"$prefix"_R1.sai $trimmeddata/"$prefix"_R2.sai $trimmeddata/"$prefix"_R1.fastq $trimmeddata/"$prefix"_R2.fastq 1> $sam/$prefix.sam 2> $log/$prefix.bwasampe.log
	# samtools view -Sb $sam/$prefix.sam > $bam/$prefix.bam	
	# $stampy -g $stampyref -h $stampyref -t4 --bamkeepgoodreads -M $bam/$prefix.bam -o $bam/$prefix.stampy.bam 2> $log/$prefix.stampy.log
	# java -jar $picardtools/CleanSam.jar INPUT=$bam/$prefix.stampy.bam OUTPUT=$bam/$prefix.clean.bam 2> $log/$prefix.cleansam.log
	# java -jar $picardtools/SortSam.jar INPUT=$bam/$prefix.clean.bam OUTPUT=$bam/$prefix.sort.bam SORT_ORDER=coordinate 2> $log/$prefix.sortsam.log
	# java -jar $picardtools/AddOrReplaceReadGroups.jar I=$bam/$prefix.sort.bam O= $bam/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True 2> $log/$prefix.addRG.log
	# rm $trimmeddata/"$prefix"_R1.sai
	# rm $trimmeddata/"$prefix"_R2.sai
	# rm $sam/$prefix.sam
	# rm $bam/$prefix.bam
	# rm $bam/$prefix.stampy.bam
	# rm $bam/$prefix.clean.bam
	# rm $bam/$prefix.sort.bam
# done < $align_samplelist


# ### CHECK IF ALIGNMENT WAS SUCCESSFUL

# # make a list of bam files
# ls lane1/whtstbk_gbs_2015_L1.bam | sed s/.sortrg.bai//g | sed s/.sortrg.bam//g | uniq  > bam_list_lane1.txt

# # compare to the input file list
# grep -v -f bam_list_lane1.txt  Samplelist.whtstbk_gbs_2015_L1.txt > rerun_samples_lane1.txt

# if [ $rerun_lines -gt 0 ]; then
	# echo "The following files failed to align:"
	# cat rerun_samples_lane1.txt
# fi

#Make bam.list for GATK
ls -d $bam/*.* | grep sortrg.bam  > $Hpath/bamlist.$plate.list

#identify local indels

java -Xmx12g -jar $Gpath/GenomeAnalysisTK_nightly.jar \
   -T RealignerTargetCreator \
   -R $ref \
   -I $Hpath/bamlist.$plate.list \
   -nt 2 \
   -log $log/$plate.RealignerTargetCreator.log \
   -o $Hpath/$plate.realign.intervals
   


#Realign around local indels
while read prefix
do
java -Xmx12g -jar $Gpath/GenomeAnalysisTK_nightly.jar \
	-T IndelRealigner \
	-R $ref \
	-I $bam/$prefix.sortrg.bam \
	-targetIntervals $Hpath/$plate.realign.intervals \
	-o $bam/$prefix.realign.bam \
	-log $log/$plate.$prefix.IndelRealigner.log 

#Call GATK HaplotypeCaller
java -Xmx12g -jar $Gpath/GenomeAnalysisTK_nightly.jar \
	-nct 1 \
	-l INFO \
	-R $ref \
	-log $log/$plate.$prefix.HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I  $bam/$prefix.realign.bam \
	--emitRefConfidence GVCF \
	--max_alternate_alleles 2 \
	-variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-o $gvcf/${prefix}.GATK.gvcf.vcf
done < $Hpath/Samplelist.${plate}.txt

#Make bam.list for GATK
ls -d $gvcf/*.GATK.gvcf.vcf > $Hpath/gvcf.$plate.list

#Genotype all gvcf together into one vcf file
java -Xmx12g -jar $Gpath/GenomeAnalysisTK_nightly.jar \
	-nt 2 \
	-l INFO \
	-R $ref \
	-log $log/$plate.GenotypeGVCFs.log \
	-T GenotypeGVCFs \
	-V gvcf.$plate.list\
	-o $Hpath/$plate.GATK.total.vcf \
	--includeNonVariantSites \
exit
