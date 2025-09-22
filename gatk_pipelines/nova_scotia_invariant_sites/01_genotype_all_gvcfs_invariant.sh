#!/bin/bash

#SBATCH --mem=200000
#SBATCH --cpus-per-task=20
#SBATCH -p noor

#jan 2016 kms
#jointly genotype combined gvcf files
bin_folder="/dscrhome/kms173/bin"
GATK="$bin_folder/GenomeAnalysisTK.jar"
java_1_8="$bin_folder/jdk1.8.0_101/bin/java"
ref='ref/revisedAssemblyMasked.fa'

set java opts
export _JAVA_OPTIONS="-Xmx180g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=20G"

#Make list of the combined gvcfs for GATK
#ls /home/samuk/gbs2012/gvcfpop/*.vcf > ns.gvcf.list 
#ls /home/samuk/gbs2015/gbs_combined/gvcfpop/*.vcf >> ns.gvcf.list 
#ls /home/samuk/north_sea/data/gvcf/*.vcf > dk.gvcf.list 
#ls /home/samuk/pac_marines/data/gvcf/*.vcf > lc.gvcf.list 
#ls /home/samuk/wheatlandii/data/gvcf/*.vcf > wheat.gvcf.list 

#Make list of the combined gvcfs for GATK
ls -d gvcf/*.vcf > ns.gvcf.list

#Genotype all gvcf for nova scotia pops
 $java_1_8 -Xmx100g -Djava.io.tmpdir=./tmp -jar $GATK \
   -nt 12 \
	 -l INFO \
	 -R $ref \
	 -log log/GenotypeGVCFs.log \
	 -T GenotypeGVCFs \
	 -V ns.gvcf.list \
 	 --includeNonVariantSites \
	 --max_alternate_alleles 4 \
	--disable_auto_index_creation_and_locking_when_reading_rods \
	 -o /dev/stdout | gzip -c > ns_invar_vcf.gz 
   	
# #Genotype all gvcf for north sea pops
#java -Xmx20g -jar $Gpath/GenomeAnalysisTK.jar \
#	 -nt 12 \
#	 -l INFO \
#	 -R $ref \
#	 -log $log/GenotypeGVCFs.log \
#  -T GenotypeGVCFs \
#	 -V dk.gvcf.list \
#	 --max_alternate_alleles 4 \
#   --includeNonVariantSites \
#	 -o $out/whtstbk_denmark_invariant.vcf
	
# Genotype all gvcf for lcr
# java -Xmx20g -jar $Gpath/GenomeAnalysisTK.jar \
	# -nt 12 \
	# -l INFO \
	# -R $ref \
	# -log $log/GenotypeGVCFs.log \
	# -T GenotypeGVCFs \
	# -V lc.gvcf.list \
	# --max_alternate_alleles 4 \
	# -o $out/whtstbk_lcr.vcf
	
#Genotype wheatlandii
#java -Xmx20g -jar $Gpath/GenomeAnalysisTK.jar \
#	-nt 10 \
#	-l INFO \
#	-R $ref \
#	-log $log/GenotypeGVCFs.log \
#	-T GenotypeGVCFs \
#	-V wheat.gvcf.list \
#	--max_alternate_alleles 2 \
#	-o $out/whtstbk_wheat.vcf
	

rm *.list
