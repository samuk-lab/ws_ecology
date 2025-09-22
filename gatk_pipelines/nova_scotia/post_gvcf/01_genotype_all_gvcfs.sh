#!/bin/bash

#jan 2016 kms
#jointly genotype combined gvcf files
Gpath='/home/samuk/bin/GATK'
ref='/home/samuk/gbs2015/ref/revisedAssemblyMasked.fa'
tmp="tmp"
log="log"
out="gvcf"

if [ ! -d "$tmp" ]; then
	mkdir $tmp
fi

if [ ! -d "$out" ]; then
	mkdir $out
fi

if [ ! -d "$log" ]; then
	mkdir $log
fi

set java opts
export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=8G"

#Make list of the combined gvcfs for GATK
ls /home/samuk/gbs2012/gvcfpop/*.vcf > ns.gvcf.list 
ls /home/samuk/gbs2015/gbs_combined/gvcfpop/*.vcf >> ns.gvcf.list 
ls /home/samuk/north_sea/data/gvcf/*.vcf > dk.gvcf.list 
ls /home/samuk/pac_marines/data/gvcf/*.vcf > lc.gvcf.list 
ls /home/samuk/wheatlandii/data/gvcf/*.vcf > wheat.gvcf.list 

# #Genotype all gvcf for nova scotia pops
# java -Xmx20g -jar $Gpath/GenomeAnalysisTK.jar \
	# -nt 12 \
	# -l INFO \
	# -R $ref \
	# -log $log/GenotypeGVCFs.log \
	# -T GenotypeGVCFs \
	# -V ns.gvcf.list\
	# --max_alternate_alleles 4 \
	# -o $out/whtstbk_nova_scotia.vcf
	
# #Genotype all gvcf for north sea pops
# java -Xmx20g -jar $Gpath/GenomeAnalysisTK.jar \
	# -nt 12 \
	# -l INFO \
	# -R $ref \
	# -log $log/GenotypeGVCFs.log \
	# -T GenotypeGVCFs \
	# -V dk.gvcf.list \
	# --max_alternate_alleles 4 \
	# -o $out/whtstbk_denmark.vcf
	
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
java -Xmx20g -jar $Gpath/GenomeAnalysisTK.jar \
	-nt 10 \
	-l INFO \
	-R $ref \
	-log $log/GenotypeGVCFs.log \
	-T GenotypeGVCFs \
	-V wheat.gvcf.list \
	--max_alternate_alleles 2 \
	-o $out/whtstbk_wheat.vcf
	

rm *.list
