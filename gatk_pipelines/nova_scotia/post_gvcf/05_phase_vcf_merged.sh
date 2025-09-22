#!/bin/bash

#jan 2016 kms
#compute population genetic parameters for each population in a merged VCF
#populations are defined by an input file formatted like <id> <pop> (no header)

ref='/home/samuk/gbs2015/ref/revisedAssemblyMasked.fa'
tmp="tmp"
log="log"
by_chromo="by_chromo"
merged_vcf="subset_vcfs/whtstbk_bial_maf_2014.vcf.gz"
bcftools="~bin/bcftools"
shapeit_maps="utils/shapeit_maps"
phased="phased_vcfs"

if [ ! -d "$by_chromo" ]; then
	mkdir $by_chromo
fi

if [ ! -d "$phased" ]; then
	mkdir $phased
fi

# # create the bial file for shape it
# if [ -f $merged_vcf ];
	# then
		# echo "Bial file exists, skipping..." 
	# else
		# echo "Creating bial file..." 
		# ~/bin/bcftools view --min-alleles 2 --max-alleles 2 --types snps --output-type v $original_vcf > $merged_vcf
# fi

#chromosome set names (for subsetting)
chr_set="chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI chrXVII chrXVIII chrXIX chrXX chrXXI" 

# split by chromosome
echo "Splitting $merged_vcf by chromosome..."
for lg in $chr_set
do
	echo "$lg..."
	if [ -f $by_chromo/$lg.merged.vcf ];
	then
		echo "File exists, skipping..." 
	else
		vcftools --recode --gzvcf $merged_vcf --chr $lg --max-missing 0.5 -c > $by_chromo/$lg.merged.vcf
	fi
	
done

# phase each split vcf

for lg in $chr_set
do
	echo "Phasing $lg for $pop..."
	~/bin/shapeit --input-vcf $by_chromo/$lg.merged.vcf \
			-O $phased/$lg.merged.phased \
			-M $shapeit_maps/"$lg".genmap \
			--output-log $log/$lg.$merged.phased.log \
			--thread 8
			
	~/bin/shapeit -convert \
		--input-haps $phased/$lg.merged.phased \
		--output-log $log/$lg.mergedphased_convert.log \
		--output-vcf $phased/$lg.merged.phased.vcf
done









