#!/bin/bash

#jan 2016 kms
#compute population genetic parameters for each population in a merged VCF
#populations are defined by an input file formatted like <id> <pop> (no header)

ref='/home/samuk/gbs2015/ref/revisedAssemblyMasked.fa'
tmp="tmp"
log="log"
gvcf="gvcf"
merged="merged"
by_pop="by_pop"
stats_out="$by_pop/stats"
vcf="merged/whtstbk_nova_scotia_filtered.vcf.gz"
merged_vcf="$tmp/whtstbk_bial_maf_2014.vcf.gz"
popfile="utils/popfile_filtered.txt"
fst="fst"

if [ ! -d "$fst" ]; then
	mkdir $fst
fi

bcftools="~/bin/bcftools"

# determine the unique populations from the pop file
populations=$(awk '{print $2}' $popfile | sort | uniq)
samples=$(cat $popfile | awk '{print $1}' | paste -s -d,)

# # # create unique popfiles
# # for pop in $populations
# # do
	# # grep $pop $popfile | awk '{print $1}' > $tmp/$pop.popfile.txt
# # done

# # # create the bial file for shape it
# # if [ -f subset_vcfs/whtstbk_bial_maf_2012.vcf.gz ];
	# # then
		# # echo "Filtered vcf files exist, skipping..." 
	# # else
		# # echo "Filtering master vcf (biallelic, maf > 0.05, snps only, samples in popfile...)"
		# # samples=$(cat $popfile | grep -v 2012 | awk '{print $1}' | paste -s -d,)
		# # ~/bin/bcftools view --min-alleles 2 --max-alleles 2 --samples $samples --min-af 0.05:minor $vcf --output-type z --types snps > subset_vcfs/whtstbk_bial_maf_2014.vcf.gz
		
		# # samples=$(cat $popfile | awk '{print $1}' | paste -s -d,)
		# # ~/bin/bcftools view --min-alleles 2 --max-alleles 2 --samples $samples --min-af 0.05:minor $vcf --output-type z --types snps > subset_vcfs/whtstbk_bial_maf_all.vcf.gz
		
		# # samples=$(cat $popfile | grep 2012 | awk '{print $1}' | paste -s -d,)
		# # ~/bin/bcftools view --min-alleles 2 --max-alleles 2 --samples $samples --min-af 0.05:minor $vcf --output-type z --types snps > subset_vcfs/whtstbk_bial_maf_2012.vcf.gz
# # fi

# # sh utils/prepare_pruned_plink_raw.sh subset_vcfs/whtstbk_bial_maf_2014.vcf.gz pruned_snp_tables/whtstbk_2014_pruned no_sex
# # sh utils/prepare_pruned_plink_raw.sh subset_vcfs/whtstbk_bial_maf_2012.vcf.gz pruned_snp_tables/whtstbk_2012_pruned no_sex
# # sh utils/prepare_pruned_plink_raw.sh subset_vcfs/whtstbk_bial_maf_all.vcf.gz pruned_snp_tables/whtstbk_all_pruned no_sex
# # sh utils/prepare_pruned_plink_raw.sh subset_vcfs/whtstbk_bial_maf_all.vcf.gz pruned_  snp_tables/whtstbk_all_pruned_sex sex

# # # also create an unpruned file WITH sex chromosome
# # # (for genetic sexing of individuals...for science)

# echo "Creating unpruned sex file..." 
# vcftools --plink \
	# --max-alleles 2 \
	# --min-alleles 2 \
	# --max-missing 0.8 \
	# --gzvcf subset_vcfs/whtstbk_bial_maf_all.vcf.gz \
	# --out $tmp/temp_plink_pre_prune

# plink --file $tmp/temp_plink_pre_prune --recodeA --out analysis_ready/snp_tables/whtstbk_all_sex

# gzip -c analysis_ready/snp_tables/whtstbk_all_sex.raw > analysis_ready/snp_tables/whtstbk_all_sex.gz

# # # also create an unpruned file WITHOUT sex chromosome
# # # (for genetic sexing of individuals...for science)

# echo "Creating unpruned NO sex file..." 
# vcftools --plink \
	# --max-alleles 2 \
	# --min-alleles 2 \
	# --max-missing 0.8 \
	# --not-chr chrXIX \
	# --gzvcf subset_vcfs/whtstbk_bial_maf_all.vcf.gz \
	# --out $tmp/temp_plink_pre_prune

# plink --file $tmp/temp_plink_pre_prune --recodeA --out analysis_ready/snp_tables/whtstbk_all_no_sex

# gzip -c analysis_ready/snp_tables/whtstbk_all_no_sex.raw > analysis_ready/snp_tables/whtstbk_all_no_sex.gz

# # also create an unpruned file WITHOUT sex chromosome (2014 only)

samples=$(cat $popfile | grep -v 2012 | awk '{print $1}' | paste -s -d,)
~/bin/bcftools view --min-alleles 2 --max-alleles 2 --samples $samples --min-af 0.05:minor $vcf --output-type z --types snps > $tmp/2014_only.vcf.gz

echo "Creating unpruned NO sex file..." 
vcftools --plink \
	--max-alleles 2 \
	--min-alleles 2 \
	--max-missing 0.8 \
	--not-chr chrXIX \
	--gzvcf $tmp/2014_only.vcf.gz \
	--out $tmp/temp_plink_pre_prune

plink --file $tmp/temp_plink_pre_prune --recodeA --out analysis_ready/snp_tables/whtstbk_2014_no_sex

gzip -c analysis_ready/snp_tables/whtstbk_2014_no_sex.raw > analysis_ready/snp_tables/whtstbk_2014_no_sex.gz

# # also create an unpruned file WITHOUT sex chromosome (2012 only)

samples=$(cat $popfile | grep 2012 | awk '{print $1}' | paste -s -d,)
~/bin/bcftools view --min-alleles 2 --max-alleles 2 --samples $samples --min-af 0.05:minor $vcf --output-type z --types snps > $tmp/2012_only.vcf.gz

echo "Creating unpruned NO sex file..." 
vcftools --plink \
	--max-alleles 2 \
	--min-alleles 2 \
	--max-missing 0.8 \
	--not-chr chrXIX \
	--gzvcf $tmp/2012_only.vcf.gz \
	--out $tmp/temp_plink_pre_prune

plink --file $tmp/temp_plink_pre_prune --recodeA --out analysis_ready/snp_tables/whtstbk_2012_no_sex

gzip -c analysis_ready/snp_tables/whtstbk_2012_no_sex.raw > analysis_ready/snp_tables/whtstbk_2012_no_sex.gz
