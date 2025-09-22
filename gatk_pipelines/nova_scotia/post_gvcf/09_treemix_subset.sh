#!/bin/bash

#jan 2016 kms
#jointly genotype combined gvcf files

lcr_vcf="merged/whtstbk_lcr_filtered.vcf.gz"
vcf="merged/whtstbk_merged_multi.vcf.gz"

# get site list the lcr vcf

echo "creating site list..."
#zcat $lcr_vcf | grep -v "#" | cut -f1,2 | sed 's/ /\t/g' > tmp/lcr_sites.txt

echo "subsetting master vcf..."
#vcftools --recode --positions tmp/lcr_sites.txt --gzvcf $vcf --stdout | gzip > treemix_snps/whtstbk_lcr_dk_join.vcf.gz

echo "pruning for ld..."
# prune for ld
#sh utils/prepare_pruned_plink_raw.sh treemix_snps/whtstbk_lcr_dk_join.vcf.gz whtstbk_lcr_dk_join_pruned no_sex treemix_snps

Rscript utils/format_raw_to_treemix.R analysis_ready/whtstbk_lcr_dk_join_pruned.raw treemix_snps
