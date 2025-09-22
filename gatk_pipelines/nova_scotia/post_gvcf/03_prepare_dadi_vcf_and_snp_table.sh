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
fst="fst"

bcftools="~/bin/bcftools"

### NO OUTGROUP

# echo "Creating unpruned sex file..." 
 # vcftools --recode \
	 # --max-alleles 2 \
	 # --min-alleles 2 \
	 # --max-missing 0.8 \
	 # --not-chr chrXIX \
	 # --remove-indels \
	 # --gzvcf $vcf \
	 # --stdout > subset_vcfs/whtstbk_bial_nomaf_nosex.vcf

# gzip -c subset_vcfs/whtstbk_bial_nomaf_nosex.vcf > subset_vcfs/whtstbk_bial_nomaf_nosex.vcf.gz
# rm subset_vcfs/whtstbk_bial_nomaf_nosex.vcf

#### OUTGROUP(S)

vcf="merged/whtstbk_dk_merged.vcf.gz"

echo "Creating unpruned sex file..." 
 vcftools --recode \
	 --max-alleles 2 \
	 --min-alleles 2 \
	 --max-missing 0.8 \
	 --not-chr chrXIX \
	 --remove-indels \
	 --gzvcf $vcf \
	 --stdout > subset_vcfs/whtstbk_bial_nomaf_nosex_outgroup.vcf

gzip -c subset_vcfs/whtstbk_bial_nomaf_nosex_outgroup.vcf > subset_vcfs/whtstbk_bial_nomaf_nosex_outgroup.vcf.gz
rm subset_vcfs/whtstbk_bial_nomaf_nosex_outgroup.vcf

sh utils/prepare_pruned_plink_raw.sh subset_vcfs/whtstbk_bial_nomaf_nosex_outgroup.vcf.gz whtstbk_bial_nomaf_nosex_outgroup_pruned no_sex

