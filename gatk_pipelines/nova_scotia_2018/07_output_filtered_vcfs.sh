#!/bin/bash

#SBATCH --mem=40G
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J vcf_output-%j
#SBATCH -o tmp/vcf_output-%j.out

# output snps tables for each groupomosome, for each vcf

vcf="data/vcf/dpse_iso_seq_hard_filtered_snps.vcf"

# split by chromosome
echo "Converting $vcf to filtered VCFs..."

vcftools --vcf $vcf \
 --chr 2 \
 --recode --out data/vcf/final/dpse_iso_seq_CHR2.vcf

vcftools --vcf $vcf \
 --chr 3 \
 --recode --out data/vcf/final/dpse_iso_seq_CHR3.vcf

vcftools --vcf $vcf \
--chr XL_group3a \
--chr XL_group3b \
--chr XL_group1e \
--chr XL_group1a \
--recode --out data/vcf/final/dpse_iso_seq_CHRXL.vcf

vcftools --vcf $vcf \
 --chr XR_group3a \
 --chr XR_group8 \
 --chr XR_group6 \
 --chr XR_group5 \
 --recode --out data/vcf/final/dpse_iso_seq_CHRXR.vcf

vcftools --vcf $vcf \
 --chr 4_group1 \
 --chr 4_group2 \
 --chr 4_group3 \
 --chr 4_group4 \
 --chr 4_group5 \
 --recode --out data/vcf/final/dpse_iso_seq_CHR4.vcf
