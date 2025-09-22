#!/bin/bash

#SBATCH --mem=200000
#SBATCH --cpus-per-task=20
#SBATCH -p noor

vcftools --recode \
	--max-missing 0.4 \
	--gzvcf ns_invar_vcf.gz \
	-c | gzip -c > ns_invar_max.missing_0.4.vcf.gz

zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrI' | gzip -c > by_chromo/chrI.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrII' | gzip -c > by_chromo/chrII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrIII' | gzip -c > by_chromo/chrIII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrIV' | gzip -c > by_chromo/chrIV.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrV' | gzip -c > by_chromo/chrV.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrVI' | gzip -c > by_chromo/chrVI.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrVII' | gzip -c > by_chromo/chrVII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrVIII' | gzip -c > by_chromo/chrVIII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrIX' | gzip -c > by_chromo/chrIX.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrX' | gzip -c > by_chromo/chrX.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXI' | gzip -c > by_chromo/chrXI.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXII' | gzip -c > by_chromo/chrXII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXIII' | gzip -c > by_chromo/chrXIII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXIV' | gzip -c > by_chromo/chrXIV.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXV' | gzip -c > by_chromo/chrXV.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXVI' | gzip -c > by_chromo/chrXVI.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXVII' | gzip -c > by_chromo/chrXVII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXVIII' | gzip -c > by_chromo/chrXVIII.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXIX' | gzip -c > by_chromo/chrXIX.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXX' | gzip -c > by_chromo/chrXX.vcf.gz
zcat ns_invar_max.missing_0.4.vcf.gz | grep -w '^#\|^#CHROM\|^chrXXI' | gzip -c > by_chromo/chrXXI.vcf.gz
