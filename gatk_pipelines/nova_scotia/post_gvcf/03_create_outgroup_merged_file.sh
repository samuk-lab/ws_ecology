#!/bin/bash

#jan 2016 kms
#jointly genotype combined gvcf files
Gpath='/home/samuk/bin/GATK'
ref='/home/samuk/gbs2015/ref/revisedAssemblyMasked.fa'
tmp="tmp"
log="log"
gvcf="gvcf"
merged="merged"

if [ ! -d "merged" ]; then
	mkdir $merged
fi

set java opts
export _JAVA_OPTIONS="-Xmx20g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=8G"

ls $gvcf/*.vcf > vcf.list

while read vcf 
do

outfile_name=`echo $vcf | sed 's/\.vcf/_filtered.vcf/g' | sed 's/gvcf/merged/g'`

echo "Applying filters for $outfile_name..."

#Max Missing filter

if [ -f $outfile_name.gz ];
	then
		echo "$outfile_name.gz exists, skipping..." 
	else
		echo "Filtering $vcf..."
		vcftools --recode \
		--max-missing 0.66 \
		--vcf $vcf \
		-c > $outfile_name
fi
	

if [ -f $outfile_name.gz ];
	then
		echo "$outfile_name.gz exists, skipping..." 
	else
		echo "Compressing and indexing $outfile_name..."
		bgzip $outfile_name
		tabix -p vcf $outfile_name.gz
fi

done < vcf.list

if [ -f $merged/whtstbk_dk_merged.vcf ];
	then
		echo "$merged/whtstbk_dk_merged.vcf exists, skipping..." 
	else
		echo "Merging filtered files..."
		~/bin/bcftools merge merged/*_filt*.vcf.gz --output-type v --merge all > merged/whtstbk_dk_merged.vcf
fi

echo "Filtering merged file..."
~/bin/bcftools view \
	--min-alleles 2 \
	--max-alleles 2 \
	--types snps \
	$merged/whtstbk_dk_merged.vcf \
	--output-type v > $merged/whtstbk_dk_merged_bial.vcf

echo "Final max missing filter..."
#final max missing filter + drop INFO
vcftools --recode --max-missing 0.75 --vcf $merged/whtstbk_dk_merged_bial.vcf -c > $merged/whtstbk_dk_merged.vcf

cd merged

echo "Final compressing and indexing..."
bgzip -c whtstbk_dk_merged.vcf > whtstbk_dk_merged.vcf.gz
tabix -p vcf whtstbk_dk_merged.vcf.gz

#rm whtstbk_dk_merged.vcf
#rm whtstbk_dk_merged_bial.vcf

# final max missing filter + drop INFO
#vcftools --recode --max-missing 0.975 --vcf whtstbk_merged.vcf -c > whtstbk_final_strict.vcf


#bgzip -c whtstbk_dk_merged_bial.vcf > whtstbk_dk_merged.vcf.gz
#tabix -p vcf whtstbk_dk_merged.vcf.gz
#rm $merged/whtstbk_merged.vcf.gz


