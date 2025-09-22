#!/bin/bash
# reformat a filtered vcf for structure, hierfstat, etc.
# requires PGDspider
# jan 2016 kieran samuk

# initials
pgd_spider='/home/samuk/whtstbk_pop_gen/vcf_processing/pgd/PGDSpider_2.0.9.1/PGDSpider2-cli.jar'
spid_files='/home/samuk/whtstbk_pop_gen/vcf_processing/pgd/spid_files'
tmp="tmp"
filtered_vcf="analysis_ready/pruned_snp_tables/whtstbk_2014_pruned.vcf"
reformat="reformat"
out_name="reformat/whtstbk_2014_pruned"

if [ ! -d "$reformat" ]; then
	mkdir $reformat
fi

# java opts
export _JAVA_OPTIONS="-Xmx50g"
export _JAVA_OPTIONS="-XX:MaxDirectMemorySize=8G"

### PGD Spider reformatting
### !!spid files are made using the PGDspider gui!!

# convert filtered vcf to faststructure format
java -Xmx10000m -Xms10000m -jar $pgd_spider \
-inputfile $filtered_vcf \
-inputformat VCF \
-outputfile $out_name.str \
-outputformat STRUCTURE \
-spid $spid_files/vcf_to_faststructure.spid

cp $out_name.str ~/vcf_processing/analysis_ready

cp $out_name.str ~/whtstbk_pop_gen/analysis_faststructure