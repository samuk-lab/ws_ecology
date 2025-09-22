#!/bin/bash

input_file='whtstbk_lcr_dk_join_pruned.treemix.gz'
output_slug='lcr_dk_outgroups'

rm -r $output_slug

if [ ! -d "$output_slug" ]; then
	mkdir $output_slug
fi

# build a ML tree

kvalues="1 5 25"
mvalues="0 1 2 3 4 5 6 7 8 9 10"
roots="LC_NA"

for k in $kvalues
do

	for m in $mvalues
	do
	
		for root in $roots
		do
			treemix -i $input_file -k $k -root $root -m $m -o $output_slug/k"$k"_m"$m"_root"$root"
		done
	
	done

done


# build a ML tree w/ outgroup
#treemix -i input file.gz -root DK -o out stem

# build a ML tree + set SNP window size (for controlling for LD) and 
#treemix -i input file.gz -k 1000 -o out stem

# build a ML tree allowing for 2 migration events
#treemix -i input file.gz -m 2 -o out stem

# 4.5 Input a previously generated tree/graph (-g)
# There are two ways to input a previously generated tree/graph. The most simple is to input from
# TreeMix format using the -g flag, which take a file of vertices and a file of edges as input. For
# example:
# >treemix -i input file.gz -m 2 -g out stem.vertices.gz out stem.edges.gz -o out stem2
# 4.6 Generate a bootstrap replicate (-bootstrap)
# For judging the confidence in a given tree topology, it is often of interest to generate a bootstrap
# replicate. Bootstrapping is done over blocks of contiguous SNPs. The following example shows
# how to generate a single bootstrap replicate by resampling blocks of 500 SNPs:
# >treemix -i input file.gz -bootstrap -k 500 -o replicate1
