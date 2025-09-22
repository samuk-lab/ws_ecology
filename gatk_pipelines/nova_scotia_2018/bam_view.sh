#!/bin/bash

coords=$1
folder=$2

prefixes="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"
sample_ids="A56 MV2-25 M14 A60 M17 A19 Flag14 VY-F16 A12 A47 A14 M27 A48 A49 A57 A30 M20 M15 A24 M13 M6 A6"
arr=($sample_ids)

#echo $sample_ids

bamdir=$(echo $1 | sed 's/\:/_/g')

mkdir data/bam_view/"$folder"
mkdir data/bam_view/"$folder"/"$bamdir"


for prefix in ${prefixes[@]}
do
    echo "$coords S$prefix ${arr[$prefix-1]}" > tmp/"$prefix"_"$coords"_tmp.txt
    samtools tview -p $coords data/bam/S"$prefix"_merged.bam data/ref_genome/dpse-all-chromosome-r3.04.fasta > tmp/"$prefix"_"$coords"_tmp2.txt
    cat tmp/"$prefix"_"$coords"_tmp.txt tmp/"$prefix"_"$coords"_tmp2.txt > data/bam_view/"$folder"/"$bamdir"/"$bamdir"_"S$prefix"_"${arr[$prefix-1]}"_bamview.txt
    rm tmp/"$prefix"_"$coords"_tmp.txt
    rm tmp/"$prefix"_"$coords"_tmp2.txt  
done

# bamdir=$(echo $1 | sed 's/\:/_/g')
# 
# mkdir data/bam_view/"$bamdir"
# 
# foo () {
#     local prefix=$1
#     echo "$coords S$prefix ${arr[$prefix-1]}" > tmp/"$prefix"_"$coords"_tmp.txt
#     samtools tview -p $coords data/bam/S"$prefix"_merged.bam data/ref_genome/dpse-all-chromosome-r3.04.fasta > tmp/"$prefix"_"$coords"_tmp2.txt
#     cat tmp/"$prefix"_"$coords"_tmp.txt tmp/"$prefix"_"$coords"_tmp2.txt > data/bam_view/"$bamdir"/"$bamdir"_"S$prefix"_"${arr[$prefix-1]}"_bamview.txt
#     rm tmp/"$prefix"_"$coords"_tmp.txt 
#     rm tmp/"$prefix"_"$coords"_tmp2.txt
# }
# 
# for prefix in ${prefixes[@]}; do foo "$prefix" & done
