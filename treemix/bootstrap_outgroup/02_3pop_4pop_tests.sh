#!/bin/bash

input_file='whtstbk_outgroup_2014_pruned.treemix.gz'
output_slug='wht_cmn_outgroups'



fourpop -i $input_file -k 10 > fourpop_k10.txt
threepop -i $input_file -k 10 > threepop_k10.txt

fourpop -i $input_file -k 5 > fourpop_k5.txt
threepop -i $input_file -k 5 > threepop_k5.txt