#!/bin/bash

# run faststructure for K1-8

for k in `seq 8`
do
	echo "K=$k"
	python fastStructure/structure.py -K $k --input=data/test --output=data/output_plink/whtstbk_structure --seed=1859	
done
