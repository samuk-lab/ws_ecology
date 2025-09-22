#!/bin/bash

# run faststructure for K1-8

for k in `seq 8`
do
	echo "K=$k"
	python fastStructure/distruct.py -K $k --input=data/output/whtstbk_structure --output=fig/distruct_${k}.svg
done


