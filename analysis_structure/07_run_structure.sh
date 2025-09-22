#!/bin/bash

# run admixture


for k in `seq 6`
do
	echo "K=$k"
	admixture --cv whtstbk_2020_nochr.ped $k | tee log${K}.out
done
