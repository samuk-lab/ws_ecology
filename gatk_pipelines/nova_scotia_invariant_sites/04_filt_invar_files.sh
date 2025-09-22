#!/bin/bash

#SBATCH --mem=200000
#SBATCH --cpus-per-task=20
#SBATCH -p noor

Rscript 07_filter_invariant_vcf.R
