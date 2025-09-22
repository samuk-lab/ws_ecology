#!/bin/bash

#SBATCH --mem=200000
#SBATCH --cpus-per-task=20
#SBATCH -p noor

Rscript 06_invariant_vcf_prep_02.R
