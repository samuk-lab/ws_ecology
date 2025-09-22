#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J copy_klk

cp /dscrhome/kms173/noor2/kms173/iso_seq/data/demulti/* /datacommons/noor/klk37/AFC_MC_Genomes/raw_fasta
