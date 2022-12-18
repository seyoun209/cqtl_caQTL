#!/bin/bash
#SBATCH -J atac_qc
#SBATCH -p general
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --mem=15g
#SBATCH -t 1-00:00:00
#SBATCH -e ./tmp/atacqc_e.%j
#SBATCH -o ./tmp/atacqc_o.%j

ml r/4.2.1

Rscript ./qc_atac.R
