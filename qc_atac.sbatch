#!/bin/bash
#SBATCH -J atac_qc
#SBATCH -p general
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=400g
#SBATCH -t 12:00:00
#SBATCH -e ./tmp/splitebam_v3_e.%j
#SBATCH -o ./tmp/splitebam_v3_o.%j

ml r/4.2.1

#Rscript qc_atac_v2.R
Rscript deseq_atac.R $1
