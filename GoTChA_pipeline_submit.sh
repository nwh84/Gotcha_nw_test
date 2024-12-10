#!/bin/bash
#BSUB -J gotcha_data
#BSUB -P acc_MDS 
#BSUB -q premium
#BSUB -n 10
#BSUB -W 24:00
#BSUB -R rusage[mem=15GB]
#BSUB -R span[hosts=1] 
#BSUB -oo %J.stdout
#BSUB -eo %J.stderr
#BSUB -R "select[osmajor!=CENT7]"
#BSUB -L /bin/bash


module load R/4.2.0
module load nlopt

Rscript GoTChA_test_full_MutCall.R

