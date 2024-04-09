#!/bin/bash
#BSUB -J 5_test_differential
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/5_test_differential_%J.stdout
#BSUB -eo logs/5_test_differential_%J.stderr

# initialize conda environment
source ~/.bashrc
conda activate hypertribe
cd $LS_SUBCWD

# R: format output
Rscript 5_test_differential.R
