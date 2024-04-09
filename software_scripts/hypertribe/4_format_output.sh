#!/bin/bash
#BSUB -J 4_format_output
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/4_format_output_%J.stdout
#BSUB -eo logs/4_format_output_%J.stderr

# initialize conda environment
source ~/.bashrc
conda activate hypertribe
cd $LS_SUBCWD

# R: format output
Rscript 4_format_output.R
