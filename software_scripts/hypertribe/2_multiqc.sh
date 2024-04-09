#!/bin/bash
#BSUB -J 2_multiqc
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/2_multiqc_%J.stdout
#BSUB -eo logs/2_multiqc_%J.stderr

# initialize conda environment
source ~/.bashrc
conda activate hypertribe
cd $LS_SUBCWD

# input directories
fastqc_folder=../../output_data/hypertribe/1_star_align

# output directories
multiqc_folder=../../output_data/hypertribe/1_star_align

# multiqc: generate files
multiqc \
${fastqc_folder} \
-o ${multiqc_folder}
