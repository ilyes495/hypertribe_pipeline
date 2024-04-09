#!/bin/bash
#BSUB -J picard_dictionary_genome
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/picard_dictionary_genome_%J.stdout
#BSUB -eo logs/picard_dictionary_genome_%J.stderr

# initialize conda environment
source ~/.bashrc
conda activate genome
cd $LS_SUBCWD

# input directories
genome_file=hg38_RBP_ADAR.fa

# output directories
genome_dict_file=hg38_RBP_ADAR.dict

# picard: create sequence dictionary of genome
echo "Run picard CreateSequenceDictionary on genome"
picard CreateSequenceDictionary \
R=${genome_file} \
O=${genome_dict_file}
