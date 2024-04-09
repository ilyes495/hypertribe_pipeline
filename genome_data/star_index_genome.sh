#!/bin/bash
#BSUB -J star_index_genome
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 16
#BSUB -C 0
#BSUB -o logs/star_index_genome_%J.stdout
#BSUB -eo logs/star_index_genome_%J.stderr

# initialize conda environment
source ~/.bashrc
cd $LS_SUBCWD
conda activate genome

# input directories
genome_file=hg38_RBP_ADAR.fa
gtf_file=hg38_RBP_ADAR.gtf

# output directories
index_folder=./STAR_INDEX_OUTPUT

# make directories
mkdir -p ${index_folder}

# star index transcriptome
echo "STAR genomeGenerate...."
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${index_folder} \
--genomeFastaFiles ${genome_file} \
--sjdbGTFfile ${gtf_file} \
--sjdbOverhang 101 # this is the length of a mate read - 1
