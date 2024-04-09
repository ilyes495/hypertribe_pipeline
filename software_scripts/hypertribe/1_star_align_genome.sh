#!/bin/bash
#BSUB -J 1_star_align_genome[1]
#BSUB -W 8:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/1_star_align_genome_%J_%I.stdout
#BSUB -eo logs/1_star_align_genome_%J_%I.stderr

# initialize conda environment
source ~/.bashrc
conda activate hypertribe
cd $LS_SUBCWD

# define samples
samples=("Sample_1" "Sample_1" "Sample_3" \
         "Sample_4" "Sample_5" "Sample_6")

new_names=("Sample_New_Name_1" "Sample_New_Name_1" "Sample_New_Name_3" \
           "Sample_New_Name_4" "Sample_New_Name_5" "Sample_New_Name_6")


i=$((LSB_JOBINDEX-1))
sample=${samples[$i]};
name=${new_names[$i]};
echo ${sample} $i ;

# input directories
index_folder=../../genome_data/STAR_INDEX_OUTPUT
sequence_folder=../../input_data/

# output directories
align_folder=../../output_data/hypertribe/1_star_align

# make directories
mkdir -p logs
mkdir -p ${align_folder}
mkdir -p ${align_folder}/${name}
mkdir -p ${align_folder}/${name}/bamqc

# star: align reads to transcriptome
echo "Run STAR alignReads on " ${sample}
STAR --runMode alignReads \
--runThreadN 16 \
--genomeDir ${index_folder} \
--readFilesIn ${sequence_folder}/${sample}/*L001_R1_001.fastq.gz  ${sequence_folder}/${sample}/*L001_R2_001.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ${align_folder}/${name}/${name}. \
--outSAMtype BAM SortedByCoordinate \
--outFilterMultimapNmax 1 \
--outFilterMultimapScoreRange 0 \
--outFilterMismatchNmax 5 \
--alignIntronMin 70 \
--alignIntronMax 100000

# rename bam files
mv \
${align_folder}/${name}/${name}.Aligned.sortedByCoord.out.bam \
${align_folder}/${name}/${name}.bam 

# samtools: index aligned reads
echo "Run samtools index on " ${name}
samtools index \
${align_folder}/${name}/${name}.bam

# qualimap: bamqc aligned reads
echo "Run qualimap bamqc on " ${name}
qualimap bamqc \
--java-mem-size=20G \
-bam ${align_folder}/${name}/${name}.bam \
-outdir ${align_folder}/${name}/bamqc/
