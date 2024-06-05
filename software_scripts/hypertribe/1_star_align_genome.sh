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

# Trim raw reads with cutadapt
echo "Run cutadapt on " ${sample}
cutadapt -O 5 \
--match-read-wildcards \
--times 2 -e 0.0 \
--quality-cutoff 6 -m 18 \
-o ${align_folder}/${name}/${name}_L001_R1.fastqTr.fq -p ${align_folder}/${name}/${name}_L001_R2.fastqTr.fq \
-b TCGTATGCCGTCTTCTGCTTG \
-b ATCTCGTATGCCGTCTTCTGCTTG \
-b CGACAGGTTCAGAGTTCTACAGTCCGACGATC \
-b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-b AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA \
-b TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
${sequence_folder}/${sample}/*L001_R1_001.fastq.gz  ${sequence_folder}/${sample}/*L001_R2_001.fastq.gz

 



# star: align reads to transcriptome
echo "Run STAR alignReads on " ${sample}
STAR --runMode alignReads \
--runThreadN 16 \
--genomeDir ${index_folder} \
--readFilesIn ${align_folder}/${name}/${name}_L001_R1.fastqTr.fq ${align_folder}/${name}/${name}_L001_R2.fastqTr.fq \
--outFileNamePrefix ${align_folder}/${name}/${name}. \
--outSAMtype BAM SortedByCoordinate \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--genomeLoad NoSharedMemory \
--outBAMcompression 10 \
--outFilterMultimapNmax 10 \
--outFilterMultimapScoreRange 1 \
--outFilterScoreMin 10 \
--outReadsUnmapped Fastx \
--outSAMattrRGline ID:${name} LB:library PL:illumina PU:machine SM:GM12878 \
--outSAMattributes All \
--outSAMmode Full \
--outSAMunmapped Within \
--outStd Log \

echo "Run STAR 2-pass mapping (improve alignmets using table of splice junctions and create a new index)" ${sample}
STAR --runMode alignReads \
--runThreadN 16 \
--genomeDir ${index_folder} \
--readFilesIn ${align_folder}/${name}/${name}_L001_R1.fastqTr.fq ${align_folder}/${name}/${name}_L001_R2.fastqTr.fq \
--outFileNamePrefix ${align_folder}/${name}/${name}. \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--sjdbFileChrStartEnd ${align_folder}/${name}/${name}.SJ.out.tab \
--outSAMtype BAM SortedByCoordinate \
--outSAMattrRGline ID:${name} LB:library PL:illumina PU:machine SM:GM12878 \
--genomeLoad NoSharedMemory \
--outBAMcompression 10 \
--outFilterMultimapNmax 10 \
--outFilterMultimapScoreRange 1 \
--outFilterScoreMin 10 \
--outFilterType BySJout \
--outReadsUnmapped Fastx \
--outSAMattributes All \
--outSAMmode Full \
--outSAMunmapped Within \
--outStd Log \

# Select only unique alignments, no multimaps
(samtools view -H ${align_folder}/${name}/${name}.Aligned.sortedByCoord.out.bam; samtools view ${align_folder}/${name}/${name}.Aligned.sortedByCoord.out.bam| grep -w 'NH:i:1') \
|samtools view -Sb - > ${align_folder}/${name}/${name}.Aligned.sortedByCoord.uniq.bam

# rename bam files
mv \
${align_folder}/${name}/${name}.Aligned.sortedByCoord.uniq.bam \
${align_folder}/${name}/${name}.bam 


# samtools: index aligned reads
echo "Run samtools index on " ${name}
samtools index \
${align_folder}/${name}/${name}.bam


# create tmp bam without RG tag for qualimap
echo "Run samtools view on " ${name}
samtools view -h ${align_folder}/${name}/${name}.bam | \
grep -v "^@RG" | \
samtools view -b > \
${align_folder}/${name}/${name}_tmp.bam


# qualimap: bamqc aligned reads
echo "Run qualimap bamqc on " ${name}
qualimap bamqc \
--java-mem-size=20G \
-bam ${align_folder}/${name}/${name}_tmp.bam \
-outdir ${align_folder}/${name}/bamqc/

# remove tmp bam
rm ${align_folder}/${name}/${name}_tmp.bam

