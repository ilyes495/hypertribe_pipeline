#!/bin/bash
#BSUB -J 3_variant_calling[1-6]
#BSUB -W 48:00
#BSUB -R rusage[mem=20]
#BSUB -n 8
#BSUB -C 0
#BSUB -o logs/3_variant_calling_%J_%I.stdout
#BSUB -eo logs/3_variant_calling_%J_%I.stderr

# initialize conda environment
source ~/.bashrc
conda activate hypertribe
cd $LS_SUBCWD

# define samples
samples=("Sample_New_Name_1" "Sample_New_Name_1" "Sample_New_Name_3" \
           "Sample_New_Name_4" "Sample_New_Name_5" "Sample_New_Name_6")

i=$((LSB_JOBINDEX-1))
sample=${samples[$i]}
echo $sample

# input directories
align_folder=../../output_data/hypertribe/bam_files
genome_file=../../genome_data/hg38_RBP_ADAR.fa
genome_dict_file=../../genome_data/hg38_RBP_ADAR.dict
dbsnp_file=../../genome_data/dbsnp/dbsnp.vcf.gz

# output directories
variant_folder=../../output_data/hypertribe/3_variant

# executable directories
GenomeAnalysisTK_jar=../../software/hypertribe/GenomeAnalysisTK.jar

# load Java 8
module load java

make directories
mkdir -p ${variant_folder}/${sample}

# picard: replace read groups of aligned reads
echo "Run picard AddOrReplaceReadGroups on samples"
picard AddOrReplaceReadGroups \
I=${align_folder}/${sample}/${sample}.bam \
O=${variant_folder}/${sample}/1_rg_added.bam \
RGID=id \
RGLB=library \
RGPL=ILLUMINA \
RGPU=machine \
RGSM=sample \

# picard: mark duplicates of aligned reads. Reserve 30G memory for each file.
echo "Run picard MarkDuplicates on samples"
picard MarkDuplicates \
I=${variant_folder}/${sample}/1_rg_added.bam \
O=${variant_folder}/${sample}/2_dedupped.bam \
M=${variant_folder}/${sample}/2_output.metrics \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true \

# picard: reorder sam files of aligned reads
echo "Run picard ReorderSam on samples"
picard ReorderSam \
I=${variant_folder}/${sample}/2_dedupped.bam \
O=${variant_folder}/${sample}/3_reordered.bam \
R=${genome_file} \
SD=${genome_dict_file} \
CREATE_INDEX=TRUE \

# GenomeAnalysisTK: split aligned reads
echo "Run GenomeAnalysisTK SplitNCigarReads on samples"
java -jar ${GenomeAnalysisTK_jar} \
-T SplitNCigarReads \
-R ${genome_file} \
-I ${variant_folder}/${sample}/3_reordered.bam \
-o ${variant_folder}/${sample}/4_split.bam \
-rf ReassignOneMappingQuality \
-RMQF 255 \
-RMQT 60 \
-U ALLOW_N_CIGAR_READS \

# GenomeAnalysisTK: call SNPs of aligned reads
echo "Run GenomeAnalysisTK HaplotypeCaller on samples"
java -jar ${GenomeAnalysisTK_jar} \
-T HaplotypeCaller \
-R ${genome_file} \
-I ${variant_folder}/${sample}/4_split.bam \
-o ${variant_folder}/${sample}/5_haplotypes.vcf \
--dbsnp ${dbsnp_file} \
-dontUseSoftClippedBases \
-stand_call_conf 10.0 \

# GenomeAnalysisTK: annotate SNPs
echo "Run GenomeAnalysisTK VariantFiltration on samples"
java -jar ${GenomeAnalysisTK_jar} \
-T VariantFiltration \
-R ${genome_file} \
-V ${variant_folder}/${sample}/5_haplotypes.vcf \
-o ${variant_folder}/${sample}/6_output.vcf \
-window 35 \
-cluster 3 \
-filterName FS \
-filter "FS > 30.0" \
-filterName QD \
-filter "QD < 2.0" \
