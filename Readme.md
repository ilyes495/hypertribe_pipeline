# Instructions to Run HyperTribe Pipeline

## 1. Installation of Conda Environment

### Prerequisites:
- Anaconda or Miniconda installed on your system.
- R
- All the scripts are written to run on LSF environment

### install the conda/mamba environment
Create the Conda environment named "hypertribe" using the provided environment.yml file:
   ```bash
   conda env create -f environment.yml 
   ```

## 2. Downloading Genome and GTF Annotation Files

### Steps:
Navigate to the "download_scripts" directory.
Run the provided scripts to download the genome and GTF annotation files.

```bash
./download_genome.sh
./download_gtf.sh
```

## 3. Place the fastq files under `input_data` folder

## 4. HyperTRIBE Psudogene Files

we provided a sample sequence and gtf annotation of a unique sequence of the HyperTRIBE construct that can be used to quantify its expression level.
modify the two files to include your specific experiment construct, such as the name of the file and pseudo gene name. also set `L` in the gtf file to the lenght of the unique sequence.

after that, concatenate the `.fa` files from the geneme and the construct pseudo gene into a single `.fa` file 

Example:
```bash
    cat hg38.fa RBP_ADAR.fa > hg38_RBP_ADAR.fa
```
do the same thing with the gtf files.

## 5. Star Indexer Scripts

### Steps:
1. Navigate to the "genome_data" folder.
2. Modify the paths of the files in the `star_index_genome.sh` and `picard_dictionary_genome.sh` scripts
3. Run the star indexer script to generate necessary index file.

```bash
bsub < star_index_genome.sh 
```

then 

```bash
bsub < picard_dictionary_genome.sh 
```


## 5. Running the Pipeline

### Steps:

1. Navigate to the "software_scripts/hypertribe" folder.
2. Execute the pipeline scripts in the given order:

### Step 1: Alignemnt
a. modify the samples' file name inside `1_star_align_genome.sh` script 

```bash
samples=("Sample_1" "Sample_1" "Sample_3" \
         "Sample_4" "Sample_5" "Sample_6")
```
also if you want to update the filenames to a more descriptive filenames, modify the following line accordingly, otherwise provide the same input as for `samples`

``` bash
new_names=("Sample_New_Name_1" "Sample_New_Name_1" "Sample_New_Name_3" \
           "Sample_New_Name_4" "Sample_New_Name_5" "Sample_New_Name_6")
```
b. Set the path to Star Index folder 
```bash
index_folder=../../genome_data/STAR_INDEX_OUTPUT
```

c. run the alignemnt step

```bash
bsub < 1_star_align_genome.sh
```

### Step 2: MultiQC
```bash
bsub < 2_multiqc.sh
```

### Step 3: Variant Calling
a. modify the following lines as in previous steps:

```bash 
samples=("Sample_New_Name_1" "Sample_New_Name_1" "Sample_New_Name_3" \
           "Sample_New_Name_4" "Sample_New_Name_5" "Sample_New_Name_6")           
```
and 

```bash
genome_file=../../genome_data/hg38_RBP_ADAR.fa
genome_dict_file=../../genome_data/hg38_RBP_ADAR.dict
dbsnp_file=../../genome_data/dbsnp/dbsnp.vcf.gz
```

and finally:

```bash
bsub < 3_variant_calling.sh
```

### Step 4: Format Output

This step aims at aggregating the results of the variant calling for each sample into a single file that can be used to exploration or downstream analysis.

This step requires R.

a. modify the following lines in the `4_format_output.R` script:

```R
sample_list <- c(
  "Sample_New_Name_1", "Sample_New_Name_1", "Sample_New_Name_3",
  "Sample_New_Name_4", "Sample_New_Name_5", "Sample_New_Name_6"
  )

gtf_path <- "../../genome_data/hg38_RBP_ADAR.gtf"
```

b. run the script

```bash 
bsub < 4_format_output.sh
```

### Step 5: Differential Analysis

This step runs differential analysis between the control and HyperTRIBE samples to identify significant edited sites. it applies the same filtering and processing steps as described in the original HyperTRIBE paper.

There are various part of the script that need to be modified depending on the number of samples in each group. the current script assumes 3 samples per group.

The following lines need to be modified accordingly:

```R
fit1 <- mle_custom_h1(ref_list[(1:3)], alt_list[(1:3)], ref_list[-c(1:3)], alt_list[-c(1:3)])
```

```R
ctrl_freq_list <- rowMeans(alt_freq_df[, c(1:3)], na.rm = TRUE)
test_freq_list <- rowMeans(alt_freq_df[, c(4:6)], na.rm = TRUE)
stats_df <- data.frame(
  diff_mean = test_freq_list - ctrl_freq_list,
  ctrl_mean = ctrl_freq_list,
  test_mean = test_freq_list,
  pval = res_df$pval,
  Control_Sample_1_freq =  alt_freq_df[, 1],
  Control_Sample_2_freq =  alt_freq_df[, 2],
  Control_Sample_3_freq =  alt_freq_df[, 3],
  Treatment_Sample_1_freq =  alt_freq_df[, 4],
  Treatment_Sample_2_freq =  alt_freq_df[, 5], 
  Treatment_Sample_3_freq =  alt_freq_df[, 6]
)
```

```R
stats_df <- stats_df[c(
  "diff_mean", "ctrl_mean", "test_mean",
  "pval", "padj",
  "Control_Sample_1_freq", 
  "Control_Sample_2_freq",
  "Control_Sample_3_freq",
  "Treatment_Sample_1_freq",
  "Treatment_Sample_2_freq",
  "Treatment_Sample_3_freq" 
)]
```

Rename the ouput filename:

```R
write.csv(stats_df,
  paste0(output_folder, "5_CELL_LINE_Control_Treatment.csv"),
  row.names = FALSE
)
```

Finally, run the script

```bash
bsub < 5_test_differential.sh
```



