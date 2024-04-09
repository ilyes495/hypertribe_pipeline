# Load libraries ---------------------------
# Import libraries
library(biomaRt)
library(BiocParallel)
library(DESeq2)
library(dplyr)
library(GenomicAlignments)
library(GenomicFeatures)
library(parallel)
library(reshape2)
library(Rsamtools)
library(stringr)
library(VariantAnnotation)

# Load data ---------------------------
# Data directories
sample_list <- c(
  "Sample_New_Name_1", "Sample_New_Name_1", "Sample_New_Name_3",
  "Sample_New_Name_4", "Sample_New_Name_5", "Sample_New_Name_6"
  )

gtf_path <- "../../genome_data/GENOME_BUILD.gtf"
output_folder <- "../../output_data/hypertribe/"

variant_folder <- "../../output_data/hypertribe/3_variant/"
bam_path_list <- paste0(variant_folder, sample_list, "/3_reordered.bam")
names(bam_path_list) <- sample_list
vcf_path_list <- paste0(variant_folder, sample_list, "/6_output.vcf")


# Read transcripts from GTF file
TxDb <- makeTxDbFromGFF(gtf_path, format = "gtf")
txbygene_grl <- transcriptsBy(TxDb, "gene")

# Get read counts per transcript ---------------------------
read_counts_gr <- summarizeOverlaps(
  features = txbygene_grl,
  reads = BamFileList(bam_path_list,
    yieldSize = 1e6
  ),
  mode = "Union",
  singleEnd = FALSE,
  ignore.strand = TRUE,
  fragments = TRUE,
  BPPARAM = SnowParam(workers = 3)
)
saveRDS(read_counts_gr, paste0(output_folder, "4a_read_counts_gr.rds"))

# Filter VCF files ---------------------------
#' Read a VCF and filter entries
#'
#' This function reads a VCF file, and only keeps A/G edits in the forward
#' strand, T/C edits in the reverse strand, as well as those not found in
#' dbSNP and are present within annotated transcripts.
#'
#' @param vcf_path The path to the VCF file.
#' @param genome The genome name.
#'
#' @return The filtered VCF table as a GenomicRanges object.
filter_vcf <- function(vcf_path, genome) {

  # Read VCF file
  vcf <- readVcf(vcf_path, genome)

  # Only keep edits not found in dbSNP
  filtered_vcf <- vcf[!info(vcf)$DB]

  # Only keep A/G edits in forward strand
  pos_mask <- rowRanges(filtered_vcf)$REF == "A" &
    sapply(rowRanges(filtered_vcf)$ALT, function(seq) {
      return("G" %in% seq)
    })
  pos_gr <- rowRanges(filtered_vcf)[pos_mask]
  strand(pos_gr) <- "+"

  # Only keep T/C edits in reverse strand
  neg_mask <- rowRanges(filtered_vcf)$REF == "T" &
    sapply(rowRanges(filtered_vcf)$ALT, function(seq) {
      return("C" %in% seq)
    })
  neg_gr <- rowRanges(filtered_vcf)[neg_mask]
  strand(neg_gr) <- "-"

  # Combine forward and reverse strands
  filtered_vcf_gr <- c(pos_gr, neg_gr)

  # Only keep edits found in transcripts
  filtered_vcf_gr <- filtered_vcf_gr[countOverlaps(
    filtered_vcf_gr,
    txbygene_grl
  ) > 0]
}

# Read and filter VCF files
filtered_vcf_gr_list <- mclapply(vcf_path_list, filter_vcf, "hg19",
  mc.cores = length(vcf_path_list)
)
names(filtered_vcf_gr_list) <- sample_list
filtered_vcf_grl <- GRangesList(filtered_vcf_gr_list)

# Save VCF dataframe
saveRDS(filtered_vcf_grl, paste0(output_folder, "4b_filtered_vcf_grl.rds"))

# Get all edit sites that only match to one transcript
unique_vcf_gr <- sort(unique(unlist(filtered_vcf_grl, use.names = FALSE)))
subset_vcf_gr <- unique_vcf_gr[which(countOverlaps(unique_vcf_gr, txbygene_grl) == 1)]

# Generate pileup dataframes ---------------------------
#' Generate pileup dataframe for all edit sites for a sample
#'
#' This function reads the BAM file for a sample and generate a pileup
#' dataframe for the union of edit sites of all samples.
#'
#' @param bam_path The path to the BAM file.
#' @param unique_vcf_gr The GenomicRanges object describing the location
#'   of the union of edit sites of all samples.
#'
#' @return The pileup dataframe.
generate_pileup_df <- function(bam_path, unique_vcf_gr) {

  # Read BAM file
  bai_path <- paste0(str_split(bam_path, "bam")[[1]][1], "bai")
  bam <- BamFile(bam_path, bai_path)

  # Generate pileup dataframe
  bam_param <- ScanBamParam(
    flag = scanBamFlag(
      hasUnmappedMate = FALSE,
      isProperPair = TRUE,
      isDuplicate = FALSE
    ),
    which = unique_vcf_gr
  )
  pileup_param <- PileupParam(
    distinguish_strands = FALSE,
    min_base_quality = 10,
    max_depth = 1e4
  )
  pileup_df <- pileup(bam, unique_vcf_gr,
    scanBamParam = bam_param, pileupParam = pileup_param
  )
}

# Generate pileup dataframes for all samples
pileup_df_list <- mclapply(bam_path_list, generate_pileup_df, unique_vcf_gr,
                           mc.cores = length(bam_path_list))
names(pileup_df_list) <- sample_list

# Save list of pileup dataframes
saveRDS(pileup_df_list, paste0(output_folder, "4c_pileup_df_list.rds"))

# Process pileup dataframes ---------------------------
#' Process pileup dataframe for a sample
#'
#' This function adds strand information to and sort a pileup dataframe
#'
#' @param pileup_df The pileup dataframe.
#' @param unique_vcf_gr The GenomicRanges object describing the location
#'   of the union of edit sites of all samples.
#'
#' @return The processed pileup dataframe.
process_pileup_df <- function(pileup_df, unique_vcf_gr) {

  # Melt table by nucleotide coordinates
  pileup_df <- dcast(pileup_df, which_label ~ nucleotide,
    value.var = "count",
    fill = 0, drop = FALSE
  )

  # Add strand information
  pileup_df$strand <- as.character(strand(unique_vcf_gr))

  # Check if there are missing rows
  snp_id_list <- sprintf(
    "%s:%d-%d",
    seqnames(unique_vcf_gr),
    start(unique_vcf_gr),
    start(unique_vcf_gr)
  )
  stopifnot(all(snp_id_list == pileup_df$which_label))
  row.names(pileup_df) <- pileup_df$which_label

  # Split and reorganize dataframe based on strand
  pileup_df <- split(pileup_df, pileup_df$strand)
  pileup_df$`+` <- data.frame(
    ref_count = pileup_df$`+`$A,
    alt_count = pileup_df$`+`$G,
    row.names = row.names(pileup_df$`+`)
  )
  pileup_df$`-` <- data.frame(
    ref_count = pileup_df$`-`$T,
    alt_count = pileup_df$`-`$C,
    row.names = row.names(pileup_df$`-`)
  )

  # Combine and sort data frames
  rbind(pileup_df$`+`, pileup_df$`-`)[snp_id_list, ]
}

# Process pileup dataframe for all samples
processed_pileup_df_list <- lapply(
  pileup_df_list, process_pileup_df,
  unique_vcf_gr
)

# Save list of processed pileup dataframes
saveRDS(
  processed_pileup_df_list,
  paste0(output_folder, "4d_processed_pileup_df_list.rds")
)

# Get genomic data ---------------------------
#' Obtain genomic information for all edit sites
#'
#' This function queries the Ensembl database and returns the gene symbol,
#' for the genes of all edit sites. It also queries for the genomic feature
#' where the edit site is at.
#'
#' @param subset_vcf_gr All edit sites that are present in only one transcript.
#' @param txbygene_grl The list of all genes and their transcripts.
#' @param gtf_path The path to the GTF file.
#'
#' @return A dataframe consisting of all genomic information of the edits.
get_anno_df <- function(subset_vcf_gr, txbygene_grl, gtf_path) {

  # Get matches between edits and transcripts
  hits <- findOverlaps(subset_vcf_gr, txbygene_grl)

  # Get names of transcripts
  ensg_id_list <- names(txbygene_grl)[subjectHits(hits)]

  # Load Ensembl database
  ensembl <- useEnsembl(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl"
  )

  # Get gene symbols for all ENSG IDs
  genomic_df <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensg_id_list,
    mart = ensembl,
    useCache = FALSE
  )

  # Remove entries with duplicated ENSG IDs
  genomic_df <- genomic_df %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)

  # Generate annotation dataframe
  anno_df <- DataFrame(ensg_id = ensg_id_list)
  anno_df <- merge(anno_df, genomic_df,
    by.x = "ensg_id", by.y = "ensembl_gene_id",
    all.x = TRUE, sort = FALSE
  )
  anno_df <- anno_df[match(ensg_id_list, anno_df$ensg_id), ]
  names(anno_df) <- c("ensg_id", "gene_symbol")

  # Determine whether edit is in intro, exon, 5'UTR or 3'UTR
  gtf <- read.table(gtf_path, header = FALSE, sep = "\t")

  exon_gtf <- gtf[gtf$V3 == "exon", c(1, 4, 5, 7)]
  names(exon_gtf) <- c("seqname", "start", "end", "strand")
  exon_gr <- makeGRangesFromDataFrame(exon_gtf)

  utr5_gtf <- gtf[gtf$V3 == "5UTR", c(1, 4, 5, 7)]
  names(utr5_gtf) <- c("seqname", "start", "end", "strand")
  utr5_gr <- makeGRangesFromDataFrame(utr5_gtf)

  utr3_gtf <- gtf[gtf$V3 == "3UTR", c(1, 4, 5, 7)]
  names(utr3_gtf) <- c("seqname", "start", "end", "strand")
  utr3_gr <- makeGRangesFromDataFrame(utr3_gtf)

  anno_df$feature <- "intron"
  anno_df$feature[countOverlaps(subset_vcf_gr, exon_gr) > 0] <- "exon"
  anno_df$feature[countOverlaps(subset_vcf_gr, utr5_gr) > 0] <- "utr5"
  anno_df$feature[countOverlaps(subset_vcf_gr, utr3_gr) > 0] <- "utr3"
  anno_df
}
anno_df <- get_anno_df(subset_vcf_gr, txbygene_grl, gtf_path)

# Combine data ---------------------------
#' Adds counts to all edit sites
#'
#' This function adds counts to the locations of all edit sites.
#'
#' @param pileup_df The pileup dataframe.
#' @param subset_vcf_gr All edit sites that are present in only one transcript.
#'
#' @return A dataframe consisting of both the counts and location of edit sites.
add_counts <- function(pileup_df, subset_vcf_gr) {
  snp_id_list <- sprintf(
    "%s:%d-%d",
    seqnames(subset_vcf_gr),
    start(subset_vcf_gr),
    start(subset_vcf_gr)
  )
  pileup_df <- pileup_df[snp_id_list, ]
  formatted_df <- cbind(as.data.frame(subset_vcf_gr), pileup_df)
  formatted_df
}

# Add genomic information
mcols(subset_vcf_gr) <- anno_df

# Add counts to all pileup dataframes
formatted_df <- lapply(processed_pileup_df_list, add_counts, subset_vcf_gr)

# Get genomic information
formatted_anno_df <- formatted_df[[1]][(1:(ncol(formatted_df[[1]]) - 2))]

# Get counts
formatted_counts_df <- do.call(
  "cbind",
  lapply(
    formatted_df,
    function(df) {
      df[-(1:(ncol(formatted_df[[1]]) - 2))]
    }
  )
)


# Combine annotations and count dataframes
formatted_df <- cbind(formatted_anno_df, formatted_counts_df)

# Format data ---------------------------
# Add number of edits per gene
event_id_list <- row.names(formatted_df)
gene_num_edits_df <- as.data.frame(table(formatted_df$ensg_id))
names(gene_num_edits_df) <- c("ensg_id", "gene_num_events")
formatted_df <- merge(formatted_df, gene_num_edits_df,
  by = "ensg_id",
  keep.x = TRUE, keep.y = FALSE, sort = FALSE
)
formatted_df$event_id <- event_id_list

# Add FPKM
colData(read_counts_gr)$condition <- factor(sample_list)
fpkm_df <- fpkm(DESeqDataSet(read_counts_gr, ~condition))
fpkm_df <- data.frame(
  ensg_id = row.names(fpkm_df),
  ctrl_fpkm = rowMeans(fpkm_df[, 1:3]),
  test_fpkm = rowMeans(fpkm_df[, 4:6])
)
formatted_df <- merge(formatted_df, fpkm_df, by = "ensg_id")

# Rename columns
old_col_list <- c(
  "ensg_id", "chr", "start", "end", "length",
  "strand", "gene_symbol", "feature",
  "C1_ref_counts", "C1_alt_counts",
  "C2_ref_counts", "C2_alt_counts",
  "C3_ref_counts", "C3_alt_counts",
  "AS_1_ref_counts", "AS_1_alt_counts",
  "AS_2_ref_counts", "AS_2_alt_counts",
  "AS_3_ref_counts", "AS_3_alt_counts",
  "gene_num_events", "event_id", "ctrl_fpkm", "test_fpkm"
)
new_col_list <- c(
  "ensg_id", "gene_symbol", "chr", "strand", "start",
  "end", "length", "feature", "event_id", "gene_num_events",
  "C1_ref_counts", "C1_alt_counts",
  "C2_ref_counts", "C2_alt_counts",
  "C3_ref_counts", "C3_alt_counts",
  "AS_1_ref_counts", "AS_1_alt_counts",
  "AS_2_ref_counts", "AS_2_alt_counts",
  "AS_3_ref_counts", "AS_3_alt_counts",
  "ctrl_fpkm", "test_fpkm"
)
names(formatted_df) <- old_col_list
formatted_df <- formatted_df %>%
  dplyr::select(new_col_list) %>%
  arrange(ensg_id)

# Save dataframe
write.csv(formatted_df,
  paste0(output_folder, "4_formatted_output.csv"),
  row.names = FALSE
)
