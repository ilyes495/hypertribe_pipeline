# Load libraries ---------------------------
# Import libraries
library(bbmle)
library(parallel)
library(VGAM)

# Load data ---------------------------
# Data directories
output_folder <- "../../output_data/hypertribe/"

# Read formatted dataframe
formatted_df <- read.csv(paste0(output_folder, "4_formatted_output.csv"))

# Extract counts of reference (A) and alternative (G) nucleotides
ref_counts_df <- formatted_df[, grep("ref_counts", colnames(formatted_df))]
alt_counts_df <- formatted_df[, grep("alt_counts", colnames(formatted_df))]

# Likelihood models ---------------------------
#' Likelihood model for the null hypothesis
#'
#' This function computes the likelihood of the data under the null hypothesis.
#'
#' @param ref The number of edits in the control samples.
#' @param alt The number of edits in the treatment samples.
#'
#' @return The likelihood.
mle_custom_h0 <- function(ref, alt, debug = FALSE) {
  x <- alt
  size <- ref + alt

  bbll <- function(prob, rho) {
    if ((prob > 0) & (rho > 0) & (prob < 1) & (rho < 1)) {
      -sum(dbetabinom(x, size, prob, rho, log = TRUE))
    } else {
      NA
    }
  }

  fit <- mle2(bbll,
    start = list(prob = .1, rho = .5),
    method = "Nelder-Mead",
    control = list(maxit = 1e5, trace = as.integer(debug))
  )
}

#' Likelihood model for the alternative hypothesis
#'
#' This function computes the likelihood of the data under the alternative
#' hypothesis.
#'
#' @param ref The number of edits in the control samples.
#' @param alt The number of edits in the treatment samples.
#'
#' @return The likelihood.
mle_custom_h1 <- function(ref1, alt1, ref2, alt2, debug = FALSE) {
  x1 <- alt1
  size1 <- ref1 + alt1
  prob1_init <- mean(x1 / size1) + 1e-3
  if (is.na(prob1_init) | prob1_init >= 1 | prob1_init <= 0) {
    prob1_init <- .05
  }
  x2 <- alt2
  size2 <- ref2 + alt2
  prob2_init <- mean(x2 / size2) + 1e-3

  if (is.na(prob2_init) | prob2_init >= 1 | prob2_init <= 0) {
    prob2_init <- .05
  }

  bbll <- function(prob1, prob2, rho) {
    if ((prob1 > 0) & (prob2 > 0) & (rho > 0) &
      (prob1 < 1) & (prob2 < 1) & (rho < 1)) {
      -(sum(dbetabinom(x1, size1, prob1, rho, log = TRUE)) +
        sum(dbetabinom(x2, size2, prob2, rho, log = TRUE)))
    } else {
      NA
    }
  }

  fit <- mle2(bbll,
    start = list(prob1 = prob1_init, prob2 = prob2_init, rho = .1),
    method = "Nelder-Mead",
    control = list(maxit = 1e5, trace = as.integer(debug))
  )
}

# Computing the likelihood of data ---------------------------
#' Compute the differential likelihood of the data to be under the alternative
#' hypothesis
#'
#' This function computes the differential likelihood of the data to be under
#' the alternative hypothesis.
#'
#' @param i Row number of the formatted dataframe.
#' @param alt The formatted dataframe.
#'
#' @return The p value.
fit_mle <- function(i, formatted_df) {

  # Get the six counts of row i
  counts_list <- as.integer(formatted_df[i, grep(
    "counts",
    colnames(formatted_df)
  )])

  # Split the six counts into three references (A) and three alternatives (G)
  ref_list <- counts_list[seq(1, length(counts_list), 2)]
  alt_list <- counts_list[seq(2, length(counts_list), 2)]

  # Null hypothesis: ctrl = test
  fit0 <- mle_custom_h0(ref_list, alt_list)

  # Tested hypothesis: ctrl < test
  fit1 <- mle_custom_h1(ref_list[(1:3)], alt_list[(1:3)], ref_list[-c(1:3)], alt_list[-c(1:3)])

  # Compute p value, fold change and convergence
  pval <- pchisq(2 * (fit0@min - fit1@min), 1, lower.tail = FALSE)
  value <- fit1@min
  c(pval = pval, value = value)
}
res_list <- mclapply(seq_len(nrow(formatted_df)), fit_mle,
  formatted_df,
  mc.cores = 16
)
res_df <- data.frame(do.call("rbind", res_list))

# Format data ---------------------------
# Compute editing frequency
alt_freq_df <- alt_counts_df / (ref_counts_df + alt_counts_df)
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

# Only retain edits with high read counts
rowmask <- rowSums(alt_counts_df != 0) > 1 &
  (rowMeans(ref_counts_df) + rowMeans(alt_counts_df) >= 5) &
  res_df$value > 0
stats_df <- stats_df[rowmask, ]

# Compute adjusted p values
stats_df$padj <- p.adjust(stats_df$pval, "BH")
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

# Combine annotations and results
stats_df <- cbind(
  formatted_df[rowmask, 1:10],
  stats_df,
  formatted_df[rowmask, -(1:10)]
)

# Save dataframe
write.csv(stats_df,
  paste0(output_folder, "5_CELL_LINE_Control_Treatment.csv"),
  row.names = FALSE
)


