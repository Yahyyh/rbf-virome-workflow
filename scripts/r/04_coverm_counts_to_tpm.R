#!/usr/bin/env Rscript

# Convert per-sample CoverM count tables to a feature-by-sample count matrix and length-normalized TPM matrix.
#
# Usage:
#   Rscript 04_coverm_counts_to_tpm.R \
#     lengths.tsv count_directory output_prefix
#
# lengths.tsv: exactly one feature-ID column and one positive length_bp column.
# count_directory: one two-column TSV/TXT per sample (feature ID, read count).
# The basename (with .counts/_counts removed) becomes sample_id.
#
# Outputs:
#   <prefix>.counts.tsv       feature_id plus sample counts (NO length column)
#   <prefix>.tpm.tsv          feature_id plus sample TPM
#   <prefix>.lengths.tsv      feature_id plus length_bp
#   <prefix>.library_qc.tsv   count totals, RPK totals, and TPM-sum check

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 04_coverm_counts_to_tpm.R lengths.tsv count_directory output_prefix",
    call. = FALSE
  )
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path[1L]), winslash = "/"))
source(file.path(script_dir, "lib", "common.R"))

length_table <- read_tsv(args[1L])
if (ncol(length_table) != 2L) {
  stop("Length table must have exactly two columns: feature_id and length_bp.", call. = FALSE)
}
names(length_table) <- c("feature_id", "length_bp")
length_table$feature_id <- trimws(as.character(length_table$feature_id))
length_table$length_bp <- suppressWarnings(as.numeric(length_table$length_bp))
assert_unique_nonmissing(length_table$feature_id, "length-table feature_id")
if (
  anyNA(length_table$length_bp) ||
  any(!is.finite(length_table$length_bp)) ||
  any(length_table$length_bp <= 0)
) {
  stop("All feature lengths must be finite, numeric, non-missing, and positive.", call. = FALSE)
}

count_directory <- args[2L]
if (!dir.exists(count_directory)) {
  stop(sprintf("Count directory does not exist: %s", count_directory), call. = FALSE)
}
count_files <- list.files(
  count_directory,
  pattern = "\\.(tsv|txt)$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (length(count_files) == 0L) {
  stop("No .tsv or .txt count files were found in the count directory.", call. = FALSE)
}

sample_name_from_path <- function(path) {
  sample_id <- tools::file_path_sans_ext(basename(path))
  sample_id <- sub("([._-]counts?)$", "", sample_id, ignore.case = TRUE)
  trimws(sample_id)
}
sample_ids <- vapply(count_files, sample_name_from_path, character(1))
assert_unique_nonmissing(sample_ids, "sample IDs derived from count filenames")

known_features <- length_table$feature_id
count_tables <- Map(function(path, sample_id) {
  input <- read_tsv(path)
  if (ncol(input) != 2L) {
    stop(sprintf("Count file must contain exactly two columns: %s", path), call. = FALSE)
  }
  names(input) <- c("feature_id", sample_id)
  input$feature_id <- trimws(as.character(input$feature_id))
  assert_unique_nonmissing(input$feature_id, sprintf("feature IDs in %s", basename(path)))
  input[[sample_id]] <- suppressWarnings(as.numeric(input[[sample_id]]))
  if (
    anyNA(input[[sample_id]]) ||
    any(!is.finite(input[[sample_id]])) ||
    any(input[[sample_id]] < 0)
  ) {
    stop(sprintf("Counts must be finite, numeric, non-missing, and non-negative: %s", path), call. = FALSE)
  }
  unknown <- setdiff(input$feature_id, known_features)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "%s contains %d feature(s) absent from the length table; first: %s",
        basename(path), length(unknown), paste(head(unknown, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  input
}, count_files, sample_ids)

counts_merged <- Reduce(
  function(x, y) merge(x, y, by = "feature_id", all = TRUE, sort = FALSE),
  count_tables
)
counts_merged[is.na(counts_merged)] <- 0
counts <- merge(
  length_table["feature_id"],
  counts_merged,
  by = "feature_id",
  all.x = TRUE,
  sort = FALSE
)
counts[is.na(counts)] <- 0
counts <- counts[match(length_table$feature_id, counts$feature_id), , drop = FALSE]
sample_columns <- setdiff(names(counts), "feature_id")
count_matrix <- as.matrix(counts[sample_columns])
storage.mode(count_matrix) <- "double"

if (any(colSums(count_matrix) <= 0)) {
  empty <- sample_columns[colSums(count_matrix) <= 0]
  stop(
    sprintf("Zero-total sample(s) cannot be converted to TPM: %s", paste(empty, collapse = ", ")),
    call. = FALSE
  )
}

rpk <- sweep(count_matrix, 1L, length_table$length_bp / 1000, "/")
rpk_totals <- colSums(rpk)
if (any(rpk_totals <= 0)) {
  stop("At least one sample has a zero RPK total and cannot be converted to TPM.", call. = FALSE)
}
tpm_matrix <- sweep(rpk, 2L, rpk_totals / 1e6, "/")
tpm <- data.frame(feature_id = length_table$feature_id, tpm_matrix, check.names = FALSE)

tpm_sums <- colSums(tpm_matrix)
if (any(abs(tpm_sums - 1e6) > 1e-4)) {
  stop("Internal TPM validation failed: non-empty sample columns do not sum to 1,000,000.", call. = FALSE)
}

qc <- data.frame(
  sample_id = sample_columns,
  total_counts = colSums(count_matrix),
  nonzero_features = colSums(count_matrix > 0),
  total_rpk = rpk_totals,
  tpm_sum = tpm_sums,
  stringsAsFactors = FALSE
)

output_prefix <- args[3L]
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
write_tsv(counts, paste0(output_prefix, ".counts.tsv"))
write_tsv(tpm, paste0(output_prefix, ".tpm.tsv"))
write_tsv(length_table, paste0(output_prefix, ".lengths.tsv"))
write_tsv(qc, paste0(output_prefix, ".library_qc.tsv"))
writeLines(
  c(
    "TPM = (count / length_kb) / sum(count / length_kb) * 1,000,000.",
    "TPM is within-library relative abundance and must not be interpreted as absolute concentration."
  ),
  paste0(output_prefix, ".README.txt")
)
writeLines(
  capture.output(utils::sessionInfo()),
  paste0(output_prefix, ".sessionInfo.txt")
)
