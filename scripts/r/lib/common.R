RBF_RANDOM_SEED <- 20260925L

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Required package '%s' is not installed.", package), call. = FALSE)
  }
}

assert_columns <- function(data, required, label = "input") {
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "%s is missing required column(s): %s",
        label,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

assert_unique_nonmissing <- function(values, label) {
  values <- as.character(values)
  if (anyNA(values) || any(trimws(values) == "")) {
    stop(sprintf("%s contains missing or empty values.", label), call. = FALSE)
  }
  if (anyDuplicated(values)) {
    stop(sprintf("%s must be unique.", label), call. = FALSE)
  }
  invisible(values)
}

read_tsv <- function(path) {
  require_namespace("data.table")
  if (!file.exists(path)) {
    stop(sprintf("Input file does not exist: %s", path), call. = FALSE)
  }
  data.table::fread(
    path,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE
  )
}

write_tsv <- function(data, path) {
  require_namespace("data.table")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data, path, sep = "\t", quote = FALSE, na = "NA")
}

read_feature_matrix <- function(path, label = "feature matrix") {
  table <- read_tsv(path)
  if (ncol(table) < 2L) {
    stop(
      sprintf("%s must contain one feature-ID column and at least one sample column.", label),
      call. = FALSE
    )
  }

  feature_ids <- as.character(table[[1L]])
  assert_unique_nonmissing(feature_ids, sprintf("%s feature IDs", label))
  assert_unique_nonmissing(names(table)[-1L], sprintf("%s sample column names", label))

  values <- table[, -1L, drop = FALSE]
  values[] <- lapply(values, function(column) suppressWarnings(as.numeric(column)))
  if (anyNA(values) || any(!is.finite(as.matrix(values)))) {
    stop(sprintf("All %s values must be finite, numeric, and non-missing.", label), call. = FALSE)
  }
  if (any(as.matrix(values) < 0)) {
    stop(sprintf("All %s values must be non-negative.", label), call. = FALSE)
  }

  matrix <- as.matrix(values)
  rownames(matrix) <- feature_ids
  storage.mode(matrix) <- "double"
  matrix
}

read_metadata <- function(path, required = c("sample_id", "source", "season")) {
  metadata <- read_tsv(path)
  assert_columns(metadata, required, "metadata")
  metadata[] <- lapply(metadata, function(column) {
    if (is.character(column)) trimws(column) else column
  })
  metadata$sample_id <- as.character(metadata$sample_id)
  assert_unique_nonmissing(metadata$sample_id, "metadata sample_id")

  for (column in intersect(c("source", "season", "pair_id", "group", "nucleic_acid"), names(metadata))) {
    metadata[[column]] <- as.character(metadata[[column]])
    if (anyNA(metadata[[column]]) || any(trimws(metadata[[column]]) == "")) {
      stop(sprintf("metadata column '%s' contains missing or empty values.", column), call. = FALSE)
    }
  }
  metadata
}

align_feature_matrix <- function(feature_matrix, metadata, minimum_samples = 3L) {
  shared <- metadata$sample_id[metadata$sample_id %in% colnames(feature_matrix)]
  if (length(shared) < minimum_samples) {
    stop(
      sprintf(
        "Only %d samples are shared by the matrix and metadata; at least %d are required.",
        length(shared), minimum_samples
      ),
      call. = FALSE
    )
  }
  missing_metadata <- setdiff(colnames(feature_matrix), metadata$sample_id)
  if (length(missing_metadata) > 0L) {
    warning(
      sprintf(
        "Dropping matrix sample(s) absent from metadata: %s",
        paste(missing_metadata, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  list(
    abundance = t(feature_matrix[, shared, drop = FALSE]),
    metadata = metadata[match(shared, metadata$sample_id), , drop = FALSE]
  )
}

align_two_feature_matrices <- function(first, second, metadata, minimum_samples = 4L) {
  shared <- metadata$sample_id[
    metadata$sample_id %in% colnames(first) & metadata$sample_id %in% colnames(second)
  ]
  if (length(shared) < minimum_samples) {
    stop(
      sprintf(
        "Only %d samples are shared by both matrices and metadata; at least %d are required.",
        length(shared), minimum_samples
      ),
      call. = FALSE
    )
  }
  list(
    first = t(first[, shared, drop = FALSE]),
    second = t(second[, shared, drop = FALSE]),
    metadata = metadata[match(shared, metadata$sample_id), , drop = FALSE]
  )
}

drop_zero_features <- function(sample_by_feature, label = "matrix") {
  keep <- colSums(sample_by_feature) > 0
  output <- sample_by_feature[, keep, drop = FALSE]
  if (ncol(output) == 0L) {
    stop(sprintf("%s contains no non-zero features.", label), call. = FALSE)
  }
  if (any(rowSums(output) <= 0)) {
    stop(sprintf("Every sample in %s must have a positive total.", label), call. = FALSE)
  }
  output
}

assert_integer_counts <- function(sample_by_feature, tolerance = 1e-8) {
  if (any(abs(sample_by_feature - round(sample_by_feature)) > tolerance)) {
    stop(
      "Count data must contain non-negative integers; do not use TPM or relative abundance for richness estimators.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

to_relative_abundance <- function(sample_by_feature) {
  totals <- rowSums(sample_by_feature)
  if (any(totals <= 0)) {
    stop("Cannot normalize samples with zero total abundance.", call. = FALSE)
  }
  sweep(sample_by_feature, 1L, totals, "/")
}

analysis_group <- function(metadata) {
  assert_columns(metadata, c("source", "season"), "metadata")
  factor(
    paste(metadata$source, metadata$season, sep = ":"),
    levels = unique(paste(metadata$source, metadata$season, sep = ":"))
  )
}

within_block_permutations <- function(blocks, nperm = 999L) {
  require_namespace("permute")
  permute::how(nperm = nperm, blocks = factor(blocks))
}

whole_block_permutations <- function(blocks, nperm = 999L) {
  require_namespace("permute")
  permute::how(
    nperm = nperm,
    within = permute::Within(type = "none"),
    plots = permute::Plots(strata = factor(blocks), type = "free")
  )
}

extract_adonis_term <- function(model, term, comparison, scheme, n_samples) {
  if (!term %in% rownames(model)) {
    stop(sprintf("Term '%s' is absent from the PERMANOVA model.", term), call. = FALSE)
  }
  data.frame(
    comparison = comparison,
    term = term,
    n_samples = n_samples,
    df = unname(model[term, "Df"]),
    sum_of_squares = unname(model[term, "SumOfSqs"]),
    pseudo_F = unname(model[term, "F"]),
    r2 = unname(model[term, "R2"]),
    p_value = unname(model[term, "Pr(>F)"]),
    permutation_scheme = scheme,
    stringsAsFactors = FALSE
  )
}

save_plot_pair <- function(plot, path_prefix, width, height) {
  require_namespace("ggplot2")
  require_namespace("svglite")
  dir.create(dirname(path_prefix), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    paste0(path_prefix, ".svg"), plot = plot, width = width, height = height,
    units = "in", device = svglite::svglite, bg = "white"
  )
  ggplot2::ggsave(
    paste0(path_prefix, ".pdf"), plot = plot, width = width, height = height,
    units = "in", device = grDevices::cairo_pdf, bg = "white"
  )
}

write_session_info <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
}

format_p_value <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  sprintf("P = %.3f", p)
}
