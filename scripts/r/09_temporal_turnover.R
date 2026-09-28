#!/usr/bin/env Rscript

# Time-decay analysis using actual sampling dates and one upper-triangle distance per sample pair. 
# Mantel permutations test association between Bray-Curtis dissimilarity and temporal distance within each source.
#
# Usage:
#   Rscript 09_temporal_turnover.R \
#     abundance.tsv metadata.tsv output_dir [analysis_prefix]
#
# metadata.tsv requires sample_id, source, and sampling_date (YYYY-MM-DD).

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop(
    paste(
      "Usage: 09_temporal_turnover.R abundance.tsv metadata.tsv",
      "output_dir [analysis_prefix]"
    ),
    call. = FALSE
  )
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path[1L]), winslash = "/"))
source(file.path(script_dir, "lib", "common.R"))

for (package in c("vegan", "ggplot2", "svglite")) {
  require_namespace(package)
}

abundance_path <- args[1L]
metadata_path <- args[2L]
output_dir <- args[3L]
prefix <- if (length(args) == 4L) args[4L] else "virome"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(RBF_RANDOM_SEED)

feature_matrix <- read_feature_matrix(abundance_path, "abundance matrix")
metadata <- read_metadata(
  metadata_path,
  required = c("sample_id", "source", "sampling_date")
)
metadata$sampling_date <- as.Date(as.character(metadata$sampling_date), format = "%Y-%m-%d")
if (anyNA(metadata$sampling_date)) {
  stop("Every sampling_date must be a valid YYYY-MM-DD date.", call. = FALSE)
}
aligned <- align_feature_matrix(feature_matrix, metadata, minimum_samples = 8L)
abundance <- drop_zero_features(aligned$abundance, "abundance matrix")
metadata <- aligned$metadata
relative <- to_relative_abundance(abundance)

sources <- unique(metadata$source)
pair_rows <- list()
test_rows <- list()
pair_index <- 1L
test_index <- 1L

for (source_value in sources) {
  keep <- metadata$source == source_value
  source_metadata <- metadata[keep, , drop = FALSE]
  source_relative <- relative[keep, , drop = FALSE]
  if (nrow(source_metadata) < 4L) {
    warning(
      sprintf("Skipping source '%s': fewer than four dated samples.", source_value),
      call. = FALSE
    )
    next
  }
  order_index <- order(source_metadata$sampling_date, source_metadata$sample_id)
  source_metadata <- source_metadata[order_index, , drop = FALSE]
  source_relative <- source_relative[order_index, , drop = FALSE]

  bray <- vegan::vegdist(source_relative, method = "bray")
  temporal <- stats::dist(as.numeric(source_metadata$sampling_date), method = "euclidean")
  bray_matrix <- as.matrix(bray)
  temporal_matrix <- as.matrix(temporal)
  upper <- which(upper.tri(bray_matrix), arr.ind = TRUE)

  pair_rows[[pair_index]] <- data.frame(
    source = source_value,
    sample_1 = source_metadata$sample_id[upper[, 1L]],
    sample_2 = source_metadata$sample_id[upper[, 2L]],
    date_1 = source_metadata$sampling_date[upper[, 1L]],
    date_2 = source_metadata$sampling_date[upper[, 2L]],
    time_lag_days = temporal_matrix[upper],
    bray_curtis = bray_matrix[upper],
    stringsAsFactors = FALSE
  )
  pair_index <- pair_index + 1L

  set.seed(RBF_RANDOM_SEED)
  mantel_test <- vegan::mantel(
    bray,
    temporal,
    method = "spearman",
    permutations = 9999L,
    na.rm = FALSE
  )
  test_rows[[test_index]] <- data.frame(
    source = source_value,
    n_samples = nrow(source_metadata),
    n_unique_pairs = nrow(pair_rows[[pair_index - 1L]]),
    method = "Spearman Mantel test",
    alternative = "greater (community dissimilarity increases with time lag)",
    statistic = unname(mantel_test$statistic),
    permutations = mantel_test$permutations,
    p_value = mantel_test$signif,
    stringsAsFactors = FALSE
  )
  test_index <- test_index + 1L
}

if (length(pair_rows) == 0L) {
  stop("No source contains enough dated samples for temporal-turnover analysis.", call. = FALSE)
}
pair_table <- do.call(rbind, pair_rows)
test_table <- do.call(rbind, test_rows)
test_table$p_adjust_bh <- stats::p.adjust(test_table$p_value, method = "BH")
write_tsv(pair_table, file.path(output_dir, paste0(prefix, "_temporal_pair_distances.tsv")))
write_tsv(test_table, file.path(output_dir, paste0(prefix, "_temporal_mantel_tests.tsv")))

caption <- paste(vapply(seq_len(nrow(test_table)), function(index) {
  sprintf(
    "%s: rho = %.3f, BH q = %.3g",
    test_table$source[index],
    test_table$statistic[index],
    test_table$p_adjust_bh[index]
  )
}, character(1)), collapse = "; ")

temporal_plot <- ggplot2::ggplot(
  pair_table,
  ggplot2::aes(x = time_lag_days, y = bray_curtis, colour = source)
) +
  ggplot2::geom_point(size = 2.1, alpha = 0.58) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  ggplot2::facet_wrap(~source) +
  ggplot2::labs(
    x = "Temporal separation (days)",
    y = "Bray-Curtis dissimilarity",
    colour = "Source",
    caption = paste(
      caption,
      "Points are unique sample pairs; the fitted line is descriptive and inference uses Mantel permutations."
    )
  ) +
  ggplot2::theme_bw(base_size = 10.5) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "none")
save_plot_pair(
  temporal_plot,
  file.path(output_dir, paste0(prefix, "_temporal_turnover")),
  8.5,
  5.2
)

write_session_info(output_dir)
