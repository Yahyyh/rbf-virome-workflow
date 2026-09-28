#!/usr/bin/env Rscript

# Community-level virus-host concordance using the result-producing method:
# sample-wise relative abundance -> Hellinger transformation -> PCA/RDA -> symmetric Procrustes and PROTEST.
# CRISPR spacer-supported virus-host associations.
#
# Usage:
#   Rscript 02_virus_host_procrustes.R \
#     virus_abundance.tsv host_abundance.tsv metadata.tsv output_dir [prefix]

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop(
    paste(
      "Usage: 02_virus_host_procrustes.R virus_abundance.tsv",
      "host_abundance.tsv metadata.tsv output_dir [prefix]"
    ),
    call. = FALSE
  )
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path[1L]), winslash = "/"))
source(file.path(script_dir, "lib", "common.R"))

for (package in c("vegan", "permute", "ggplot2", "svglite")) {
  require_namespace(package)
}

virus_path <- args[1L]
host_path <- args[2L]
metadata_path <- args[3L]
output_dir <- args[4L]
prefix <- if (length(args) == 5L) args[5L] else "virus_host"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(RBF_RANDOM_SEED)

virus_feature <- read_feature_matrix(virus_path, "virus abundance matrix")
host_feature <- read_feature_matrix(host_path, "host abundance matrix")
metadata <- read_metadata(metadata_path, required = c("sample_id", "source", "season"))

aligned <- align_two_feature_matrices(
  virus_feature,
  host_feature,
  metadata,
  minimum_samples = 6L
)
virus <- drop_zero_features(aligned$first, "virus abundance matrix")
host <- drop_zero_features(aligned$second, "host abundance matrix")
metadata <- aligned$metadata
metadata$group <- analysis_group(metadata)

virus_hellinger <- vegan::decostand(virus, method = "hellinger", MARGIN = 1L)
host_hellinger <- vegan::decostand(host, method = "hellinger", MARGIN = 1L)

virus_pca <- vegan::rda(virus_hellinger, scale = FALSE)
host_pca <- vegan::rda(host_hellinger, scale = FALSE)

fit <- vegan::procrustes(
  X = host_pca,
  Y = virus_pca,
  symmetric = TRUE,
  scaling = 1
)
set.seed(RBF_RANDOM_SEED)
test_unrestricted <- vegan::protest(
  X = host_pca,
  Y = virus_pca,
  permutations = 999L,
  symmetric = TRUE,
  scaling = 1
)

group_control <- within_block_permutations(metadata$group, 999L)
set.seed(RBF_RANDOM_SEED)
test_within_group <- vegan::protest(
  X = host_pca,
  Y = virus_pca,
  permutations = group_control,
  symmetric = TRUE,
  scaling = 1
)

if (ncol(fit$X) < 2L || ncol(fit$Yrot) < 2L) {
  stop("The matched ordinations contain fewer than two usable dimensions.", call. = FALSE)
}

scores <- data.frame(
  sample_id = rownames(fit$X),
  source = metadata$source[match(rownames(fit$X), metadata$sample_id)],
  season = metadata$season[match(rownames(fit$X), metadata$sample_id)],
  group = metadata$group[match(rownames(fit$X), metadata$sample_id)],
  host_axis_1 = fit$X[, 1L],
  host_axis_2 = fit$X[, 2L],
  virus_axis_1_rotated = fit$Yrot[, 1L],
  virus_axis_2_rotated = fit$Yrot[, 2L],
  pointwise_residual = stats::residuals(fit),
  stringsAsFactors = FALSE
)
write_tsv(scores, file.path(output_dir, paste0(prefix, "_procrustes_scores.tsv")))

summary_table <- data.frame(
  analysis = c("unrestricted PROTEST", "within source-by-season sensitivity"),
  n_matched_samples = nrow(metadata),
  dimensions_used = ncol(fit$X),
  procrustes_correlation = c(unname(test_unrestricted$t0), unname(test_within_group$t0)),
  m2 = c(unname(test_unrestricted$ss), unname(test_within_group$ss)),
  permutations = c(test_unrestricted$permutations, test_within_group$permutations),
  p_value = c(unname(test_unrestricted$signif), unname(test_within_group$signif)),
  interpretation = c(
    "overall community concordance",
    "concordance beyond shared source-by-season membership"
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary_table, file.path(output_dir, paste0(prefix, "_protest_summary.tsv")))

plot_label <- sprintf(
  "M² = %.4f; %s",
  unname(test_unrestricted$ss),
  format_p_value(unname(test_unrestricted$signif))
)

procrustes_plot <- ggplot2::ggplot(scores) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey78", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey78", linewidth = 0.35) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = host_axis_1,
      y = host_axis_2,
      xend = virus_axis_1_rotated,
      yend = virus_axis_2_rotated,
      colour = group
    ),
    linewidth = 0.55,
    alpha = 0.68,
    arrow = grid::arrow(length = grid::unit(0.12, "cm"))
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = host_axis_1, y = host_axis_2, colour = group),
    shape = 17,
    size = 3.0,
    alpha = 0.85
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = virus_axis_1_rotated, y = virus_axis_2_rotated, colour = group),
    shape = 16,
    size = 3.0,
    alpha = 0.85
  ) +
  ggplot2::annotate(
    "text", x = -Inf, y = Inf, label = plot_label,
    hjust = -0.05, vjust = 1.25, size = 3.6
  ) +
  ggplot2::labs(
    x = "Procrustes dimension 1",
    y = "Procrustes dimension 2",
    colour = "Source:season",
    caption = "Triangles: host community; circles: rotated viral community"
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

save_plot_pair(
  procrustes_plot,
  file.path(output_dir, paste0(prefix, "_procrustes")),
  7.0,
  5.6
)

write_session_info(output_dir)
