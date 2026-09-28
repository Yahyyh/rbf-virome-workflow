#!/usr/bin/env Rscript

# Reproducible group-level feature-presence intersections (an UpSet-style summary) without inferring metadata from sample-name substrings.
#
# Usage:
#   Rscript 08_group_intersections.R \
#     abundance.tsv metadata.tsv output_dir [prefix] [detection_threshold]
#

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:5) {
  stop(
    paste(
      "Usage: 08_group_intersections.R abundance.tsv metadata.tsv output_dir",
      "[prefix] [detection_threshold]"
    ),
    call. = FALSE
  )
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path[1L]), winslash = "/"))
source(file.path(script_dir, "lib", "common.R"))

for (package in c("ggplot2", "patchwork", "svglite")) {
  require_namespace(package)
}

abundance_path <- args[1L]
metadata_path <- args[2L]
output_dir <- args[3L]
prefix <- if (length(args) >= 4L) args[4L] else "features"
detection_threshold <- if (length(args) == 5L) {
  suppressWarnings(as.numeric(args[5L]))
} else 0
if (is.na(detection_threshold) || !is.finite(detection_threshold) || detection_threshold < 0) {
  stop("detection_threshold must be a finite non-negative number.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

feature_matrix <- read_feature_matrix(abundance_path, "abundance matrix")
metadata <- read_metadata(metadata_path, required = c("sample_id", "source", "season"))
aligned <- align_feature_matrix(feature_matrix, metadata, minimum_samples = 4L)
abundance <- aligned$abundance
metadata <- aligned$metadata
metadata$group <- factor(
  paste(metadata$source, metadata$season, sep = ":"),
  levels = unique(paste(metadata$source, metadata$season, sep = ":"))
)
group_levels <- levels(metadata$group)
if (length(group_levels) < 2L) {
  stop("At least two source:season groups are required.", call. = FALSE)
}

presence <- vapply(group_levels, function(group_value) {
  colSums(abundance[metadata$group == group_value, , drop = FALSE] > detection_threshold) > 0
}, logical(ncol(abundance)))
if (is.null(dim(presence))) {
  presence <- matrix(presence, ncol = 1L, dimnames = list(colnames(abundance), group_levels))
}
rownames(presence) <- colnames(abundance)
colnames(presence) <- group_levels
present_anywhere <- rowSums(presence) > 0
presence <- presence[present_anywhere, , drop = FALSE]
if (nrow(presence) == 0L) {
  stop("No feature exceeds the detection threshold in any group.", call. = FALSE)
}

intersection_label <- apply(presence, 1L, function(row) {
  paste(group_levels[row], collapse = " & ")
})
assignments <- data.frame(
  feature_id = rownames(presence),
  as.data.frame(presence, check.names = FALSE),
  intersection = intersection_label,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(assignments, file.path(output_dir, paste0(prefix, "_intersection_membership.tsv")))

summary_table <- as.data.frame(table(intersection_label), stringsAsFactors = FALSE)
names(summary_table) <- c("intersection", "n_features")
summary_table <- summary_table[order(summary_table$n_features, decreasing = TRUE), , drop = FALSE]
summary_table$rank <- seq_len(nrow(summary_table))
summary_table$pattern_id <- sprintf("I%02d", summary_table$rank)
summary_table$detection_rule <- sprintf("abundance > %s in at least one group sample", detection_threshold)
write_tsv(summary_table, file.path(output_dir, paste0(prefix, "_intersection_summary.tsv")))

top_summary <- head(summary_table, 20L)
top_summary$pattern_id <- factor(top_summary$pattern_id, levels = top_summary$pattern_id)
bar_plot <- ggplot2::ggplot(
  top_summary,
  ggplot2::aes(x = pattern_id, y = n_features)
) +
  ggplot2::geom_col(fill = "#4C78A8", width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = n_features),
    vjust = -0.25,
    size = 3
  ) +
  ggplot2::labs(x = NULL, y = "Number of features") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank()
  )

matrix_rows <- list()
row_index <- 1L
for (summary_index in seq_len(nrow(top_summary))) {
  active_groups <- strsplit(
    as.character(top_summary$intersection[summary_index]),
    " & ",
    fixed = TRUE
  )[[1L]]
  matrix_rows[[row_index]] <- data.frame(
    pattern_id = top_summary$pattern_id[summary_index],
    group = factor(group_levels, levels = rev(group_levels)),
    present = group_levels %in% active_groups,
    stringsAsFactors = FALSE
  )
  row_index <- row_index + 1L
}
matrix_data <- do.call(rbind, matrix_rows)
matrix_data$pattern_id <- factor(matrix_data$pattern_id, levels = levels(top_summary$pattern_id))

matrix_plot <- ggplot2::ggplot(
  matrix_data,
  ggplot2::aes(x = pattern_id, y = group)
) +
  ggplot2::geom_point(colour = "grey82", size = 2.8) +
  ggplot2::geom_point(
    data = matrix_data[matrix_data$present, , drop = FALSE],
    colour = "#1F1F1F",
    size = 3.2
  ) +
  ggplot2::labs(x = "Intersection pattern", y = NULL) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

combined_plot <- patchwork::wrap_plots(
  bar_plot,
  matrix_plot,
  ncol = 1,
  heights = c(2.2, 1.2)
)
save_plot_pair(
  combined_plot,
  file.path(output_dir, paste0(prefix, "_group_intersections")),
  max(8.0, 0.45 * nrow(top_summary) + 3.0),
  7.0
)

write_session_info(output_dir)
