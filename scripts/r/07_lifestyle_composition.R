#!/usr/bin/env Rscript

# Descriptive abundance and within-lifestyle diversity for predicted virulent and temperate vOTUs.
#
# Usage:
#   Rscript 07_lifestyle_composition.R \
#     abundance.tsv virus_annotations.tsv metadata.tsv output_dir [prefix]
#
# virus_annotations.tsv requires virus_id, lifestyle, and taxonomy.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop(
    paste(
      "Usage: 07_lifestyle_composition.R abundance.tsv virus_annotations.tsv",
      "metadata.tsv output_dir [prefix]"
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
annotation_path <- args[2L]
metadata_path <- args[3L]
output_dir <- args[4L]
prefix <- if (length(args) == 5L) args[5L] else "lifestyle"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

feature_matrix <- read_feature_matrix(abundance_path, "viral abundance matrix")
metadata <- read_metadata(
  metadata_path,
  required = c("sample_id", "source", "season", "pair_id")
)
aligned <- align_feature_matrix(feature_matrix, metadata, minimum_samples = 8L)
abundance <- drop_zero_features(aligned$abundance, "viral abundance matrix")
metadata <- aligned$metadata
relative <- to_relative_abundance(abundance)

annotations <- read_tsv(annotation_path)
assert_columns(annotations, c("virus_id", "lifestyle", "taxonomy"), "virus annotation table")
annotations <- annotations[c("virus_id", "lifestyle", "taxonomy")]
annotations[] <- lapply(annotations, function(column) trimws(as.character(column)))
assert_unique_nonmissing(annotations$virus_id, "virus annotation virus_id")
annotations$taxonomy[is.na(annotations$taxonomy) | annotations$taxonomy == ""] <- "Unclassified"

lifestyle_key <- tolower(annotations$lifestyle)
lifestyle_key[lifestyle_key %in% c("lytic", "virulent")] <- "virulent"
lifestyle_key[lifestyle_key %in% c("lysogenic", "temperate")] <- "temperate"
annotations$lifestyle <- lifestyle_key
annotations <- annotations[
  annotations$lifestyle %in% c("virulent", "temperate"),
  ,
  drop = FALSE
]
if (!all(c("virulent", "temperate") %in% annotations$lifestyle)) {
  stop("Both virulent and temperate predicted lifestyle classes are required.", call. = FALSE)
}
shared_features <- intersect(colnames(relative), annotations$virus_id)
if (length(shared_features) < 2L) {
  stop("Fewer than two vOTUs have both abundance and retained lifestyle annotation.", call. = FALSE)
}
excluded_votus <- setdiff(colnames(relative), annotations$virus_id)
write_tsv(
  data.frame(
    virus_id = excluded_votus,
    exclusion_reason = rep(
      "missing retained virulent/temperate lifestyle annotation",
      length(excluded_votus)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, paste0(prefix, "_excluded_votus.tsv"))
)
relative <- relative[, shared_features, drop = FALSE]
relative <- to_relative_abundance(relative)
annotations <- annotations[match(shared_features, annotations$virus_id), , drop = FALSE]
if (!all(c("virulent", "temperate") %in% annotations$lifestyle)) {
  stop("Both lifestyle classes must have at least one vOTU in the abundance matrix.", call. = FALSE)
}

long <- as.data.frame(as.table(relative), stringsAsFactors = FALSE)
names(long) <- c("sample_id", "virus_id", "relative_abundance")
long <- merge(long, annotations, by = "virus_id", all.x = TRUE, sort = FALSE)
long <- merge(
  long,
  metadata[c("sample_id", "source", "season", "pair_id")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
long$lifestyle <- factor(long$lifestyle, levels = c("virulent", "temperate"))

lifestyle_abundance <- aggregate(
  relative_abundance ~ sample_id + source + season + pair_id + lifestyle,
  data = long,
  FUN = sum
)
sample_grid <- expand.grid(
  sample_id = metadata$sample_id,
  lifestyle = c("virulent", "temperate"),
  stringsAsFactors = FALSE
)
sample_grid <- merge(
  sample_grid,
  metadata[c("sample_id", "source", "season", "pair_id")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
lifestyle_abundance <- merge(
  sample_grid,
  lifestyle_abundance,
  by = c("sample_id", "source", "season", "pair_id", "lifestyle"),
  all.x = TRUE,
  sort = FALSE
)
lifestyle_abundance$relative_abundance[is.na(lifestyle_abundance$relative_abundance)] <- 0
lifestyle_abundance$lifestyle <- factor(
  lifestyle_abundance$lifestyle,
  levels = c("virulent", "temperate")
)
lifestyle_abundance$percent <- 100 * lifestyle_abundance$relative_abundance
write_tsv(
  lifestyle_abundance,
  file.path(output_dir, paste0(prefix, "_sample_relative_abundance.tsv"))
)

sample_order <- metadata$sample_id
lifestyle_abundance$sample_id <- factor(lifestyle_abundance$sample_id, levels = sample_order)
lifestyle_plot <- ggplot2::ggplot(
  lifestyle_abundance,
  ggplot2::aes(x = sample_id, y = percent, fill = lifestyle)
) +
  ggplot2::geom_col(width = 0.82) +
  ggplot2::facet_grid(~source + season, scales = "free_x", space = "free_x") +
  ggplot2::labs(
    x = "Sample", y = "Annotated viral community (%)",
    fill = "Predicted lifestyle",
    caption = "Percentages are calculated after restricting to vOTUs with a retained lifestyle prediction."
  ) +
  ggplot2::theme_bw(base_size = 9.5) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )
save_plot_pair(
  lifestyle_plot,
  file.path(output_dir, paste0(prefix, "_relative_abundance")),
  10.5,
  4.8
)

taxonomy_abundance <- aggregate(
  relative_abundance ~ sample_id + source + season + lifestyle + taxonomy,
  data = long,
  FUN = sum
)
lifestyle_total <- aggregate(
  relative_abundance ~ sample_id + lifestyle,
  data = taxonomy_abundance,
  FUN = sum
)
names(lifestyle_total)[3L] <- "lifestyle_total"
taxonomy_abundance <- merge(
  taxonomy_abundance,
  lifestyle_total,
  by = c("sample_id", "lifestyle"),
  all.x = TRUE,
  sort = FALSE
)
taxonomy_abundance$within_lifestyle_percent <- ifelse(
  taxonomy_abundance$lifestyle_total > 0,
  100 * taxonomy_abundance$relative_abundance / taxonomy_abundance$lifestyle_total,
  NA_real_
)
write_tsv(
  taxonomy_abundance,
  file.path(output_dir, paste0(prefix, "_taxonomy_composition.tsv"))
)

taxonomy_totals <- aggregate(relative_abundance ~ taxonomy, taxonomy_abundance, sum)
taxonomy_totals <- taxonomy_totals[order(taxonomy_totals$relative_abundance, decreasing = TRUE), ]
top_taxa <- head(taxonomy_totals$taxonomy, 12L)
taxonomy_plot_data <- taxonomy_abundance
taxonomy_plot_data$taxonomy_plot <- ifelse(
  taxonomy_plot_data$taxonomy %in% top_taxa,
  taxonomy_plot_data$taxonomy,
  "Other"
)
taxonomy_plot_data <- aggregate(
  within_lifestyle_percent ~ sample_id + source + season + lifestyle + taxonomy_plot,
  data = taxonomy_plot_data,
  FUN = sum
)
taxonomy_plot_data$sample_id <- factor(taxonomy_plot_data$sample_id, levels = sample_order)
taxonomy_plot <- ggplot2::ggplot(
  taxonomy_plot_data,
  ggplot2::aes(x = sample_id, y = within_lifestyle_percent, fill = taxonomy_plot)
) +
  ggplot2::geom_col(width = 0.82) +
  ggplot2::facet_grid(lifestyle ~ source + season, scales = "free_x", space = "free_x") +
  ggplot2::labs(
    x = "Sample", y = "Within-lifestyle composition (%)", fill = "Taxonomy"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )
save_plot_pair(
  taxonomy_plot,
  file.path(output_dir, paste0(prefix, "_taxonomy_composition")),
  11.5,
  6.8
)

# Diversity is computed across vOTUs within each predicted lifestyle. A sample with no detected member of a lifestyle receives NA rather than a fabricated 0.
diversity_rows <- list()
row_index <- 1L
for (sample_id in rownames(relative)) {
  for (lifestyle_value in c("virulent", "temperate")) {
    feature_ids <- annotations$virus_id[annotations$lifestyle == lifestyle_value]
    values <- relative[sample_id, feature_ids, drop = TRUE]
    observed <- sum(values > 0)
    total <- sum(values)
    shannon <- if (total > 0) vegan::diversity(values, index = "shannon") else NA_real_
    pielou <- if (!is.na(shannon) && observed > 1L) shannon / log(observed) else NA_real_
    diversity_rows[[row_index]] <- data.frame(
      sample_id = sample_id,
      lifestyle = lifestyle_value,
      observed_votus = observed,
      shannon = shannon,
      pielou = pielou,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
diversity <- do.call(rbind, diversity_rows)
diversity <- merge(
  diversity,
  metadata[c("sample_id", "source", "season", "pair_id")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
write_tsv(diversity, file.path(output_dir, paste0(prefix, "_within_class_diversity.tsv")))

diversity_tests <- list()
test_index <- 1L
for (lifestyle_value in c("virulent", "temperate")) {
  lifestyle_data <- diversity[diversity$lifestyle == lifestyle_value & !is.na(diversity$shannon), ]
  for (season_value in unique(lifestyle_data$season)) {
    subset_data <- lifestyle_data[lifestyle_data$season == season_value, ]
    first <- subset_data[
      subset_data$source == unique(metadata$source)[1L],
      c("pair_id", "shannon")
    ]
    second <- subset_data[
      subset_data$source == unique(metadata$source)[2L],
      c("pair_id", "shannon")
    ]
    names(first)[2L] <- "first"
    names(second)[2L] <- "second"
    paired <- merge(first, second, by = "pair_id")
    test <- if (nrow(paired) >= 2L && any(paired$first != paired$second)) {
      stats::wilcox.test(paired$first, paired$second, paired = TRUE, exact = FALSE)
    } else if (nrow(paired) >= 2L) {
      list(statistic = 0, p.value = 1)
    } else NULL
    diversity_tests[[test_index]] <- data.frame(
      lifestyle = lifestyle_value,
      comparison = sprintf("source within %s", season_value),
      test = "paired Wilcoxon signed-rank",
      n_first = nrow(paired),
      n_second = nrow(paired),
      median_difference = if (nrow(paired)) stats::median(paired$first - paired$second) else NA_real_,
      statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
      p_value = if (is.null(test)) NA_real_ else test$p.value,
      stringsAsFactors = FALSE
    )
    test_index <- test_index + 1L
  }

  for (source_value in unique(lifestyle_data$source)) {
    first <- lifestyle_data$shannon[
      lifestyle_data$source == source_value & lifestyle_data$season == unique(metadata$season)[1L]
    ]
    second <- lifestyle_data$shannon[
      lifestyle_data$source == source_value & lifestyle_data$season == unique(metadata$season)[2L]
    ]
    test <- if (length(first) >= 2L && length(second) >= 2L && length(unique(c(first, second))) > 1L) {
      stats::wilcox.test(first, second, paired = FALSE, exact = FALSE)
    } else if (length(first) >= 2L && length(second) >= 2L) {
      list(statistic = 0, p.value = 1)
    } else NULL
    diversity_tests[[test_index]] <- data.frame(
      lifestyle = lifestyle_value,
      comparison = sprintf("season within %s", source_value),
      test = "Wilcoxon rank-sum",
      n_first = length(first),
      n_second = length(second),
      median_difference = if (length(first) && length(second)) {
        stats::median(first) - stats::median(second)
      } else NA_real_,
      statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
      p_value = if (is.null(test)) NA_real_ else test$p.value,
      stringsAsFactors = FALSE
    )
    test_index <- test_index + 1L
  }
}
diversity_tests <- do.call(rbind, diversity_tests)
diversity_tests$p_adjust_bh <- stats::p.adjust(diversity_tests$p_value, method = "BH")
write_tsv(
  diversity_tests,
  file.path(output_dir, paste0(prefix, "_within_class_shannon_tests.tsv"))
)

diversity$group <- factor(
  paste(diversity$source, diversity$season, sep = ":"),
  levels = unique(paste(metadata$source, metadata$season, sep = ":"))
)
diversity$lifestyle <- factor(diversity$lifestyle, levels = c("virulent", "temperate"))
diversity_plot <- ggplot2::ggplot(
  diversity,
  ggplot2::aes(x = group, y = shannon, fill = group)
) +
  ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  ggplot2::geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.72) +
  ggplot2::facet_wrap(~lifestyle, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "Within-lifestyle Shannon diversity") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    panel.grid = ggplot2::element_blank(),
    legend.position = "none"
  )
save_plot_pair(
  diversity_plot,
  file.path(output_dir, paste0(prefix, "_within_class_shannon")),
  8.0,
  5.0
)

write_session_info(output_dir)
