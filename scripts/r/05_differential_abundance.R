#!/usr/bin/env Rscript

# Planned feature-wise differential-abundance screening for the paired RBF design.
#
# Usage:
#   Rscript 05_differential_abundance.R \
#     abundance.tsv metadata.tsv output_dir [analysis_prefix] [min_prevalence]
#
# abundance.tsv uses feature IDs in the first column and sample IDs thereafter.
# metadata.tsv requires sample_id, source, season, and pair_id. Source contrasts are paired within sampling date and season contrasts are unpaired across dates.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:5) {
  stop(
    paste(
      "Usage: 05_differential_abundance.R abundance.tsv metadata.tsv",
      "output_dir [analysis_prefix] [min_prevalence]"
    ),
    call. = FALSE
  )
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path[1L]), winslash = "/"))
source(file.path(script_dir, "lib", "common.R"))

for (package in c("ggplot2", "svglite")) {
  require_namespace(package)
}

abundance_path <- args[1L]
metadata_path <- args[2L]
output_dir <- args[3L]
prefix <- if (length(args) >= 4L) args[4L] else "features"
min_prevalence <- if (length(args) == 5L) suppressWarnings(as.numeric(args[5L])) else 0.20
if (is.na(min_prevalence) || min_prevalence < 0 || min_prevalence > 1) {
  stop("min_prevalence must be a number from 0 to 1.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

feature_matrix <- read_feature_matrix(abundance_path, "abundance matrix")
metadata <- read_metadata(
  metadata_path,
  required = c("sample_id", "source", "season", "pair_id")
)
aligned <- align_feature_matrix(feature_matrix, metadata, minimum_samples = 8L)
abundance <- drop_zero_features(aligned$abundance, "abundance matrix")
metadata <- aligned$metadata
relative <- to_relative_abundance(abundance)

metadata$source <- factor(metadata$source, levels = unique(metadata$source))
metadata$season <- factor(metadata$season, levels = unique(metadata$season))
metadata$pair_id <- factor(metadata$pair_id)
if (nlevels(metadata$source) != 2L || nlevels(metadata$season) != 2L) {
  stop("This study-specific script requires exactly two source and two season levels.", call. = FALSE)
}
if (any(table(metadata$pair_id) != 2L)) {
  stop("Each pair_id must contain exactly one sample from each source.", call. = FALSE)
}
if (any(vapply(split(metadata$source, metadata$pair_id), function(x) length(unique(x)) != 2L, logical(1)))) {
  stop("Each pair_id must contain both source levels.", call. = FALSE)
}

safe_wilcox <- function(first, second, paired) {
  if (paired && all(first - second == 0)) {
    return(c(statistic = 0, p_value = 1))
  }
  if (!paired && length(unique(c(first, second))) == 1L) {
    return(c(statistic = 0, p_value = 1))
  }
  test <- tryCatch(
    stats::wilcox.test(first, second, paired = paired, exact = FALSE),
    error = function(error) NULL
  )
  if (is.null(test)) return(c(statistic = NA_real_, p_value = NA_real_))
  c(statistic = unname(test$statistic), p_value = test$p.value)
}

analyse_contrast <- function(first, second, comparison, paired, first_label, second_label) {
  if (nrow(first) < 2L || nrow(second) < 2L) {
    stop(sprintf("Fewer than two observations per level for: %s", comparison), call. = FALSE)
  }
  if (paired && nrow(first) != nrow(second)) {
    stop(sprintf("Paired contrast has unequal group sizes: %s", comparison), call. = FALSE)
  }
  pooled <- rbind(first, second)
  prevalence <- colMeans(pooled > 0)
  keep <- prevalence >= min_prevalence
  if (!any(keep)) {
    stop(sprintf("No features pass the prevalence filter for: %s", comparison), call. = FALSE)
  }
  first <- first[, keep, drop = FALSE]
  second <- second[, keep, drop = FALSE]
  feature_ids <- colnames(first)

  rows <- lapply(seq_along(feature_ids), function(index) {
    first_values <- first[, index]
    second_values <- second[, index]
    positive <- c(first_values[first_values > 0], second_values[second_values > 0])
    pseudocount <- if (length(positive) > 0L) min(positive) / 2 else 1e-12
    effect <- if (paired) {
      stats::median(log2((first_values + pseudocount) / (second_values + pseudocount)))
    } else {
      log2(
        (stats::median(first_values) + pseudocount) /
          (stats::median(second_values) + pseudocount)
      )
    }
    test <- safe_wilcox(first_values, second_values, paired)
    data.frame(
      comparison = comparison,
      first_level = first_label,
      second_level = second_label,
      test = if (paired) "paired Wilcoxon signed-rank" else "Wilcoxon rank-sum",
      feature_id = feature_ids[index],
      n_first = length(first_values),
      n_second = length(second_values),
      prevalence_first = mean(first_values > 0),
      prevalence_second = mean(second_values > 0),
      median_first = stats::median(first_values),
      median_second = stats::median(second_values),
      effect_log2_ratio_first_over_second = effect,
      pseudocount = pseudocount,
      statistic = test[["statistic"]],
      p_value = test[["p_value"]],
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$p_adjust_bh <- stats::p.adjust(result$p_value, method = "BH")
  result$significant_bh_0_05 <- !is.na(result$p_adjust_bh) & result$p_adjust_bh < 0.05
  result
}

results <- list()
result_index <- 1L
source_levels <- levels(metadata$source)
season_levels <- levels(metadata$season)

for (season_value in season_levels) {
  subset_metadata <- metadata[metadata$season == season_value, , drop = FALSE]
  first_map <- subset_metadata[
    subset_metadata$source == source_levels[1L],
    c("pair_id", "sample_id")
  ]
  second_map <- subset_metadata[
    subset_metadata$source == source_levels[2L],
    c("pair_id", "sample_id")
  ]
  names(first_map)[2L] <- "first_sample"
  names(second_map)[2L] <- "second_sample"
  paired_map <- merge(first_map, second_map, by = "pair_id", sort = TRUE)
  first <- relative[match(paired_map$first_sample, rownames(relative)), , drop = FALSE]
  second <- relative[match(paired_map$second_sample, rownames(relative)), , drop = FALSE]
  results[[result_index]] <- analyse_contrast(
    first,
    second,
    sprintf("%s vs %s within %s", source_levels[1L], source_levels[2L], season_value),
    TRUE,
    source_levels[1L],
    source_levels[2L]
  )
  result_index <- result_index + 1L
}

for (source_value in source_levels) {
  first_ids <- metadata$sample_id[
    metadata$source == source_value & metadata$season == season_levels[1L]
  ]
  second_ids <- metadata$sample_id[
    metadata$source == source_value & metadata$season == season_levels[2L]
  ]
  first <- relative[match(first_ids, rownames(relative)), , drop = FALSE]
  second <- relative[match(second_ids, rownames(relative)), , drop = FALSE]
  results[[result_index]] <- analyse_contrast(
    first,
    second,
    sprintf("%s vs %s within %s", season_levels[1L], season_levels[2L], source_value),
    FALSE,
    season_levels[1L],
    season_levels[2L]
  )
  result_index <- result_index + 1L
}

all_results <- do.call(rbind, results)
write_tsv(
  all_results,
  file.path(output_dir, paste0(prefix, "_planned_differential_abundance.tsv"))
)

plot_rows <- do.call(rbind, lapply(split(all_results, all_results$comparison), function(table) {
  significant <- table[table$significant_bh_0_05, , drop = FALSE]
  selected <- if (nrow(significant) > 0L) significant else table
  selected <- selected[order(selected$p_adjust_bh, -abs(selected$effect_log2_ratio_first_over_second)), , drop = FALSE]
  head(selected, 20L)
}))
plot_rows$feature_label <- paste(plot_rows$comparison, plot_rows$feature_id, sep = " | ")
plot_rows$feature_label <- factor(
  plot_rows$feature_label,
  levels = plot_rows$feature_label[order(plot_rows$effect_log2_ratio_first_over_second)]
)
plot_rows$minus_log10_q <- -log10(pmax(plot_rows$p_adjust_bh, .Machine$double.xmin))
plot_rows$max_prevalence <- pmax(plot_rows$prevalence_first, plot_rows$prevalence_second)

differential_plot <- ggplot2::ggplot(
  plot_rows,
  ggplot2::aes(
    x = effect_log2_ratio_first_over_second,
    y = feature_label,
    size = max_prevalence,
    fill = minus_log10_q
  )
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.45) +
  ggplot2::geom_point(shape = 21, colour = "grey20", stroke = 0.25, alpha = 0.88) +
  ggplot2::scale_fill_viridis_c(option = "B") +
  ggplot2::scale_size_continuous(range = c(2, 7), limits = c(0, 1)) +
  ggplot2::labs(
    x = "Median log2 ratio (first / second)",
    y = NULL,
    size = "Maximum prevalence",
    fill = "-log10(BH q)",
    caption = paste(
      "BH correction is applied within each planned contrast.",
      "If no feature passes q < 0.05, the 20 lowest-q features are shown."
    )
  ) +
  ggplot2::theme_bw(base_size = 9.5) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
save_plot_pair(
  differential_plot,
  file.path(output_dir, paste0(prefix, "_differential_abundance")),
  10.0,
  max(6.0, 0.22 * nrow(plot_rows) + 2.0)
)

write_session_info(output_dir)
