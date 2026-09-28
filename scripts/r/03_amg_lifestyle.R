#!/usr/bin/env Rscript

# AMG abundance and genomic-category association by predicted viral lifestyle.
#
# Usage:
#   Rscript 03_amg_lifestyle.R \
#     amg_abundance_long.tsv metadata.tsv output_dir [analysis_prefix]
#
# amg_abundance_long.tsv requires:
#   virus_id       stable vOTU identifier
#   amg_id         unique AMG/ORF identifier (not merely a repeated KO label)
#   amg_category   functional category
#   lifestyle      predicted lifestyle; virulent/lytic or temperate/lysogenic
#   sample_id      sample identifier
#   tpm            non-negative abundance assigned to the AMG-bearing vOTU
#
# Each virus_id/amg_id/sample_id combination must occur once.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop(
    paste(
      "Usage: 03_amg_lifestyle.R amg_abundance_long.tsv",
      "metadata.tsv output_dir [analysis_prefix]"
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

amg_path <- args[1L]
metadata_path <- args[2L]
output_dir <- args[3L]
prefix <- if (length(args) == 4L) args[4L] else "amg"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(RBF_RANDOM_SEED)

amg <- read_tsv(amg_path)
assert_columns(
  amg,
  c("virus_id", "amg_id", "amg_category", "lifestyle", "sample_id", "tpm"),
  "AMG abundance table"
)

text_columns <- c("virus_id", "amg_id", "amg_category", "lifestyle", "sample_id")
amg[text_columns] <- lapply(amg[text_columns], function(column) trimws(as.character(column)))
if (anyNA(amg[text_columns]) || any(as.matrix(amg[text_columns]) == "")) {
  stop("AMG identifier, category, lifestyle, and sample columns cannot be missing or empty.", call. = FALSE)
}
amg$tpm <- suppressWarnings(as.numeric(amg$tpm))
if (anyNA(amg$tpm) || any(!is.finite(amg$tpm)) || any(amg$tpm < 0)) {
  stop("AMG TPM values must be finite, numeric, non-missing, and non-negative.", call. = FALSE)
}

lifestyle_key <- tolower(amg$lifestyle)
lifestyle_key[lifestyle_key %in% c("lytic", "virulent")] <- "virulent"
lifestyle_key[lifestyle_key %in% c("lysogenic", "temperate")] <- "temperate"
keep_lifestyle <- lifestyle_key %in% c("virulent", "temperate")
if (any(!keep_lifestyle)) {
  warning(
    sprintf(
      "Excluded %d row(s) whose lifestyle was not virulent/lytic or temperate/lysogenic.",
      sum(!keep_lifestyle)
    ),
    call. = FALSE
  )
}
amg <- amg[keep_lifestyle, , drop = FALSE]
amg$lifestyle <- factor(lifestyle_key[keep_lifestyle], levels = c("virulent", "temperate"))
if (nrow(amg) == 0L || any(table(amg$lifestyle) == 0L)) {
  stop("Both virulent and temperate predicted lifestyle classes are required.", call. = FALSE)
}

key <- paste(amg$virus_id, amg$amg_id, amg$sample_id, sep = "\r")
if (anyDuplicated(key)) {
  stop(
    "AMG table contains duplicate virus_id/amg_id/sample_id rows; resolve the join before analysis.",
    call. = FALSE
  )
}

virus_lifestyle <- unique(amg[c("virus_id", "lifestyle")])
if (anyDuplicated(virus_lifestyle$virus_id)) {
  stop("At least one virus_id has conflicting lifestyle predictions.", call. = FALSE)
}
amg_definition <- unique(amg[c("virus_id", "amg_id", "amg_category")])
if (anyDuplicated(amg_definition[c("virus_id", "amg_id")])) {
  stop("At least one virus_id/amg_id maps to multiple AMG categories.", call. = FALSE)
}

metadata <- read_metadata(
  metadata_path,
  required = c("sample_id", "source", "season")
)
unknown_samples <- setdiff(unique(amg$sample_id), metadata$sample_id)
if (length(unknown_samples) > 0L) {
  stop(
    sprintf(
      "AMG table contains sample(s) absent from metadata: %s",
      paste(unknown_samples, collapse = ", ")
    ),
    call. = FALSE
  )
}
missing_samples <- setdiff(metadata$sample_id, unique(amg$sample_id))
if (length(missing_samples) > 0L) {
  stop(
    paste(
      "Metadata contains sample(s) absent from the AMG table.",
      "Pass metadata filtered to this analysis and include explicit zero-TPM rows",
      "for samples with no detected AMG-bearing vOTU. First missing sample(s):",
      paste(head(missing_samples, 5L), collapse = ", ")
    ),
    call. = FALSE
  )
}
metadata$source <- factor(metadata$source, levels = unique(metadata$source))
metadata$season <- factor(metadata$season, levels = unique(metadata$season))
amg <- merge(
  amg,
  metadata[c("sample_id", "source", "season")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

# Genomic load is counted once per vOTU/AMG, never once per sample. With this
# input schema, viruses carrying zero retained AMGs are not represented, so the
# output is explicitly labelled as an AMG-positive-virus analysis.
amg_catalog <- unique(amg[c("virus_id", "amg_id", "amg_category", "lifestyle")])
virus_load <- aggregate(
  amg_id ~ virus_id + lifestyle,
  data = amg_catalog,
  FUN = function(x) length(unique(x))
)
names(virus_load)[names(virus_load) == "amg_id"] <- "n_distinct_amgs"
write_tsv(
  virus_load,
  file.path(output_dir, paste0(prefix, "_amg_positive_virus_load.tsv"))
)

load_split <- split(virus_load$n_distinct_amgs, virus_load$lifestyle)
if (length(load_split) != 2L || any(lengths(load_split) < 2L)) {
  stop("At least two AMG-positive viruses per lifestyle are required for the load test.", call. = FALSE)
}
if (length(unique(unlist(load_split, use.names = FALSE))) == 1L) {
  load_test <- list(statistic = 0, p.value = 1)
} else {
  load_test <- stats::wilcox.test(
    load_split[["virulent"]],
    load_split[["temperate"]],
    paired = FALSE,
    exact = FALSE
  )
}
load_summary <- data.frame(
  comparison = "virulent versus temperate among AMG-positive vOTUs",
  test = "Wilcoxon rank-sum",
  n_virulent = length(load_split[["virulent"]]),
  n_temperate = length(load_split[["temperate"]]),
  median_virulent = stats::median(load_split[["virulent"]]),
  median_temperate = stats::median(load_split[["temperate"]]),
  median_difference = stats::median(load_split[["virulent"]]) -
    stats::median(load_split[["temperate"]]),
  statistic = unname(load_test$statistic),
  p_value = load_test$p.value,
  stringsAsFactors = FALSE
)
write_tsv(load_summary, file.path(output_dir, paste0(prefix, "_virus_load_test.tsv")))

load_plot <- ggplot2::ggplot(
  virus_load,
  ggplot2::aes(x = lifestyle, y = n_distinct_amgs, fill = lifestyle)
) +
  ggplot2::geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.78) +
  ggplot2::geom_jitter(width = 0.08, height = 0, size = 2.1, alpha = 0.72) +
  ggplot2::annotate(
    "text", x = 1.5, y = Inf,
    label = format_p_value(load_test$p.value),
    vjust = 1.35, size = 3.5
  ) +
  ggplot2::labs(
    x = "Predicted lifestyle",
    y = "Distinct retained AMGs per AMG-positive vOTU",
    caption = "Zero-AMG vOTUs require a separate all-vOTU catalog and are not represented here."
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none")
save_plot_pair(load_plot, file.path(output_dir, paste0(prefix, "_virus_load")), 6.4, 5.2)

# Complete zero-valued sample/category combinations before averaging; otherwise
# absent categories would be omitted and group means would be biased upward.
sample_totals <- aggregate(
  tpm ~ sample_id + source + season + lifestyle + amg_category,
  data = amg,
  FUN = sum
)
complete_grid <- expand.grid(
  sample_id = metadata$sample_id,
  lifestyle = levels(amg$lifestyle),
  amg_category = sort(unique(amg$amg_category)),
  stringsAsFactors = FALSE
)
complete_grid <- merge(
  complete_grid,
  metadata[c("sample_id", "source", "season")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
sample_totals <- merge(
  complete_grid,
  sample_totals,
  by = c("sample_id", "source", "season", "lifestyle", "amg_category"),
  all.x = TRUE,
  sort = FALSE
)
sample_totals$tpm[is.na(sample_totals$tpm)] <- 0
sample_totals$lifestyle <- factor(sample_totals$lifestyle, levels = c("virulent", "temperate"))

lifestyle_totals <- aggregate(
  tpm ~ sample_id + lifestyle,
  data = sample_totals,
  FUN = sum
)
names(lifestyle_totals)[3L] <- "lifestyle_tpm"
sample_totals <- merge(
  sample_totals,
  lifestyle_totals,
  by = c("sample_id", "lifestyle"),
  all.x = TRUE,
  sort = FALSE
)
sample_totals$percent_within_lifestyle <- ifelse(
  sample_totals$lifestyle_tpm > 0,
  100 * sample_totals$tpm / sample_totals$lifestyle_tpm,
  NA_real_
)
write_tsv(sample_totals, file.path(output_dir, paste0(prefix, "_sample_category_abundance.tsv")))

group_keys <- c("source", "season", "lifestyle", "amg_category")
group_mean <- aggregate(sample_totals$tpm, sample_totals[group_keys], mean)
group_sd <- aggregate(sample_totals$tpm, sample_totals[group_keys], stats::sd)
group_n <- aggregate(sample_totals$tpm, sample_totals[group_keys], length)
names(group_mean)[5L] <- "mean_tpm"
names(group_sd)[5L] <- "sd_tpm"
names(group_n)[5L] <- "n_samples"
group_summary <- Reduce(
  function(x, y) merge(x, y, by = group_keys, all = TRUE),
  list(group_mean, group_sd, group_n)
)
group_summary$se_tpm <- group_summary$sd_tpm / sqrt(group_summary$n_samples)
write_tsv(group_summary, file.path(output_dir, paste0(prefix, "_group_category_summary.tsv")))

category_totals <- aggregate(tpm ~ amg_category, data = sample_totals, sum)
category_totals <- category_totals[order(category_totals$tpm, decreasing = TRUE), , drop = FALSE]
top_categories <- head(category_totals$amg_category, 12L)
dot_data <- sample_totals
dot_data$category_plot <- ifelse(
  dot_data$amg_category %in% top_categories,
  dot_data$amg_category,
  "Other"
)
dot_summary <- aggregate(
  tpm ~ source + season + lifestyle + category_plot,
  data = dot_data,
  mean
)
dot_plot <- ggplot2::ggplot(
  dot_summary,
  ggplot2::aes(x = season, y = category_plot, size = tpm, colour = tpm)
) +
  ggplot2::geom_point(alpha = 0.82) +
  ggplot2::facet_grid(lifestyle ~ source, scales = "free_x", space = "free_x") +
  ggplot2::scale_size_continuous(range = c(0.5, 8)) +
  ggplot2::scale_colour_viridis_c(option = "C", trans = "sqrt") +
  ggplot2::labs(
    x = "Season", y = "AMG category",
    size = "Mean TPM", colour = "Mean TPM"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
    panel.grid.minor = ggplot2::element_blank()
  )
save_plot_pair(dot_plot, file.path(output_dir, paste0(prefix, "_category_abundance_dotplot")), 9.0, 7.2)

composition_data <- sample_totals
composition_data$category_plot <- ifelse(
  composition_data$amg_category %in% top_categories,
  composition_data$amg_category,
  "Other"
)
composition_data <- aggregate(
  composition_data["percent_within_lifestyle"],
  by = composition_data[c("sample_id", "source", "season", "lifestyle", "category_plot")],
  FUN = function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
)
composition_mean <- aggregate(
  composition_data["percent_within_lifestyle"],
  by = composition_data[c("source", "season", "lifestyle", "category_plot")],
  FUN = function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
)
write_tsv(
  composition_mean,
  file.path(output_dir, paste0(prefix, "_group_mean_category_percent.tsv"))
)
composition_plot <- ggplot2::ggplot(
  composition_mean,
  ggplot2::aes(x = lifestyle, y = percent_within_lifestyle, fill = category_plot)
) +
  ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  ggplot2::facet_grid(source ~ season) +
  ggplot2::labs(
    x = "Predicted lifestyle",
    y = "Mean within-lifestyle AMG composition (%)",
    fill = "AMG category"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot_pair(composition_plot, file.path(output_dir, paste0(prefix, "_functional_composition")), 9.2, 6.2)

# Category enrichment is tested on unique virus-category presence, not repeated
# sample rows. Shuffling lifestyle labels at the virus level preserves each virus's multi-category AMG profile and avoids pseudoreplication.
presence <- unique(amg_catalog[c("virus_id", "lifestyle", "amg_category")])
lifestyle_levels <- c("virulent", "temperate")
category_levels <- sort(unique(presence$amg_category))
if (length(category_levels) < 2L) {
  stop("At least two AMG categories are required for category association.", call. = FALSE)
}
observed <- table(
  factor(presence$lifestyle, levels = lifestyle_levels),
  factor(presence$amg_category, levels = category_levels)
)

chisq_components <- function(tab) {
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  statistic <- sum(ifelse(expected > 0, (tab - expected)^2 / expected, 0))
  row_fraction <- rowSums(tab) / sum(tab)
  column_fraction <- colSums(tab) / sum(tab)
  denominator <- sqrt(
    expected *
      outer(1 - row_fraction, rep(1, ncol(tab))) *
      outer(rep(1, nrow(tab)), 1 - column_fraction)
  )
  standardized <- (tab - expected) / denominator
  list(statistic = statistic, expected = expected, standardized = standardized)
}

observed_components <- chisq_components(observed)
virus_labels <- unique(presence[c("virus_id", "lifestyle")])
n_permutations <- 9999L
permuted_statistics <- numeric(n_permutations)
for (iteration in seq_len(n_permutations)) {
  shuffled <- sample(as.character(virus_labels$lifestyle), replace = FALSE)
  names(shuffled) <- virus_labels$virus_id
  permuted_table <- table(
    factor(shuffled[presence$virus_id], levels = lifestyle_levels),
    factor(presence$amg_category, levels = category_levels)
  )
  permuted_statistics[iteration] <- chisq_components(permuted_table)$statistic
}
permutation_p <- (1 + sum(permuted_statistics >= observed_components$statistic)) /
  (n_permutations + 1)
degrees_of_freedom <- (nrow(observed) - 1L) * (ncol(observed) - 1L)

association_summary <- data.frame(
  analysis_unit = "unique vOTU-category presence",
  statistic = observed_components$statistic,
  degrees_of_freedom = degrees_of_freedom,
  asymptotic_p_value = stats::pchisq(
    observed_components$statistic,
    df = degrees_of_freedom,
    lower.tail = FALSE
  ),
  virus_label_permutation_p_value = permutation_p,
  permutations = n_permutations,
  minimum_expected_count = min(observed_components$expected),
  expected_cells_below_5 = sum(observed_components$expected < 5),
  primary_p_value = "virus_label_permutation_p_value",
  stringsAsFactors = FALSE
)
write_tsv(
  association_summary,
  file.path(output_dir, paste0(prefix, "_lifestyle_category_association.tsv"))
)

association_cells <- expand.grid(
  lifestyle = lifestyle_levels,
  amg_category = category_levels,
  stringsAsFactors = FALSE
)
association_cells$observed <- as.vector(observed)
association_cells$expected <- as.vector(observed_components$expected)
association_cells$standardized_residual <- as.vector(observed_components$standardized)
write_tsv(
  association_cells,
  file.path(output_dir, paste0(prefix, "_lifestyle_category_cells.tsv"))
)

residual_plot <- ggplot2::ggplot(
  association_cells,
  ggplot2::aes(
    x = lifestyle,
    y = amg_category,
    size = abs(standardized_residual),
    fill = standardized_residual
  )
) +
  ggplot2::geom_point(shape = 21, colour = "grey25", stroke = 0.25) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0
  ) +
  ggplot2::scale_size_continuous(range = c(1.5, 9)) +
  ggplot2::labs(
    x = "Predicted lifestyle", y = "AMG category",
    size = "|standardized residual|",
    fill = "Standardized residual",
    caption = sprintf(
      "Virus-label permutation test: %s (%d permutations)",
      format_p_value(permutation_p), n_permutations
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())
save_plot_pair(residual_plot, file.path(output_dir, paste0(prefix, "_category_residuals")), 7.6, 7.0)

write_session_info(output_dir)
