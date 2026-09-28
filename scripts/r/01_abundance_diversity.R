#!/usr/bin/env Rscript
#
# Usage:
#   Rscript 01_abundance_diversity.R \
#     counts.tsv composition.tsv metadata.tsv output_dir [analysis_prefix]
#
# counts.tsv contains genuine integer read counts and is used for richness estimators. 
# composition.tsv contains TPM or another non-negative abundance measure and is converted to sample-wise relative abundance for Bray-Curtis.
# Both matrices use feature IDs in the first column and sample IDs thereafter.
# metadata.tsv requires sample_id, source, season, and pair_id.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop(
    paste(
      "Usage: 01_abundance_diversity.R counts.tsv composition.tsv",
      "metadata.tsv output_dir [analysis_prefix]"
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

counts_path <- args[1L]
composition_path <- args[2L]
metadata_path <- args[3L]
output_dir <- args[4L]
prefix <- if (length(args) == 5L) args[5L] else "virome"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(RBF_RANDOM_SEED)

counts_feature <- read_feature_matrix(counts_path, "count matrix")
composition_feature <- read_feature_matrix(composition_path, "composition matrix")
if (!setequal(rownames(counts_feature), rownames(composition_feature))) {
  stop("Count and composition matrices must contain the same feature IDs.", call. = FALSE)
}
composition_feature <- composition_feature[rownames(counts_feature), , drop = FALSE]

metadata <- read_metadata(
  metadata_path,
  required = c("sample_id", "source", "season", "pair_id")
)
aligned <- align_two_feature_matrices(
  counts_feature,
  composition_feature,
  metadata,
  minimum_samples = 8L
)
counts <- drop_zero_features(aligned$first, "count matrix")
composition <- aligned$second[, colnames(counts), drop = FALSE]
metadata <- aligned$metadata
assert_integer_counts(counts)

metadata$source <- factor(metadata$source, levels = unique(metadata$source))
metadata$season <- factor(metadata$season, levels = unique(metadata$season))
metadata$pair_id <- factor(metadata$pair_id)
sample_order <- order(metadata$pair_id, metadata$source)
metadata <- metadata[sample_order, , drop = FALSE]
counts <- counts[sample_order, , drop = FALSE]
composition <- composition[sample_order, , drop = FALSE]
metadata$group <- analysis_group(metadata)

if (nlevels(metadata$source) != 2L || nlevels(metadata$season) != 2L) {
  stop("This study-specific script requires exactly two source levels and two season levels.", call. = FALSE)
}
if (any(table(metadata$pair_id) != 2L)) {
  stop("Each pair_id must contain exactly one sample from each source.", call. = FALSE)
}
if (any(vapply(split(metadata$source, metadata$pair_id), function(x) length(unique(x)) != 2L, logical(1)))) {
  stop("Each pair_id must contain both source levels.", call. = FALSE)
}

richness <- vegan::estimateR(counts)
observed <- as.numeric(richness["S.obs", ])
chao1 <- as.numeric(richness["S.chao1", ])
ace <- as.numeric(richness["S.ACE", ])
shannon <- vegan::diversity(counts, index = "shannon")
simpson <- vegan::diversity(counts, index = "simpson")
pielou <- ifelse(observed > 1, shannon / log(observed), NA_real_)
goods_coverage <- 1 - rowSums(counts == 1) / rowSums(counts)

alpha <- data.frame(
  sample_id = metadata$sample_id,
  source = metadata$source,
  season = metadata$season,
  pair_id = metadata$pair_id,
  group = metadata$group,
  observed_features = observed,
  chao1 = chao1,
  ace = ace,
  shannon = shannon,
  simpson = simpson,
  pielou = pielou,
  goods_coverage = goods_coverage,
  stringsAsFactors = FALSE
)
write_tsv(alpha, file.path(output_dir, paste0(prefix, "_alpha_diversity.tsv")))

planned_univariate_tests <- function(data, value_column) {
  source_levels <- levels(data$source)
  season_levels <- levels(data$season)
  results <- list()
  index <- 1L

  for (season_value in season_levels) {
    subset_data <- data[data$season == season_value, , drop = FALSE]
    first <- subset_data[subset_data$source == source_levels[1L], c("pair_id", value_column)]
    second <- subset_data[subset_data$source == source_levels[2L], c("pair_id", value_column)]
    names(first)[2L] <- "first"
    names(second)[2L] <- "second"
    paired <- merge(first, second, by = "pair_id")
    if (nrow(paired) < 2L) {
      stop(sprintf("Fewer than two complete pairs in season '%s'.", season_value), call. = FALSE)
    }
    if (all(paired$first - paired$second == 0)) {
      statistic <- 0
      p_value <- 1
    } else {
      test <- stats::wilcox.test(
        paired$first,
        paired$second,
        paired = TRUE,
        exact = FALSE
      )
      statistic <- unname(test$statistic)
      p_value <- test$p.value
    }
    results[[index]] <- data.frame(
      comparison = sprintf("%s vs %s within %s", source_levels[1L], source_levels[2L], season_value),
      test = "paired Wilcoxon signed-rank",
      n_first = nrow(paired),
      n_second = nrow(paired),
      effect_median_difference = stats::median(paired$first - paired$second),
      statistic = statistic,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }

  for (source_value in source_levels) {
    subset_data <- data[data$source == source_value, , drop = FALSE]
    first <- subset_data[subset_data$season == season_levels[1L], value_column]
    second <- subset_data[subset_data$season == season_levels[2L], value_column]
    if (length(first) < 2L || length(second) < 2L) {
      stop(sprintf("Fewer than two dates per season for source '%s'.", source_value), call. = FALSE)
    }
    if (length(unique(c(first, second))) == 1L) {
      statistic <- 0
      p_value <- 1
    } else {
      test <- stats::wilcox.test(first, second, paired = FALSE, exact = FALSE)
      statistic <- unname(test$statistic)
      p_value <- test$p.value
    }
    results[[index]] <- data.frame(
      comparison = sprintf("%s vs %s within %s", season_levels[1L], season_levels[2L], source_value),
      test = "Wilcoxon rank-sum",
      n_first = length(first),
      n_second = length(second),
      effect_median_difference = stats::median(first) - stats::median(second),
      statistic = statistic,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }

  output <- do.call(rbind, results)
  output$p_adjust_bh <- stats::p.adjust(output$p_value, method = "BH")
  output
}

alpha_tests <- planned_univariate_tests(alpha, "shannon")
write_tsv(alpha_tests, file.path(output_dir, paste0(prefix, "_shannon_planned_tests.tsv")))

alpha_plot <- ggplot2::ggplot(alpha, ggplot2::aes(x = group, y = shannon, fill = group)) +
  ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.78) +
  ggplot2::geom_jitter(width = 0.08, height = 0, size = 2.2, alpha = 0.75) +
  ggplot2::labs(x = NULL, y = "Shannon diversity", fill = "Source:season") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    panel.grid = ggplot2::element_blank(),
    legend.position = "none"
  )
save_plot_pair(alpha_plot, file.path(output_dir, paste0(prefix, "_shannon")), 7.0, 5.0)

# Beta diversity
composition <- drop_zero_features(composition, "composition matrix")
relative <- to_relative_abundance(composition)
bray <- vegan::vegdist(relative, method = "bray")
pcoa <- vegan::wcmdscale(bray, k = 2L, eig = TRUE, add = "lingoes")
positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
if (ncol(pcoa$points) < 2L || length(positive_eigenvalues) < 2L) {
  stop("Bray-Curtis PCoA produced fewer than two usable positive axes.", call. = FALSE)
}
explained <- 100 * positive_eigenvalues / sum(positive_eigenvalues)

pcoa_scores <- data.frame(
  sample_id = rownames(pcoa$points),
  source = metadata$source[match(rownames(pcoa$points), metadata$sample_id)],
  season = metadata$season[match(rownames(pcoa$points), metadata$sample_id)],
  group = metadata$group[match(rownames(pcoa$points), metadata$sample_id)],
  axis_1 = pcoa$points[, 1L],
  axis_2 = pcoa$points[, 2L],
  stringsAsFactors = FALSE
)
write_tsv(pcoa_scores, file.path(output_dir, paste0(prefix, "_bray_pcoa_scores.tsv")))
write_tsv(
  data.frame(
    axis = seq_along(pcoa$eig),
    eigenvalue = pcoa$eig,
    positive_eigenvalue_percent = ifelse(
      pcoa$eig > 0,
      100 * pcoa$eig / sum(positive_eigenvalues),
      NA_real_
    ),
    correction = "Lingoes",
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, paste0(prefix, "_bray_pcoa_eigenvalues.tsv"))
)

pcoa_plot <- ggplot2::ggplot(
  pcoa_scores,
  ggplot2::aes(x = axis_1, y = axis_2, colour = group, shape = source)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.35) +
  ggplot2::geom_point(size = 3.4, alpha = 0.82) +
  ggplot2::labs(
    x = sprintf("PCoA1 (%.2f%%)", explained[1L]),
    y = sprintf("PCoA2 (%.2f%%)", explained[2L]),
    colour = "Source:season",
    shape = "Source"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

if (all(table(pcoa_scores$group) >= 3L)) {
  pcoa_plot <- pcoa_plot + ggplot2::stat_ellipse(
    ggplot2::aes(group = group, colour = group),
    level = 0.90,
    linewidth = 0.6,
    linetype = 2,
    show.legend = FALSE
  )
}
save_plot_pair(pcoa_plot, file.path(output_dir, paste0(prefix, "_bray_pcoa")), 7.2, 5.4)

# Planned PERMANOVA effects respect the split design: source is permuted within sampling dates; season permutes whole date pairs.
source_control <- within_block_permutations(metadata$pair_id, 999L)
season_control <- whole_block_permutations(metadata$pair_id, 999L)

set.seed(RBF_RANDOM_SEED)
source_model <- vegan::adonis2(
  bray ~ season + source,
  data = metadata,
  permutations = source_control,
  by = "margin"
)
set.seed(RBF_RANDOM_SEED)
season_model <- vegan::adonis2(
  bray ~ source + season,
  data = metadata,
  permutations = season_control,
  by = "margin"
)
set.seed(RBF_RANDOM_SEED)
interaction_model <- vegan::adonis2(
  bray ~ source * season,
  data = metadata,
  permutations = season_control,
  by = "margin"
)

planned_effects <- rbind(
  extract_adonis_term(
    source_model, "source", "overall source effect adjusted for season",
    "within-date-pair permutations", nrow(metadata)
  ),
  extract_adonis_term(
    season_model, "season", "overall season effect adjusted for source",
    "whole-date-pair permutations", nrow(metadata)
  ),
  extract_adonis_term(
    interaction_model, "source:season", "source-by-season interaction",
    "whole-date-pair permutations", nrow(metadata)
  )
)
planned_effects$p_adjust_bh <- stats::p.adjust(planned_effects$p_value, method = "BH")
write_tsv(planned_effects, file.path(output_dir, paste0(prefix, "_permanova_planned_effects.tsv")))

pairwise_permanova <- function(relative_abundance, metadata) {
  source_levels <- levels(metadata$source)
  season_levels <- levels(metadata$season)
  results <- list()
  index <- 1L

  for (season_value in season_levels) {
    keep <- metadata$season == season_value
    sub_metadata <- droplevels(metadata[keep, , drop = FALSE])
    sub_distance <- vegan::vegdist(relative_abundance[keep, , drop = FALSE], method = "bray")
    control <- within_block_permutations(sub_metadata$pair_id, 999L)
    set.seed(RBF_RANDOM_SEED)
    model <- vegan::adonis2(
      sub_distance ~ source,
      data = sub_metadata,
      permutations = control
    )
    results[[index]] <- extract_adonis_term(
      model,
      "source",
      sprintf("%s vs %s within %s", source_levels[1L], source_levels[2L], season_value),
      "within-date-pair permutations",
      nrow(sub_metadata)
    )
    index <- index + 1L
  }

  for (source_value in source_levels) {
    keep <- metadata$source == source_value
    sub_metadata <- droplevels(metadata[keep, , drop = FALSE])
    sub_distance <- vegan::vegdist(relative_abundance[keep, , drop = FALSE], method = "bray")
    set.seed(RBF_RANDOM_SEED)
    model <- vegan::adonis2(
      sub_distance ~ season,
      data = sub_metadata,
      permutations = 999L
    )
    results[[index]] <- extract_adonis_term(
      model,
      "season",
      sprintf("%s vs %s within %s", season_levels[1L], season_levels[2L], source_value),
      "free permutations of sampling dates within source",
      nrow(sub_metadata)
    )
    index <- index + 1L
  }

  output <- do.call(rbind, results)
  output$p_adjust_bh <- stats::p.adjust(output$p_value, method = "BH")
  output
}

pairwise_results <- pairwise_permanova(relative, metadata)
write_tsv(
  pairwise_results,
  file.path(output_dir, paste0(prefix, "_permanova_pairwise_planned.tsv"))
)

dispersion <- vegan::betadisper(
  bray,
  metadata$group,
  type = "median",
  bias.adjust = TRUE,
  add = "lingoes"
)
set.seed(RBF_RANDOM_SEED)
dispersion_test <- vegan::permutest(dispersion, permutations = 999L)
dispersion_table <- data.frame(
  term = rownames(dispersion_test$tab),
  as.data.frame(dispersion_test$tab, check.names = FALSE),
  row.names = NULL,
  check.names = FALSE
)
write_tsv(
  dispersion_table,
  file.path(output_dir, paste0(prefix, "_betadisper_global.tsv"))
)

write_session_info(output_dir)
