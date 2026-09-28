#!/usr/bin/env Rscript
# Usage:
#   Rscript 06_virus_host_covariation.R \
#     virus_abundance.tsv host_abundance.tsv links.tsv metadata.tsv output_dir \
#     [analysis_prefix] [min_prevalence]
#
# links.tsv requires virus_id and host_id; evidence_type is optional.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 5:7) {
  stop(
    paste(
      "Usage: 06_virus_host_covariation.R virus_abundance.tsv",
      "host_abundance.tsv links.tsv metadata.tsv output_dir",
      "[analysis_prefix] [min_prevalence]"
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

virus_path <- args[1L]
host_path <- args[2L]
links_path <- args[3L]
metadata_path <- args[4L]
output_dir <- args[5L]
prefix <- if (length(args) >= 6L) args[6L] else "virus_host"
min_prevalence <- if (length(args) == 7L) suppressWarnings(as.numeric(args[7L])) else 0.20
if (is.na(min_prevalence) || min_prevalence < 0 || min_prevalence > 1) {
  stop("min_prevalence must be a number from 0 to 1.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

virus_feature <- read_feature_matrix(virus_path, "virus abundance matrix")
host_feature <- read_feature_matrix(host_path, "host abundance matrix")
metadata <- read_metadata(metadata_path, required = c("sample_id", "source", "season"))
aligned <- align_two_feature_matrices(
  virus_feature,
  host_feature,
  metadata,
  minimum_samples = 8L
)
virus <- to_relative_abundance(drop_zero_features(aligned$first, "virus abundance matrix"))
host <- to_relative_abundance(drop_zero_features(aligned$second, "host abundance matrix"))
metadata <- aligned$metadata

links <- read_tsv(links_path)
assert_columns(links, c("virus_id", "host_id"), "putative-link table")
links$virus_id <- trimws(as.character(links$virus_id))
links$host_id <- trimws(as.character(links$host_id))
if (anyNA(links[c("virus_id", "host_id")]) || any(as.matrix(links[c("virus_id", "host_id")]) == "")) {
  stop("virus_id and host_id in the putative-link table cannot be missing or empty.", call. = FALSE)
}
if (!"evidence_type" %in% names(links)) links$evidence_type <- "not_provided"
links$evidence_type <- trimws(as.character(links$evidence_type))
links$evidence_type[is.na(links$evidence_type) | links$evidence_type == ""] <- "not_provided"
links <- aggregate(
  evidence_type ~ virus_id + host_id,
  data = links,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)

links$virus_found <- links$virus_id %in% colnames(virus)
links$host_found <- links$host_id %in% colnames(host)
write_tsv(
  links[!links$virus_found | !links$host_found, , drop = FALSE],
  file.path(output_dir, paste0(prefix, "_unmatched_predicted_links.tsv"))
)
links <- links[links$virus_found & links$host_found, , drop = FALSE]
if (nrow(links) == 0L) {
  stop("No putative links match both abundance matrices.", call. = FALSE)
}

strata <- list(all_samples = rep(TRUE, nrow(metadata)))
for (value in unique(metadata$source)) {
  strata[[paste0("source_", value)]] <- metadata$source == value
}
for (value in unique(metadata$season)) {
  strata[[paste0("season_", value)]] <- metadata$season == value
}
strata <- strata[vapply(strata, function(x) sum(x), integer(1)) >= 6L]

correlation_rows <- list()
row_index <- 1L
for (stratum_name in names(strata)) {
  keep <- strata[[stratum_name]]
  for (link_index in seq_len(nrow(links))) {
    virus_values <- virus[keep, links$virus_id[link_index]]
    host_values <- host[keep, links$host_id[link_index]]
    virus_prevalence <- mean(virus_values > 0)
    host_prevalence <- mean(host_values > 0)
    status <- "tested"
    rho <- NA_real_
    statistic <- NA_real_
    p_value <- NA_real_

    if (virus_prevalence < min_prevalence || host_prevalence < min_prevalence) {
      status <- "below_prevalence_filter"
    } else if (length(unique(virus_values)) < 2L || length(unique(host_values)) < 2L) {
      status <- "constant_abundance"
    } else {
      test <- stats::cor.test(
        virus_values,
        host_values,
        method = "spearman",
        exact = FALSE,
        alternative = "two.sided"
      )
      rho <- unname(test$estimate)
      statistic <- unname(test$statistic)
      p_value <- test$p.value
    }

    correlation_rows[[row_index]] <- data.frame(
      stratum = stratum_name,
      virus_id = links$virus_id[link_index],
      host_id = links$host_id[link_index],
      evidence_type = links$evidence_type[link_index],
      n_samples = sum(keep),
      virus_prevalence = virus_prevalence,
      host_prevalence = host_prevalence,
      spearman_rho = rho,
      statistic = statistic,
      p_value = p_value,
      status = status,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}

correlations <- do.call(rbind, correlation_rows)
correlations$p_adjust_bh <- NA_real_
for (stratum_name in unique(correlations$stratum)) {
  index <- correlations$stratum == stratum_name & correlations$status == "tested"
  correlations$p_adjust_bh[index] <- stats::p.adjust(
    correlations$p_value[index],
    method = "BH"
  )
}
correlations$selected_for_plot <-
  correlations$status == "tested" &
  !is.na(correlations$p_adjust_bh) &
  correlations$p_adjust_bh < 0.05 &
  abs(correlations$spearman_rho) >= 0.60
write_tsv(
  correlations,
  file.path(output_dir, paste0(prefix, "_predicted_link_covariation.tsv"))
)

selected <- correlations[correlations$selected_for_plot, , drop = FALSE]
if (nrow(selected) > 0L) {
  selected <- do.call(rbind, lapply(split(selected, selected$stratum), function(table) {
    table <- table[order(table$p_adjust_bh, -abs(table$spearman_rho)), , drop = FALSE]
    head(table, 25L)
  }))

  edge_rows <- list()
  node_rows <- list()
  edge_index <- 1L
  node_index <- 1L
  for (stratum_name in unique(selected$stratum)) {
    table <- selected[selected$stratum == stratum_name, , drop = FALSE]
    virus_nodes <- sort(unique(table$virus_id))
    host_nodes <- sort(unique(table$host_id))
    virus_y <- setNames(seq(1, 0, length.out = length(virus_nodes)), virus_nodes)
    host_y <- setNames(seq(1, 0, length.out = length(host_nodes)), host_nodes)
    for (index in seq_len(nrow(table))) {
      edge_rows[[edge_index]] <- data.frame(
        stratum = stratum_name,
        x = 0,
        y = virus_y[[table$virus_id[index]]],
        xend = 1,
        yend = host_y[[table$host_id[index]]],
        spearman_rho = table$spearman_rho[index],
        q_value = table$p_adjust_bh[index],
        stringsAsFactors = FALSE
      )
      edge_index <- edge_index + 1L
    }
    node_rows[[node_index]] <- rbind(
      data.frame(
        stratum = stratum_name, node_id = virus_nodes, node_type = "virus",
        x = 0, y = unname(virus_y), stringsAsFactors = FALSE
      ),
      data.frame(
        stratum = stratum_name, node_id = host_nodes, node_type = "host",
        x = 1, y = unname(host_y), stringsAsFactors = FALSE
      )
    )
    node_index <- node_index + 1L
  }
  edge_plot_data <- do.call(rbind, edge_rows)
  edge_plot_data$plot_significance <- -log10(
    pmax(edge_plot_data$q_value, .Machine$double.xmin)
  )
  node_plot_data <- do.call(rbind, node_rows)

  network_plot <- ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = edge_plot_data,
      ggplot2::aes(
        x = x, y = y, xend = xend, yend = yend,
        colour = spearman_rho, linewidth = plot_significance
      ),
      curvature = 0.12, alpha = 0.62
    ) +
    ggplot2::geom_point(
      data = node_plot_data,
      ggplot2::aes(x = x, y = y, fill = node_type),
      shape = 21, size = 3.2, colour = "grey20"
    ) +
    ggplot2::geom_text(
      data = node_plot_data[node_plot_data$node_type == "virus", , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = node_id),
      hjust = 1.08, size = 2.5
    ) +
    ggplot2::geom_text(
      data = node_plot_data[node_plot_data$node_type == "host", , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = node_id),
      hjust = -0.08, size = 2.5
    ) +
    ggplot2::facet_wrap(~stratum, scales = "free_y") +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC", mid = "grey88", high = "#B2182B", midpoint = 0,
      limits = c(-1, 1)
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 1.5)) +
    ggplot2::scale_fill_manual(values = c(virus = "#80B1D3", host = "#FDB462")) +
    ggplot2::coord_cartesian(xlim = c(-0.35, 1.35), clip = "off") +
    ggplot2::labs(
      x = NULL, y = NULL, colour = "Spearman rho", linewidth = "-log10(BH q)",
      fill = "Node type",
      caption = paste(
        "Only a priori predicted links with |rho| >= 0.60 and BH q < 0.05 are drawn.",
        "Covariation is not evidence of infection."
      )
    ) +
    ggplot2::theme_void(base_size = 9.5) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(5.5, 70, 5.5, 70)
    )
  save_plot_pair(
    network_plot,
    file.path(output_dir, paste0(prefix, "_covariation_network")),
    12.0,
    max(6.0, 3.5 * ceiling(length(unique(selected$stratum)) / 2))
  )
} else {
  empty_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0,
      label = "No predicted virus-host pair passed |rho| >= 0.60 and BH q < 0.05.",
      size = 4
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()
  save_plot_pair(
    empty_plot,
    file.path(output_dir, paste0(prefix, "_covariation_network")),
    8.0,
    3.5
  )
}

write_session_info(output_dir)
