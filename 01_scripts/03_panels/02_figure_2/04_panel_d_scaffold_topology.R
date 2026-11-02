#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(
  sys.frames(),
  function(x) {
    z <- x$ofile
    if (is.null(z) || !length(z)) "" else as.character(z[[1]])
  },
  character(1)
)
script_path <- if (any(nzchar(source_files))) {
  normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path("01_scripts", "03_panels", "02_figure_2", "04_panel_d_scaffold_topology.R"),
    winslash = "/", mustWork = TRUE
  )
}
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
p <- function(...) file.path(root, ...)

data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
if (nzchar(data_override) && !dir.exists(data_override)) {
  stop("FIGURE_DATA_DIR does not exist: ", data_override)
}
data_candidates <- c(
  data_override,
  p("04_plot_data")
)
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
data_dir <- data_candidates[[1]]
read_fig <- function(file) {
  read_tsv(
    file.path(data_dir, "02_figure_2", file),
    show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN")
  )
}

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) {
  output_override
} else {
  p("06_panel_output", "02_figure_2")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) {
  if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x))
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), x,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), x,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_png, bg = "white", limitsize = FALSE
  )
}

ink <- "#30363B"
border <- "#DCE3E7"
blue_dark <- "#B5C8D1"
blue <- "#D2DEE4"
blue_pale <- "#F1F5F7"
red_dark <- "#CDB0AD"
red <- "#E4D0CD"
red_pale <- "#FAF1F0"
core_fill <- red
core_dark <- red_dark
scaffold_fill <- blue
scaffold_dark <- blue_dark
fdr_sig <- red_dark
fdr_ns <- "#7D868E"

as_flag <- function(x) {
  if (is.null(x)) return(logical())
  if (is.logical(x)) replace(x, is.na(x), FALSE) else {
    tolower(as.character(dplyr::coalesce(x, ""))) %in% c("true", "t", "1", "yes", "y")
  }
}

network_nodes <- read_fig("05_panel_g_pathogen_network_nodes.tsv")
pin1_mask <- grepl("Q4UG71", dplyr::coalesce(as.character(network_nodes$node_id), ""), ignore.case = TRUE) |
  grepl("(^|[^A-Za-z0-9])(Ta)?Pin1([^A-Za-z0-9]|$)", paste(
    dplyr::coalesce(as.character(network_nodes$display_label), ""),
    dplyr::coalesce(as.character(network_nodes$canonical_gene), ""),
    dplyr::coalesce(as.character(network_nodes$gene_symbol), ""),
    dplyr::coalesce(as.character(network_nodes$preferred_name), "")
  ), ignore.case = TRUE)

if ("local_core" %in% names(network_nodes)) {
  network_nodes$local_core <- as_flag(network_nodes$local_core)
}
if ("accessible_in_network" %in% names(network_nodes)) {
  network_nodes$accessible_in_network <- as_flag(network_nodes$accessible_in_network)
}
# Pin1 is never local-core. Do not promote recurrent_core to Core.
if ("local_core" %in% names(network_nodes)) {
  network_nodes$local_core[pin1_mask] <- FALSE
}

if (all(c("local_core", "accessible_in_network") %in% names(network_nodes))) {
  network_nodes$network_role <- dplyr::case_when(
    network_nodes$local_core ~ "Core",
    network_nodes$accessible_in_network ~ "Accessible scaffold",
    TRUE ~ "Network context"
  )
} else {
  if (!"network_role" %in% names(network_nodes) && "scaffold_tier" %in% names(network_nodes)) {
    network_nodes$network_role <- dplyr::case_when(
      network_nodes$scaffold_tier %in% c("strict", "Core") ~ "Core",
      network_nodes$scaffold_tier %in% c("relaxed", "accessible_scaffold", "Accessible scaffold") ~ "Accessible scaffold",
      TRUE ~ "Network context"
    )
  }
  if ("network_role" %in% names(network_nodes)) {
    network_nodes$network_role <- dplyr::case_when(
      pin1_mask ~ "Network context",
      network_nodes$network_role %in% c("Core", "strict") ~ "Core",
      network_nodes$network_role %in% c("Accessible scaffold", "relaxed", "accessible_scaffold") ~ "Accessible scaffold",
      TRUE ~ network_nodes$network_role
    )
  }
}
network_nodes$network_role[pin1_mask] <- "Network context"

if (!"module_status" %in% names(network_nodes)) {
  network_nodes$module_status <- ifelse(
    !is.na(network_nodes$module_id) & nzchar(as.character(network_nodes$module_id)),
    "network_module",
    "unassigned"
  )
}
required_columns <- c(
  "module_status", "network_role", "degree", "strength",
  "betweenness", "eigenvector"
)
missing_columns <- setdiff(required_columns, names(network_nodes))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (any(pin1_mask & network_nodes$network_role == "Core")) {
  stop("Pin1 is never local-core.")
}

# Only connected scaffold proteins are compared. Network-context nodes and
# isolated scaffold proteins are deliberately outside this panel's question.
# Core = local-core; non-core scaffold = accessible-in-network minus core.
analysis_nodes <- network_nodes %>%
  filter(
    module_status == "network_module",
    network_role %in% c("Core", "Accessible scaffold")
  ) %>%
  transmute(
    group = if_else(
      network_role == "Core",
      "Connected core",
      "Connected non-core scaffold"
    ),
    degree = as.numeric(degree),
    strength = as.numeric(strength),
    betweenness = as.numeric(betweenness),
    eigenvector = as.numeric(eigenvector)
  ) %>%
  filter(if_all(c(degree, strength, betweenness, eigenvector), is.finite)) %>%
  mutate(
    group = factor(
      group,
      levels = c("Connected core", "Connected non-core scaffold")
    )
  )

group_counts <- count(analysis_nodes, group, name = "n")
n_core <- group_counts$n[group_counts$group == "Connected core"]
n_scaffold <- group_counts$n[group_counts$group == "Connected non-core scaffold"]
if (length(n_core) != 1 || length(n_scaffold) != 1 || min(n_core, n_scaffold) < 2) {
  stop("Both connected-core and connected non-core scaffold groups require at least two proteins.")
}

metric_definitions <- tibble(
  metric_raw = c("degree", "strength", "betweenness", "eigenvector"),
  metric = c("Degree", "Strength", "Betweenness", "Eigenvector")
) %>%
  mutate(metric = factor(metric, levels = metric))

cliffs_delta <- function(x, y) {
  pooled_ranks <- rank(c(x, y), ties.method = "average")
  u_statistic <- sum(pooled_ranks[seq_along(x)]) - length(x) * (length(x) + 1) / 2
  2 * u_statistic / (length(x) * length(y)) - 1
}

bootstrap_delta <- function(x, y, replicates = 2500L) {
  replicate(
    replicates,
    cliffs_delta(
      sample(x, length(x), replace = TRUE),
      sample(y, length(y), replace = TRUE)
    )
  )
}

set.seed(20260713)
statistics <- lapply(seq_len(nrow(metric_definitions)), function(i) {
  metric_name <- metric_definitions$metric_raw[[i]]
  x <- analysis_nodes[[metric_name]][analysis_nodes$group == "Connected core"]
  y <- analysis_nodes[[metric_name]][analysis_nodes$group == "Connected non-core scaffold"]
  wilcox_result <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE, correct = FALSE)
  delta_bootstrap <- bootstrap_delta(x, y)
  tibble(
    metric_raw = metric_name,
    metric = metric_definitions$metric[[i]],
    p_value = wilcox_result$p.value,
    cliffs_delta = cliffs_delta(x, y),
    ci_low = unname(quantile(delta_bootstrap, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(delta_bootstrap, 0.975, na.rm = TRUE))
  )
}) %>%
  bind_rows() %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significant = p_value < 0.05,
    metric = factor(as.character(metric), levels = levels(metric_definitions$metric))
  )

format_probability <- function(x) {
  ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
}
statistics <- statistics %>%
  mutate(
    statistic_label = paste("P", format_probability(p_value))
  )

distribution_data <- analysis_nodes %>%
  pivot_longer(
    cols = c(degree, strength, betweenness, eigenvector),
    names_to = "metric_raw", values_to = "raw_value"
  ) %>%
  left_join(metric_definitions, by = "metric_raw") %>%
  mutate(
    metric = factor(as.character(metric), levels = levels(metric_definitions$metric)),
    plot_value = case_when(
      metric_raw %in% c("degree", "strength") ~ log1p(raw_value),
      metric_raw == "betweenness" ~ log10(raw_value + 1e-8),
      metric_raw == "eigenvector" ~ log10(raw_value + 1e-10),
      TRUE ~ raw_value
    )
  )

annotation_data <- distribution_data %>%
  distinct(metric) %>%
  left_join(
    statistics %>%
      select(metric, statistic_label, significant),
    by = "metric"
  )

group_colours <- c(
  `Connected core` = core_dark,
  `Connected non-core scaffold` = scaffold_dark
)
group_fills <- c(
  `Connected core` = core_fill,
  `Connected non-core scaffold` = scaffold_fill
)

set.seed(20260713)
panel <- ggplot(distribution_data, aes(x = group, y = plot_value)) +
  geom_violin(
    aes(fill = group, colour = group),
    width = 0.88, alpha = 0.34, linewidth = 0.38,
    trim = FALSE, scale = "width", na.rm = TRUE
  ) +
  geom_jitter(
    aes(colour = group),
    width = 0.13, height = 0,
    size = 1.32, alpha = 0.82, stroke = 0,
    show.legend = FALSE
  ) +
  geom_boxplot(
    aes(fill = group, colour = group),
    width = 0.19, alpha = 0.52, linewidth = 0.42,
    outlier.shape = NA, median.linewidth = 0.56,
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(annotation_data, !significant),
    aes(x = 1.5, y = Inf, label = statistic_label),
    inherit.aes = FALSE,
    colour = fdr_ns, size = 2.12,
    hjust = 0.5, vjust = 1.35,
    show.legend = FALSE
  ) +
  geom_text(
    data = filter(annotation_data, significant),
    aes(x = 1.5, y = Inf, label = statistic_label),
    inherit.aes = FALSE,
    colour = fdr_sig, fontface = "bold", size = 2.12,
    hjust = 0.5, vjust = 1.35,
    show.legend = FALSE
  ) +
  facet_wrap(
    vars(metric), nrow = 1, scales = "free_y",
    axes = "all_y", axis.labels = "all_y"
  ) +
  scale_x_discrete(
    breaks = levels(analysis_nodes$group), labels = NULL,
    expand = expansion(add = c(0.36, 0.36))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.12, 0.13)),
    labels = scales::label_number(accuracy = 0.1, trim = TRUE)
  ) +
  scale_colour_manual(values = group_colours, name = NULL) +
  labs(x = NULL, y = "Transformed topology value") +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = border, linewidth = 0.38),
    axis.ticks = element_line(colour = ink, linewidth = 0.30),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6.4, colour = ink),
    axis.title.y = element_text(size = 7.2, colour = ink, margin = margin(r = 5)),
    strip.background = element_rect(
      fill = "#F8F9FA", colour = border, linewidth = 0.48
    ),
    strip.text = element_text(
      size = 7.2, colour = ink, face = "bold",
      lineheight = 0.98, margin = margin(3.8, 2, 3.8, 2)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.key = element_blank(),
    legend.key.width = unit(4.1, "mm"),
    legend.key.height = unit(2.8, "mm"),
    legend.text = element_text(size = 6.4),
    legend.spacing.x = unit(1.5, "mm"),
    legend.box.spacing = unit(2.0, "mm"),
    legend.margin = margin(t = 1.5),
    panel.spacing.x = unit(3.2, "mm"),
    plot.margin = margin(4, 8, 3, 7),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  scale_fill_manual(
    values = group_fills,
    breaks = c("Connected core", "Connected non-core scaffold"),
    labels = c(
      sprintf("Core (n=%d)", n_core),
      sprintf("Non-core scaffold (n=%d)", n_scaffold)
    ),
    name = NULL
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 0.82, linewidth = 0.45)
    ),
    colour = "none"
  ) +
  theme(
    legend.justification = "center",
    legend.box.just = "center"
  )

save_panel(panel, "Figure_2_panel_d_scaffold_topology", 184, 126)

cat("Connected core n = ", n_core, "; connected non-core scaffold n = ", n_scaffold, ".\n", sep = "")
print(
  statistics %>%
    transmute(metric = as.character(metric), p_value, fdr, cliffs_delta, ci_low, ci_high)
)
cat("Rendered Figure_2_panel_d_scaffold_topology.\n")
