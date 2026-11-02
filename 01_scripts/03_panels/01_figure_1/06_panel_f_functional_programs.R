#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
  library(igraph); library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) { z <- x$ofile; if (is.null(z) || !length(z)) "" else as.character(z[[1]]) }, character(1))
script_path <- if (any(nzchar(source_files))) normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE) else if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else normalizePath(file.path("01_scripts", "03_panels", "01_figure_1", "06_panel_f_functional_programs.R"), winslash = "/", mustWork = TRUE)
data_env <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
root <- if (nzchar(data_env)) normalizePath(file.path(data_env, "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE); p <- function(...) file.path(root, ...)
data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = ""); if (nzchar(data_override) && !dir.exists(data_override)) stop("FIGURE_DATA_DIR does not exist: ", data_override)
data_candidates <- c(data_override, p("04_plot_data")); data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]; if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first."); data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(file.path(data_dir, "01_figure_1", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = ""); output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "01_figure_1"); dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) { if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x)); ggsave(file.path(output_dir, paste0(stem, ".pdf")), x, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE); ggsave(file.path(output_dir, paste0(stem, ".png")), x, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE) }

ink <- "#30363B"; blue_dark <- "#B5C8D1"; blue <- "#D2DEE4"; blue_pale <- "#F1F5F7"; red_dark <- "#CDB0AD"; red <- "#E4D0CD"; red_pale <- "#FAF1F0"
theme_panel <- theme_classic(base_size = 8, base_family = "sans") + theme(panel.grid = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill = NA, colour = "#DCE3E7", linewidth = 0.35), axis.ticks = element_line(linewidth = 0.3, colour = ink), axis.text = element_text(colour = ink), axis.title = element_text(colour = ink), legend.position = "bottom", legend.justification = "center", legend.box.just = "center", legend.direction = "horizontal", legend.box = "horizontal", legend.key.height = unit(2.2, "mm"), legend.key.width = unit(3.8, "mm"), plot.margin = margin(6, 8, 6, 8), panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA))

program_name <- function(x) {
  z <- coalesce(x, "")
  case_when(
    str_detect(z, "host_signaling_broad|signal_transduction") ~ "Signaling",
    str_detect(z, "host_rna_processing|rna_processing") ~ "RNA processing",
    str_detect(z, "host_ubiquitin|proteostasis|H_NET_13") ~ "Proteostasis",
    str_detect(z, "host_metabolism|H_NET_2") ~ "Metabolism",
    str_detect(z, "H_NET_9") ~ "Translation",
    str_detect(z, "host_microtubule|cytoskeleton|microtubule") ~ "Microtubule",
    str_detect(z, "H_NET_3") ~ "Trafficking",
    str_detect(z, "H_NET_6") ~ "Cell cycle",
    TRUE ~ "Other"
  )
}

as_flag <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) replace(x, is.na(x), FALSE) else str_to_lower(coalesce(as.character(x), "")) %in% c("true", "t", "1", "yes", "y")
}
composition <- read_fig("05_panel_e_host_functions.tsv")
if (!"module_biological_label" %in% names(composition) && "display_name" %in% names(composition)) {
  composition$module_biological_label <- composition$display_name
}
if (!"module_biological_label" %in% names(composition) && "module_label" %in% names(composition)) {
  composition$module_biological_label <- composition$module_label
}
if (!"module_biological_label" %in% names(composition)) composition$module_biological_label <- "Other / unresolved"
composition <- composition %>%
  mutate(program = program_name(paste(module_id, module_label, module_biological_label))) %>%
  group_by(program) %>%
  summarise(scaffold_n = sum(scaffold_n), core_n = sum(core_n), .groups = "drop") %>%
  filter(scaffold_n > 0)

nodes <- read_fig("06_panel_g_host_network_nodes.tsv")
if (!"module_biological_label" %in% names(nodes) && "display_name" %in% names(nodes)) nodes$module_biological_label <- nodes$display_name
if (!"module_biological_label" %in% names(nodes) && "module_label" %in% names(nodes)) nodes$module_biological_label <- nodes$module_label
if (!"module_biological_label" %in% names(nodes)) nodes$module_biological_label <- "Other / unresolved"
nodes <- nodes %>%
  mutate(
    relaxed_member = as_flag(relaxed_member),
    strict_member = as_flag(strict_member),
    program = program_name(paste(coalesce(as.character(module_id), ""), module_biological_label))
  )
edges <- read_fig("07_panel_g_host_network_edges.tsv")
if (!"weight" %in% names(edges)) edges$weight <- 1
graph <- graph_from_data_frame(edges %>% select(source, target), directed = FALSE, vertices = nodes %>% transmute(name = node_id))
scaffold_hubs <- nodes %>%
  filter(relaxed_member | strict_member, degree > 0) %>%
  arrange(desc(strength), desc(degree)) %>%
  slice_head(n = 8) %>%
  pull(node_id)
hub_neighborhood_ids <- unique(c(
  scaffold_hubs,
  unlist(ego(graph, order = 1, nodes = scaffold_hubs, mode = "all") %>% lapply(as_ids), use.names = FALSE)
))
annotated <- nodes %>% mutate(hub_neighborhood = node_id %in% hub_neighborhood_ids)
network_n <- nrow(annotated)
neighborhood_n <- sum(annotated$hub_neighborhood)
enrichment <- annotated %>%
  group_by(program) %>%
  summarise(background_n = n(), neighborhood_genes = sum(hub_neighborhood), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    p_value = fisher.test(matrix(c(
      neighborhood_genes, neighborhood_n - neighborhood_genes,
      background_n - neighborhood_genes,
      network_n - neighborhood_n - background_n + neighborhood_genes
    ), nrow = 2))$p.value,
    odds_ratio = ((neighborhood_genes + 0.5) * (network_n - neighborhood_n - background_n + neighborhood_genes + 0.5)) /
      ((neighborhood_n - neighborhood_genes + 0.5) * (background_n - neighborhood_genes + 0.5))
  ) %>%
  ungroup() %>%
  mutate(fdr_bh = p.adjust(p_value, "BH"), log2_odds_ratio = log2(odds_ratio), significant = fdr_bh < 0.05)

program_order <- composition %>% arrange(scaffold_n) %>% pull(program)
composition <- composition %>% mutate(label = factor(program, levels = program_order))
enrichment <- enrichment %>% filter(program %in% program_order) %>% mutate(label = factor(program, levels = program_order))

p_composition <- ggplot(composition, aes(scaffold_n, label)) +
  geom_col(fill = blue_pale, colour = blue_dark, linewidth = 0.35, width = 0.60) +
  geom_point(aes(size = core_n), shape = 21, fill = red, colour = red_dark, stroke = 0.38) +
  scale_size_continuous(range = c(1.5, 4.6), breaks = pretty_breaks(n = 3), guide = guide_legend(nrow = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, colour = red_dark, stroke = 0.70))) +
  scale_y_discrete(labels = function(x) str_wrap(x, 24)) +
  labs(x = "Recurrent scaffold genes", y = NULL, size = "Core genes") +
  theme_panel +
  theme(axis.text.y = element_text(size = 6.4, lineheight = 0.88))

any_fdr <- isTRUE(any(enrichment$significant, na.rm = TRUE))
fdr_guide <- if (!any_fdr) {
  "none"
} else if (n_distinct(enrichment$significant) > 1) {
  guide_legend(nrow = 1, override.aes = list(shape = c(1, 21), fill = NA, colour = c("#6F767D", red_dark), size = 2.5, stroke = c(0.60, 0.78)))
} else {
  guide_legend(nrow = 1, override.aes = list(shape = 21, fill = NA, colour = red_dark, size = 2.5, stroke = 0.78))
}
p_enrichment <- ggplot(enrichment, aes(log2_odds_ratio, label)) +
  geom_vline(xintercept = 0, colour = "#E4E9EC", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = log2_odds_ratio, yend = label, colour = log2_odds_ratio), linewidth = 2.7, lineend = "round", alpha = 0.55) +
  geom_point(aes(size = neighborhood_genes, fill = log2_odds_ratio, shape = significant), colour = ink, stroke = 0.38) +
  scale_colour_gradient2(low = blue_dark, mid = "#DCE3E7", high = red_dark, midpoint = 0, guide = "none") +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, guide = "none") +
  scale_shape_manual(
    values = c(`TRUE` = 21, `FALSE` = 1),
    breaks = c(FALSE, TRUE),
    drop = FALSE,
    labels = c(`FALSE` = "FDR >= 0.05", `TRUE` = "FDR < 0.05"),
    guide = fdr_guide
  ) +
  scale_size_continuous(range = c(1.8, 4.8), breaks = pretty_breaks(n = 3), guide = guide_legend(nrow = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.65))) +
  labs(x = "Hub-neighborhood enrichment (log2 odds ratio)", y = NULL, size = "Neighborhood genes", shape = NULL) +
  theme_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

panel <- p_composition + p_enrichment +
  plot_layout(widths = c(0.95, 1.15), guides = "collect") &
  theme(
    legend.position = "bottom", legend.justification = "center", legend.box.just = "center",
    legend.direction = "horizontal", legend.box = "horizontal", legend.key = element_blank(),
    legend.title = element_text(size = 5.6, margin = margin(r = 1.5)), legend.text = element_text(size = 5.3),
    legend.key.height = unit(2.2, "mm"), legend.key.width = unit(3.1, "mm"),
    legend.spacing.x = unit(1.4, "mm"), legend.spacing.y = unit(0, "mm"),
    legend.box.spacing = unit(0.6, "mm"), legend.margin = margin(t = 0.4)
  )

save_panel(panel, "Figure_1_panel_f_functional_programs", 183, 92)
cat("Rendered Figure_1_panel_f_functional_programs.\n")
