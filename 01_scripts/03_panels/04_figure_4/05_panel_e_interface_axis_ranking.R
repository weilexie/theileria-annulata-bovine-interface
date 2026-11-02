#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(grid)
})

ofiles <- vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else as.character(x$ofile[[1]]), character(1))
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- if (any(nzchar(ofiles))) normalizePath(tail(ofiles[nzchar(ofiles)], 1), winslash = "/") else if (length(arg)) normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/") else normalizePath("01_scripts/03_panels/04_figure_4/05_panel_e_interface_axis_ranking.R", winslash = "/")
root <- normalizePath(file.path(dirname(script), "..", "..", ".."), winslash = "/")
p <- function(...) file.path(root, ...)
lock_cands <- unique(c(
  file.path(getwd(), "01_scripts", "03_panels", "00_locks.R"),
  file.path(dirname(if (exists("script_path")) script_path else script), "..", "00_locks.R"),
  file.path(dirname(if (exists("script_path")) script_path else script), "..", "00_locks.R"),
  file.path(dirname(if (exists("script_path")) script_path else script), "..", "00_locks.R")
))
lock_file <- lock_cands[file.exists(lock_cands)][1]
if (is.na(lock_file) || !file.exists(lock_file)) stop("Cannot find 00_locks.R")
source(normalizePath(lock_file, winslash = "/"), local = TRUE)
root <- normalizePath(file.path(dirname(lock_file), "..", ".."), winslash = "/")
dd <- Sys.getenv("FIGURE_DATA_DIR", "")
cand <- c(dd, p("04_plot_data"))
data_dir <- cand[nzchar(cand) & dir.exists(cand)][1]
if (is.na(data_dir)) stop("Run 05_build_figure_data.R first.")
read_fig <- function(x) read_tsv(file.path(data_dir, "04_figure_4", x), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
od <- Sys.getenv("FIGURE_OUTPUT_DIR", "")
out <- if (nzchar(od)) od else p("06_panel_output", "04_figure_4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Gold uses a grey point stroke; fills/rings stay in the locked interface purple family.
outline <- "#748087"
lock_order <- c(core_axis_ids, extended_axis_ids)

x <- read_fig("07_panel_d_axis_descriptor.tsv")
for (col in c(
  "biological_axis_id", "axis_id", "interaction_name", "pathogen_gene", "host_gene",
  "parasite_factor_name", "host_target_name", "final_tier", "pathogen_program",
  "host_program", "proximity_pathogen_module", "proximity_host_module",
  "direct_physical_edge", "interface_score", "persistence", "rank",
  "aggregate_robustness_score", "stability_score", "threshold_consensus_frequency",
  "pcsf_edge_frequency", "pcsf_endpoint_min_frequency", "retained", "context_only"
)) {
  if (!col %in% names(x)) x[[col]] <- NA
}
x <- x %>%
  mutate(
    axis_id_key = axis_key(coalesce(biological_axis_id, axis_id), interaction_name, pathogen_gene, host_gene),
    final_tier = factor(axis_tier(axis_id_key, final_tier), levels = c("Core", "Extended")),
    context_only = flag_col(context_only),
    direct_physical_edge = flag_col(direct_physical_edge),
    pathogen_gene = endpoint_display(sanitize_label(coalesce(pathogen_gene, parasite_factor_name, interaction_name))),
    host_gene = endpoint_display(sanitize_label(coalesce(host_gene, host_target_name, interaction_name))),
    pathogen_program_short = pretty_module(coalesce(pathogen_program, proximity_pathogen_module)),
    host_program_short = pretty_module(coalesce(host_program, proximity_host_module)),
    # TSV interface_score is saturated RWR min percentile; rank on 108-run consensus frequency.
    interface_score = num_col(coalesce(threshold_consensus_frequency, stability_score)),
    persistence = num_col(coalesce(persistence, threshold_consensus_frequency, 0))
  ) %>%
  filter(
    axis_id_key %in% c(core_axis_ids, extended_axis_ids),
    final_tier %in% c("Core", "Extended"),
    !context_only
  ) %>%
  arrange(
    final_tier,
    desc(interface_score),
    match(axis_id_key, lock_order),
    axis_id_key
  ) %>%
  mutate(
    rank = row_number(),
    axis_text = sprintf(
      "%02d  %s\u2013%s\n%s\n%s",
      rank, pathogen_gene, host_gene, pathogen_program_short, host_program_short
    ),
    rank_block = final_tier,
    axis = factor(axis_text, levels = rev(unique(axis_text)))
  )

if (nrow(x) != 11L) {
  stop("Expected 11 retained biological axes; found ", nrow(x), ".")
}
if (sum(x$final_tier == "Core") != 3L || sum(x$final_tier == "Extended") != 8L) {
  stop("Expected 3 core + 8 extended axes.")
}
if (any(x$axis_id_key == "p104--CLASP1")) {
  stop("Context-only CLASP1 reached the ranking panel.")
}
if ((max(x$interface_score) - min(x$interface_score)) < 0.05) {
  stop("Ranking score is saturated across retained axes.")
}

x_min <- max(0, floor(min(x$interface_score, na.rm = TRUE) * 20) / 20 - 0.02)
x_max <- ceiling(max(x$interface_score, na.rm = TRUE) * 20) / 20 + 0.04

plot <- ggplot(x, aes(interface_score, axis)) +
  geom_segment(aes(x = x_min, xend = interface_score, yend = axis), colour = line, linewidth = 1.05, lineend = "round") +
  geom_point(aes(size = persistence, fill = final_tier), shape = 21, colour = outline, stroke = 0.38) +
  geom_point(data = filter(x, direct_physical_edge), shape = 21, fill = NA, colour = interface_dark, size = 5.0, stroke = 0.72) +
  geom_text(aes(label = sprintf("%.2f", interface_score)), hjust = -0.40, size = 2.0, colour = ink) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.01, 0.12))
  ) +
  scale_y_discrete(drop = TRUE, expand = expansion(add = 0.72)) +
  scale_fill_manual(
    values = c(Core = interface_fill, Extended = interface_pale),
    breaks = c("Core", "Extended"),
    name = "Interface tier"
  ) +
  scale_size_continuous(
    range = c(2.5, 4.25),
    limits = c(0, 1),
    breaks = c(0.5, 0.75, 1),
    labels = scales::percent,
    name = "Persistence"
  ) +
  facet_wrap(vars(rank_block), nrow = 1, scales = "free_y") +
  labs(x = "Integrated interface score", y = NULL) +
  guides(
    fill = guide_legend(
      order = 1, title.position = "left", nrow = 1,
      override.aes = list(
        shape = 21,
        fill = c(interface_fill, interface_pale),
        colour = c(interface_dark, interface_mid),
        size = 3.0, stroke = 0.55
      )
    ),
    size = guide_legend(
      order = 2, title.position = "left", nrow = 1,
      override.aes = list(shape = 21, fill = NA, colour = outline, stroke = 0.55)
    )
  ) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.35),
    panel.spacing.x = unit(5.5, "mm"),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 5.8, colour = ink),
    axis.text.y = element_text(size = 5.0, colour = ink, lineheight = 0.88),
    axis.title.x = element_text(size = 6.7, margin = margin(t = 4)),
    strip.background = element_rect(fill = interface_pale, colour = frame, linewidth = 0.35),
    strip.text = element_text(size = 6.1, face = "bold", colour = ink, margin = margin(2.0, 0, 2.0, 0)),
    legend.position = "bottom",
    legend.location = "plot",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.key = element_blank(),
    legend.title = element_text(size = 5.6, margin = margin(r = 3.2)),
    legend.text = element_text(size = 5.3),
    legend.key.width = unit(4.0, "mm"),
    legend.key.height = unit(2.8, "mm"),
    legend.spacing.x = unit(5.5, "mm"),
    legend.spacing.y = unit(1.2, "mm"),
    legend.box.spacing = unit(6.5, "mm"),
    legend.margin = margin(t = 1.2, r = 5.5, b = 0.4, l = 5.5),
    plot.margin = margin(5, 8, 5, 6)
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(file.path(out, "Figure_4_panel_e_interface_axis_ranking.pdf"), plot, width = 205, height = 112, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure_4_panel_e_interface_axis_ranking.png"), plot, width = 205, height = 112, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
}
cat(
  "Rendered Figure_4_panel_e_interface_axis_ranking with ", nrow(x),
  " retained biological axes (3 core + 8 extended); x = threshold_consensus_frequency.\n",
  sep = ""
)
cat(
  paste0(
    sprintf(
      "%02d  %s  score=%.4f  persistence=%.2f  %s\n",
      x$rank, x$axis_id_key, x$interface_score, x$persistence, as.character(x$final_tier)
    ),
    collapse = ""
  )
)
