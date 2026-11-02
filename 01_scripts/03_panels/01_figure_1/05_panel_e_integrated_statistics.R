#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) { z <- x$ofile; if (is.null(z) || !length(z)) "" else as.character(z[[1]]) }, character(1))
script_path <- if (any(nzchar(source_files))) normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE) else if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else normalizePath(file.path("01_scripts", "03_panels", "01_figure_1", "05_panel_e_integrated_statistics.R"), winslash = "/", mustWork = TRUE)
data_env <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
root <- if (nzchar(data_env)) normalizePath(file.path(data_env, "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE); p <- function(...) file.path(root, ...)
data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = ""); if (nzchar(data_override) && !dir.exists(data_override)) stop("FIGURE_DATA_DIR does not exist: ", data_override)
data_candidates <- c(data_override, p("04_plot_data")); data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]; if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first."); data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(file.path(data_dir, "01_figure_1", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = ""); output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "01_figure_1"); dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) { if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x)); ggsave(file.path(output_dir, paste0(stem, ".pdf")), x, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE); ggsave(file.path(output_dir, paste0(stem, ".png")), x, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE) }
helper_candidates <- unique(c(
  file.path(dirname(script_path), "00_display_core.R"),
  p("01_scripts", "03_panels", "01_figure_1", "00_display_core.R"),
  file.path(getwd(), "01_scripts", "03_panels", "01_figure_1", "00_display_core.R")
))
helper_candidates <- helper_candidates[file.exists(helper_candidates)]
if (!length(helper_candidates)) stop("Missing 00_display_core.R")
source(helper_candidates[[1]], local = TRUE)

ink <- "#30363B"; blue_dark <- "#D2DEE4"; red_dark <- "#E4D0CD"
theme_panel <- theme_classic(base_size = 8, base_family = "sans") + theme(
  panel.grid = element_blank(), axis.line = element_blank(),
  panel.border = element_rect(fill = NA, colour = "#DCE3E7", linewidth = 0.35),
  axis.ticks = element_line(linewidth = 0.3, colour = ink),
  axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
  legend.position = "bottom", legend.justification = "center", legend.box.just = "center",
  legend.direction = "horizontal", legend.box = "horizontal", legend.key = element_blank(),
  legend.title = element_text(size = 6.0, margin = margin(b = 2.2)), legend.text = element_text(size = 5.7),
  legend.box.spacing = unit(3.0, "mm"), legend.margin = margin(t = 1.5),
  plot.margin = margin(6, 8, 4, 8),
  panel.background = element_rect(fill = "white", colour = NA),
  plot.background = element_rect(fill = "white", colour = NA)
)

core <- display_core_catalog(read_fig, p)
statistics_long <- display_core_statistics_long(core)
if (n_distinct(statistics_long$gene_label) != nrow(core)) stop("F1e gene order drifted from the shared display-core catalog.")

panel <- ggplot(statistics_long, aes(metric, gene_label, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = blue_dark, mid = "#FAFAFA", high = red_dark, midpoint = 0.5, limits = c(0, 1),
    labels = percent,
    guide = guide_colourbar(
      direction = "horizontal", title.position = "top", title.hjust = 0.5,
      barwidth = unit(36, "mm"), barheight = unit(2.2, "mm")
    )
  ) +
  labs(x = NULL, y = NULL, fill = "Integrated statistic") +
  theme_panel +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.2), axis.text.y = element_text(size = 4.4))

save_panel(panel, "Figure_1_panel_e_integrated_statistics", 183, 210)
cat("Rendered Figure_1_panel_e_integrated_statistics.\n")
