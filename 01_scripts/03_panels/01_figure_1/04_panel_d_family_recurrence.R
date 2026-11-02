#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) { z <- x$ofile; if (is.null(z) || !length(z)) "" else as.character(z[[1]]) }, character(1))
script_path <- if (any(nzchar(source_files))) normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE) else if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else normalizePath(file.path("01_scripts", "03_panels", "01_figure_1", "04_panel_d_family_recurrence.R"), winslash = "/", mustWork = TRUE)
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

ink <- "#30363B"; up <- "#E4D0CD"; down <- "#D2DEE4"; none <- "#F4F6F7"
d_x_labels <- c(
  pxd010713_ev = "TBL20\nEV",
  pxd056906_ev = "Serum\nEV",
  dual_phospho = "Dual\nphospho",
  schizont_copur = "Schizont\ncopur",
  mc4404 = "eIF5A\npulldown",
  bioid = "CLASP1\nBioID",
  GSE250088_Cheeseman = "Cheeseman",
  GSE250088_Villares = "Villares",
  GSE119138 = "GSE119138"
)
theme_panel <- theme_classic(base_size = 8, base_family = "sans") + theme(panel.grid = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill = NA, colour = "#DCE3E7", linewidth = 0.35), axis.ticks = element_line(linewidth = 0.3, colour = ink), axis.text = element_text(colour = ink), axis.title = element_text(colour = ink), legend.position = "bottom", legend.justification = "center", legend.box.just = "center", legend.direction = "horizontal", legend.box = "horizontal", legend.key = element_blank(), legend.key.height = unit(2.8, "mm"), legend.key.width = unit(4.1, "mm"), legend.spacing.x = unit(1.4, "mm"), legend.box.spacing = unit(2.0, "mm"), legend.margin = margin(t = 1.5), plot.margin = margin(6, 8, 5, 8), panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA))

core <- display_core_catalog(read_fig, p)
recurrence <- display_core_recurrence_gold(core, p)
if (n_distinct(recurrence$gene_label) != nrow(core)) stop("F1d gene order drifted from the shared display-core catalog.")
if (!n_distinct(recurrence$vote) || !"Up" %in% recurrence$vote) stop("F1d recurrence matrix has no Up votes.")
dots <- recurrence %>% filter(vote %in% c("Up", "Down"), is.finite(sig), sig > 0)
fill_present <- intersect(c("Up", "Down"), as.character(unique(dots$vote)))
fill_cols <- c(Up = up, Down = down)[fill_present]
fill_edge <- c(Up = up, Down = down)[fill_present]

panel <- ggplot(recurrence, aes(study_family, gene_label)) +
  geom_point(data = dots, aes(fill = vote, size = sig), shape = 21, colour = "white", stroke = 0.28) +
  scale_fill_manual(values = fill_cols, breaks = fill_present, guide = guide_legend(nrow = 1, order = 1, override.aes = list(shape = 21, fill = NA, colour = unname(fill_edge), size = 2.4, stroke = 0.75))) +
  scale_size_continuous(name = expression(-log[10](FDR)), range = c(0.85, 2.15), limits = c(0, 6), breaks = c(2, 4, 6), guide = guide_legend(nrow = 1, order = 2, title.position = "left", override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.65))) +
  scale_x_discrete(labels = d_x_labels, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, fill = NULL) + theme_panel + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 5.15, lineheight = 0.86), axis.text.y = element_text(size = 4.4), legend.title = element_text(size = 5.9, margin = margin(r = 3)), legend.text = element_text(size = 6.0))

save_panel(panel, "Figure_1_panel_d_family_recurrence", 183, 165)
cat("Rendered Figure_1_panel_d_family_recurrence.\n")
