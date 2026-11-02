#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) { z <- x$ofile; if (is.null(z) || !length(z)) "" else as.character(z[[1]]) }, character(1))
script_path <- if (any(nzchar(source_files))) normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE) else if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else normalizePath(file.path("01_scripts", "03_panels", "01_figure_1", "02_panel_b_integration_funnel.R"), winslash = "/", mustWork = TRUE)
data_env <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
root <- if (nzchar(data_env)) normalizePath(file.path(data_env, "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE); p <- function(...) file.path(root, ...)
data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = ""); if (nzchar(data_override) && !dir.exists(data_override)) stop("FIGURE_DATA_DIR does not exist: ", data_override)
data_candidates <- c(data_override, p("04_plot_data")); data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]; if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first."); data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(file.path(data_dir, "01_figure_1", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = ""); output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "01_figure_1"); dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) { if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x)); ggsave(file.path(output_dir, paste0(stem, ".pdf")), x, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE); ggsave(file.path(output_dir, paste0(stem, ".png")), x, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE) }

ink <- "#30363B"; line <- "#E4E9EC"; blue <- "#D2DEE4"; blue_dark <- "#B5C8D1"; blue_pale <- "#F1F5F7"; red <- "#E4D0CD"; red_dark <- "#CDB0AD"; red_pale <- "#FAF1F0"
theme_panel <- theme_classic(base_size = 8, base_family = "sans") + theme(panel.grid = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill = NA, colour = "#DCE3E7", linewidth = 0.35), axis.ticks = element_line(linewidth = 0.3, colour = ink), axis.text = element_text(colour = ink), axis.title = element_text(colour = ink), legend.position = "none", plot.margin = margin(6, 8, 6, 8), panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA))

stage_labels <- c("Mapped\nproteins", "Independent\nvote", "Independent\nscaffold", "Display\ncore")
attrition <- read_fig("02_panel_b_host_attrition.tsv") %>%
  mutate(
    stage = factor(stage, levels = unique(stage)),
    x = row_number(),
    stage_label = if (n() == length(stage_labels)) stage_labels else str_replace_all(as.character(stage), " ", "\n"),
    fill_group = factor(
      c("Input", "Integrated", "Scaffold", "Core")[pmin(seq_len(n()), 4)],
      levels = c("Input", "Integrated", "Scaffold", "Core")
    ),
    node_label = comma(value)
  )
y_floor <- 10
if (is.finite(min(attrition$value, na.rm = TRUE)) && min(attrition$value, na.rm = TRUE) < 10) {
  y_floor <- 10^floor(log10(max(min(attrition$value, na.rm = TRUE), 1)))
}
steps <- attrition %>% arrange(x) %>% mutate(next_x = lead(x), next_value = lead(value), x_start = x + 0.18, x_end = next_x - 0.18, label_x = (x + next_x) / 2, label_y = sqrt(pmax(value, 1) * pmax(next_value, 1)) * 1.35, retained = next_value / value, step_label = paste0(percent(retained, accuracy = 0.1), " retained")) %>% filter(!is.na(next_x))
panel <- ggplot(attrition, aes(x, value)) +
  geom_segment(aes(xend = x, y = y_floor, yend = value), colour = line, linewidth = 0.55) +
  geom_curve(data = steps, aes(x = x_start, xend = x_end, y = value, yend = next_value), inherit.aes = FALSE, curvature = -0.08, colour = "#6F767D", linewidth = 0.35, arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_label(data = steps, aes(x = label_x, y = label_y, label = step_label), inherit.aes = FALSE, size = 2.0, linewidth = 0.12, label.padding = unit(0.1, "lines"), fill = "white", colour = "#6F767D") +
  geom_point(aes(size = value, fill = fill_group), shape = 21, colour = ink, stroke = 0.3) +
  geom_text(aes(label = node_label, colour = fill_group), vjust = -1.1, size = 2.5) +
  scale_x_continuous(breaks = attrition$x, labels = attrition$stage_label) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale()), limits = c(y_floor, max(attrition$value) * 3.8), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c(Input = "#F4F6F7", Integrated = blue_pale, Scaffold = red_pale, Core = red), guide = "none") +
  scale_colour_manual(values = c(Input = "#6F767D", Integrated = blue_dark, Scaffold = red_dark, Core = red_dark), guide = "none") +
  scale_size_continuous(range = c(2.1, 7.2), guide = "none") + coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Proteins / effects (log scale)") + theme_panel + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.4, lineheight = 0.88), plot.margin = margin(12, 8, 6, 8))

save_panel(panel, "Figure_1_panel_b_integration_funnel", 125, 82)
cat("Rendered Figure_1_panel_b_integration_funnel.\n")
