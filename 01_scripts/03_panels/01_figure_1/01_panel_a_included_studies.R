#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) { z <- x$ofile; if (is.null(z) || !length(z)) "" else as.character(z[[1]]) }, character(1))
script_path <- if (any(nzchar(source_files))) normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE) else if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else normalizePath(file.path("01_scripts", "03_panels", "01_figure_1", "01_panel_a_included_studies.R"), winslash = "/", mustWork = TRUE)
data_env <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
root <- if (nzchar(data_env)) normalizePath(file.path(data_env, "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
p <- function(...) file.path(root, ...)
data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
if (nzchar(data_override) && !dir.exists(data_override)) stop("FIGURE_DATA_DIR does not exist: ", data_override)
data_candidates <- c(data_override, p("04_plot_data"))
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(file.path(data_dir, "01_figure_1", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "01_figure_1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) {
  if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x))
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), x, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(stem, ".png")), x, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)
}
source(p("01_scripts", "03_panels", "00_studies_helper.R"), local = TRUE)

pretty_host_family <- function(x) {
  case_when(
    str_detect(x, regex("pxd010713", ignore_case = TRUE)) ~ "TBL20 vs BL20 EV",
    str_detect(x, regex("pxd056906", ignore_case = TRUE)) ~ "Infected-cattle serum EV",
    str_detect(x, regex("dual_phospho", ignore_case = TRUE)) ~ "BL3 vs TA phospho",
    str_detect(x, regex("schizont_copur", ignore_case = TRUE)) ~ "Schizont co-purification",
    str_detect(x, regex("mc4404", ignore_case = TRUE)) ~ "eIF5A pulldown",
    str_detect(x, regex("bioid", ignore_case = TRUE)) ~ "CLASP1 BioID",
    TRUE ~ str_replace_all(x, "_", " ")
  )
}

group_levels <- c("Independent proteome", "Phosphoproteome", "Interface-proximal")
studies <- read_fig("01_panel_a_host_studies.tsv") %>%
  mutate(
    primary_family = if ("primary_family" %in% names(.)) as.logical(primary_family) else TRUE,
    group = case_when(
      str_detect(coalesce(omics_type, ""), regex("phospho", ignore_case = TRUE)) ~ "Phosphoproteome",
      coalesce(primary_family, FALSE) | integration_role == "independent_host_state" ~ "Independent proteome",
      TRUE ~ "Interface-proximal"
    ),
    group = factor(group, levels = intersect(group_levels, unique(group))),
    study_label = pretty_host_family(study_family)
  ) %>%
  arrange(group == "Interface-proximal", effect_rows) %>%
  mutate(study_label = factor(study_label, levels = unique(study_label)))

panel <- draw_studies_lollipop(
  studies,
  y = "study_label",
  group = "group",
  x_lab = "Detected host proteins",
  y_labels = identity
)
save_panel(panel, "Figure_1_panel_a_included_studies", 125, 82)
cat("Rendered Figure_1_panel_a_included_studies.\n")
