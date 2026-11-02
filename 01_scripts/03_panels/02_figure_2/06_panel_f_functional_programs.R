#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) { z <- x$ofile; if (is.null(z) || !length(z)) "" else as.character(z[[1]]) }, character(1))
script_path <- if (any(nzchar(source_files))) normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE) else if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE) else normalizePath(file.path("01_scripts", "03_panels", "02_figure_2", "06_panel_f_functional_programs.R"), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE); p <- function(...) file.path(root, ...)
data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = ""); if (nzchar(data_override) && !dir.exists(data_override)) stop("FIGURE_DATA_DIR does not exist: ", data_override)
data_candidates <- c(data_override, p("04_plot_data")); data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]; if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first."); data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(file.path(data_dir, "02_figure_2", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = ""); output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "02_figure_2"); dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) { if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x)); ggsave(file.path(output_dir, paste0(stem, ".pdf")), x, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE); ggsave(file.path(output_dir, paste0(stem, ".png")), x, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE) }

ink <- "#30363B"; blue_dark <- "#B5C8D1"; blue <- "#D2DEE4"; blue_pale <- "#F1F5F7"; red_dark <- "#CDB0AD"; red <- "#E4D0CD"; red_pale <- "#FAF1F0"
theme_panel <- theme_classic(base_size = 8, base_family = "sans") + theme(panel.grid = element_blank(), axis.line = element_blank(), panel.border = element_rect(fill = NA, colour = "#DCE3E7", linewidth = 0.35), axis.ticks = element_line(linewidth = 0.3, colour = ink), axis.text = element_text(colour = ink), axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title = element_text(colour = ink), legend.position = "bottom", legend.justification = "center", legend.box.just = "center", legend.direction = "horizontal", legend.box = "horizontal", legend.key.height = unit(2.4, "mm"), legend.key.width = unit(4, "mm"), plot.margin = margin(7, 9, 7, 9), panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA))

is_network_id <- function(x) {
  z <- str_replace_all(coalesce(as.character(x), ""), "\\s+", "_")
  str_detect(z, "^[HP]_NET(_|$)") | str_detect(z, "^[HP]_FUN(_|$)")
}

program_name <- function(x) {
  z <- coalesce(as.character(x), "")
  case_when(
    !nzchar(str_squish(z)) ~ NA_character_,
    is_network_id(z) ~ NA_character_,
    str_detect(z, "surface_export|traffick") ~ "Surface / export",
    str_detect(z, "membrane_bioenergetics|bioenergetic") ~ "Membrane bioenergetics",
    str_detect(z, "genome_expression") ~ "Genome-expression",
    str_detect(z, "cytoskeletal") ~ "Cytoskeletal",
    str_detect(z, "protein_folding") ~ "Protein folding",
    str_detect(z, "central_metabolism|cytoplasmic_metabolism") ~ "Central metabolism",
    str_detect(z, "translation_ribosome") ~ "Translation",
    str_detect(z, "rna_processing") ~ "RNA processing",
    str_detect(z, "proteostasis") ~ "Proteostasis",
    str_detect(z, "Unassigned|hypothetical|unresolved|other_function|cellular_process_broad|other_accessible") ~ "Other / unresolved",
    TRUE ~ NA_character_
  )
}

functions <- read_fig("04_panel_e_pathogen_functions.tsv")
if (!"module_biological_label" %in% names(functions) && "display_name" %in% names(functions)) {
  functions$module_biological_label <- functions$display_name
}
if (!"module_biological_label" %in% names(functions) && "module_label" %in% names(functions)) {
  functions$module_biological_label <- functions$module_label
}
if (!"module_biological_label" %in% names(functions)) functions$module_biological_label <- "Unassigned"
if (!"module_label" %in% names(functions)) functions$module_label <- NA_character_
if (!"dominant_function_class" %in% names(functions)) functions$dominant_function_class <- NA_character_
if (!"module_function_signature" %in% names(functions)) functions$module_function_signature <- NA_character_
functions <- functions %>%
  mutate(
    module_biological_label = coalesce(module_biological_label, "Unassigned"),
    program = coalesce(
      program_name(module_biological_label),
      program_name(module_label),
      program_name(dominant_function_class),
      program_name(module_function_signature),
      "Other / unresolved"
    )
  ) %>%
  group_by(program) %>%
  summarise(scaffold_n = sum(scaffold_n), core_n = sum(core_n), .groups = "drop") %>%
  filter(scaffold_n > 0)

total_scaffold <- sum(functions$scaffold_n)
total_core <- sum(functions$core_n)
enrichment <- functions %>%
  rowwise() %>%
  mutate(
    a = core_n,
    b = scaffold_n - core_n,
    c = total_core - core_n,
    d = total_scaffold - total_core - b,
    p_value = fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value,
    odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)),
    log_or_se = sqrt(1 / (a + 0.5) + 1 / (b + 0.5) + 1 / (c + 0.5) + 1 / (d + 0.5)),
    log2_odds_ratio = log2(odds_ratio),
    ci_low = log2_odds_ratio - 1.96 * log_or_se / log(2),
    ci_high = log2_odds_ratio + 1.96 * log_or_se / log(2)
  ) %>%
  ungroup() %>%
  mutate(
    fdr_bh = p.adjust(p_value, method = "BH"),
    significant = fdr_bh < 0.05,
    fdr_label = paste0("FDR ", format.pval(fdr_bh, digits = 2, eps = 0.001))
  )

program_order <- functions %>% arrange(scaffold_n) %>% pull(program)
composition <- functions %>% mutate(label = factor(program, levels = program_order))
enrichment <- enrichment %>% mutate(label = factor(program, levels = program_order))
any_fdr <- isTRUE(any(enrichment$significant, na.rm = TRUE))

p_composition <- ggplot(composition, aes(scaffold_n, label)) +
  geom_col(fill = red_pale, colour = red_dark, linewidth = 0.35, width = 0.60) +
  geom_point(aes(size = core_n), shape = 21, fill = red, colour = ink, stroke = 0.38) +
  scale_size_continuous(range = c(1.7, 5.0), breaks = pretty_breaks(n = 3), guide = guide_legend(nrow = 1, order = 1, title.position = "left", override.aes = list(shape = 21, fill = NA, colour = red_dark, stroke = 0.70))) +
  scale_y_discrete(labels = function(x) str_wrap(x, 24)) +
  labs(x = "Accessible scaffold proteins", y = NULL, size = "Core proteins") +
  theme_panel +
  theme(axis.text.y = element_text(size = 6.4, lineheight = 0.88))

p_enrichment <- ggplot(enrichment, aes(log2_odds_ratio, label)) +
  geom_vline(xintercept = 0, colour = "#E4E9EC", linewidth = 0.35) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y", width = 0.15, colour = blue_dark, linewidth = 0.48
  ) +
  geom_point(aes(size = scaffold_n), shape = 21, fill = if (any_fdr) red else "white", colour = ink, stroke = 0.42) +
  scale_size_continuous(range = c(2.1, 5.4), breaks = pretty_breaks(n = 3), guide = guide_legend(nrow = 1, order = 2, title.position = "left", override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.65))) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(x = "Core enrichment (log2 odds ratio)", y = NULL, size = "Scaffold proteins") +
  theme_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

panel <- p_composition + p_enrichment +
  plot_layout(widths = c(0.96, 1.18), guides = "collect") &
  theme(
    legend.position = "bottom", legend.justification = "center", legend.box.just = "center",
    legend.direction = "horizontal", legend.box = "horizontal", legend.key = element_blank(),
    legend.title = element_text(size = 5.6, margin = margin(r = 3.2)), legend.text = element_text(size = 5.3),
    legend.key.height = unit(2.2, "mm"), legend.key.width = unit(3.1, "mm"),
    legend.spacing.x = unit(3.2, "mm"), legend.spacing.y = unit(1.2, "mm"),
    legend.box.spacing = unit(6.5, "mm"), legend.margin = margin(t = 1.2, r = 4, b = 0.4, l = 4)
  )

save_panel(panel, "Figure_2_panel_f_functional_programs", 183, 92)
cat("Programs:", paste(composition$program[order(-composition$scaffold_n)], collapse = "; "), "\n")
cat("FDR range:", paste(range(enrichment$fdr_bh), collapse = " - "), "; any FDR<0.05:", any_fdr, "; FDR legend:", if (any_fdr) "shown" else "dropped", "\n")
cat("Rendered Figure_2_panel_f_functional_programs.\n")
