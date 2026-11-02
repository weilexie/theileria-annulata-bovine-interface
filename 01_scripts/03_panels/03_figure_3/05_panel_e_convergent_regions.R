#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
})

source_files <- vapply(sys.frames(), function(x) {
  z <- x$ofile
  if (is.null(z) || !length(z)) "" else as.character(z[[1]])
}, character(1))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (any(nzchar(source_files))) {
  normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path("01_scripts", "03_panels", "03_figure_3", "05_panel_e_convergent_regions.R"), winslash = "/", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
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

data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
data_candidates <- c(data_override, p("04_plot_data"))
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(
  file.path(data_dir, "03_figure_3", file),
  show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN")
)

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "03_figure_3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(plot, stem, width_mm = 210, height_mm = 156) {
  if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(plot))
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_png, bg = "white", limitsize = FALSE
  )
}

pretty_module_lock <- pretty_module
pretty_module <- function(z, wrap = FALSE) {
  raw <- sanitize_label(coalesce(as.character(z), "Unassigned"))
  aliased <- recode(
    raw,
    host_glycolysis_metabolism = "H_FUN_glycolysis_metabolism",
    host_mediator_transcription = "H_FUN_mediator_transcription",
    host_ubiquitin_cullin = "H_FUN_ubiquitin_cullin",
    host_chaperone_folding = "H_FUN_chaperone_folding",
    host_endocytosis_phagosome = "H_FUN_endocytosis_phagosome",
    host_proteasome = "H_FUN_proteasome",
    host_ribosome_translation = "H_FUN_ribosome_translation",
    host_spliceosome = "H_FUN_spliceosome",
    host_vesicle_COPII = "H_FUN_vesicle_COPII",
    host_microtubule_cytoskeleton = "H_FUN_cytoskeleton_microtubule",
    host_cytoskeleton_microtubule = "H_FUN_cytoskeleton_microtubule",
    host_signaling_broad = "H_NET_1",
    host_rna_processing = "H_NET_4",
    host_metabolism = "H_NET_5",
    host_chromatin_organization = "H_NET_7",
    host_ubiquitin_proteostasis = "H_NET_8",
    host_transcription_machinery = "H_NET_12",
    parasite_cellular_process_broad = "P_NET_1",
    parasite_cytoplasmic_metabolism = "P_NET_3",
    parasite_rna_processing = "P_NET_4",
    parasite_surface_family = "P_FUN_surface_family",
    parasite_PIase_PPI = "P_FUN_PIase_PPI",
    parasite_PPIase_PPI = "P_FUN_PIase_PPI",
    parasite_translation_folding = "P_FUN_translation_folding",
    parasite_DNA_RNA_binding = "P_FUN_DNA_RNA_binding",
    parasite_kinase_phosphatase = "P_FUN_kinase_phosphatase",
    parasite_membrane_transport = "P_FUN_membrane_transport",
    parasite_metabolism_ATP = "P_FUN_metabolism_ATP",
    TA08425 = "p104",
    `p104=TA08425` = "p104",
    TA_15795 = "TaPin1",
    `TA_15795?` = "TaPin1",
    .default = raw
  )
  out <- pretty_module_lock(aliased, wrap = wrap)
  out <- str_replace_all(out, "PIase", "PPIase")
  out <- str_replace_all(out, "Plase", "PPIase")
  out <- str_replace(out, "^host ", "")
  out <- str_replace(out, "^parasite ", "")
  out <- str_replace(out, "^FUN ", "")
  sanitize_label(out)
}

proximity <- read_fig("04_panel_e_module_proximity.tsv")
if (!"host_reachable_fraction" %in% names(proximity) && all(c("host_reachable_n", "host_size") %in% names(proximity))) {
  proximity <- proximity %>%
    mutate(
      host_reachable_fraction = num_col(host_reachable_n) / pmax(num_col(host_size), 1),
      pathogen_reachable_fraction = num_col(pathogen_reachable_n) / pmax(num_col(pathogen_size), 1)
    )
}
proximity <- proximity %>%
  mutate(
    host_size = num_col(host_size),
    pathogen_size = num_col(pathogen_size),
    host_reachable_fraction = num_col(host_reachable_fraction),
    pathogen_reachable_fraction = num_col(pathogen_reachable_fraction),
    proximity_fdr = num_col(proximity_fdr, default = 1)
  )
significant_pairs <- proximity %>% filter(is.finite(proximity_fdr), proximity_fdr < 0.05)
if (!nrow(significant_pairs)) stop("No FDR-significant functional module pairs were found.")

used_groups <- FALSE
group_file <- file.path(data_dir, "03_figure_3", "panel_e_program_groups.tsv")
if (file.exists(group_file)) {
  group_map <- read_tsv(group_file, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
  join_by <- intersect(names(significant_pairs), names(group_map))
  if (length(join_by)) {
    significant_pairs <- significant_pairs %>%
      left_join(group_map, by = join_by, suffix = c("", "_grp"))
    used_groups <- TRUE
  }
}

significant_pairs <- ensure_cols(
  significant_pairs,
  c("host_display_name", "pathogen_display_name", "host_module", "pathogen_module")
)
significant_pairs <- significant_pairs %>%
  mutate(
    host_pretty = pretty_module(coalesce(host_module, host_display_name)),
    pathogen_pretty = pretty_module(coalesce(pathogen_module, pathogen_display_name))
  )

host_regions <- significant_pairs %>%
  group_by(module_id = host_pretty) %>%
  summarise(
    significant_partner_n = n_distinct(pathogen_pretty),
    module_size = max(host_size, na.rm = TRUE),
    reachable_fraction = host_reachable_fraction[which.max(host_size)][[1]],
    .groups = "drop"
  ) %>%
  mutate(side = "Host")
pathogen_regions <- significant_pairs %>%
  group_by(module_id = pathogen_pretty) %>%
  summarise(
    significant_partner_n = n_distinct(host_pretty),
    module_size = max(pathogen_size, na.rm = TRUE),
    reachable_fraction = pathogen_reachable_fraction[which.max(pathogen_size)][[1]],
    .groups = "drop"
  ) %>%
  mutate(side = "Pathogen")

region_summary <- bind_rows(host_regions, pathogen_regions) %>%
  mutate(
    side = factor(side, levels = c("Host", "Pathogen")),
    module_label = paste0(as.character(side), "  ·  ", module_id)
  ) %>%
  arrange(side, significant_partner_n, module_size, module_id) %>%
  mutate(module_label = factor(module_label, levels = rev(unique(module_label))))

base_theme <- theme_classic(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(), axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.30),
    axis.title = element_text(size = 6.2, colour = ink),
    axis.text = element_text(size = 5.7, colour = ink),
    axis.ticks = element_line(colour = frame, linewidth = 0.25),
    legend.position = "none", legend.title = element_blank(),
    legend.direction = "horizontal", legend.key.width = unit(3.5, "mm"),
    legend.key.height = unit(2.2, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(4, 6, 4, 4)
  )

side_scale <- scale_colour_manual(values = c(Host = host_dark, Pathogen = pathogen_dark), guide = "none")
reach_min <- max(0, floor((min(region_summary$reachable_fraction, na.rm = TRUE) - 0.04) * 50) / 50)
reach_max <- 1.01
reach_breaks <- seq(0, 1, by = 0.02)
reach_breaks <- reach_breaks[reach_breaks >= reach_min - 1e-8 & reach_breaks <= 1 + 1e-8]
size_from <- 10
size_min <- min(region_summary$module_size, na.rm = TRUE)
if (is.finite(size_min) && size_min < 10) {
  size_from <- 10^floor(log10(pmax(size_min, 1e-6)))
}
partner_max <- max(region_summary$significant_partner_n, na.rm = TRUE)
partner_hi <- max(2 * ceiling(partner_max / 2), 2)
partner_breaks <- seq(0, partner_hi, by = 2)

p_partner <- ggplot(region_summary, aes(significant_partner_n, module_label, colour = side)) +
  geom_segment(aes(x = 0, xend = significant_partner_n, yend = module_label), colour = line, linewidth = 0.45) +
  geom_point(size = 2.15) +
  side_scale +
  scale_x_continuous(
    limits = c(0, partner_hi),
    breaks = partner_breaks,
    labels = function(x) format(round(x), trim = TRUE),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Significant partner programs", y = NULL) +
  base_theme +
  theme(
    axis.text.y = element_text(size = 5.55, colour = ink),
    plot.margin = margin(4, 4, 4, 6),
    legend.position = "none"
  )

p_reachable <- ggplot(region_summary, aes(reachable_fraction, module_label, colour = side)) +
  geom_segment(aes(x = reach_min, xend = reachable_fraction, yend = module_label), colour = line, linewidth = 0.45) +
  geom_point(size = 2.15) +
  side_scale +
  scale_x_continuous(
    limits = c(reach_min, reach_max),
    breaks = reach_breaks,
    labels = percent,
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Reachable proteins", y = NULL) +
  base_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

p_size <- ggplot(region_summary, aes(module_size, module_label, colour = side)) +
  geom_segment(aes(x = size_from, xend = module_size, yend = module_label), colour = line, linewidth = 0.45) +
  geom_point(size = 2.15) +
  side_scale +
  scale_x_log10(
    labels = comma,
    breaks = c(10, 100, 1000),
    expand = expansion(mult = c(0.04, 0.14))
  ) +
  coord_cartesian(xlim = c(size_from, NA), clip = "off") +
  labs(x = "Proteins in program", y = NULL) +
  base_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

panel <- (p_partner | p_reachable | p_size) +
  plot_layout(widths = c(1.65, 1.0, 1.0)) &
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

save_panel(panel, "Figure_3_panel_e_convergent_regions")
row_report <- region_summary %>%
  transmute(
    row = as.character(module_label),
    partners = significant_partner_n,
    size = module_size,
    reachable = round(100 * reachable_fraction, 2)
  )
cat(
  "Rendered Figure_3_panel_e_convergent_regions: ",
  nrow(significant_pairs), " significant pairs collapsed to ",
  nrow(region_summary), " convergence-hub programs; groups TSV ",
  if (used_groups) "used" else "absent",
  ".\n",
  sep = ""
)
cat(sprintf("%s | partners=%s | size=%s | reachable=%s%%\n", row_report$row, row_report$partners, row_report$size, row_report$reachable), sep = "")
