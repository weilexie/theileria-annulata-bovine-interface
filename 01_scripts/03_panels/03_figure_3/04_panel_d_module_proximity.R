#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
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
  normalizePath(file.path("01_scripts", "03_panels", "03_figure_3", "04_panel_d_module_proximity.R"), winslash = "/", mustWork = TRUE)
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
save_panel <- function(plot, stem, width_mm = 144, height_mm = 156) {
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

ink <- "#30363B"
frame <- "#DCE3E7"
interface_pale <- "#F6F2F6"
interface_dark <- "#B7A6B8"

# pretty_module names only; leftover NET k -> Other; PIase reads as Plase.
axis_label <- function(z) {
  lab <- pretty_module(z, wrap = FALSE)
  lab <- if_else(str_detect(lab, "^NET[[:space:]]+[0-9]+$"), "Other", lab)
  str_replace_all(lab, "PIase", "PPIase")
}

# Gold two-line ticks: wrap at " / ", then str_wrap; do not hyphenate mid-word.
wrap_axis <- function(z) {
  z <- as.character(z)
  z <- str_replace(z, " / ", " /\n")
  single <- !str_detect(z, "\n")
  z[single] <- str_wrap(z[single], width = 12)
  still <- !str_detect(z, "\n")
  z[still] <- str_replace(z[still], "-broad", "-\nbroad")
  z
}

proximity <- read_fig("04_panel_e_module_proximity.tsv")
required <- c("host_module", "pathogen_module", "proximity_z", "empirical_p", "proximity_fdr")
missing <- setdiff(required, names(proximity))
if (length(missing)) stop("Module-proximity data are missing: ", paste(missing, collapse = ", "))

universe_n <- nrow(proximity)
universe_fdr <- sum(proximity$proximity_fdr < 0.05, na.rm = TRUE)

plot_data <- proximity %>%
  mutate(
    host_lab = axis_label(host_module),
    pathogen_lab = axis_label(pathogen_module),
    attraction = pmax(-proximity_z, 0),
    neg_log10_p = -log10(pmax(empirical_p, .Machine$double.xmin)),
    significant = proximity_fdr < 0.05
  ) %>%
  group_by(host_lab, pathogen_lab) %>%
  arrange(proximity_fdr, desc(attraction), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

host_order <- plot_data %>%
  group_by(host_lab) %>%
  summarise(
    significant_n = sum(significant, na.rm = TRUE),
    best_fdr = min(proximity_fdr, na.rm = TRUE),
    max_attraction = max(attraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(significant_n), best_fdr, desc(max_attraction)) %>%
  pull(host_lab)

pathogen_order <- plot_data %>%
  group_by(pathogen_lab) %>%
  summarise(
    significant_n = sum(significant, na.rm = TRUE),
    best_fdr = min(proximity_fdr, na.rm = TRUE),
    max_attraction = max(attraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(significant_n), best_fdr, desc(max_attraction)) %>%
  pull(pathogen_lab)

plot_data <- plot_data %>%
  mutate(
    host_plot = factor(host_lab, levels = rev(host_order)),
    pathogen_plot = factor(pathogen_lab, levels = pathogen_order)
  )

any_fdr <- isTRUE(any(plot_data$significant, na.rm = TRUE))

panel <- ggplot(plot_data, aes(pathogen_plot, host_plot)) +
  geom_point(
    aes(size = attraction, fill = neg_log10_p),
    shape = 21, colour = "#B0B9BF", stroke = 0.18, alpha = 0.84
  )

if (any_fdr) {
  panel <- panel +
    geom_point(
      data = filter(plot_data, significant),
      aes(size = attraction, shape = "FDR < 0.05"),
      fill = NA, colour = ink, stroke = 0.72
    ) +
    scale_shape_manual(values = c(`FDR < 0.05` = 21), name = NULL)
}

legend_guides <- list(
  size = guide_legend(
    order = 1, title.position = "top", title.hjust = 0.5, nrow = 1,
    keywidth = unit(6.2, "mm"), keyheight = unit(5.2, "mm"),
    override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.60)
  ),
  fill = guide_colourbar(
    order = 2, direction = "horizontal", title.position = "top", title.hjust = 0.5,
    barwidth = unit(32, "mm"), barheight = unit(2.6, "mm")
  )
)
if (any_fdr) {
  legend_guides$shape <- guide_legend(
    order = 3, nrow = 1, title.position = "top",
    keywidth = unit(6.2, "mm"), keyheight = unit(5.2, "mm"),
    override.aes = list(size = 3.2, fill = NA, colour = ink, stroke = 0.72)
  )
}

panel <- panel +
  scale_size_continuous(
    name = "Attraction (-z)", range = c(0.35, 5.0),
    breaks = c(2, 4, 6),
    limits = c(0, max(6, max(plot_data$attraction, na.rm = TRUE)))
  ) +
  scale_fill_gradient(
    name = expression(-log[10]("empirical P")),
    low = interface_pale, high = interface_dark,
    breaks = pretty_breaks(n = 3)
  ) +
  scale_x_discrete(labels = wrap_axis, drop = FALSE) +
  scale_y_discrete(labels = wrap_axis, drop = FALSE) +
  labs(x = "Pathogen program", y = "Host program") +
  do.call(guides, legend_guides) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.30),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.text.x = element_text(size = 5.35, lineheight = 0.86, hjust = 0.5, angle = 0, vjust = 1),
    axis.text.y = element_text(size = 5.35, lineheight = 0.86),
    axis.title = element_text(size = 6.6),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.key = element_blank(),
    legend.title = element_text(size = 5.7, margin = margin(b = 1.6)),
    legend.text = element_text(size = 5.5),
    legend.key.width = unit(5.2, "mm"),
    legend.key.height = unit(4.4, "mm"),
    legend.spacing.x = unit(10, "mm"),
    legend.spacing.y = unit(2.0, "mm"),
    legend.box.spacing = unit(3.2, "mm"),
    legend.box.margin = margin(t = 1.5, r = 2, b = 0.5, l = 2),
    legend.margin = margin(t = 1.2, r = 6, b = 0, l = 6),
    plot.margin = margin(6, 10, 7, 8)
  )

save_panel(panel, "Figure_3_panel_d_module_proximity")
cat(
  "Rendered Figure_3_panel_d_module_proximity as a dot matrix: ",
  nrow(plot_data), " display cells from ", universe_n, " module pairs; ",
  sum(plot_data$significant), " FDR<0.05 cells (universe ",
  universe_fdr, "/", universe_n, ").\n",
  "Host programs (", nlevels(plot_data$host_plot), "): ",
  paste(rev(levels(plot_data$host_plot)), collapse = " | "), "\n",
  "Pathogen programs (", nlevels(plot_data$pathogen_plot), "): ",
  paste(levels(plot_data$pathogen_plot), collapse = " | "), "\n",
  sep = ""
)
