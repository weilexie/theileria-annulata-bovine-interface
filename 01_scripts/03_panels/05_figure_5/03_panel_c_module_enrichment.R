#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
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
  normalizePath(file.path("01_scripts", "03_panels", "05_figure_5", "03_panel_c_module_enrichment.R"), winslash = "/", mustWork = TRUE)
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
data_dir <- data_candidates[
  nzchar(data_candidates) & dir.exists(data_candidates) &
    file.exists(file.path(data_candidates, "05_figure_5", "01_panel_a_final_nodes.tsv"))
][1]
if (is.na(data_dir)) stop("Figure-data directory not found.")
fig_dir <- file.path(data_dir, "figure_5")
read_fig <- function(file) read_tsv(
  file.path(fig_dir, file), show_col_types = FALSE, progress = FALSE,
  na = c("", "NA", "NaN")
)

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "05_figure_5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Named via pretty_module / table10. Leftover NET k is Other / unresolved, not a GO name.
is_leftover_module <- function(pretty, raw) {
  pretty <- str_squish(coalesce(as.character(pretty), ""))
  raw <- str_squish(coalesce(as.character(raw), ""))
  !nzchar(pretty) |
    pretty %in% c("Unassigned", "unassigned") |
    pretty == raw |
    str_detect(pretty, "^[HP]_NET(_|\\s*)?[0-9]+$") |
    str_detect(pretty, "^[HP] NET [0-9]+$") |
    str_detect(pretty, "^NET[ _]?[0-9]+$")
}

pretty_enrich <- function(x) {
  raw <- coalesce(as.character(x), "")
  lab <- pretty_module(raw, wrap = FALSE)
  lab <- if_else(is_leftover_module(lab, raw), "Other / unresolved", lab)
  lab <- gsub(" / ", " /\n", lab, fixed = TRUE)
  lab <- gsub("-broad", "-\nbroad", lab, fixed = TRUE)
  lab <- gsub("Cytoplasmic metabolism", "Cytoplasmic\nmetabolism", lab, fixed = TRUE)
  lab
}

finite_log2 <- function(z) {
  z <- suppressWarnings(as.numeric(z))
  out <- ifelse(is.finite(z) & z > 0, log2(z), NA_real_)
  out
}

prepare_enrichment <- function(file) {
  res <- read_fig(file) %>%
    filter(selected_n > 0) %>%
    mutate(
      module_name = pretty_enrich(module_label),
      leftover = is_leftover_module(pretty_module(module_label, wrap = FALSE), module_label)
    ) %>%
    # Named forest only. Unnamed NET clusters are the 5b Other / unresolved
    # sector; merging their OR/P would invent a statistic the TSV does not have.
    filter(!leftover) %>%
    arrange(desc(odds_ratio), desc(selected_n), module_label) %>%
    mutate(
      y_label = paste0(module_name, "  (n=", selected_n, ")"),
      y_id = sprintf("r%02d", row_number()),
      or_log = finite_log2(odds_ratio),
      lo_log = finite_log2(ci_low),
      hi_log = finite_log2(ci_high)
    )
  stopifnot(nrow(res) > 0)
  stopifnot(!any(str_detect(res$y_label, "H_NET|P_NET|NET [0-9]|Other")))

  y_labels <- setNames(res$y_label, res$y_id)
  res$y_id <- factor(res$y_id, levels = rev(res$y_id))

  finite_x <- c(res$lo_log, res$hi_log, res$or_log, 0)
  xlim <- range(finite_x, finite = TRUE, na.rm = TRUE)
  xlim <- c(min(xlim[[1]], -0.5) - 0.35, max(xlim[[2]], 0.5) + 0.45)
  p_max <- max(-log10(pmax(res$p_value, .Machine$double.xmin)), na.rm = TRUE)
  p_max <- max(1, ceiling(p_max))
  res <- res %>% mutate(
    or_plot = if_else(is.finite(or_log), pmin(pmax(or_log, xlim[[1]]), xlim[[2]]), NA_real_),
    ci_low_plot = if_else(is.finite(lo_log), pmax(lo_log, xlim[[1]]), NA_real_),
    ci_high_plot = if_else(is.finite(hi_log), pmin(hi_log, xlim[[2]]), NA_real_),
    p_score = -log10(pmax(p_value, .Machine$double.xmin)),
    p_plot = xlim[[1]] + p_score / p_max * diff(xlim),
    fdr_hit = is.finite(fdr) & fdr < 0.05
  )
  list(data = res, xlim = xlim, p_max = p_max, y_labels = y_labels)
}

ink <- "#30363B"
border <- "#DCE3E7"
host <- "#D0DEE5"
host_dark <- "#8099AA"
pathogen <- "#E2C2BF"
pathogen_dark <- "#AE7772"
core <- "#927D96"

host_res <- prepare_enrichment("05_panel_c_host_module_enrichment.tsv")
pathogen_res <- prepare_enrichment("06_panel_c_pathogen_module_enrichment.tsv")
size_max <- max(c(host_res$data$selected_n, pathogen_res$data$selected_n), 1, na.rm = TRUE)
size_breaks <- c(5, 10, 20)
size_breaks <- size_breaks[size_breaks <= size_max]
if (length(size_breaks) < 2) {
  size_breaks <- sort(unique(c(
    min(c(host_res$data$selected_n, pathogen_res$data$selected_n), na.rm = TRUE),
    size_max
  )))
}
any_fdr <- isTRUE(any(c(host_res$data$fdr_hit, pathogen_res$data$fdr_hit), na.rm = TRUE))

make_plot <- function(x, side, fill_colour, line_colour, show_legend = TRUE) {
  fdr_rows <- filter(x$data, fdr_hit, is.finite(or_plot))
  p <- ggplot(x$data, aes(y = y_id)) +
    geom_vline(xintercept = 0, colour = border, linewidth = 0.42) +
    geom_errorbar(
      aes(xmin = ci_low_plot, xmax = ci_high_plot),
      orientation = "y", width = 0.16, linewidth = 0.55,
      colour = line_colour, na.rm = TRUE
    ) +
    geom_point(
      aes(x = or_plot, size = selected_n, shape = "log2 OR (95% CI)"),
      fill = fill_colour, colour = ink, stroke = 0.35, na.rm = TRUE
    ) +
    geom_point(
      aes(x = p_plot, shape = "-log10(P)"),
      size = 1.65, fill = NA, colour = "#6F767D", stroke = 0.55, na.rm = TRUE
    )
  if (nrow(fdr_rows) > 0) {
    p <- p + geom_point(
      data = fdr_rows,
      aes(x = or_plot, size = selected_n, colour = "BH-FDR < 0.05"),
      shape = 21, fill = NA, stroke = 1.05, na.rm = TRUE
    ) + scale_colour_manual(values = c(`BH-FDR < 0.05` = core), name = NULL)
  }
  p +
    scale_shape_manual(
      values = c(`log2 OR (95% CI)` = 21, `-log10(P)` = 23),
      name = "Statistic"
    ) +
    scale_size_continuous(
      limits = c(1, size_max), range = c(2.2, 5.5),
      name = "Scaffold nodes", guide = "none"
    ) +
    scale_x_continuous(
      limits = x$xlim, breaks = pretty_breaks(5),
      labels = function(z) number(z, accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.04)),
      sec.axis = sec_axis(
        transform = ~ (. - x$xlim[[1]]) / diff(x$xlim) * x$p_max,
        breaks = pretty_breaks(if (x$p_max <= 4) 3 else 4),
        labels = function(z) number(z, accuracy = if (x$p_max <= 4) 0.5 else 1),
        name = "-log10(P)"
      )
    ) +
    scale_y_discrete(labels = x$y_labels) +
    labs(title = side, x = "log2 odds ratio (95% CI)", y = NULL) +
    guides(
      shape = guide_legend(
        order = 1, nrow = 1, title.position = "left",
        override.aes = list(
          size = c(2.0, 2.4), fill = NA,
          colour = c("#6F767D", ink), stroke = c(0.55, 0.35)
        )
      ),
      size = guide_legend(
        order = 2, nrow = 1, title.position = "left",
        override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.35)
      ),
      colour = if (any_fdr) {
        guide_legend(
          order = 3, nrow = 1, title.position = "left",
          override.aes = list(shape = 21, fill = NA, size = 3.2, stroke = 1.05)
        )
      } else {
        "none"
      }
    ) +
    theme_minimal(base_size = 7, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = border, linewidth = 0.32),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.24, colour = ink),
      axis.text = element_text(colour = ink),
      axis.text.x.top = element_text(colour = "#6F767D", size = 5.3),
      axis.text.y = element_text(size = 5.25, lineheight = 0.9),
      axis.title = element_text(colour = ink, size = 6.2),
      axis.title.x.top = element_text(colour = "#6F767D", size = 5.8),
      plot.title = element_text(colour = line_colour, size = 7, face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.title = element_text(size = 5.6, face = "plain"),
      legend.text = element_text(size = 5.2),
      legend.key.height = unit(2.4, "mm"),
      legend.key.width = unit(4.4, "mm"),
      legend.spacing.x = unit(0.7, "mm"),
      legend.box.spacing = unit(1.0, "mm"),
      plot.margin = margin(4, 5, 3, 5),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

host_plot <- make_plot(host_res, "Host modules", host, host_dark, show_legend = FALSE)
pathogen_plot <- make_plot(pathogen_res, "Pathogen modules", pathogen, pathogen_dark, show_legend = FALSE)

# A single drawn legend keeps all three groups on one baseline with fixed spacing.
legend_plot <- ggplot() +
  annotate("text", x = 0.00, y = 0.50, label = "Statistic", hjust = 0, size = 2.15, colour = ink, family = "sans") +
  annotate("point", x = 0.82, y = 0.50, shape = 23, size = 2.1, fill = NA, colour = "#6F767D", stroke = 0.55) +
  annotate("text", x = 1.08, y = 0.50, label = "-log10(P)", hjust = 0, size = 2.05, colour = ink, family = "sans") +
  annotate("point", x = 2.23, y = 0.50, shape = 21, size = 2.35, fill = NA, colour = ink, stroke = 0.35) +
  annotate("text", x = 2.52, y = 0.50, label = "log2 OR (95% CI)", hjust = 0, size = 2.05, colour = ink, family = "sans") +
  annotate("text", x = 4.35, y = 0.50, label = "Scaffold nodes", hjust = 0, size = 2.15, colour = ink, family = "sans") +
  annotate("point", x = 5.82, y = 0.50, shape = 21, size = 2.45, fill = NA, colour = ink, stroke = 0.35) +
  annotate("text", x = 6.12, y = 0.50, label = as.character(size_breaks[[1]]), hjust = 0, size = 2.05, colour = ink, family = "sans") +
  annotate("point", x = 6.72, y = 0.50, shape = 21, size = 3.55, fill = NA, colour = ink, stroke = 0.35) +
  annotate("text", x = 7.10, y = 0.50, label = as.character(size_breaks[[min(2, length(size_breaks))]]), hjust = 0, size = 2.05, colour = ink, family = "sans") +
  annotate("point", x = 7.72, y = 0.50, shape = 21, size = 4.65, fill = NA, colour = ink, stroke = 0.35) +
  annotate("text", x = 8.20, y = 0.50, label = as.character(size_breaks[[length(size_breaks)]]), hjust = 0, size = 2.05, colour = ink, family = "sans")
if (any_fdr) {
  legend_plot <- legend_plot +
    annotate("point", x = 9.18, y = 0.50, shape = 21, size = 3.15, fill = NA, colour = core, stroke = 1.05) +
    annotate("text", x = 9.58, y = 0.50, label = "BH-FDR < 0.05", hjust = 0, size = 2.15, colour = ink, family = "sans")
}
legend_plot <- legend_plot +
  coord_cartesian(xlim = c(0, if (any_fdr) 11.2 else 9.2), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(0, 4, 0, 4), plot.background = element_rect(fill = "white", colour = NA))

# Align the full gtables first, so unequal y-label widths do not change panel widths.
aligned_plots <- cowplot::align_plots(
  host_plot, pathogen_plot,
  align = "hv", axis = "tblr"
)
plots_row <- cowplot::plot_grid(
  plotlist = aligned_plots, nrow = 1,
  rel_widths = c(1, 1), align = "h", axis = "tb"
)
combined <- cowplot::plot_grid(
  plots_row, legend_plot, ncol = 1,
  rel_heights = c(1, 0.075)
)

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", unset = ""), "1")) {
  stem <- file.path(output_dir, "Figure_5_panel_c_module_enrichment")
  ggsave(paste0(stem, ".pdf"), combined, width = 183, height = 112, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), combined, width = 183, height = 112, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)
}

cat(
  "Rendered Figure_5_panel_c_module_enrichment with paired host/pathogen module enrichment.\n",
  "FDR < 0.05 rows: ",
  sum(c(host_res$data$fdr_hit, pathogen_res$data$fdr_hit), na.rm = TRUE),
  "; FDR legend: ", if (any_fdr) "shown" else "dropped", "\n",
  sep = ""
)
