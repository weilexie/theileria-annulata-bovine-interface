#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

source_files <- vapply(
  sys.frames(),
  function(x) if (is.null(x$ofile)) "" else as.character(x$ofile[[1]]),
  character(1)
)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (any(nzchar(source_files))) {
  normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/")
} else if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
} else {
  normalizePath(
    "01_scripts/03_panels/03_figure_3/02_panel_b_rwr_convergence.R",
    winslash = "/"
  )
}

root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/")
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

figure_data_override <- Sys.getenv("FIGURE_DATA_DIR", "")
data_candidates <- c(
  figure_data_override,
  p("04_plot_data")
)
data_dir <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)][1]
if (is.na(data_dir)) stop("Run stage 09 before rendering Figure 3b.")

input_file <- file.path(data_dir, "03_figure_3", "03_panel_b_rwr_landscape.tsv")
if (!file.exists(input_file)) stop("Missing Figure 3b data: ", input_file)

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", "")
output_dir <- if (nzchar(output_override)) {
  output_override
} else {
  p("06_panel_output", "03_figure_3")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rw <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA", "NaN")
)
required <- c(
  "node_id", "layer", "direct_evidence", "opposite_signal", "recruited"
)
missing_cols <- setdiff(required, names(rw))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

rw <- rw %>%
  filter(layer %in% c("host", "pathogen")) %>%
  mutate(
    layer = factor(layer, levels = c("host", "pathogen")),
    recruited = as.logical(recruited),
    gate_member = direct_evidence < 0.50 & opposite_signal >= 0.95
  )

facet_key <- rw %>%
  group_by(layer) %>%
  summarise(
    node_n = n(),
    recruited_n = sum(recruited),
    .groups = "drop"
  ) %>%
  mutate(
    facet_label = sprintf(
      "%s  |  %s nodes  |  %s recruited",
      if_else(layer == "host", "Host", "Pathogen"),
      comma(node_n), comma(recruited_n)
    )
  )
facet_labels <- setNames(facet_key$facet_label, as.character(facet_key$layer))

recruited_nodes <- rw %>%
  filter(recruited) %>%
  arrange(node_id) %>%
  mutate(
    deterministic_index = row_number(),
    jitter_x = ((((deterministic_index * 73) %% 1009) / 1008) - 0.5) * 0.008,
    jitter_y = ((((deterministic_index * 193) %% 1013) / 1012) - 0.5) * 0.004,
    direct_plot = case_when(
      direct_evidence <= 0 ~ abs(jitter_x),
      direct_evidence >= 1 ~ 1 - abs(jitter_x),
      TRUE ~ pmin(pmax(direct_evidence + jitter_x, 0), 1)
    ),
    opposite_plot = case_when(
      opposite_signal <= 0 ~ abs(jitter_y),
      opposite_signal >= 1 ~ 1 - abs(jitter_y),
      TRUE ~ pmin(pmax(opposite_signal + jitter_y, 0), 1)
    )
  )

# Sit the two-line label on the gate line, in a gap that does not cover
# the recruited cluster (host ~32%, pathogen ~25%).
gate_labels <- tibble(
  layer = factor(c("host", "pathogen"), levels = levels(rw$layer)),
  x = c(0.16, 0.39),
  y = c(0.918, 0.925),
  label = "RWR recruitment\ngate"
)

gate_fill <- pathogen_pale
gate_line <- pathogen_dark

plot <- ggplot(rw, aes(direct_evidence, opposite_signal)) +
  annotate(
    "rect",
    xmin = 0, xmax = 0.50,
    ymin = 0.95, ymax = 1,
    fill = gate_fill,
    colour = NA,
    alpha = 0.72
  ) +
  stat_bin_2d(
    aes(
      fill = layer,
      alpha = after_stat(count / max(count))
    ),
    bins = c(58, 48),
    colour = NA,
    show.legend = FALSE
  ) +
  annotate(
    "segment",
    x = 0, xend = 0.50,
    y = 0.95, yend = 0.95,
    colour = gate_line,
    linewidth = 0.42,
    linetype = "22"
  ) +
  annotate(
    "segment",
    x = 0.50, xend = 0.50,
    y = 0.95, yend = 1,
    colour = gate_line,
    linewidth = 0.42,
    linetype = "22"
  ) +
  geom_label(
    data = gate_labels,
    aes(x, y, label = label),
    inherit.aes = FALSE,
    colour = gate_line,
    fill = scales::alpha("white", 0.90),
    linewidth = 0.18,
    label.padding = unit(0.7, "mm"),
    size = 1.62,
    fontface = "bold"
  ) +
  geom_point(
    data = recruited_nodes,
    aes(direct_plot, opposite_plot, colour = layer),
    inherit.aes = FALSE,
    shape = 16,
    size = 1.22,
    stroke = 0,
    alpha = 0.96,
    show.legend = FALSE
  ) +
  facet_wrap(
    vars(layer),
    nrow = 1,
    labeller = labeller(layer = facet_labels)
  ) +
  scale_fill_manual(
    values = c(host = host_fill, pathogen = pathogen_fill),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(host = host_dark, pathogen = pathogen_dark),
    guide = "none"
  ) +
  scale_alpha_continuous(
    range = c(0.32, 0.84),
    limits = c(0, 1),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = 0)
  ) +
  coord_fixed(ratio = 1, clip = "on") +
  labs(
    x = "Direct evidence percentile",
    y = "Opposite-side RWR percentile"
  ) +
  theme_classic(base_size = 8.0, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.38),
    axis.ticks = element_line(linewidth = 0.32, colour = ink),
    axis.text = element_text(colour = ink, size = 7.1),
    axis.title = element_text(colour = ink, size = 7.7),
    strip.background = element_rect(fill = "#F8F9FA", colour = frame, linewidth = 0.38),
    strip.text = element_text(colour = ink, size = 7.5, face = "bold", margin = margin(3, 2, 3, 2)),
    panel.spacing.x = unit(8, "mm"),
    legend.position = "none",
    plot.margin = margin(5, 7, 5, 7),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

stem <- file.path(output_dir, "Figure_3_panel_b_rwr_convergence")
if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(
  paste0(stem, ".pdf"), plot,
  width = 183, height = 93, units = "mm",
  device = cairo_pdf, bg = "white", limitsize = FALSE
)
ggsave(
  paste0(stem, ".png"), plot,
  width = 183, height = 93, units = "mm", dpi = 600,
  device = ragg::agg_png, bg = "white", limitsize = FALSE
)
}

gate_summary <- rw %>%
  group_by(layer) %>%
  summarise(
    nodes = n(),
    recruited = sum(recruited),
    recruited_in_gate = sum(recruited & gate_member),
    .groups = "drop"
  )
cat("Rendered Figure 3b using all ", comma(nrow(rw)), " nodes in the 2D bins.\n", sep = "")
print(gate_summary)
print(gate_labels %>% select(layer, x, y))
