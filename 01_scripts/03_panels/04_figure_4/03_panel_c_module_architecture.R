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
script <- if (any(nzchar(ofiles))) normalizePath(tail(ofiles[nzchar(ofiles)], 1), winslash = "/") else if (length(arg)) normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/") else normalizePath("01_scripts/03_panels/04_figure_4/03_panel_c_module_architecture.R", winslash = "/")
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

# Named programs only; unnamed NET IDs collapse to Other / unresolved (Figure 2f honesty).
is_raw_net <- function(z) {
  z <- str_squish(coalesce(as.character(z), ""))
  !nzchar(z) |
    z %in% c("Unassigned", "unassigned") |
    str_detect(z, "^[HP]_NET(_|\\s*)?[0-9]+$") |
    str_detect(z, "^[HP] NET [0-9]+$") |
    str_detect(z, "^NET[ _]?[0-9]+$")
}

pretty_display <- function(id, label) {
  id <- coalesce(as.character(id), "")
  label <- coalesce(as.character(label), "")
  from_id <- pretty_module(id, wrap = FALSE)
  from_label <- pretty_module(label, wrap = FALSE)
  leftover <- function(z) {
    z <- str_squish(coalesce(z, ""))
    !nzchar(z) | z %in% c("Unassigned", id, label) | is_raw_net(z)
  }
  dplyr::case_when(
    !leftover(from_id) ~ from_id,
    !leftover(from_label) ~ from_label,
    TRUE ~ "Other / unresolved"
  )
}

x <- ensure_cols(
  read_fig("06_panel_c_pcsf_module_burden.tsv"),
  c("module_id", "module_label", "selected_nodes", "mean_selection_frequency", "endpoint_n", "hub_n", "total_strength")
) %>%
  mutate(
    layer = factor(layer, levels = c("host", "pathogen"), labels = c("Host", "Pathogen")),
    selected_nodes = num_col(selected_nodes),
    mean_selection_frequency = num_col(mean_selection_frequency),
    endpoint_n = num_col(endpoint_n),
    hub_n = num_col(hub_n),
    total_strength = num_col(total_strength),
    module_display = pretty_display(module_id, module_label),
    w_sel = pmax(selected_nodes, 1e-9)
  ) %>%
  group_by(layer, module_display) %>%
  summarise(
    mean_selection_frequency = weighted.mean(mean_selection_frequency, w = w_sel),
    selected_nodes = sum(selected_nodes),
    endpoint_n = sum(endpoint_n),
    hub_n = sum(hub_n),
    total_strength = sum(total_strength),
    .groups = "drop"
  ) %>%
  arrange(layer, selected_nodes, mean_selection_frequency, module_display) %>%
  mutate(
    y_key = paste(layer, module_display, sep = "|"),
    y_key = factor(y_key, levels = unique(y_key)),
    point_mm = scales::rescale(selected_nodes, to = c(3.0, 7.1))
  )

x_hi <- max(x$selected_nodes, na.rm = TRUE)
x_by <- if (x_hi >= 40) 10 else 5
size_breaks <- seq(x_by, x_by * floor(x_hi / x_by), by = x_by)
if (length(size_breaks) > 3) {
  size_breaks <- size_breaks[unique(round(seq(1, length(size_breaks), length.out = 3)))]
}

endpoint_modules <- x %>% filter(endpoint_n > 0)
hub_modules <- x %>% filter(hub_n > 0)

plot <- ggplot(x, aes(selected_nodes, y_key)) +
  geom_segment(aes(x = 0, xend = selected_nodes, yend = y_key), colour = line, linewidth = 1.0, lineend = "round") +
  geom_point(aes(size = selected_nodes, fill = mean_selection_frequency), shape = 21, colour = ink, stroke = 0.45) +
  geom_point(
    data = endpoint_modules,
    aes(selected_nodes, y_key, colour = "Interface endpoint"), inherit.aes = FALSE,
    shape = 21, fill = NA,
    size = endpoint_modules$point_mm + 1.25, stroke = 0.85
  ) +
  geom_point(
    data = hub_modules,
    aes(selected_nodes, y_key, colour = "Hub module"), inherit.aes = FALSE,
    shape = 21, fill = NA,
    size = hub_modules$point_mm + 0.55, stroke = 0.62
  ) +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(z) str_remove(z, "^[^|]+\\|"), expand = expansion(add = 0.55)) +
  scale_x_continuous(breaks = seq(0, x_hi + x_by, by = x_by), expand = expansion(mult = c(0.02, 0.16))) +
  scale_size_continuous(range = c(3.0, 7.1), breaks = size_breaks, name = "Selected proteins") +
  scale_fill_gradient(low = interface_pale, high = interface_mid, limits = c(0.5, 1), breaks = c(0.5, 0.75, 1), labels = scales::percent, name = "Mean persistence") +
  scale_colour_manual(values = c("Interface endpoint" = interface_dark, "Hub module" = ink), name = "Outer ring") +
  labs(x = "PCSF-selected proteins", y = NULL) +
  guides(
    size = guide_legend(
      order = 1, title.position = "top", title.hjust = 0.5, nrow = 1,
      keywidth = unit(6.2, "mm"), keyheight = unit(5.2, "mm"),
      override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.60)
    ),
    fill = guide_colourbar(
      order = 2, title.position = "top", title.hjust = 0.5,
      barwidth = unit(26, "mm"), barheight = unit(2.6, "mm")
    ),
    colour = guide_legend(
      order = 3, title.position = "top", title.hjust = 0.5, nrow = 1,
      override.aes = list(shape = 21, fill = NA, size = 3.1, stroke = 0.75)
    )
  ) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.35),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(colour = ink),
    axis.text.y = element_text(size = 5.8, lineheight = 0.90),
    axis.title.x = element_text(size = 6.6, margin = margin(t = 4)),
    strip.background = element_rect(fill = "#F8F9FA", colour = frame, linewidth = 0.35),
    strip.text = element_text(size = 6.8, face = "bold", colour = ink),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.justification = "center",
    legend.key = element_blank(),
    legend.key.width = unit(5.0, "mm"),
    legend.key.height = unit(4.4, "mm"),
    legend.spacing.x = unit(6.5, "mm"),
    legend.spacing.y = unit(1.5, "mm"),
    legend.box.spacing = unit(2.4, "mm"),
    legend.box.margin = margin(t = 1.0, r = 0, b = 0, l = 0),
    legend.title = element_text(size = 5.8, margin = margin(b = 1.6)),
    legend.text = element_text(size = 5.5),
    legend.margin = margin(t = 1.2, r = 1.5, b = 0, l = 1.5),
    plot.margin = margin(5, 6, 4, 5)
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(file.path(out, "Figure_4_panel_c_module_architecture.pdf"), plot, width = 140, height = 124, units = "mm", device = cairo_pdf, bg = "white")
ggsave(file.path(out, "Figure_4_panel_c_module_architecture.png"), plot, width = 140, height = 124, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
}
cat("Rendered Figure_4_panel_c_module_architecture.\n")
