#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
  library(cowplot)
  library(patchwork)
  library(readr)
  library(ragg)
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
  normalizePath(
    file.path("01_scripts", "03_panels", "05_figure_5", "05_panel_e_bpq_reversal.R"),
    winslash = "/", mustWork = TRUE
  )
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
data_candidates <- data_candidates[nzchar(data_candidates)]
data_dir <- data_candidates[file.exists(file.path(data_candidates, "05_figure_5", "01_panel_a_final_nodes.tsv"))][1]
if (is.na(data_dir)) stop("Run 01_scripts/02_build_figure_data.R first.")
read_fig <- function(file) {
  fp <- file.path(data_dir, "05_figure_5", file)
  if (!file.exists(fp)) return(tibble())
  read_tsv(fp, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
}

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "05_figure_5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(plot, stem, width_mm = 205, height_mm = 108) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
    units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot, width = width_mm, height = height_mm,
    units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE
  )
}

looks_like_gene <- function(x) {
  z <- str_squish(coalesce(as.character(x), ""))
  nzchar(z) &
    !str_detect(z, "^TA[0-9]{4,}$") &
    !str_detect(z, "^TA_") &
    !str_detect(z, "^(9913|5874)\\.") &
    !str_detect(z, "^ENS") &
    !str_detect(z, "^A0A") &
    !str_detect(z, "^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
}

tidy_gene <- function(x) {
  z <- endpoint_display(sanitize_label(x))
  z <- str_replace(z, "P68402$", "")
  stripped <- str_replace(z, "[OPQ][0-9][A-Z0-9]{3}[0-9]$", "")
  z <- if_else(nchar(stripped) >= 3L, stripped, z)
  str_squish(z)
}

is_ribo_symbol <- function(x) str_detect(coalesce(as.character(x), ""), "^(RPS|RPL)[A-Z0-9P]")

cd2ap_id <- "9913.ENSBTAP00000072401"
pkm_id <- "9913.ENSBTAP00000044619"
rbm39_id <- "9913.ENSBTAP00000060007"

layer_path <- p("05_extra_tables", "02_bpq_reversal_layers.tsv")
protein_path <- p("03_algorithm_tables", "09_bpq_reversal_site_protein_table.tsv")
if (!file.exists(layer_path)) stop("Missing ", layer_path)
if (!file.exists(protein_path)) stop("Missing ", protein_path)

layers <- read_tsv(layer_path, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
proteins <- read_tsv(protein_path, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
nodes <- read_fig("01_panel_a_final_nodes.tsv")
if (!nrow(nodes)) stop("Missing figure_5/panel_a_final_nodes.tsv for gene labels.")

for (col in c("display_label", "canonical_gene", "gene_symbol", "node_id")) {
  if (!col %in% names(nodes)) nodes[[col]] <- NA_character_
}
for (col in c(
  "protein_id", "string_id", "bpq_delta", "trans_delta", "reversal_consistent",
  "in_rwr_target_set", "in_axis_endpoint", "in_strict_endpoint", "in_persistent_pcsf"
)) {
  if (!col %in% names(proteins)) stop("Protein table missing column: ", col)
}

proteins <- proteins %>%
  mutate(
    string_id = as.character(string_id),
    protein_id = as.character(protein_id),
    bpq_delta = suppressWarnings(as.numeric(bpq_delta)),
    trans_delta = suppressWarnings(as.numeric(trans_delta)),
    in_rwr_target_set = flag_col(in_rwr_target_set),
    in_axis_endpoint = flag_col(in_axis_endpoint),
    in_strict_endpoint = flag_col(in_strict_endpoint),
    in_persistent_pcsf = flag_col(in_persistent_pcsf),
    reversal_consistent = flag_col(reversal_consistent)
  )

both <- proteins %>% filter(is.finite(bpq_delta), is.finite(trans_delta))
computed_rev <- both$trans_delta * both$bpq_delta < 0
if (!identical(as.logical(both$reversal_consistent), computed_rev)) {
  stop("reversal_consistent does not match opposite-sign trans_delta * bpq_delta < 0.")
}

layer_lookup <- layers %>%
  mutate(
    n_with_transformation = as.integer(n_with_transformation),
    n_reversal_consistent = as.integer(n_reversal_consistent),
    reversal_fraction = as.numeric(reversal_fraction)
  )
need_layers <- c("axis_endpoint", "persistent_pcsf", "rwr_target_set", "all_bpq")
if (!all(need_layers %in% layer_lookup$layer)) {
  stop("table05 is missing a plotted layer: ", paste(setdiff(need_layers, layer_lookup$layer), collapse = ", "))
}

count_layer <- function(mask) {
  d <- both %>% filter(.data[[mask]])
  c(n_rev = sum(d$reversal_consistent), n_tot = nrow(d))
}
observed <- list(
  all_bpq = c(n_rev = sum(both$reversal_consistent), n_tot = nrow(both)),
  rwr_target_set = count_layer("in_rwr_target_set"),
  persistent_pcsf = count_layer("in_persistent_pcsf"),
  axis_endpoint = count_layer("in_axis_endpoint")
)
for (key in names(observed)) {
  row <- layer_lookup %>% filter(layer == key)
  if (observed[[key]][["n_tot"]] != row$n_with_transformation[[1]] ||
      observed[[key]][["n_rev"]] != row$n_reversal_consistent[[1]]) {
    stop(sprintf(
      "Protein table %s is %d/%d but table05 is %d/%d",
      key, observed[[key]][["n_rev"]], observed[[key]][["n_tot"]],
      row$n_reversal_consistent[[1]], row$n_with_transformation[[1]]
    ))
  }
}

label_map <- nodes %>%
  transmute(
    string_id = as.character(node_id),
    display_label = tidy_gene(coalesce(display_label, canonical_gene, gene_symbol))
  ) %>%
  filter(nzchar(string_id)) %>%
  group_by(string_id) %>%
  slice_head(n = 1) %>%
  ungroup()

proteins <- proteins %>%
  left_join(label_map, by = "string_id") %>%
  mutate(
    display_label = tidy_gene(coalesce(display_label, NA_character_)),
    display_label = if_else(looks_like_gene(display_label), display_label, NA_character_)
  )

both <- proteins %>% filter(is.finite(bpq_delta), is.finite(trans_delta))
pcsf_both <- both %>% filter(in_persistent_pcsf)
if (nrow(both) != 262L) stop("Expected 262 proteins with both contrasts, got ", nrow(both))
if (nrow(pcsf_both) != 22L) stop("Expected 22 PCSF hosts with both contrasts, got ", nrow(pcsf_both))
if (sum(pcsf_both$reversal_consistent) != 12L) {
  stop("Expected 12 reversal-consistent PCSF hosts, got ", sum(pcsf_both$reversal_consistent))
}

cd2ap <- both %>% filter(string_id == cd2ap_id)
if (nrow(cd2ap) != 1L) stop("CD2AP is missing from the BPQ table.")
if (!isTRUE(cd2ap$reversal_consistent[[1]]) || !isTRUE(cd2ap$in_persistent_pcsf[[1]])) {
  stop("CD2AP must be a reversal-consistent PCSF host in the BPQ table.")
}
if (!identical(cd2ap$display_label[[1]], "CD2AP")) {
  proteins$display_label[proteins$string_id == cd2ap_id] <- "CD2AP"
  both$display_label[both$string_id == cd2ap_id] <- "CD2AP"
  pcsf_both$display_label[pcsf_both$string_id == cd2ap_id] <- "CD2AP"
  cd2ap$display_label[[1]] <- "CD2AP"
}

bpq_only <- proteins %>%
  filter(string_id %in% c(pkm_id, rbm39_id), is.finite(bpq_delta), !is.finite(trans_delta)) %>%
  mutate(display_label = recode(string_id, !!pkm_id := "PKM", !!rbm39_id := "RBM39", .default = display_label)) %>%
  arrange(match(string_id, c(pkm_id, rbm39_id)))
if (nrow(bpq_only) != 2L || !setequal(bpq_only$display_label, c("PKM", "RBM39"))) {
  stop("PKM and RBM39 must appear as BPQ-only (no transformation contrast).")
}

label_pool <- pcsf_both %>%
  filter(looks_like_gene(display_label), !is_ribo_symbol(display_label)) %>%
  mutate(mag = pmax(abs(trans_delta), abs(bpq_delta))) %>%
  arrange(desc(string_id == cd2ap_id), desc(mag), display_label) %>%
  distinct(string_id, .keep_all = TRUE) %>%
  slice_head(n = 8)
if (!"CD2AP" %in% label_pool$display_label) stop("CD2AP was dropped from scatter labels.")
if (any(str_detect(label_pool$display_label, "TA_15795")) ||
    any(str_detect(label_pool$display_label, "[\u4e00-\u9fff]"))) {
  stop("Forbidden label reached the BPQ scatter.")
}

in_node_label <- function(z) {
  z <- as.character(z)
  if_else(str_detect(z, "^TA[0-9]{5}$"), str_replace(z, "^(TA)([0-9]+)$", "\\1\n\\2"), z)
}

# Two retained core axes that appear in the BPQ table. The table05 2/2
# layer also contains context-only p104–CLASP1; that pair stays off Figure 5.
axis_pairs <- tibble(
  axis_id_key = c("TA20980--CD2AP", "TaPin1--PKM"),
  interaction = c("TA20980-CD2AP", "TaPin1-PKM"),
  host_target = c("CD2AP", "PKM"),
  pathogen_factor = c("TA20980", "TaPin1"),
  mechanism = c("Trafficking / signaling", "Metabolic / proteostasis"),
  note = c("reversal-consistent", "BPQ only"),
  cx = c(0.25, 0.75),
  cy = 0.50,
  angle = c(8, -8)
) %>%
  mutate(
    half_length = 0.12,
    angle_rad = angle * pi / 180,
    x_host = cx - half_length * cos(angle_rad),
    x_pathogen = cx + half_length * cos(angle_rad),
    y_host = cy - half_length * sin(angle_rad) * 0.83,
    y_pathogen = cy + half_length * sin(angle_rad) * 0.83,
    label_x = cx,
    label_y = cy + 0.155,
    note_y = cy - 0.175,
    host_dot = in_node_label(host_target),
    pathogen_dot = in_node_label(pathogen_factor),
    mechanism_colour = unname(mechanism_cols[mechanism])
  )
if (any(str_detect(paste(axis_pairs$interaction, axis_pairs$host_target), "CLASP1|TA_15795"))) {
  stop("Forbidden axis reached the BPQ pair diagrams.")
}

pair_plot <- ggplot(axis_pairs) +
  annotate(
    "text", x = 0.50, y = 0.92, label = "BPQ-measured axes  |  n=2",
    hjust = 0.5, vjust = 1, family = "sans", size = 2.02,
    fontface = "bold", colour = ink
  ) +
  geom_segment(
    aes(x = x_host, xend = x_pathogen, y = y_host, yend = y_pathogen, colour = mechanism_colour),
    linewidth = 1.20, lineend = "round", show.legend = FALSE
  ) +
  geom_point(
    aes(x = x_host, y = y_host), shape = 21, size = 6.4, stroke = 0.50,
    fill = host_fill, colour = "#8099AA"
  ) +
  geom_point(
    aes(x = x_pathogen, y = y_pathogen), shape = 21, size = 6.4, stroke = 0.50,
    fill = pathogen_fill, colour = "#AE7772"
  ) +
  geom_text(
    aes(x = x_host, y = y_host, label = host_dot),
    family = "sans", colour = "#8099AA", fontface = "bold", size = 1.34, lineheight = 0.78
  ) +
  geom_text(
    aes(x = x_pathogen, y = y_pathogen, label = pathogen_dot),
    family = "sans", colour = "#AE7772", fontface = "bold", size = 1.22, lineheight = 0.78
  ) +
  geom_text(
    aes(x = label_x, y = label_y, label = interaction, angle = angle),
    family = "sans", colour = ink, size = 1.55, fontface = "bold"
  ) +
  geom_text(
    aes(x = cx, y = note_y, label = note),
    family = "sans", colour = host_dark, size = 1.50
  ) +
  scale_colour_identity() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.12, 0.96), clip = "off") +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    plot.margin = margin(4, 6, 0, 4),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

gauge_meta <- tibble(
  layer = c("persistent_pcsf", "rwr_target_set", "all_bpq"),
  name = c("PCSF\npersistent host", "RWR\ntarget set", "All BPQ-\nmeasured proteins"),
  fill = c(interface_fill, host_dark, host_fill)
)
gauges <- gauge_meta %>%
  left_join(
    layer_lookup %>% select(layer, n_with_transformation, n_reversal_consistent, reversal_fraction),
    by = "layer"
  ) %>%
  mutate(
    count_lab = paste(n_reversal_consistent, "/", n_with_transformation),
    name = factor(name, levels = rev(name))
  )
if (any(!is.finite(gauges$reversal_fraction))) stop("Gauge fractions missing from table05.")

bg_frac <- gauges$reversal_fraction[gauges$layer == "all_bpq"][[1]]
bg_lab <- sprintf("%.1f%%", 100 * bg_frac)

bar_plot <- ggplot(gauges, aes(y = name)) +
  geom_col(aes(x = 1), fill = line, width = 0.62, colour = NA) +
  geom_col(aes(x = reversal_fraction, fill = fill), width = 0.62, colour = NA, show.legend = FALSE) +
  geom_vline(xintercept = bg_frac, colour = ink, linewidth = 0.32, linetype = "22", alpha = 0.58) +
  annotate(
    "text",
    x = bg_frac - 0.012, y = 3.58, label = paste0(bg_lab, " background"),
    colour = ink, size = 1.55, family = "sans", hjust = 1, vjust = 0.5
  ) +
  geom_text(
    aes(x = 0.985, label = count_lab),
    hjust = 1, vjust = 0.5, family = "sans", size = 2.05, colour = ink
  ) +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0%", "50%", "100%"),
    expand = c(0, 0),
    name = "Reversal-consistent fraction"
  ) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.42, 0.70))) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  labs(y = NULL) +
  theme_classic(base_size = 7.0, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.34),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = ink, linewidth = 0.24),
    axis.text.y = element_text(colour = ink, size = 5.50, hjust = 1, lineheight = 0.84, margin = margin(r = 3.2)),
    axis.text.x = element_text(colour = ink, size = 5.4, margin = margin(t = 1.6)),
    axis.title.x = element_text(colour = ink, size = 6.0, margin = margin(t = 2.2)),
    axis.ticks.length.x = unit(1.2, "mm"),
    plot.margin = margin(0, 6, 0, 4),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

e_lim <- 4.0
lim <- e_lim
x_lim <- c(-e_lim, e_lim)
y_lim <- c(-e_lim, e_lim)
quads <- tibble(
  xmin = c(x_lim[1], 0),
  xmax = c(0, x_lim[2]),
  ymin = c(0, y_lim[1]),
  ymax = c(y_lim[2], 0)
)
quad_lab <- tibble(
  trans_delta = c(-0.90 * e_lim, 0.90 * e_lim),
  bpq_delta = c(0.86 * e_lim, -0.86 * e_lim),
  lab = "reversal",
  hjust = c(0, 1)
)
legend_df <- tibble(
  trans_delta = NA_real_,
  bpq_delta = NA_real_,
  cls = factor(
    c("Background", "PCSF", "Reversal-consistent PCSF"),
    levels = c("Background", "PCSF", "Reversal-consistent PCSF")
  )
)
right_plot <- ggplot() +
  geom_rect(
    data = quads,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = interface_pale, alpha = 0.55, colour = NA, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, colour = frame, linewidth = 0.32) +
  geom_vline(xintercept = 0, colour = frame, linewidth = 0.32) +
  geom_text(
    data = quad_lab,
    aes(trans_delta, bpq_delta, label = lab, hjust = hjust),
    vjust = 0.5, family = "sans", size = 1.68, colour = interface_dark, alpha = 0.85
  ) +
  geom_point(
    data = both,
    aes(trans_delta, bpq_delta),
    colour = frame, size = 0.70, alpha = 0.45, stroke = 0, show.legend = FALSE
  ) +
  geom_point(
    data = pcsf_both,
    aes(trans_delta, bpq_delta),
    shape = 21, fill = NA, colour = interface_dark, size = 2.10, stroke = 0.55, show.legend = FALSE
  ) +
  geom_point(
    data = filter(pcsf_both, reversal_consistent),
    aes(trans_delta, bpq_delta),
    shape = 21, fill = interface_fill, colour = "white", size = 2.15, stroke = 0.55, show.legend = FALSE
  ) +
  geom_point(
    data = legend_df,
    aes(trans_delta, bpq_delta, shape = cls, fill = cls, colour = cls, size = cls),
    stroke = 0.55, na.rm = TRUE
  ) +
  geom_text_repel(
    data = label_pool,
    aes(trans_delta, bpq_delta, label = display_label),
    colour = ink, size = 1.68, family = "sans", seed = 20260904,
    box.padding = 0.22, point.padding = 0.18, min.segment.length = 0.12,
    segment.colour = ink, segment.size = 0.22, segment.alpha = 0.50,
    max.overlaps = Inf, force = 1.0, max.time = 2,
    bg.color = "white", bg.r = 0.08,
    xlim = x_lim, ylim = y_lim, show.legend = FALSE
  ) +
  scale_shape_manual(
    values = c(Background = 16, PCSF = 21, `Reversal-consistent PCSF` = 21),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(Background = frame, PCSF = "white", `Reversal-consistent PCSF` = interface_fill),
    name = NULL
  ) +
  scale_colour_manual(
    values = c(Background = frame, PCSF = interface_dark, `Reversal-consistent PCSF` = "white"),
    name = NULL
  ) +
  scale_size_manual(
    values = c(Background = 1.15, PCSF = 2.10, `Reversal-consistent PCSF` = 2.15),
    name = NULL
  ) +
  guides(
    shape = guide_legend(nrow = 1, override.aes = list(
      size = c(1.55, 2.15, 2.15),
      fill = c(frame, "white", interface_fill),
      colour = c(host_dark, interface_dark, "white"),
      stroke = c(0.35, 0.55, 0.55)
    )),
    fill = "none",
    colour = "none",
    size = "none"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(5), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(5), expand = c(0, 0)) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, clip = "on") +
  labs(
    x = "Transformation  (TA vs BL3), log2",
    y = "BPQ  (treated vs infected), log2"
  ) +
  theme_classic(base_size = 7.0, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.34),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.22, colour = ink),
    axis.text = element_text(colour = ink, size = 5.4, margin = margin(t = 1.6)),
    axis.title.x = element_text(colour = ink, size = 6.0, margin = margin(t = 2.2)),
    axis.title.y = element_text(colour = ink, size = 6.0),
    axis.ticks.length.x = unit(1.2, "mm"),
    legend.position = c(0.50, 0.012),
    legend.justification = c(0.50, 0),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(colour = ink, size = 4.70),
    legend.key = element_blank(),
    legend.key.width = unit(3.4, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.spacing.x = unit(1.1, "mm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    legend.background = element_rect(fill = alpha("white", 0.88), colour = NA),
    plot.margin = margin(4, 8, 0, 4),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# Original layout: pairs over bars on the left, scatter spanning the right.
# Shared patchwork row so the bar panel.border bottom and scatter
# panel.border bottom sit on the same baseline.
combined <- patchwork::wrap_plots(
  pair_plot, right_plot, bar_plot,
  design = c(
    patchwork::area(1, 1),
    patchwork::area(1, 2, 2, 2),
    patchwork::area(2, 1)
  ),
  widths = c(0.80, 1.02),
  heights = c(0.26, 0.74)
) +
  patchwork::plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", unset = ""), "1")) {
  save_panel(combined, "Figure_5_panel_e_bpq_reversal")
}
cat(sprintf(
  paste0(
    "Rendered Figure_5_panel_e_bpq_reversal: gauges %s; ",
    "scatter n=%d background / %d PCSF / %d reversal PCSF; ",
    "labels %s; BPQ-only %s.\n"
  ),
  paste(sprintf("%s %s (%.1f%%)", gauges$layer, gauges$count_lab, 100 * gauges$reversal_fraction), collapse = "; "),
  nrow(both), nrow(pcsf_both), sum(pcsf_both$reversal_consistent),
  paste(label_pool$display_label, collapse = ", "),
  paste(axis_pairs$interaction, collapse = ", ")
))
