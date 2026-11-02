#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

source_files <- vapply(
  sys.frames(),
  function(x) if (is.null(x$ofile)) "" else as.character(x$ofile[[1]]),
  character(1)
)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (any(nzchar(source_files))) {
  normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path("01_scripts", "03_panels", "03_figure_3", "01_panel_a_module_bridge_map.R"),
    winslash = "/", mustWork = TRUE
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

ink <- "#30363B"
line <- "#E4E9EC"
host_dark <- "#B5C8D1"
host_fill <- "#D2DEE4"
host_pale <- "#F1F5F7"
pathogen_dark <- "#CDB0AD"
pathogen_fill <- "#E4D0CD"
pathogen_pale <- "#FAF1F0"

data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
data_candidates <- c(
  data_override,
  p("04_plot_data")
)
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
data_dir <- data_candidates[[1]]
read_fig <- function(file) {
  read_tsv(
    file.path(data_dir, "03_figure_3", file),
    show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN")
  )
}

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) {
  output_override
} else {
  p("06_panel_output", "03_figure_3")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(plot, stem, width_mm = 142, height_mm = 104) {
  if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(plot))
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_png, bg = "white", limitsize = FALSE
  )
}

module_nodes <- ensure_cols(
  read_fig("01_panel_a_module_nodes.tsv"),
  c("module_label", "module_id", "side", "scaffold_n", "core_n")
) %>%
  mutate(
    module_id = na_if(str_squish(as.character(coalesce(module_id, ""))), ""),
    module_label = na_if(str_squish(as.character(coalesce(module_label, ""))), ""),
    module_key = coalesce(module_id, module_label, "__unassigned__"),
    side = str_to_lower(as.character(side)),
    scaffold_n = num_col(scaffold_n),
    core_n = num_col(core_n),
    pretty = pretty_module(coalesce(module_id, module_label), wrap = FALSE)
  )
module_bridges <- ensure_cols(
  read_fig("02_panel_a_module_bridges.tsv"),
  c("host_network_module", "pathogen_network_module", "host_module", "pathogen_module", "bridge_weight", "axis_n")
) %>%
  mutate(
    host_key = coalesce(host_network_module, host_module, "__unassigned__"),
    pathogen_key = coalesce(pathogen_network_module, pathogen_module, "__unassigned__"),
    host_key = if_else(nzchar(as.character(host_key)), as.character(host_key), "__unassigned__"),
    pathogen_key = if_else(nzchar(as.character(pathogen_key)), as.character(pathogen_key), "__unassigned__"),
    bridge_weight = num_col(bridge_weight, 1),
    axis_n = num_col(axis_n, 1)
  )

if (!nrow(module_nodes)) stop("Figure 3a requires module nodes.")
if (!nrow(module_bridges)) stop("Figure 3a requires module bridges.")
if (any(!is.finite(module_nodes$scaffold_n)) || any(module_nodes$scaffold_n <= 0)) {
  stop("Every module requires a positive scaffold_n.")
}

n_modules_full <- nrow(module_nodes)
n_bridges_full <- nrow(module_bridges)
host_total <- sum(module_nodes$scaffold_n[module_nodes$side == "host"])
pathogen_total <- sum(module_nodes$scaffold_n[module_nodes$side == "pathogen"])

# Display collapse: keep named programs with real bridges; remainder -> Other / unassigned.
# Same pretty name (e.g. microtubule FUN + NET) merges. Do not copy gold 7/8 node counts.
host_keep_pretty <- c(
  "Signaling-broad",
  "Microtubule / cytoskeleton",
  "RNA processing",
  "Ubiquitin / proteostasis",
  "Metabolism"
)
pathogen_keep_pretty <- c(
  "Cellular process-broad",
  "Cytoplasmic metabolism",
  "RNA processing",
  "Translation / folding",
  "Surface family",
  "PIase / PPI"
)

module_nodes <- module_nodes %>%
  mutate(
    display_program = case_when(
      side == "host" & pretty %in% host_keep_pretty ~ pretty,
      side == "pathogen" & pretty %in% pathogen_keep_pretty ~ pretty,
      TRUE ~ "Other / unassigned"
    )
  )

alias_lookup <- bind_rows(
  module_nodes %>% transmute(alias = module_id, side, display_program),
  module_nodes %>% transmute(alias = module_label, side, display_program),
  module_nodes %>% transmute(alias = module_key, side, display_program)
) %>%
  filter(!is.na(alias), nzchar(alias)) %>%
  distinct(alias, side, .keep_all = TRUE)

host_y <- c(
  "Signaling-broad" = 0.80,
  "Microtubule / cytoskeleton" = 0.42,
  "RNA processing" = 0.04,
  "Ubiquitin / proteostasis" = -0.34,
  "Metabolism" = -0.68,
  "Other / unassigned" = -1.00
)
pathogen_y <- c(
  "Cellular process-broad" = 0.80,
  "Cytoplasmic metabolism" = 0.42,
  "RNA processing" = 0.04,
  "Translation / folding" = -0.34,
  "Surface family" = -0.68,
  "PIase / PPI" = -1.00
)

host_fills <- c(
  "Signaling-broad" = host_pale,
  "Microtubule / cytoskeleton" = "#E3ECF1",
  "RNA processing" = host_fill,
  "Ubiquitin / proteostasis" = "#D7E1E7",
  "Metabolism" = "#CDDBE3",
  "Other / unassigned" = "#ECEFF1"
)
pathogen_fills <- c(
  "Cellular process-broad" = pathogen_pale,
  "Cytoplasmic metabolism" = "#F6EEEC",
  "RNA processing" = "#F3E8E6",
  "Translation / folding" = pathogen_fill,
  "Surface family" = "#E8D6D3",
  "PIase / PPI" = "#E0C9C6"
)

place_side <- function(data, side_name, y_map) {
  data %>%
    filter(side == side_name) %>%
    mutate(
      y = unname(y_map[display_program]),
      xmin = if_else(side_name == "host", -0.92, 0.70),
      xmax = if_else(side_name == "host", -0.70, 0.92),
      label_x = (xmin + xmax) / 2,
      y_min = y - 0.04,
      y_max = y + 0.04
    )
}

host_modules <- module_nodes %>%
  group_by(side, display_program) %>%
  summarise(
    scaffold_n = sum(scaffold_n),
    core_n = sum(core_n),
    .groups = "drop"
  ) %>%
  place_side("host", host_y)
pathogen_modules <- module_nodes %>%
  group_by(side, display_program) %>%
  summarise(
    scaffold_n = sum(scaffold_n),
    core_n = sum(core_n),
    .groups = "drop"
  ) %>%
  place_side("pathogen", pathogen_y)

if (anyNA(host_modules$y)) {
  stop("Host display y is incomplete: ", paste(host_modules$display_program[is.na(host_modules$y)], collapse = ", "))
}
if (anyNA(pathogen_modules$y)) {
  stop("Pathogen display y is incomplete: ", paste(pathogen_modules$display_program[is.na(pathogen_modules$y)], collapse = ", "))
}

module_colours <- c(
  setNames(unname(host_fills[host_modules$display_program]), host_modules$display_program),
  setNames(unname(pathogen_fills[pathogen_modules$display_program]), paste0("pathogen::", pathogen_modules$display_program))
)

bubble_name <- function(program) {
  dplyr::case_when(
    program == "Other / unassigned" ~ "Other /\nunassigned",
    program == "PIase / PPI" ~ "PPIase",
    program == "Signaling-broad" ~ "Signaling-\nbroad",
    program == "RNA processing" ~ "RNA\nprocessing",
    program == "Ubiquitin / proteostasis" ~ "Ubiquitin /\nproteostasis",
    program == "Microtubule / cytoskeleton" ~ "Microtubule /\ncytoskeleton",
    program == "Cytoplasmic metabolism" ~ "Cytoplasmic\nmetabolism",
    program == "Cellular process-broad" ~ "Cellular process-\nbroad",
    program == "Translation / folding" ~ "Translation /\nfolding",
    program == "Surface family" ~ "Surface\nfamily",
    TRUE ~ program
  )
}

size_min <- 17
size_max <- 28
panel_width_mm <- 175
x_span <- 2.04
n_min <- min(c(host_modules$scaffold_n, pathogen_modules$scaffold_n))
n_max <- max(c(host_modules$scaffold_n, pathogen_modules$scaffold_n))
bubble_size <- function(n) {
  size_min + (n - n_min) / pmax(n_max - n_min, 1) * (size_max - size_min)
}
bubble_diameter <- function(n) bubble_size(n) * x_span / panel_width_mm

module_layout <- bind_rows(host_modules, pathogen_modules) %>%
  mutate(
    colour_key = if_else(side == "host", display_program, paste0("pathogen::", display_program)),
    block_fill = unname(module_colours[colour_key]),
    ggplot_size = bubble_size(scaffold_n),
    diameter_data = bubble_diameter(scaffold_n),
    display_label = paste0(bubble_name(display_program), "\nn = ", comma(scaffold_n))
  )
if (any(is.na(module_layout$block_fill) | !nzchar(module_layout$block_fill))) {
  stop("Every display program requires a pale fill.")
}

bridge_layout <- module_bridges %>%
  left_join(
    alias_lookup %>%
      filter(side == "host") %>%
      select(host_key = alias, host_program = display_program),
    by = "host_key"
  ) %>%
  left_join(
    alias_lookup %>%
      filter(side == "pathogen") %>%
      select(pathogen_key = alias, pathogen_program = display_program),
    by = "pathogen_key"
  ) %>%
  mutate(
    host_program = coalesce(host_program, "Other / unassigned"),
    pathogen_program = coalesce(pathogen_program, "Other / unassigned")
  ) %>%
  group_by(host_program, pathogen_program) %>%
  summarise(
    bridge_weight = max(bridge_weight),
    axis_n = sum(axis_n),
    .groups = "drop"
  ) %>%
  mutate(bridge_id = row_number()) %>%
  left_join(
    host_modules %>%
      transmute(
        host_program = display_program, host_y = y,
        host_ymin = y_min, host_ymax = y_max,
        host_diameter = bubble_diameter(scaffold_n)
      ),
    by = "host_program"
  ) %>%
  left_join(
    pathogen_modules %>%
      transmute(
        pathogen_program = display_program, pathogen_y = y,
        pathogen_ymin = y_min, pathogen_ymax = y_max,
        pathogen_diameter = bubble_diameter(scaffold_n)
      ),
    by = "pathogen_program"
  )
if (anyNA(bridge_layout$host_y) || anyNA(bridge_layout$pathogen_y)) {
  stop("At least one bridge endpoint did not match a display program.")
}

ribbon_gap <- 0.004
bridge_layout <- bridge_layout %>%
  mutate(
    ribbon_width = rescale(bridge_weight, to = c(0.010, 0.028)),
    host_ribbon_width = ribbon_width,
    pathogen_ribbon_width = ribbon_width
  ) %>%
  group_by(host_program) %>%
  mutate(
    host_cap = 0.42 * first(host_diameter),
    host_raw_total = sum(host_ribbon_width) + ribbon_gap * (n() - 1),
    host_scale = pmin(1, host_cap / pmax(host_raw_total, 1e-6)),
    host_ribbon_width = host_ribbon_width * host_scale
  ) %>%
  ungroup() %>%
  group_by(pathogen_program) %>%
  mutate(
    pathogen_cap = 0.42 * first(pathogen_diameter),
    pathogen_raw_total = sum(pathogen_ribbon_width) + ribbon_gap * (n() - 1),
    pathogen_scale = pmin(1, pathogen_cap / pmax(pathogen_raw_total, 1e-6)),
    pathogen_ribbon_width = pathogen_ribbon_width * pathogen_scale
  ) %>%
  ungroup() %>%
  group_by(host_program) %>%
  arrange(desc(pathogen_y), .by_group = TRUE) %>%
  mutate(
    host_band_total = sum(host_ribbon_width) + ribbon_gap * (n() - 1),
    host_used_before = lag(
      cumsum(host_ribbon_width) + row_number() * ribbon_gap,
      default = 0
    ),
    y_start = host_y + host_band_total / 2 - host_used_before - host_ribbon_width / 2
  ) %>%
  ungroup() %>%
  group_by(pathogen_program) %>%
  arrange(desc(host_y), .by_group = TRUE) %>%
  mutate(
    pathogen_band_total = sum(pathogen_ribbon_width) + ribbon_gap * (n() - 1),
    pathogen_used_before = lag(
      cumsum(pathogen_ribbon_width) + row_number() * ribbon_gap,
      default = 0
    ),
    y_end = pathogen_y + pathogen_band_total / 2 -
      pathogen_used_before - pathogen_ribbon_width / 2
  ) %>%
  ungroup() %>%
  mutate(
    x_start = -0.81,
    x_end = 0.81,
    host_colour_key = host_program,
    pathogen_colour_key = paste0("pathogen::", pathogen_program)
  )

bezier_ribbon <- function(
    id, x0, y0, x3, y3, start_width, end_width,
    bridge_weight, host_key, pathogen_key,
    start_colour, end_colour, n = 120L) {
  t <- seq(0, 1, length.out = n)
  x1 <- -0.31
  x2 <- 0.31
  x <- (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * x1 +
    3 * (1 - t) * t^2 * x2 + t^3 * x3
  smooth_y <- function(start, end) {
    (1 - t)^3 * start + 3 * (1 - t)^2 * t * start +
      3 * (1 - t) * t^2 * end + t^3 * end
  }
  upper <- smooth_y(y0 + start_width / 2, y3 + end_width / 2)
  lower <- smooth_y(y0 - start_width / 2, y3 - end_width / 2)
  segment_colours <- grDevices::colorRampPalette(
    c(start_colour, end_colour)
  )(n - 1L)

  bind_rows(lapply(seq_len(n - 1L), function(j) {
    tibble(
      bridge_id = id,
      segment_id = j,
      x = c(x[j], x[j + 1L], x[j + 1L], x[j]),
      y = c(upper[j], upper[j + 1L], lower[j + 1L], lower[j]),
      bridge_weight = bridge_weight,
      host_key = host_key,
      pathogen_key = pathogen_key,
      segment_fill = segment_colours[j]
    )
  }))
}

bridge_ribbons <- bind_rows(lapply(seq_len(nrow(bridge_layout)), function(i) {
  row <- bridge_layout[i, ]
  bezier_ribbon(
    row$bridge_id, row$x_start, row$y_start, row$x_end, row$y_end,
    row$host_ribbon_width, row$pathogen_ribbon_width, row$bridge_weight,
    row$host_program, row$pathogen_program,
    unname(module_colours[row$host_colour_key]),
    unname(module_colours[row$pathogen_colour_key])
  )
})) %>%
  arrange(desc(bridge_weight), bridge_id, segment_id)

weight_legend <- tibble(
  x = -5, xend = -4, y = -5, yend = -5,
  bridge_weight = c(1, 3, 5)
)

panel <- ggplot() +
  geom_polygon(
    data = bridge_ribbons,
    aes(x, y, group = interaction(bridge_id, segment_id), fill = segment_fill),
    colour = NA, alpha = 0.62
  ) +
  geom_segment(
    data = weight_legend,
    aes(x = x, xend = xend, y = y, yend = yend, linewidth = bridge_weight),
    colour = "#6F767D", lineend = "butt", show.legend = TRUE
  ) +
  geom_point(
    data = module_layout,
    aes(x = label_x, y = y, size = ggplot_size, fill = block_fill, colour = side),
    shape = 21, stroke = 0.72
  ) +
  geom_text(
    data = module_layout,
    aes(x = label_x, y = y, label = display_label),
    colour = ink, size = 1.85, lineheight = 0.84
  ) +
  annotate(
    "text", x = -0.81, y = 1.10,
    label = paste0("HOST PROGRAMS\n", comma(host_total), " proteins"),
    colour = host_dark, fontface = "bold", size = 2.65, lineheight = 0.88
  ) +
  annotate(
    "text", x = 0.81, y = 1.10,
    label = paste0("PATHOGEN PROGRAMS\n", comma(pathogen_total), " proteins"),
    colour = pathogen_dark, fontface = "bold", size = 2.65, lineheight = 0.88
  ) +
  scale_fill_identity() +
  scale_colour_manual(
    values = c(host = host_dark, pathogen = pathogen_dark),
    guide = "none"
  ) +
  scale_size_identity(guide = "none") +
  scale_linewidth_continuous(
    range = c(0.55, 3.0),
    breaks = c(1, 3, 5),
    labels = label_number(accuracy = 1),
    name = "Bridge weight",
    guide = guide_legend(
      title.position = "left", title.hjust = 0.5,
      nrow = 1, keywidth = grid::unit(11, "mm")
    )
  ) +
  coord_cartesian(xlim = c(-1.02, 1.02), ylim = c(-1.16, 1.22), clip = "off") +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 5.7, face = "bold", colour = ink),
    legend.text = element_text(size = 5.3, colour = ink),
    legend.key.height = grid::unit(3, "mm"),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_rect(fill = "white", colour = NA)
  )

save_panel(panel, "Figure_3_panel_a_module_bridge_map", width_mm = 175, height_mm = 148)

cat(
  "Rendered Figure_3_panel_a_module_bridge_map with ", n_modules_full,
  " modules and ", n_bridges_full, " bridges (display ",
  nrow(module_layout), " programs, ", nrow(bridge_layout), " ribbons).\n",
  sep = ""
)
