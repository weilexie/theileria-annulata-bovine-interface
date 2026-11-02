#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggforce)
  library(ggplot2)
  library(igraph)
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
  normalizePath(file.path("01_scripts", "03_panels", "05_figure_5", "01_panel_a_interface_atlas.R"), winslash = "/", mustWork = TRUE)
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
read_fig <- function(file) read_tsv(file.path(data_dir, "05_figure_5", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "05_figure_5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(plot, stem, width_mm = 183, height_mm = 112) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)
}

host <- host_fill
pathogen <- pathogen_fill

nodes_raw <- read_fig("01_panel_a_final_nodes.tsv")
for (col in c("display_label", "canonical_gene", "pathogen_name", "core_axis_endpoint", "hub", "module_label")) {
  if (!col %in% names(nodes_raw)) nodes_raw[[col]] <- NA
}
nodes_raw <- nodes_raw %>%
  mutate(
    core_axis_endpoint = flag_col(core_axis_endpoint),
    hub = flag_col(hub),
    module_label = coalesce(module_label, "Unassigned"),
    module_key = paste(layer, module_label, sep = "::"),
    label = endpoint_display(sanitize_label(coalesce(display_label, canonical_gene, pathogen_name, node_id)))
  )

edges_raw <- read_fig("02_panel_a_final_edges.tsv") %>%
  select(-any_of(c("x_from", "y_from", "x_to", "y_to")))
for (col in c(
  "known_cross_layer_anchor", "mechanism_class", "axis_family", "interaction_name",
  "axis_id", "biological_axis_id", "pathogen_gene", "host_gene", "final_tier",
  "species_adjusted_confidence", "edge_scope"
)) {
  if (!col %in% names(edges_raw)) edges_raw[[col]] <- NA
}

axis_file_candidates <- c(
  "08_panel_d_core_mechanisms.tsv",
  "panel_d_retained_axes.tsv"
)
axis_descriptor <- bind_rows(lapply(axis_file_candidates, function(f) {
  fp <- file.path(data_dir, "05_figure_5", f)
  if (!file.exists(fp)) return(NULL)
  read_tsv(fp, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
}))
axis_alt <- file.path(data_dir, "04_figure_4", "07_panel_d_axis_descriptor.tsv")
if (file.exists(axis_alt)) {
  axis_descriptor <- bind_rows(
    axis_descriptor,
    read_tsv(axis_alt, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
  )
}
if (is.null(axis_descriptor) || !nrow(axis_descriptor)) {
  axis_descriptor <- edges_raw %>%
    filter(flag_col(known_cross_layer_anchor) | edge_scope == "cross") %>%
    transmute(
      axis_id = coalesce(axis_id, biological_axis_id, interaction_name),
      interaction_name, pathogen_gene, host_gene,
      mechanism_class, axis_family, final_tier
    )
}
for (col in c("axis_id", "biological_axis_id", "interaction_name", "pathogen_gene", "host_gene", "mechanism_class", "axis_family", "final_tier")) {
  if (!col %in% names(axis_descriptor)) axis_descriptor[[col]] <- NA
}
axis_descriptor <- axis_descriptor %>%
  mutate(
    axis_id_key = axis_key(coalesce(biological_axis_id, axis_id), interaction_name, pathogen_gene, host_gene),
    axis_tier = axis_tier(axis_id_key, final_tier),
    mechanism = coalesce(
      unname(mechanism_cols[as.character(mechanism_class)]),
      mechanism_from_family(coalesce(axis_family, mechanism_class)),
      mechanism_from_family(axis_id_key)
    )
  ) %>%
  filter(axis_tier %in% c("Core", "Extended"), !axis_id_key %in% context_axis_ids) %>%
  distinct(axis_id_key, .keep_all = TRUE)
axis_descriptor$mechanism[!axis_descriptor$mechanism %in% names(mechanism_cols)] <- NA
axis_descriptor <- axis_descriptor %>%
  mutate(
    mechanism = coalesce(
      mechanism,
      recode(
        axis_id_key,
        `p104--MAPRE1` = "Microtubule docking",
        `TaSP--tubulin` = "Microtubule docking",
        `TaPin1--PKM` = "Metabolic / proteostasis",
        `TaPin1--FBXW7` = "Metabolic / proteostasis",
        `TaPHB1--RUVBL1` = "Metabolic / proteostasis",
        `TA20980--CD2AP` = "Trafficking / signaling",
        `TA11965--TFG` = "Trafficking / signaling",
        `TaCyp1--MED21` = "Trafficking / signaling",
        `TA05560--RBM39` = "Trafficking / signaling",
        `TaCP--PPP4C` = "Protease / cullin",
        `TaCP--CRBN` = "Protease / cullin",
        .default = NA_character_
      )
    )
  )

edges_raw <- edges_raw %>%
  mutate(
    known_cross_layer_anchor = flag_col(known_cross_layer_anchor),
    axis_id_key = axis_key(coalesce(biological_axis_id, axis_id), interaction_name, pathogen_gene, host_gene),
    axis_tier = axis_tier(axis_id_key, final_tier)
  ) %>%
  left_join(
    axis_descriptor %>% select(axis_id_key, axis_tier_lock = axis_tier, mechanism_lock = mechanism),
    by = "axis_id_key"
  ) %>%
  mutate(
    axis_tier = coalesce(axis_tier_lock, axis_tier),
    mechanism = coalesce(mechanism_lock, as.character(mechanism_class), mechanism_from_family(axis_family)),
    retained_axis = axis_tier %in% c("Core", "Extended") & !axis_id_key %in% context_axis_ids,
    core_axis = axis_tier == "Core",
    extended_axis = axis_tier == "Extended"
  )

# Never draw context-only p104–CLASP1 on the core atlas.
edges_raw <- edges_raw %>%
  mutate(
    retained_axis = retained_axis & !str_detect(coalesce(interaction_name, axis_id_key, ""), "CLASP1"),
    core_axis = core_axis & retained_axis,
    extended_axis = extended_axis & retained_axis
  )

side_layout <- function(side_name) {
  side_nodes <- nodes_raw %>% filter(layer == side_name) %>% arrange(module_label, node_id)
  intra_name <- paste0("intra_", side_name)
  side_edges <- edges_raw %>%
    filter(coalesce(edge_scope, "") == intra_name | (is.na(edge_scope) & FALSE))
  if (!nrow(side_edges)) {
    side_edges <- edges_raw %>%
      filter(source %in% side_nodes$node_id, target %in% side_nodes$node_id)
  }
  if (!nrow(side_edges)) {
    side_edges <- tibble(source = side_nodes$node_id[1], target = side_nodes$node_id[1], weight = 0.01)
  }
  graph <- graph_from_data_frame(
    side_edges %>% transmute(source, target, weight = pmax(num_col(weight, 0.01), 0.01)),
    directed = FALSE,
    vertices = side_nodes %>% transmute(name = node_id)
  )
  modules <- sort(unique(side_nodes$module_label))
  module_angle <- setNames(seq(0, 2 * pi, length.out = length(modules) + 1)[-1], modules)
  set.seed(if (side_name == "host") 20260717 else 20260718)
  initial <- cbind(
    cos(module_angle[side_nodes$module_label]) + rnorm(nrow(side_nodes), 0, 0.10),
    sin(module_angle[side_nodes$module_label]) + rnorm(nrow(side_nodes), 0, 0.10)
  )
  xy <- layout_with_fr(
    graph, coords = initial, weights = pmax(num_col(E(graph)$weight), 0.01),
    niter = 2200, grid = "nogrid"
  )
  xy[, 1] <- xy[, 1] - median(xy[, 1])
  xy[, 2] <- xy[, 2] - median(xy[, 2])
  sx <- max(quantile(abs(xy[, 1]), 0.95), 1e-6)
  sy <- max(quantile(abs(xy[, 2]), 0.95), 1e-6)
  ux <- xy[, 1] / sx
  uy <- xy[, 2] / sy
  rr <- sqrt(ux^2 + uy^2)
  theta <- atan2(uy, ux)
  radial <- tanh(rr / 1.35)
  radial <- radial / max(quantile(radial, 0.99), 1e-6)
  radial <- pmin(radial, 1.04)
  rx <- if (side_name == "host") 0.50 else 0.35
  ry <- if (side_name == "host") 0.76 else 0.60
  centre_x <- if (side_name == "host") -0.76 else 0.735
  x_ellipse <- centre_x + rx * radial * cos(theta)
  y_ellipse <- ry * radial * sin(theta)
  fan <- 0.045 * (abs(y_ellipse) / ry)^1.55
  x_ellipse <- x_ellipse + if (side_name == "host") -fan else fan
  tibble(
    node_id = V(graph)$name,
    x_side = x_ellipse,
    y_side = y_ellipse
  )
}
side_xy <- bind_rows(side_layout("host"), side_layout("pathogen"))

retained_keys <- c(core_axis_ids, extended_axis_ids)
host_ids <- nodes_raw$node_id[nodes_raw$layer == "host"]
axis_edges <- edges_raw %>%
  filter(retained_axis, axis_id_key %in% retained_keys) %>%
  left_join(
    axis_descriptor %>%
      select(axis_id_key, pathogen_gene_lock = pathogen_gene, host_gene_lock = host_gene),
    by = "axis_id_key"
  ) %>%
  mutate(
    mechanism = factor(mechanism, levels = names(mechanism_cols)),
    host_end = endpoint_display(sanitize_label(coalesce(
      host_gene_lock, host_gene, str_split_fixed(axis_id_key, "--", 2)[, 2]
    ))),
    pathogen_end = endpoint_display(sanitize_label(coalesce(
      pathogen_gene_lock, pathogen_gene, str_split_fixed(axis_id_key, "--", 2)[, 1]
    ))),
    host_node = if_else(source %in% host_ids, source, target),
    pathogen_node = if_else(source %in% host_ids, target, source)
  ) %>%
  group_by(axis_id_key) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(mechanism, desc(core_axis), desc(num_col(species_adjusted_confidence)), interaction_name)
if (!nrow(axis_edges)) stop("Figure 5a requires retained interface axes.")
missing_axes <- setdiff(retained_keys, axis_edges$axis_id_key)
if (length(missing_axes)) {
  stop("Figure 5a missing retained axes: ", paste(missing_axes, collapse = ", "))
}
axis_edges <- axis_edges %>% mutate(axis_y = seq(0.68, -0.68, length.out = n()))

endpoint_meta <- bind_rows(
  axis_edges %>% transmute(
    node_id = host_node, axis_y, core_axis,
    endpoint_label = host_end, mechanism = as.character(mechanism)
  ),
  axis_edges %>% transmute(
    node_id = pathogen_node, axis_y, core_axis,
    endpoint_label = pathogen_end, mechanism = as.character(mechanism)
  )
) %>%
  group_by(node_id) %>%
  summarise(
    mean_y = mean(axis_y),
    core_axis = any(core_axis),
    endpoint_label = first(endpoint_label),
    mechanism = first(mechanism),
    .groups = "drop"
  ) %>%
  left_join(nodes_raw %>% select(node_id, layer), by = "node_id") %>%
  group_by(layer) %>%
  arrange(desc(mean_y), .by_group = TRUE) %>%
  mutate(interface_y = seq(0.68, -0.68, length.out = n())) %>%
  ungroup()
if (any(is.na(endpoint_meta$layer) | is.na(endpoint_meta$interface_y))) {
  stop("Figure 5a retained endpoints must exist in panel_a_final_nodes.tsv.")
}

# Gold grammar: each protein is drawn once. Interface endpoints move into
# the middle column; intra-edges become the thin spokes. Named effectors
# stay on the interface, not as duplicate cloud hubs. Ribosomal proteins
# and raw TA loci are not the labeled cloud hubs.
is_ribo_lab <- function(z) str_detect(coalesce(as.character(z), ""), "^(RP[SL]|RPLP|EEF)")
pick_cloud_hubs <- function(layer_name, preferred = character()) {
  d <- nodes_raw %>%
    mutate(lab = endpoint_display(sanitize_label(coalesce(display_label, canonical_gene, pathogen_name, node_id)))) %>%
    filter(
      layer == layer_name,
      !node_id %in% endpoint_meta$node_id,
      !is_ribo_lab(lab),
      !str_detect(lab, "TA_15795|^TA[0-9]")
    )
  pref_ids <- d$node_id[match(preferred, d$lab)]
  pref_ids <- pref_ids[!is.na(pref_ids)]
  extra <- d %>%
    filter(!node_id %in% pref_ids) %>%
    arrange(desc(num_col(degree)), desc(num_col(strength))) %>%
    pull(node_id)
  unique(c(pref_ids, extra))[seq_len(min(3L, length(unique(c(pref_ids, extra)))))]
}
host_hub_ids <- pick_cloud_hubs("host", c("CDK1", "STAT3", "HSPA5"))
path_hub_ids <- pick_cloud_hubs("pathogen")

nodes <- nodes_raw %>%
  left_join(side_xy, by = "node_id") %>%
  left_join(
    endpoint_meta %>% select(node_id, interface_y, core_axis, endpoint_label, mechanism),
    by = "node_id"
  ) %>%
  mutate(
    side_sign = if_else(layer == "host", -1, 1),
    retained_endpoint = node_id %in% endpoint_meta$node_id,
    core_end = retained_endpoint & coalesce(core_axis, FALSE),
    global_hub = node_id %in% c(host_hub_ids, path_hub_ids) & !retained_endpoint,
    x_plot = if_else(retained_endpoint, 0.045 + side_sign * 0.17, x_side),
    y_plot = if_else(retained_endpoint, interface_y, y_side),
    node_size = rescale(log1p(pmax(num_col(degree), 0)), to = c(0.55, 2.15)),
    plot_label = coalesce(endpoint_label, endpoint_display(label))
  )

edges <- edges_raw %>%
  left_join(nodes %>% select(source = node_id, x_from = x_plot, y_from = y_plot), by = "source") %>%
  left_join(nodes %>% select(target = node_id, x_to = x_plot, y_to = y_plot), by = "target")

interface_edges <- axis_edges %>%
  left_join(nodes %>% select(host_node = node_id, x_from = x_plot, y_from = y_plot), by = "host_node") %>%
  left_join(nodes %>% select(pathogen_node = node_id, x_to = x_plot, y_to = y_plot), by = "pathogen_node")
core_edges <- filter(interface_edges, core_axis)
extended_edges <- filter(interface_edges, extended_axis)

hull_modules <- nodes %>%
  filter(!retained_endpoint) %>%
  add_count(module_key) %>%
  filter(n >= 3)

interface_nodes <- nodes %>% filter(retained_endpoint)
endpoint_labels <- interface_nodes %>%
  mutate(
    label_x = 0.045 + side_sign * 0.205,
    label_hjust = if_else(layer == "host", 1, 0)
  )
hub_labels <- nodes %>%
  filter(global_hub) %>%
  mutate(
    hub_label_x = x_plot + side_sign * if_else(plot_label == "HSPA5", 0.105, 0.088),
    hub_label_hjust = if_else(layer == "host", 1, 0)
  )

plot <- ggplot() +
  ggforce::geom_mark_ellipse(
    data = hull_modules,
    aes(x_plot, y_plot, group = module_key, fill = layer),
    expand = grid::unit(1.05, "mm"), radius = grid::unit(1.6, "mm"),
    alpha = 0.055, colour = NA, show.legend = FALSE
  ) +
  geom_curve(
    data = filter(edges, coalesce(edge_scope, "") == "intra_host" |
      (source %in% nodes$node_id[nodes$layer == "host"] &
         target %in% nodes$node_id[nodes$layer == "host"] &
         !coalesce(retained_axis, FALSE))),
    aes(x_from, y_from, xend = x_to, yend = y_to),
    curvature = 0.035, colour = host_dark, alpha = 0.13, linewidth = 0.18, lineend = "round"
  ) +
  geom_curve(
    data = filter(edges, coalesce(edge_scope, "") == "intra_pathogen" |
      (source %in% nodes$node_id[nodes$layer == "pathogen"] &
         target %in% nodes$node_id[nodes$layer == "pathogen"] &
         !coalesce(retained_axis, FALSE))),
    aes(x_from, y_from, xend = x_to, yend = y_to),
    curvature = -0.035, colour = pathogen_dark, alpha = 0.14, linewidth = 0.18, lineend = "round"
  ) +
  geom_curve(
    data = extended_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, colour = mechanism),
    curvature = 0.045, alpha = 0.58, linewidth = 0.62, lineend = "round"
  ) +
  geom_curve(
    data = core_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to),
    curvature = 0.045, colour = "white", linewidth = 2.45, lineend = "round"
  ) +
  geom_curve(
    data = core_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, colour = mechanism),
    curvature = 0.045, alpha = 0.78, linewidth = 1.35, lineend = "round"
  ) +
  geom_point(
    data = filter(nodes, !retained_endpoint),
    aes(x_plot, y_plot, fill = layer, size = node_size),
    shape = 21, colour = alpha(ink, 0.72), stroke = 0.17, alpha = 0.86
  ) +
  geom_point(
    data = filter(nodes, global_hub),
    aes(x_plot, y_plot), shape = 21, fill = NA, colour = "white", size = 4.5, stroke = 1.8
  ) +
  geom_point(
    data = filter(nodes, global_hub),
    aes(x_plot, y_plot, colour = layer), shape = 21, fill = NA, size = 4.1, stroke = 0.75
  ) +
  geom_point(
    data = filter(interface_nodes, !core_end),
    aes(x_plot, y_plot, colour = mechanism), shape = 21, fill = "white", size = 3.4, stroke = 0.85, alpha = 0.30
  ) +
  geom_point(
    data = filter(interface_nodes, !core_end),
    aes(x_plot, y_plot, fill = layer, colour = mechanism), shape = 21, size = 2.15, stroke = 0.62
  ) +
  geom_point(
    data = filter(interface_nodes, core_end),
    aes(x_plot, y_plot, colour = mechanism), shape = 21, fill = "white", size = 4.3, stroke = 1.05, alpha = 0.30
  ) +
  geom_point(
    data = filter(interface_nodes, core_end),
    aes(x_plot, y_plot, fill = layer, colour = mechanism), shape = 21, size = 2.7, stroke = 0.75
  ) +
  geom_text(
    data = filter(endpoint_labels, !core_end),
    aes(label_x, y_plot, label = plot_label, hjust = label_hjust, colour = layer),
    size = 1.62, fontface = "bold", family = "sans"
  ) +
  geom_text(
    data = filter(endpoint_labels, core_end),
    aes(label_x, y_plot, label = plot_label, hjust = label_hjust, colour = layer),
    size = 1.82, fontface = "bold", family = "sans"
  ) +
  geom_text(
    data = hub_labels,
    aes(hub_label_x, y_plot, label = plot_label, hjust = hub_label_hjust),
    colour = "white", size = 2.08, fontface = "bold", family = "sans"
  ) +
  geom_text(
    data = hub_labels,
    aes(hub_label_x, y_plot, label = plot_label, hjust = hub_label_hjust, colour = layer),
    size = 1.72, fontface = "bold", family = "sans", show.legend = FALSE
  ) +
  annotate("text", x = -0.73, y = 0.91, label = "HOST NETWORK", colour = host_dark, fontface = "bold", size = 2.55) +
  annotate("text", x = 0.045, y = 0.91, label = "CORE INTERFACE", colour = "#927D96", fontface = "bold", size = 2.55) +
  annotate("text", x = 0.73, y = 0.91, label = "PATHOGEN NETWORK", colour = pathogen_dark, fontface = "bold", size = 2.55) +
  scale_fill_manual(values = c(host = host, pathogen = pathogen), guide = "none") +
  scale_colour_manual(values = c(host = host_dark, pathogen = pathogen_dark, mechanism_cols), guide = "none") +
  scale_size_identity() +
  coord_fixed(ratio = 0.77, xlim = c(-1.31, 1.27), ylim = c(-0.86, 0.94), clip = "off") +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "none",
    plot.margin = margin(4, 8, 2, 8)
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", unset = ""), "1")) {
  save_panel(plot, "Figure_5_panel_a_interface_atlas")
}
cat(
  "Rendered Figure_5_panel_a_interface_atlas: ",
  nrow(nodes), " nodes; ", nrow(core_edges), " core axes thick; ",
  nrow(extended_edges), " extended axes thin; axes ",
  paste(axis_edges$axis_id_key, collapse = ", "), "; interface labels ",
  paste(sort(unique(endpoint_labels$plot_label)), collapse = ", "),
  "; hubs ", paste(sort(unique(hub_labels$plot_label)), collapse = ", "), ".\n",
  sep = ""
)
