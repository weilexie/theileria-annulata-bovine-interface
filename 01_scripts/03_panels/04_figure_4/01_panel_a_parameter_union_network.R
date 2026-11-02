#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(grid)
  library(igraph)
  library(graphlayouts)
  library(ggrepel)
})

ofiles <- vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else as.character(x$ofile[[1]]), character(1))
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- if (any(nzchar(ofiles))) normalizePath(tail(ofiles[nzchar(ofiles)], 1), winslash = "/") else if (length(arg)) normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/") else normalizePath("01_scripts/03_panels/04_figure_4/01_panel_a_parameter_union_network.R", winslash = "/")
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

data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
data_candidates <- c(data_override, p("04_plot_data"))
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
data_dir <- data_candidates[[1]]
read_fig <- function(x) read_tsv(
  file.path(data_dir, "04_figure_4", x), show_col_types = FALSE, progress = FALSE,
  na = c("", "NA", "NaN")
)

od <- Sys.getenv("FIGURE_OUTPUT_DIR", "")
out <- if (nzchar(od)) od else p("06_panel_output", "04_figure_4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

status <- read_fig("01_panel_a_pcsf_parameters.tsv") %>%
  mutate(success = flag_col(success))
nodes_raw <- ensure_cols(
  read_fig("02_panel_a_pcsf_union_nodes.tsv"),
  c("display_label", "pathogen_name", "canonical_gene")
) %>%
  mutate(
    display_label = endpoint_display(sanitize_label(coalesce(display_label, pathogen_name, canonical_gene, node_id)))
  )
edges_raw <- read_fig("03_panel_a_pcsf_union_edges.tsv")

n_persistent_nodes <- sum(nodes_raw$selection_frequency >= 0.5)
n_persistent_edges <- sum(edges_raw$selection_frequency >= 0.5)
stopifnot(
  nrow(status) == 18L,
  all(status$success),
  n_persistent_nodes == 104L,
  n_persistent_edges == 95L
)

nodes <- nodes_raw %>%
  mutate(
    layer = factor(layer, levels = c("host", "pathogen")),
    node_role = case_when(
      selected_runs == 18 ~ "Pan-run core",
      selection_frequency >= 0.5 ~ "Persistent backbone",
      TRUE ~ "Conditional envelope"
    )
  )

edges <- edges_raw %>%
  left_join(nodes %>% select(source = node_id, source_layer = layer), by = "source") %>%
  left_join(nodes %>% select(target = node_id, target_layer = layer), by = "target") %>%
  mutate(
    edge_layer = case_when(
      source_layer == "host" & target_layer == "host" ~ "host",
      source_layer == "pathogen" & target_layer == "pathogen" ~ "pathogen",
      TRUE ~ "cross"
    ),
    persistent_edge = selection_frequency >= 0.5,
    pan_run_edge = selected_runs == 18
  )

graph_union <- graph_from_data_frame(
  edges %>% select(source, target, selection_frequency),
  directed = FALSE,
  vertices = nodes %>% transmute(name = node_id)
)
if (components(graph_union)$no > 1L) {
  warning("PCSF union has ", components(graph_union)$no, " components; layout uses the full union.")
}

# One continuous PCSF union, following Figure 2a. Repeatedly retained edges are
# treated as shorter distances, which pulls the stable backbone into focus
# without manually separating species or arranging bridges as a diagram.
set.seed(20260714)
# Theileria has 158 pan-run nodes (gold 42) and a 60-node / 198-edge pathogen
# module. A length floor plus a shallow frequency exponent keep high-frequency
# edges short without letting same-frequency cliques collapse to a point.
freq_w <- pmax(E(graph_union)$selection_frequency, 1 / nrow(status))
edge_length <- 0.42 + 1 / freq_w^0.18
# Mild length jitter on envelope edges so hub-and-spoke cores do not form
# regular flowers; persistent distances stay intact.
edge_length <- edge_length * ifelse(freq_w >= 0.5, 1, runif(length(edge_length), 0.82, 1.32))
xy <- graphlayouts::layout_with_stress(
  graph_union,
  weights = edge_length,
  iter = 1600,
  bbox = 95
)
xy <- sweep(xy, 2, apply(xy, 2, median), FUN = "-")

layout_nodes <- tibble(node_id = V(graph_union)$name, x0 = xy[, 1], y0 = xy[, 2]) %>%
  left_join(nodes, by = "node_id")

# Rotate species centroids onto a horizontal interface while preserving the
# continuous topology inferred above.
centroids <- layout_nodes %>%
  group_by(layer) %>%
  summarise(x = mean(x0), y = mean(y0), .groups = "drop")
host_centroid <- centroids %>% filter(layer == "host")
pathogen_centroid <- centroids %>% filter(layer == "pathogen")
theta <- atan2(pathogen_centroid$y - host_centroid$y, pathogen_centroid$x - host_centroid$x)

layout_nodes <- layout_nodes %>%
  mutate(
    x_rot = x0 * cos(theta) + y0 * sin(theta),
    y_rot = -x0 * sin(theta) + y0 * cos(theta)
  )
if (mean(layout_nodes$x_rot[layout_nodes$layer == "pathogen"]) < mean(layout_nodes$x_rot[layout_nodes$layer == "host"])) {
  layout_nodes$x_rot <- -layout_nodes$x_rot
}

# Break coincident cores so later radial opening has a direction.
set.seed(20260714)
layout_nodes$x_rot <- layout_nodes$x_rot + rnorm(nrow(layout_nodes), sd = 0.012)
layout_nodes$y_rot <- layout_nodes$y_rot + rnorm(nrow(layout_nodes), sd = 0.012)

# Open packed cores around each layer centroid without making two islands.
expand_layer <- function(x, y, power = 0.72, min_frac = 0.10) {
  cx <- stats::median(x, na.rm = TRUE)
  cy <- stats::median(y, na.rm = TRUE)
  dx <- x - cx
  dy <- y - cy
  r <- sqrt(dx * dx + dy * dy)
  r0 <- max(r, na.rm = TRUE)
  tiny <- which(r < 1e-9)
  if (length(tiny)) {
    ang <- 2 * pi * seq_along(tiny) / length(tiny)
    dx[tiny] <- cos(ang)
    dy[tiny] <- sin(ang)
    r[tiny] <- 1e-9
  }
  r_new <- pmax(r0 * (r / r0)^power, min_frac * r0)
  list(x = cx + dx / r * r_new, y = cy + dy / r * r_new)
}
for (ly in c("host", "pathogen")) {
  idx <- which(layout_nodes$layer == ly)
  opened <- expand_layer(
    layout_nodes$x_rot[idx],
    layout_nodes$y_rot[idx],
    power = if (ly == "pathogen") 0.48 else 0.46,
    min_frac = if (ly == "pathogen") 0.05 else 0.07
  )
  layout_nodes$x_rot[idx] <- opened$x
  layout_nodes$y_rot[idx] <- opened$y
}

# Pathogen is a dense module inside a larger host cloud; scale it up around its
# own centroid so the lace is readable after the shared canvas fit.
for (ly in c("host", "pathogen")) {
  idx <- which(layout_nodes$layer == ly)
  cx <- stats::median(layout_nodes$x_rot[idx])
  cy <- stats::median(layout_nodes$y_rot[idx])
  s <- if (ly == "pathogen") 1.85 else 1.16
  layout_nodes$x_rot[idx] <- cx + (layout_nodes$x_rot[idx] - cx) * s
  layout_nodes$y_rot[idx] <- cy + (layout_nodes$y_rot[idx] - cy) * s
}

# Mild horizontal scale-apart of species centroids; the single p104–MAPRE1
# bridge stays a continuous link rather than a bipartite gap.
host_idx <- layout_nodes$layer == "host"
path_idx <- layout_nodes$layer == "pathogen"
mid_x <- mean(c(
  mean(layout_nodes$x_rot[host_idx]),
  mean(layout_nodes$x_rot[path_idx])
))
apart <- 1.06
layout_nodes$x_rot[host_idx] <- mid_x + (layout_nodes$x_rot[host_idx] - mid_x) * apart
layout_nodes$x_rot[path_idx] <- mid_x + (layout_nodes$x_rot[path_idx] - mid_x) * apart

# Soft species field along x so bridge ends meet at the interface instead of
# swapping deep into the opposite cloud.
for (ly in c("host", "pathogen")) {
  idx <- which(layout_nodes$layer == ly)
  cx <- stats::median(layout_nodes$x_rot[idx])
  layout_nodes$x_rot[idx] <- 0.88 * layout_nodes$x_rot[idx] + 0.12 * cx
}

# Gold 4a uses a slightly wider landscape (1.04 x, 0.93 y) so the lace of
# edges remains visible. Do not fill the canvas into a square.
xy_range <- apply(cbind(layout_nodes$x_rot, layout_nodes$y_rot), 2, range, na.rm = TRUE)
xy_scale <- min(2.12 / diff(xy_range[, 1]), 1.68 / diff(xy_range[, 2]))

plot_nodes <- layout_nodes %>%
  mutate(
    x_plot = x_rot * xy_scale * 1.10,
    y_plot = y_rot * xy_scale * 0.90,
    x_plot = x_plot - mean(range(x_plot)) - 0.06,
    y_plot = y_plot - mean(range(y_plot)) - 0.04,
    union_degree = as.numeric(degree(graph_union)[node_id])
  )

# Open locally overcrowded neighbourhoods (hub-and-spoke flowers, dense
# pathogen core) before the pairwise collision pass.
open_local <- function(xy, radius = 0.11) {
  n <- nrow(xy)
  if (n < 5L) return(xy)
  d <- as.matrix(stats::dist(xy))
  out <- xy
  for (i in seq_len(n)) {
    nb <- which(d[i, ] <= radius)
    if (length(nb) < 5L) next
    cx <- mean(xy[nb, 1])
    cy <- mean(xy[nb, 2])
    s <- min(1.58, 1 + 0.042 * (length(nb) - 4))
    out[i, 1] <- cx + (xy[i, 1] - cx) * s
    out[i, 2] <- cy + (xy[i, 2] - cy) * s
  }
  out
}
xy_push <- cbind(plot_nodes$x_plot, plot_nodes$y_plot)
deg_all <- plot_nodes$union_degree
for (ly in c("host", "pathogen")) {
  idx <- which(plot_nodes$layer == ly)
  xy_push[idx, ] <- open_local(
    xy_push[idx, , drop = FALSE],
    radius = if (ly == "pathogen") 0.13 else 0.10
  )
}
# Lengthen short spokes around high-degree hubs so stars read as lace, not flowers.
spread_spokes <- function(xy, degree, radius = 0.16, factor = 1.42, deg_cut = 6) {
  hubs <- which(degree >= deg_cut)
  if (!length(hubs) || nrow(xy) < 3L) return(xy)
  d <- as.matrix(stats::dist(xy))
  out <- xy
  for (h in hubs) {
    nb <- which(d[h, ] > 1e-9 & d[h, ] < radius)
    if (!length(nb)) next
    dx <- xy[nb, 1] - xy[h, 1]
    dy <- xy[nb, 2] - xy[h, 2]
    out[nb, 1] <- xy[h, 1] + dx * factor
    out[nb, 2] <- xy[h, 2] + dy * factor
  }
  out
}
for (ly in c("host", "pathogen")) {
  idx <- which(plot_nodes$layer == ly)
  xy_push[idx, ] <- spread_spokes(
    xy_push[idx, , drop = FALSE],
    deg_all[idx],
    radius = if (ly == "pathogen") 0.18 else 0.15,
    factor = if (ly == "pathogen") 1.38 else 1.58,
    deg_cut = if (ly == "pathogen") 8 else 5
  )
}

# Within-layer collision first, then a weaker global pass, so overlapping
# dots do not hide edges and the two clouds stay one continuous union.
separate_xy <- function(xy, min_sep, n_iter = 90, strength = 0.55) {
  n <- nrow(xy)
  if (n < 2L) return(xy)
  for (iter in seq_len(n_iter)) {
    d1 <- outer(xy[, 1], xy[, 1], "-")
    d2 <- outer(xy[, 2], xy[, 2], "-")
    dist <- sqrt(d1 * d1 + d2 * d2)
    diag(dist) <- Inf
    hit <- dist < min_sep
    if (!any(hit)) break
    dist[dist < 1e-9] <- 1e-9
    force <- (min_sep - dist) / 2
    force[!hit] <- 0
    fx <- strength * rowSums(force * d1 / dist)
    fy <- strength * rowSums(force * d2 / dist)
    step <- sqrt(fx * fx + fy * fy)
    cap <- min_sep * 0.85
    damp <- pmin(1, cap / pmax(step, 1e-9))
    xy[, 1] <- xy[, 1] + fx * damp
    xy[, 2] <- xy[, 2] + fy * damp
  }
  xy
}
for (ly in c("host", "pathogen")) {
  idx <- which(plot_nodes$layer == ly)
  xy_push[idx, ] <- separate_xy(
    xy_push[idx, , drop = FALSE],
    min_sep = if (ly == "pathogen") 0.082 else 0.088,
    n_iter = 260,
    strength = 0.62
  )
}
xy_push <- separate_xy(xy_push, min_sep = 0.034, n_iter = 50, strength = 0.35)
# Pathogen cloud is a dense module; a mild y-flatten keeps it tighter than host
# without collapsing back to a ball. Separate again afterwards.
path_idx <- which(plot_nodes$layer == "pathogen")
cy_p <- stats::median(xy_push[path_idx, 2])
xy_push[path_idx, 2] <- cy_p + (xy_push[path_idx, 2] - cy_p) * 0.90
hard_sep <- function(xy, min_sep, n_pass = 10) {
  if (nrow(xy) < 2L) return(xy)
  for (pass in seq_len(n_pass)) {
    d1 <- outer(xy[, 1], xy[, 1], "-")
    d2 <- outer(xy[, 2], xy[, 2], "-")
    dist <- sqrt(d1 * d1 + d2 * d2)
    diag(dist) <- Inf
    hit <- dist < min_sep
    if (!any(hit)) break
    dist[dist < 1e-9] <- 1e-9
    push <- (min_sep - dist) / 2
    push[!hit] <- 0
    xy[, 1] <- xy[, 1] + rowSums(push * d1 / dist)
    xy[, 2] <- xy[, 2] + rowSums(push * d2 / dist)
  }
  xy
}
for (ly in c("host", "pathogen")) {
  idx <- which(plot_nodes$layer == ly)
  xy_push[idx, ] <- hard_sep(
    xy_push[idx, , drop = FALSE],
    min_sep = if (ly == "pathogen") 0.082 else 0.088
  )
}
plot_nodes$x_plot <- xy_push[, 1]
plot_nodes$y_plot <- xy_push[, 2]
plot_nodes$x_plot <- plot_nodes$x_plot - mean(range(plot_nodes$x_plot)) - 0.04
plot_nodes$y_plot <- plot_nodes$y_plot - mean(range(plot_nodes$y_plot)) - 0.04

# HOST / PATHOGEN sit over each cloud's centroid, not the bounding-box
# centre or the pathogen right half (those read as shifted right).
title_pad <- 0.12
title_y <- max(plot_nodes$y_plot, na.rm = TRUE) + title_pad
host_xs <- plot_nodes$x_plot[plot_nodes$layer == "host"]
pathogen_xs <- plot_nodes$x_plot[plot_nodes$layer == "pathogen"]
# Host mass sits in the left-centre stars, not the interface tail, so use a
# lower quantile. Pathogen cloud is compact; the mean is already over it.
host_title_x <- as.numeric(stats::quantile(host_xs, 0.32, na.rm = TRUE))
pathogen_title_x <- mean(pathogen_xs, na.rm = TRUE)
ylim_top <- max(0.82, title_y + 0.12)
ylim_bottom <- -1.02

bridge_specs <- edges %>%
  filter(edge_layer == "cross", persistent_edge) %>%
  mutate(
    host_node = if_else(source_layer == "host", source, target),
    pathogen_node = if_else(source_layer == "pathogen", source, target),
    host_label = nodes$display_label[match(host_node, nodes$node_id)],
    pathogen_label = nodes$display_label[match(pathogen_node, nodes$node_id)],
    bridge_label = "p104–MAPRE1"
  )
if (nrow(bridge_specs) != 1L) {
  warning(
    "Expected a single persistent cross-layer edge (p104–MAPRE1); found ",
    nrow(bridge_specs), "."
  )
  if (!nrow(bridge_specs)) {
    bridge_specs <- edges %>%
      filter(edge_layer == "cross") %>%
      arrange(desc(selection_frequency)) %>%
      slice_head(n = 1) %>%
      mutate(
        host_node = if_else(source_layer == "host", source, target),
        pathogen_node = if_else(source_layer == "pathogen", source, target),
        bridge_label = "p104–MAPRE1"
      )
  }
}
bridge_ids <- unique(c(bridge_specs$host_node, bridge_specs$pathogen_node))

edge_plot <- edges %>%
  left_join(plot_nodes %>% select(source = node_id, x_from = x_plot, y_from = y_plot), by = "source") %>%
  left_join(plot_nodes %>% select(target = node_id, x_to = x_plot, y_to = y_plot), by = "target")

bridge_edges <- edge_plot %>%
  filter(edge_layer == "cross") %>%
  left_join(bridge_specs %>% select(source, target, bridge_label), by = c("source", "target")) %>%
  mutate(
    mid_x = (x_from + x_to) / 2,
    mid_y = (y_from + y_to) / 2 + 0.07
  )

bridge_nodes <- plot_nodes %>% filter(node_id %in% bridge_ids)
conditional_nodes <- plot_nodes %>% filter(node_role == "Conditional envelope")
persistent_nodes <- plot_nodes %>% filter(node_role == "Persistent backbone")
core_nodes <- plot_nodes %>% filter(node_role == "Pan-run core")
context_edges <- edge_plot %>% filter(edge_layer != "cross")
persistent_edges <- edge_plot %>% filter(edge_layer != "cross", persistent_edge)
core_edges <- edge_plot %>% filter(edge_layer != "cross", pan_run_edge)

# Gold 4a labels a few hubs plus interface ends. Raw TA loci stay unlabeled
# unless they are retained-axis parasite ends (none of TA20980/TA05560/TA11965
# are in this union).
is_raw_ta <- function(x) {
  grepl("^TA[_-]?[0-9]+\\??$", coalesce(as.character(x), "")) &
    !(coalesce(as.character(x), "") %in% c("TA20980", "TA05560", "TA11965"))
}
hub_nodes <- plot_nodes %>%
  filter(
    node_role != "Conditional envelope",
    !node_id %in% bridge_ids,
    !is.na(display_label),
    nzchar(display_label),
    !is_raw_ta(display_label)
  ) %>%
  group_by(layer) %>%
  slice_max(union_degree * selection_frequency, n = 4, with_ties = FALSE) %>%
  ungroup()
force_node_labels <- unique(c(
  core_endpoint_labels,
  "EEF2", "PSMD2", "HNRNPA1",
  hub_nodes$display_label
))
label_nodes <- plot_nodes %>%
  filter(
    display_label %in% force_node_labels,
    !node_id %in% bridge_ids,
    !is.na(display_label),
    nzchar(display_label),
    !is_raw_ta(display_label)
  )
hub_nodes <- hub_nodes %>% filter(display_label %in% force_node_labels)

network_labels <- bind_rows(
  bridge_edges %>%
    filter(!is.na(bridge_label), nzchar(bridge_label)) %>%
    transmute(label_x = mid_x, label_y = mid_y, label = bridge_label, label_group = "interface"),
  label_nodes %>%
    transmute(label_x = x_plot, label_y = y_plot, label = display_label, label_group = as.character(layer))
)
if (any(grepl("TA_15795\\?|[\u4e00-\u9fff]", network_labels$label))) {
  stop("Labels still contain TA_15795? or Chinese.")
}
if (!any(network_labels$label == "p104–MAPRE1")) {
  stop("p104–MAPRE1 must be force-labeled as the persistent cross-layer axis.")
}

marker_legend <- tibble(
  marker = factor(
    c("Pan-run core", "Hub", "Interface endpoint"),
    levels = c("Pan-run core", "Hub", "Interface endpoint")
  )
)

panel <- ggplot() +
  geom_segment(
    data = context_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, colour = edge_layer),
    alpha = 0.13, linewidth = 0.080, lineend = "round"
  ) +
  geom_segment(
    data = persistent_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, colour = edge_layer),
    alpha = 0.56, linewidth = 0.24, lineend = "round"
  ) +
  geom_segment(
    data = core_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, colour = edge_layer),
    alpha = 0.72, linewidth = 0.34, lineend = "round"
  ) +
  geom_segment(
    data = bridge_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to),
    colour = interface_mid, alpha = 0.70, linewidth = 0.70, lineend = "round"
  ) +
  geom_point(
    data = hub_nodes,
    aes(x_plot, y_plot, fill = layer),
    shape = 21, colour = "transparent", size = 3.2, alpha = 0.14
  ) +
  geom_point(
    data = conditional_nodes,
    aes(x_plot, y_plot, size = selected_runs, colour = layer),
    alpha = 0.46
  ) +
  geom_point(
    data = persistent_nodes,
    aes(x_plot, y_plot, size = selected_runs, fill = layer, colour = layer),
    shape = 21, stroke = 0.22, alpha = 0.84
  ) +
  geom_point(
    data = core_nodes,
    aes(x_plot, y_plot, colour = layer),
    shape = 21, fill = NA, size = 2.45, stroke = 0.36, alpha = 0.36
  ) +
  geom_point(
    data = core_nodes,
    aes(x_plot, y_plot, size = selected_runs, fill = layer, colour = layer),
    shape = 21, stroke = 0, alpha = 0.88
  ) +
  geom_point(
    data = bridge_nodes,
    aes(x_plot, y_plot),
    shape = 21, fill = NA, colour = interface_dark,
    size = 3.35, stroke = 0.72
  ) +
  geom_point(
    data = marker_legend,
    aes(x = Inf, y = Inf, shape = marker),
    alpha = 0, show.legend = TRUE, inherit.aes = FALSE
  ) +
  ggrepel::geom_label_repel(
    data = filter(network_labels, label_group == "host"),
    aes(label_x, label_y, label = label),
    colour = host_dark, fill = "white",
    linewidth = 0.16, label.r = unit(1.0, "mm"), label.padding = unit(1.15, "mm"),
    size = 2.00, fontface = "bold", seed = 20260714,
    box.padding = unit(2.4, "mm"), point.padding = unit(1.8, "mm"),
    min.segment.length = 0, max.overlaps = Inf, force = 2.0,
    nudge_x = -0.05, direction = "both",
    segment.colour = host_dark, segment.alpha = 0.70, segment.size = 0.26
  ) +
  ggrepel::geom_label_repel(
    data = filter(network_labels, label_group == "pathogen"),
    aes(label_x, label_y, label = label),
    colour = pathogen_dark, fill = "white",
    linewidth = 0.16, label.r = unit(1.0, "mm"), label.padding = unit(1.15, "mm"),
    size = 2.00, fontface = "bold", seed = 20260714,
    box.padding = unit(2.4, "mm"), point.padding = unit(1.8, "mm"),
    min.segment.length = 0, max.overlaps = Inf, force = 2.0,
    nudge_x = 0.06, direction = "both",
    segment.colour = pathogen_dark, segment.alpha = 0.70, segment.size = 0.26
  ) +
  ggrepel::geom_label_repel(
    data = filter(network_labels, label_group == "interface"),
    aes(label_x, label_y, label = label),
    colour = interface_dark, fill = "white",
    linewidth = 0.16, label.r = unit(1.0, "mm"), label.padding = unit(1.15, "mm"),
    size = 2.00, fontface = "bold", seed = 20260714,
    box.padding = unit(2.4, "mm"), point.padding = unit(1.8, "mm"),
    min.segment.length = 0, max.overlaps = Inf, force = 2.0,
    nudge_y = 0.07, direction = "both",
    segment.colour = interface_dark, segment.alpha = 0.70, segment.size = 0.26
  ) +
  annotate("text", x = host_title_x, y = title_y, label = "HOST", colour = host_dark, size = 2.38, fontface = "bold") +
  annotate("text", x = pathogen_title_x, y = title_y, label = "PATHOGEN", colour = pathogen_dark, size = 2.38, fontface = "bold") +
  scale_size_continuous(
    range = c(0.32, 2.20), breaks = c(6, 9, 18),
    name = "Retained runs"
  ) +
  scale_shape_manual(
    values = c("Pan-run core" = 21, "Hub" = 21, "Interface endpoint" = 21),
    name = "Node emphasis"
  ) +
  scale_colour_manual(values = c(host = host_dark, pathogen = pathogen_dark, interface = interface_dark), guide = "none") +
  scale_fill_manual(values = c(host = host_fill, pathogen = pathogen_fill), guide = "none") +
  guides(
    size = guide_legend(
      order = 1, nrow = 1, title.position = "left",
      override.aes = list(shape = 21, colour = ink, fill = NA, stroke = 0.70, alpha = 1)
    ),
    shape = guide_legend(
      order = 2, nrow = 1, title.position = "left",
      theme = theme(legend.title = element_text(margin = margin(l = 5, r = 1))),
      override.aes = list(
        alpha = 1,
        colour = c(ink, host_dark, interface_dark),
        fill = NA,
        size = c(2.6, 3.4, 3.0),
        stroke = c(0.78, 0.70, 0.90)
      )
    )
  ) +
  coord_equal(xlim = c(-1.20, 1.20), ylim = c(ylim_bottom, ylim_top), clip = "off") +
  theme_void(base_size = 8, base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.box.spacing = unit(1.8, "mm"),
    legend.margin = margin(t = 1.2, r = 0, b = 0, l = 0),
    legend.title = element_text(size = 6.1, colour = ink, face = "bold"),
    legend.text = element_text(size = 5.8, colour = ink),
    legend.key = element_blank(),
    legend.key.width = unit(3.4, "mm"),
    legend.key.height = unit(2.8, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    plot.margin = margin(7, 5, 4, 5),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(
  file.path(out, "Figure_4_panel_a_parameter_union_network.pdf"), panel,
  width = 190, height = 122, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE
)
ggsave(
  file.path(out, "Figure_4_panel_a_parameter_union_network.png"), panel,
  width = 190, height = 122, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE
)
}
cat(
  "Rendered Figure_4_panel_a_parameter_union_network: ", nrow(plot_nodes), " proteins; ",
  nrow(edge_plot), " edges; ", nrow(persistent_nodes) + nrow(core_nodes),
  " persistent nodes (104/95 at >=50%); cross-layer persistent edge is p104–MAPRE1 only.\n",
  "title_y=", round(title_y, 4),
  " host_x=", round(host_title_x, 4),
  " pathogen_x=", round(pathogen_title_x, 4),
  " y_plot_max=", round(max(plot_nodes$y_plot), 4),
  " host_x_range=", paste(round(range(host_xs), 3), collapse = ":"),
  " pathogen_x_range=", paste(round(range(pathogen_xs), 3), collapse = ":"),
  " ylim=c(", ylim_bottom, ", ", round(ylim_top, 4), ")",
  " min_d_host=", round(min(dist(plot_nodes[plot_nodes$layer == "host", c("x_plot", "y_plot")])), 4),
  " min_d_path=", round(min(dist(plot_nodes[plot_nodes$layer == "pathogen", c("x_plot", "y_plot")])), 4),
  "\n",
  sep = ""
)
