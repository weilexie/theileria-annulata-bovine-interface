#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(grid)
  library(patchwork)
})

ofiles <- vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else as.character(x$ofile[[1]]), character(1))
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- if (any(nzchar(ofiles))) normalizePath(tail(ofiles[nzchar(ofiles)], 1), winslash = "/") else if (length(arg)) normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/") else normalizePath("01_scripts/03_panels/04_figure_4/02_panel_b_backbone_convergence.R", winslash = "/")
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
if (is.na(data_dir)) stop("Run 09 first.")
read_fig <- function(x) read_tsv(file.path(data_dir, "04_figure_4", x), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
od <- Sys.getenv("FIGURE_OUTPUT_DIR", "")
out <- if (nzchar(od)) od else p("06_panel_output", "04_figure_4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

looks_like_gene <- function(x) {
  z <- str_squish(coalesce(as.character(x), ""))
  nzchar(z) &
    !str_detect(z, "^TA[0-9]{4,}$") &
    !str_detect(z, "^TA_") &
    !str_detect(z, "H_NET") &
    !str_detect(z, "^(9913|5874)\\.") &
    !str_detect(z, "^ENS") &
    !str_detect(z, "^A0A") &
    !str_detect(z, "^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
}

locked_keep <- c("p104", "MAPRE1", "CD2AP")

# Open packed FR cores without changing the bounding box.
expand_core <- function(x, y, power) {
  cx <- stats::median(x, na.rm = TRUE)
  cy <- stats::median(y, na.rm = TRUE)
  dx <- x - cx
  dy <- y - cy
  r <- sqrt(dx * dx + dy * dy)
  r0 <- max(r, na.rm = TRUE)
  r_use <- pmax(r, 1e-9)
  r_new <- r0 * (r_use / pmax(r0, 1e-9))^power
  list(x = cx + dx / r_use * r_new, y = cy + dy / r_use * r_new)
}

# Pull extreme FR islands onto the main hull so one cloud fills the box.
compact_outliers <- function(x, y, q = 0.90) {
  cx <- stats::median(x, na.rm = TRUE)
  cy <- stats::median(y, na.rm = TRUE)
  dx <- x - cx
  dy <- y - cy
  r <- sqrt(dx * dx + dy * dy)
  r_cap <- as.numeric(stats::quantile(r, q, na.rm = TRUE))
  s <- pmin(r, r_cap) / pmax(r, 1e-9)
  list(x = cx + dx * s, y = cy + dy * s)
}

# Same-layer collision push; clamp to HOST-left / PATHOGEN-right boxes.
collide_layers <- function(xp, yp, layer, min_sep, niter, xlim_host, xlim_path, ylim) {
  xy <- cbind(as.numeric(xp), as.numeric(yp))
  n <- nrow(xy)
  ly <- as.character(layer)
  set.seed(14)
  xy <- xy + matrix(stats::rnorm(n * 2L, 0, 0.001), n, 2L)
  for (iter in seq_len(niter)) {
    moved <- FALSE
    for (i in seq_len(n - 1L)) {
      j <- (i + 1L):n
      same <- ly[j] == ly[i]
      if (!any(same)) next
      dxy <- sweep(xy[j, , drop = FALSE], 2, xy[i, ], "-")
      dist <- sqrt(rowSums(dxy * dxy))
      hit <- which(same & dist < min_sep)
      if (!length(hit)) next
      moved <- TRUE
      dist_hit <- pmax(dist[hit], 1e-6)
      ux <- dxy[hit, 1] / dist_hit
      uy <- dxy[hit, 2] / dist_hit
      push <- (min_sep - dist_hit) / 2
      xy[i, ] <- xy[i, ] - c(mean(push * ux), mean(push * uy))
      jj <- j[hit]
      xy[jj, 1] <- xy[jj, 1] + push * ux
      xy[jj, 2] <- xy[jj, 2] + push * uy
    }
    host <- ly == "host"
    xy[host, 1] <- pmin(pmax(xy[host, 1], xlim_host[[1]]), xlim_host[[2]])
    xy[!host, 1] <- pmin(pmax(xy[!host, 1], xlim_path[[1]]), xlim_path[[2]])
    xy[, 2] <- pmin(pmax(xy[, 2], ylim[[1]]), ylim[[2]])
    if (!moved) break
  }
  xy
}

# Keep a spatially thinned subset so in-dot text does not collide.
thin_inner <- function(d, min_sep) {
  if (!nrow(d)) return(d)
  d <- d %>% arrange(desc(hub_score), nchar(display_label), display_label)
  keep <- rep(FALSE, nrow(d))
  for (i in seq_len(nrow(d))) {
    if (!any(keep)) {
      keep[[i]] <- TRUE
      next
    }
    dx <- d$xp[[i]] - d$xp[keep]
    dy <- d$yp[[i]] - d$yp[keep]
    if (min(dx * dx + dy * dy) >= min_sep * min_sep) keep[[i]] <- TRUE
  }
  d[keep, , drop = FALSE]
}

raw_nodes <- ensure_cols(
  read_fig("04_panel_b_pcsf_nodes.tsv"),
  c("display_label", "pathogen_name", "canonical_gene", "gene_symbol", "x0", "x", "y")
) %>%
  mutate(
    gene_symbol = str_squish(coalesce(as.character(gene_symbol), "")),
    display_raw = str_squish(coalesce(as.character(display_label), "")),
    pathogen_raw = str_squish(coalesce(as.character(pathogen_name), "")),
    canonical_raw = str_squish(coalesce(as.character(canonical_gene), "")),
    display_label = endpoint_display(sanitize_label(case_when(
      looks_like_gene(gene_symbol) ~ gene_symbol,
      looks_like_gene(display_raw) ~ display_raw,
      looks_like_gene(pathogen_raw) ~ pathogen_raw,
      looks_like_gene(canonical_raw) ~ canonical_raw,
      nzchar(gene_symbol) ~ gene_symbol,
      nzchar(display_raw) ~ display_raw,
      nzchar(pathogen_raw) ~ pathogen_raw,
      nzchar(canonical_raw) ~ canonical_raw,
      TRUE ~ str_remove(node_id, "^(9913|5874)\\.")
    )))
  )
raw_edges <- read_fig("05_panel_b_pcsf_edges.tsv") %>% select(source, target, selection_frequency)
union_nodes <- read_fig("02_panel_a_pcsf_union_nodes.tsv")
union_edges <- read_fig("03_panel_a_pcsf_union_edges.tsv")
run_n <- as.integer(round(max(raw_nodes$selected_runs / raw_nodes$selection_frequency, na.rm = TRUE)))

stopifnot(
  run_n == 18L,
  sum(union_nodes$selection_frequency >= 0.5) == 104L,
  sum(union_edges$selection_frequency >= 0.5) == 95L
)

convergence <- bind_rows(
  tibble(
    min_runs = seq_len(run_n),
    retained = vapply(seq_len(run_n), function(k) sum(union_nodes$selected_runs >= k), integer(1)),
    element = "Proteins"
  ),
  tibble(
    min_runs = seq_len(run_n),
    retained = vapply(seq_len(run_n), function(k) sum(union_edges$selected_runs >= k), integer(1)),
    element = "Edges"
  )
)
convergence_key <- convergence %>%
  filter(min_runs %in% c(1L, 9L, 18L))

incidence <- bind_rows(
  raw_edges %>% transmute(node_id = source, edge_weight = selection_frequency),
  raw_edges %>% transmute(node_id = target, edge_weight = selection_frequency)
) %>%
  group_by(node_id) %>%
  summarise(inc_degree = dplyr::n(), weighted_degree = sum(edge_weight), .groups = "drop")

nodes <- raw_nodes %>%
  left_join(incidence, by = "node_id") %>%
  mutate(degree = coalesce(inc_degree, 0L), weighted_degree = coalesce(weighted_degree, 0)) %>%
  group_by(layer) %>%
  mutate(
    hub_score = 0.50 * percent_rank(degree) + 0.35 * percent_rank(weighted_degree) + 0.15 * percent_rank(selection_frequency)
  ) %>%
  ungroup()

# Folded `x` stacks mirrored FR nodes. Open cores, drop far islands, then
# min-max place into HOST-left / PATHOGEN-right boxes (gold two-cloud grammar).
# Only unstack true overlaps; a hard min_sep hex-packs the cloud into a disk.
host_xlim <- c(-0.98, -0.11)
path_xlim <- c(0.11, 0.98)
layer_ylim <- c(-0.92, 0.88)
nodes$xp <- NA_real_
nodes$yp <- NA_real_
for (ly in c("host", "pathogen")) {
  idx <- which(nodes$layer == ly)
  opened <- expand_core(
    nodes$x0[idx], nodes$y[idx],
    power = if (ly == "pathogen") 0.32 else 0.44
  )
  compact <- compact_outliers(opened$x, opened$y, q = if (ly == "pathogen") 0.74 else 0.92)
  nodes$xp[idx] <- rescale(compact$x, to = if (ly == "host") host_xlim else path_xlim)
  nodes$yp[idx] <- rescale(compact$y, to = layer_ylim)
}
xy <- collide_layers(
  nodes$xp, nodes$yp, nodes$layer,
  min_sep = 0.078, niter = 48,
  xlim_host = c(-1.02, -0.07),
  xlim_path = c(0.07, 1.02),
  ylim = c(-0.95, 0.91)
)
nodes$xp <- xy[, 1]
nodes$yp <- xy[, 2]
# Gold rescales both clouds onto the same y-span so HOST/PATHOGEN are equal height.
for (ly in c("host", "pathogen")) {
  idx <- which(nodes$layer == ly)
  nodes$yp[idx] <- rescale(nodes$yp[idx], to = c(-0.89, 0.88))
}
nn_mat <- as.matrix(dist(cbind(nodes$xp, nodes$yp)))
diag(nn_mat) <- Inf
nodes$nn <- apply(nn_mat, 1, min)

# Gold 4b: 4 hubs/layer get geom_label_repel beside the node (inside the
# panel). Other high-frequency genes are written inside the dot. Locked
# interface ends are boxed if they are not already hubs.
hubs <- nodes %>%
  group_by(layer) %>%
  slice_max(hub_score, n = 4, with_ties = FALSE) %>%
  ungroup()
hubs <- bind_rows(hubs, nodes %>% filter(display_label %in% locked_keep)) %>%
  distinct(node_id, .keep_all = TRUE)
hub_ids <- hubs$node_id

edges <- raw_edges %>%
  left_join(nodes %>% select(source = node_id, x_from = xp, y_from = yp, edge_layer = layer), by = "source") %>%
  left_join(nodes %>% select(target = node_id, x_to = xp, y_to = yp), by = "target")
stable_edges <- edges %>% filter(selection_frequency >= 0.944)
hub_edges <- edges %>% filter(source %in% hub_ids | target %in% hub_ids)
high <- nodes %>% filter(selection_frequency >= 0.999)
# Host gene symbols fit on one line. Pathogen loci are TA##### — wrap to
# two lines so the black text stays inside the circle.
fits_in_dot <- function(z) {
  nchar(z) <= 6L & !str_detect(z, "^TA[0-9]{4,}$") & !str_detect(z, "^TA_")
}
compact_ta <- function(z) {
  ta <- str_match(z, "^TA([0-9]{4,})$")[, 2]
  ifelse(!is.na(ta), paste0("TA\n", ta), z)
}
pathogen_inner <- thin_inner(
  high %>%
    filter(layer == "pathogen", !node_id %in% hub_ids, nzchar(display_label)) %>%
    mutate(inner_label = compact_ta(display_label)),
  min_sep = 0.072
)
host_inner <- thin_inner(
  nodes %>%
    filter(
      layer == "host",
      selected_runs >= run_n - 1L,
      !node_id %in% hub_ids,
      fits_in_dot(display_label)
    ) %>%
    mutate(inner_label = display_label),
  min_sep = 0.11
)
heads <- tibble(layer = c("host", "pathogen"), x = c(-0.53, 0.53))

edge_ink <- "#8B949E"
g <- ggplot() +
  geom_segment(
    data = edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, alpha = selection_frequency, linewidth = selection_frequency),
    colour = edge_ink, lineend = "round"
  ) +
  geom_segment(
    data = stable_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, colour = edge_layer, linewidth = selection_frequency),
    alpha = 0.42, lineend = "round"
  ) +
  geom_segment(
    data = hub_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to, colour = edge_layer, linewidth = selection_frequency),
    alpha = 0.78, lineend = "round"
  ) +
  geom_point(data = hubs, aes(xp, yp, fill = layer), shape = 21, colour = "transparent", size = I(4.6), alpha = 0.22) +
  geom_point(data = high, aes(xp, yp, colour = layer), shape = 21, fill = NA, size = I(3.6), stroke = 0.55, alpha = 0.50) +
  geom_point(data = nodes, aes(xp, yp, size = selection_frequency, fill = layer), shape = 21, colour = ink, stroke = 0.28) +
  geom_point(data = hubs, aes(xp, yp, colour = layer), shape = 21, fill = NA, size = I(4.0), stroke = 0.85) +
  geom_text(
    data = filter(pathogen_inner, str_detect(inner_label, "\n")),
    aes(xp, yp, label = inner_label),
    colour = ink, size = 0.62, fontface = "bold", lineheight = 0.72
  ) +
  geom_text(
    data = filter(pathogen_inner, !str_detect(inner_label, "\n")),
    aes(xp, yp, label = inner_label),
    colour = ink, size = 0.78, fontface = "bold", lineheight = 0.82
  ) +
  geom_text(
    data = host_inner,
    aes(xp, yp, label = inner_label),
    colour = ink, size = 0.55, fontface = "bold", lineheight = 0.78
  ) +
  geom_label_repel(
    data = hubs,
    aes(xp, yp, label = display_label, colour = layer),
    fill = scales::alpha("white", 0.92), linewidth = 0.16, label.r = unit(0.8, "mm"),
    size = 2.2, fontface = "bold", seed = 14, max.overlaps = Inf,
    min.segment.length = 0, box.padding = 0.32, point.padding = 0.28,
    force = 1.35, max.time = 3,
    xlim = c(-0.98, 0.98), ylim = c(-0.90, 0.86),
    segment.size = 0.32, segment.alpha = 0.80, segment.colour = ink
  ) +
  geom_text(data = heads, aes(x, 0.99, label = toupper(layer), colour = layer), fontface = "bold", size = 3) +
  scale_fill_manual(values = c(host = host_fill, pathogen = pathogen_fill), guide = "none") +
  scale_colour_manual(values = c(host = host_dark, pathogen = pathogen_dark), guide = "none") +
  scale_size_continuous(
    range = c(1.7, 2.9), breaks = c(0.5, 0.75, 1), labels = scales::label_percent(accuracy = 1),
    name = "Selection frequency"
  ) +
  scale_alpha_continuous(range = c(0.28, 0.58), guide = "none") +
  scale_linewidth_continuous(range = c(0.22, 0.80), guide = "none") +
  labs(
    title = "PCSF-derived backbone",
    x = "Minimum runs retained (of 18)"
  ) +
  guides(
    size = guide_legend(
      nrow = 1, title.position = "left",
      override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.62)
    )
  ) +
  coord_cartesian(xlim = c(-1.04, 1.04), ylim = c(-0.98, 1.04), clip = "off") +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.3),
    plot.title.position = "panel",
    plot.title = element_text(
      size = 7.0, colour = ink, face = "bold", hjust = 0.5, margin = margin(b = 3)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.key = element_blank(),
    legend.title = element_text(size = 6, margin = margin(r = 2.5)),
    legend.text = element_text(size = 5.6),
    legend.key.width = unit(4.0, "mm"),
    legend.key.height = unit(2.8, "mm"),
    legend.box.spacing = unit(1.8, "mm"),
    legend.margin = margin(t = 1.2),
    axis.title.x = element_text(
      size = 6.1, colour = "white", margin = margin(t = 4)
    ),
    axis.text.x = element_text(size = 5.6, colour = "white"),
    axis.ticks.length = unit(0.9, "mm"),
    plot.margin = margin(4, 4, 2, 4),
    plot.background = element_rect(fill = "white", colour = NA)
  )

convergence_plot <- ggplot(convergence, aes(min_runs, retained, colour = element, group = element)) +
  annotate(
    "rect", xmin = 8.5, xmax = 18.5, ymin = -Inf, ymax = Inf,
    fill = interface_pale, alpha = 0.58
  ) +
  geom_step(linewidth = 0.52, direction = "hv", lineend = "round") +
  geom_point(
    data = convergence_key,
    aes(fill = element), shape = 21, size = 2.0, stroke = 0.38
  ) +
  geom_text(
    data = convergence_key %>%
      filter(element == "Proteins") %>%
      mutate(
        label_x = if_else(min_runs == run_n, min_runs - 0.28, as.numeric(min_runs)),
        label_y = retained + if_else(min_runs == run_n, 8, 6)
      ),
    aes(label_x, label_y, label = retained), colour = host_dark, size = 1.82,
    inherit.aes = FALSE, fontface = "bold"
  ) +
  geom_text(
    data = convergence_key %>%
      filter(element == "Edges") %>%
      mutate(
        label_x = if_else(min_runs == run_n, min_runs - 0.28, as.numeric(min_runs)),
        label_y = retained - if_else(min_runs == run_n, 8, 6)
      ),
    aes(label_x, label_y, label = retained), colour = interface_dark, size = 1.82,
    inherit.aes = FALSE, fontface = "bold"
  ) +
  scale_colour_manual(
    values = c(Proteins = host_dark, Edges = interface_dark),
    breaks = c("Proteins", "Edges"), name = NULL,
    guide = guide_legend(
      nrow = 1,
      override.aes = list(shape = 21, fill = NA, size = 2.0, stroke = 0.62)
    )
  ) +
  scale_fill_manual(
    values = c(Proteins = host_pale, Edges = interface_pale),
    breaks = c("Proteins", "Edges"), guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(1, 9, 18), labels = c("1", ">=9", "18"),
    limits = c(0.55, 18.75), expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(4),
    limits = c(0, max(convergence$retained, na.rm = TRUE) * 1.12),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Cross-parameter convergence",
    x = "Minimum runs retained (of 18)", y = "Retained count"
  ) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    axis.line = element_line(colour = ink, linewidth = 0.30),
    axis.ticks = element_line(colour = ink, linewidth = 0.30),
    axis.ticks.length = unit(0.9, "mm"),
    axis.title = element_text(colour = ink, size = 6.1),
    axis.text = element_text(colour = ink, size = 5.6),
    plot.title.position = "panel",
    plot.title = element_text(
      colour = ink, size = 7.0, face = "bold", hjust = 0.5, margin = margin(b = 3)
    ),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.3),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.text = element_text(size = 5.6, colour = ink),
    legend.key.width = unit(4.0, "mm"),
    legend.key.height = unit(2.8, "mm"),
    legend.spacing.x = unit(1.0, "mm"),
    legend.box.spacing = unit(1.8, "mm"),
    legend.margin = margin(t = 1.2),
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 2, 4),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

composite <- convergence_plot + g +
  patchwork::plot_layout(widths = c(1, 1.55), heights = 1) +
  patchwork::plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", colour = NA))
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(
  file.path(out, "Figure_4_panel_b_backbone_convergence.pdf"), composite,
  width = 185, height = 108, units = "mm", device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(out, "Figure_4_panel_b_backbone_convergence.png"), composite,
  width = 185, height = 108, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white"
)
}
