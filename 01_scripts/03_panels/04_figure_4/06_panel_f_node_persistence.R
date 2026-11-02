#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(grid)
})

ofiles <- vapply(
  sys.frames(),
  function(x) if (is.null(x$ofile)) "" else as.character(x$ofile[[1]]),
  character(1)
)
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- if (any(nzchar(ofiles))) {
  normalizePath(tail(ofiles[nzchar(ofiles)], 1), winslash = "/")
} else if (length(arg)) {
  normalizePath(sub("^--file=", "", arg[[1]]), winslash = "/")
} else {
  normalizePath(
    "01_scripts/03_panels/04_figure_4/06_panel_f_node_persistence.R",
    winslash = "/"
  )
}
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

data_override <- Sys.getenv("FIGURE_DATA_DIR", "")
data_candidates <- c(
  data_override,
  p("04_plot_data")
)
data_dir <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)][1]
if (is.na(data_dir)) stop("Run 05_build_figure_data.R first.")
read_fig <- function(file) {
  read_tsv(
    file.path(data_dir, "04_figure_4", file),
    show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN")
  )
}

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", "")
output_dir <- if (nzchar(output_override)) {
  output_override
} else {
  p("06_panel_output", "04_figure_4")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

interface_ring <- "#927D96"
reference_line <- line
display_n <- 40L

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
  # UniProt tail only when a gene-like prefix remains (PAFAH1B2P68402, not HSP90AA1).
  stripped <- str_replace(z, "[OPQ][0-9][A-Z0-9]{3}[0-9]$", "")
  z <- if_else(nchar(stripped) >= 3L, stripped, z)
  str_squish(z)
}

pick_gene <- function(gene_symbol, display_raw, pathogen_raw, node_id) {
  gs <- tidy_gene(gene_symbol)
  dr <- tidy_gene(display_raw)
  pr <- tidy_gene(pathogen_raw)
  fallback <- str_remove(as.character(node_id), "^(9913|5874)\\.")
  dplyr::case_when(
    looks_like_gene(gs) ~ gs,
    looks_like_gene(dr) ~ dr,
    looks_like_gene(pr) ~ pr,
    nzchar(gs) & !str_detect(gs, "^TA[0-9]") ~ gs,
    nzchar(dr) & !str_detect(dr, "^TA[0-9]") ~ dr,
    nzchar(pr) & !str_detect(pr, "^TA[0-9]") ~ pr,
    nzchar(gs) ~ gs,
    nzchar(dr) ~ dr,
    nzchar(pr) ~ pr,
    TRUE ~ fallback
  )
}

# Leftover NET n from pretty_module: Figure 1 aliases plus named programs for 10/11/16.
leftover_module <- c(
  H_NET_3 = "Trafficking",
  H_NET_6 = "Cell cycle",
  H_NET_9 = "Translation",
  H_NET_10 = "Nucleolus / rRNA",
  H_NET_11 = "Secreted / plasma",
  H_NET_13 = "Chaperone",
  H_NET_16 = "Nuclear transport",
  P_NET_2 = "Ribosome / translation"
)

module_name_of <- function(id, label) {
  id <- str_squish(coalesce(as.character(id), ""))
  label <- str_squish(coalesce(as.character(label), ""))
  from_label <- pretty_module(label)
  from_id <- pretty_module(id)
  leftover <- function(z, raw) {
    z <- str_squish(coalesce(z, ""))
    !nzchar(z) | z %in% c("Unassigned", raw, id, label) | str_detect(z, "^NET [0-9]+$")
  }
  dplyr::case_when(
    !leftover(from_label, label) ~ from_label,
    !leftover(from_id, id) ~ from_id,
    id %in% names(leftover_module) ~ unname(leftover_module[id]),
    label %in% names(leftover_module) ~ unname(leftover_module[label]),
    nzchar(from_label) ~ from_label,
    TRUE ~ from_id
  )
}

freq_bin <- function(x) {
  dplyr::case_when(
    x <= 0.5 + 1e-8 ~ "50%",
    x < 0.75 ~ "50-75%",
    x < 1 - 1e-8 ~ "75-99%",
    TRUE ~ "100%"
  )
}

retained_axis_genes <- unique(c(
  core_endpoint_labels,
  unlist(strsplit(c(core_axis_ids, extended_axis_ids, context_axis_ids), "--", fixed = TRUE))
))

raw_nodes <- ensure_cols(
  read_fig("08_panel_f_pcsf_node_persistence.tsv"),
  c(
    "display_label", "pathogen_name", "gene_symbol", "module_id", "module_label",
    "forest_degree", "network_degree", "endpoint", "hub"
  )
)

n_full <- nrow(raw_nodes)
n_host_full <- sum(raw_nodes$layer == "host", na.rm = TRUE)
n_path_full <- sum(raw_nodes$layer == "pathogen", na.rm = TRUE)
if (n_full != 104L) {
  warning("Expected 104 persistent PCSF nodes; found ", n_full, ".")
}
stopifnot(n_full >= 50L)

nodes_all <- raw_nodes %>%
  mutate(
    layer = factor(
      layer,
      levels = c("host", "pathogen"),
      labels = c("Host", "Pathogen")
    ),
    endpoint = flag_col(endpoint),
    hub = flag_col(hub),
    forest_degree = num_col(forest_degree),
    network_degree = num_col(network_degree),
    selection_frequency = num_col(selection_frequency),
    gene_symbol = str_squish(coalesce(as.character(gene_symbol), "")),
    display_raw = str_squish(coalesce(as.character(display_label), "")),
    pathogen_raw = str_squish(coalesce(as.character(pathogen_name), "")),
    display_label = pick_gene(gene_symbol, display_raw, pathogen_raw, node_id),
    module_key = coalesce(na_if(as.character(module_id), ""), module_label),
    module_name = module_name_of(module_id, module_label),
    node_label = paste0(display_label, "  ·  ", module_name),
    axis_keep = endpoint |
      display_label %in% retained_axis_genes |
      tidy_gene(gene_symbol) %in% retained_axis_genes |
      tidy_gene(display_raw) %in% retained_axis_genes |
      tidy_gene(pathogen_raw) %in% retained_axis_genes,
    ring_role = case_when(
      endpoint | axis_keep ~ "Interface endpoint",
      hub ~ "Hub",
      TRUE ~ "Backbone node"
    ),
    freq_bin = freq_bin(selection_frequency),
    must_keep = axis_keep | hub | selection_frequency < (1 - 1e-8)
  )

full_bins <- nodes_all %>%
  count(layer, freq_bin, name = "n") %>%
  arrange(layer, freq_bin)

nodes <- nodes_all %>%
  group_by(layer) %>%
  arrange(
    desc(must_keep),
    desc(hub),
    desc(forest_degree),
    desc(network_degree),
    display_label,
    .by_group = TRUE
  ) %>%
  mutate(
    n_must = sum(must_keep),
    rest_rank = if_else(must_keep, 0L, row_number() - n_must),
    keep = must_keep | rest_rank <= pmax(0L, display_n - n_must)
  ) %>%
  filter(keep) %>%
  arrange(selection_frequency, forest_degree, display_label, .by_group = TRUE) %>%
  mutate(y_key = paste(layer, sprintf("%03d", row_number()), node_label, sep = "|")) %>%
  ungroup() %>%
  mutate(y_key = factor(y_key, levels = unique(y_key)))

stopifnot(!any(str_detect(nodes$node_label, "H_NET|P_NET|NET [0-9]|TA_15795")))
stopifnot(!any(str_detect(nodes$display_label, "^HS$|P68402|TA_15795")))
needed_endpoints <- intersect(
  c("p104", "MAPRE1", "CD2AP", "PKM", "TA20980"),
  unique(c(nodes_all$display_label, tidy_gene(nodes_all$gene_symbol)))
)
stopifnot(all(needed_endpoints %in% nodes$display_label))

n_host_draw <- sum(nodes$layer == "Host")
n_path_draw <- sum(nodes$layer == "Pathogen")
drawn_bins <- nodes %>%
  count(layer, freq_bin, name = "n") %>%
  arrange(layer, freq_bin)
fd_lo <- min(nodes$forest_degree)
fd_hi <- max(nodes$forest_degree)
size_breaks <- unique(as.integer(c(
  if (fd_hi >= 1) 1 else 0,
  if (fd_hi >= 8) 8 else as.integer(round(stats::median(nodes$forest_degree[nodes$forest_degree > 0]))),
  fd_hi
)))
size_breaks <- sort(size_breaks[size_breaks > 0 & size_breaks <= fd_hi])
if (!length(size_breaks)) size_breaks <- max(1L, as.integer(fd_hi))
drawn_endpoints <- nodes %>%
  filter(axis_keep) %>%
  arrange(layer, display_label) %>%
  distinct(layer, display_label) %>%
  transmute(txt = paste0(display_label, " (", as.character(layer), ")")) %>%
  pull(txt)
drawn_hubs <- nodes %>%
  filter(hub) %>%
  arrange(layer, display_label) %>%
  distinct(layer, display_label) %>%
  transmute(txt = paste0(display_label, " (", as.character(layer), ")")) %>%
  pull(txt)

host_modules <- sort(unique(nodes$module_key[nodes$layer == "Host"]))
pathogen_modules <- sort(unique(nodes$module_key[nodes$layer == "Pathogen"]))
module_colours <- c(
  setNames(
    grDevices::colorRampPalette(c(blue_pale, blue_dark))(length(host_modules)),
    host_modules
  ),
  setNames(
    grDevices::colorRampPalette(c(red_pale, red_dark))(length(pathogen_modules)),
    pathogen_modules
  )
)

plot <- ggplot(nodes, aes(selection_frequency, y_key)) +
  geom_vline(
    xintercept = c(0.5, 0.75, 1),
    colour = reference_line,
    linewidth = 0.42
  ) +
  geom_segment(
    aes(x = 0.5, xend = selection_frequency, yend = y_key),
    colour = reference_line,
    linewidth = 0.72,
    lineend = "round"
  ) +
  geom_point(
    aes(size = forest_degree, fill = module_key, colour = ring_role),
    shape = 21,
    stroke = 0.55
  ) +
  facet_wrap(~layer, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(0.46, 1.06),
    breaks = c(0.5, 0.75, 1),
    labels = percent,
    expand = expansion(mult = 0)
  ) +
  scale_y_discrete(
    labels = function(x) str_remove(x, "^[^|]+\\|[0-9]+\\|"),
    expand = expansion(add = 0.55)
  ) +
  scale_fill_manual(values = module_colours, guide = "none") +
  scale_colour_manual(
    values = c(
      "Backbone node" = frame,
      "Hub" = ink,
      "Interface endpoint" = interface_ring
    ),
    breaks = c("Backbone node", "Hub", "Interface endpoint"),
    name = "Outer ring"
  ) +
  scale_size_continuous(
    name = "Forest degree",
    range = c(1.7, 3.6),
    limits = c(0, max(1, fd_hi)),
    breaks = size_breaks
  ) +
  labs(
    x = "Selection frequency across PCSF parameter runs",
    y = NULL
  ) +
  guides(
    colour = guide_legend(
      order = 1,
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(shape = 21, fill = NA, size = 3.6, stroke = 0.9)
    ),
    size = guide_legend(
      order = 2,
      title.position = "left",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(6.2, "mm"),
      keyheight = unit(5.2, "mm"),
      override.aes = list(shape = 21, fill = NA, colour = ink, stroke = 0.60)
    )
  ) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.35),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 5.9, colour = ink),
    axis.text.y = element_text(
      size = 5.0, colour = ink, lineheight = 0.88,
      margin = margin(r = 2.5)
    ),
    axis.title.x = element_text(size = 6.6, margin = margin(t = 4)),
    strip.background = element_rect(
      fill = "#F8F9FA", colour = frame, linewidth = 0.35
    ),
    strip.text = element_text(size = 6.8, face = "bold", colour = ink),
    panel.spacing.x = unit(6, "mm"),
    legend.position = "bottom",
    legend.location = "plot",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.justification = "center",
    legend.key = element_blank(),
    legend.spacing.x = unit(6.5, "mm"),
    legend.spacing.y = unit(1.2, "mm"),
    legend.box.spacing = unit(2.4, "mm"),
    legend.margin = margin(t = 1.5, r = 2, b = 0, l = 2),
    legend.title = element_text(size = 5.7, face = "bold", margin = margin(r = 4)),
    legend.text = element_text(size = 5.4),
    legend.key.width = unit(5.2, "mm"),
    legend.key.height = unit(4.4, "mm"),
    plot.margin = margin(5, 10, 4, 12)
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(
  file.path(output_dir, "Figure_4_panel_f_node_persistence.pdf"),
  plot,
  width = 340, height = 198, units = "mm",
  device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure_4_panel_f_node_persistence.png"),
  plot,
  width = 340, height = 198, units = "mm", dpi = 600,
  device = ragg::agg_png, bg = "white"
)
}

fmt_bins <- function(d) {
  paste(
    sprintf("%s %s n=%d", as.character(d$layer), d$freq_bin, d$n),
    collapse = "; "
  )
}

cat(
  "Rendered Figure_4_panel_f_node_persistence: display ",
  n_host_draw, " host + ", n_path_draw, " pathogen of ",
  n_full, " persistent (>=50%) PCSF nodes (",
  n_host_full, " host / ", n_path_full, " pathogen).",
  " Full bins: ", fmt_bins(full_bins),
  ". Drawn bins: ", fmt_bins(drawn_bins),
  ". Endpoints: ",
  if (length(drawn_endpoints)) paste(drawn_endpoints, collapse = ", ") else "none in persistent set",
  ". Hubs: ",
  if (length(drawn_hubs)) paste(drawn_hubs, collapse = ", ") else "none",
  ".\n",
  sep = ""
)
