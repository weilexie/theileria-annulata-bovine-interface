#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
  library(igraph)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) {
  z <- x$ofile
  if (is.null(z) || !length(z)) "" else as.character(z[[1]])
}, character(1))
script_path <- if (any(nzchar(source_files))) {
  normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path("01_scripts", "03_panels", "01_figure_1", "03_panel_c_host_network.R"), winslash = "/", mustWork = TRUE)
}
data_env <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
root <- if (nzchar(data_env)) normalizePath(file.path(data_env, "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
p <- function(...) file.path(root, ...)

data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
if (nzchar(data_override) && !dir.exists(data_override)) stop("FIGURE_DATA_DIR does not exist: ", data_override)
data_candidates <- c(data_override, p("04_plot_data"))
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(file.path(data_dir, "01_figure_1", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "01_figure_1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) {
  if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x))
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), x, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(stem, ".png")), x, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)
}

ink <- "#30363B"; line <- "#E4E9EC"; paper <- "#FFFFFF"
blue_dark <- "#B5C8D1"; blue <- "#D2DEE4"; blue_pale <- "#F1F5F7"
red_dark <- "#CDB0AD"; red <- "#E4D0CD"; red_pale <- "#FAF1F0"
theme_network <- theme_void(base_size = 8, base_family = "sans") +
  theme(legend.position = "none", plot.margin = margin(4, 8, 4, 8), plot.background = element_rect(fill = paper, colour = NA))

as_flag <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) replace(x, is.na(x), FALSE) else str_to_lower(coalesce(as.character(x), "")) %in% c("true", "t", "1", "yes", "y")
}
pretty_host_label <- function(display_label, protein_name = NA_character_) {
  dl <- coalesce(as.character(display_label), "")
  known <- c("MAPRE1", "PKM", "CD2AP", "RUVBL1", "TFG", "TUBA1A", "TUBG1", "PPP4C", "RBM39", "MED21", "CRBN", "FBXW7", "CLASP1")
  out <- dl
  keep <- dl %in% known | str_detect(dl, "^[A-Z][A-Za-z0-9-]{1,14}$") & !str_detect(dl, "^[A-Z0-9]{6,10}$")
  pn <- coalesce(as.character(protein_name), "")
  rp <- str_match(pn, regex("(40S|60S)(?: acidic)? ribosomal protein ([A-Z]?[SPL][0-9A-Z]*)", ignore_case = TRUE))
  ribo <- !is.na(rp[, 3])
  subu <- toupper(rp[, 3])
  ribo_lab <- case_when(
    ribo & str_starts(subu, "S") ~ paste0("RPS", str_remove(subu, "^S")),
    ribo & str_starts(subu, "L") ~ paste0("RPL", str_remove(subu, "^L")),
    ribo & str_starts(subu, "P") ~ paste0("RPLP", str_remove(subu, "^P")),
    TRUE ~ NA_character_
  )
  short <- str_trim(str_match(pn, "^([^,;]+)")[, 2])
  short_ok <- !is.na(short) & nchar(short) > 0 & nchar(short) <= 22 & !str_detect(short, "^[A-Z0-9]{6,}$")
  dplyr::case_when(
    keep ~ dl,
    !is.na(ribo_lab) ~ ribo_lab,
    short_ok ~ short,
    TRUE ~ dl
  )
}

nodes <- read_fig("06_panel_g_host_network_nodes.tsv")
edges <- read_fig("07_panel_g_host_network_edges.tsv")
required_node <- c("node_id", "strict_member", "relaxed_member", "degree", "strength", "module_id")
missing_node <- setdiff(required_node, names(nodes))
if (length(missing_node)) stop("F1c host network nodes missing columns: ", paste(missing_node, collapse = ", "))
if (!"display_label" %in% names(nodes)) {
  label_src <- intersect(c("gene_symbol", "preferred_name"), names(nodes))
  nodes$display_label <- if (length(label_src)) dplyr::coalesce(!!!nodes[label_src], nodes$node_id) else nodes$node_id
}
if (!"module_label" %in% names(nodes)) nodes$module_label <- as.character(nodes$module_id)
if (!"module_biological_label" %in% names(nodes)) nodes$module_biological_label <- nodes$module_label
if (!"weight" %in% names(edges)) edges$weight <- 1
if ("gene_symbol" %in% names(nodes)) {
  nodes <- nodes %>% mutate(display_label = coalesce(na_if(display_label, ""), gene_symbol, node_id))
}
if (!"protein_name" %in% names(nodes)) nodes$protein_name <- NA_character_
nodes$display_label <- pretty_host_label(nodes$display_label, nodes$protein_name)
if (!"network_role" %in% names(nodes)) nodes$network_role <- NA_character_
nodes <- nodes %>%
  mutate(
    strict_member = as_flag(strict_member),
    relaxed_member = as_flag(relaxed_member) | strict_member,
    degree = as.numeric(degree),
    strength = as.numeric(strength),
    network_role = case_when(
      strict_member | network_role %in% c("Core", "strict") ~ "Core",
      relaxed_member | network_role %in% c("Recurrent scaffold", "relaxed") ~ "Recurrent scaffold",
      TRUE ~ coalesce(na_if(as.character(network_role), ""), "Network context")
    )
  )
graph_all <- graph_from_data_frame(edges %>% select(source, target), directed = FALSE, vertices = nodes %>% transmute(name = node_id))
component_info <- components(graph_all)
main_component <- which.max(component_info$csize)
main_ids <- names(component_info$membership)[component_info$membership == main_component]
first_known <- function(x, fallback) {
  known <- unique(na.omit(x))
  if (length(known)) as.character(known[[1]]) else fallback
}
program_label <- function(x) {
  z <- coalesce(x, "")
  case_when(
    str_detect(z, "host_signaling_broad|signal_transduction") ~ "Signaling",
    str_detect(z, "host_chromatin|chromatin") ~ "Chromatin",
    str_detect(z, "host_transcription") ~ "Transcription",
    str_detect(z, "host_rna_processing|rna_processing") ~ "RNA processing",
    str_detect(z, "host_ubiquitin|proteostasis") ~ "Proteostasis",
    str_detect(z, "host_microtubule|cytoskeleton|microtubule") ~ "Microtubule",
    str_detect(z, "host_metabolism|mitochondrial_energy") ~ "Metabolism",
    str_detect(z, "H_NET_9") ~ "Translation",
    str_detect(z, "H_NET_13") ~ "Chaperone",
    str_detect(z, "H_NET_6") ~ "Cell cycle",
    str_detect(z, "H_NET_3") ~ "Trafficking",
    TRUE ~ "Other"
  )
}

main_nodes <- nodes %>%
  filter(node_id %in% main_ids) %>%
  mutate(module_key = if_else(is.na(module_id), "Other", as.character(module_id)))
endpoint_ids <- main_nodes %>%
  filter(
    (relaxed_member | strict_member),
    degree > 0,
    str_to_upper(coalesce(display_label, "")) %in% c("MAPRE1", "PKM", "CD2AP")
  ) %>%
  pull(node_id)
non_ribo_hubs <- main_nodes %>%
  filter(relaxed_member | strict_member, degree > 0) %>%
  filter(!str_detect(coalesce(display_label, ""), "^RP[SLP]")) %>%
  arrange(desc(strict_member), desc(strength), desc(degree)) %>%
  slice_head(n = 8) %>%
  pull(node_id)
scaffold_hubs <- unique(c(endpoint_ids, non_ribo_hubs))
main_nodes <- main_nodes %>% mutate(scaffold_hub = node_id %in% scaffold_hubs)

module_summary <- main_nodes %>%
  group_by(module_key) %>%
  summarise(
    node_n = n(), recurrent_n = sum(relaxed_member), core_n = sum(strict_member),
    module_label = first_known(module_label, if_else(first(module_key) == "Other", "Other", paste0("H", first(module_key)))),
    module_biological_label = first_known(module_biological_label, "Other program"),
    .groups = "drop"
  ) %>%
  arrange(desc(node_n)) %>%
  mutate(layout_rank = row_number())
n_inner <- min(6L, max(nrow(module_summary) - 1L, 0L))
n_outer <- max(nrow(module_summary) - 1L - n_inner, 0L)
module_summary <- module_summary %>%
  mutate(
    ring = case_when(layout_rank == 1 ~ "centre", layout_rank <= 1 + n_inner ~ "inner", TRUE ~ "outer"),
    ring_index = case_when(ring == "inner" ~ layout_rank - 1L, ring == "outer" ~ layout_rank - 1L - n_inner, TRUE ~ 1L),
    centre_x = case_when(
      ring == "centre" ~ 0,
      ring == "inner" ~ 0.64 * cos(2 * pi * (ring_index - 1) / pmax(n_inner, 1) + pi / 6),
      TRUE ~ 1.04 * cos(2 * pi * (ring_index - 1) / pmax(n_outer, 1) + pi / 18)
    ),
    centre_y = case_when(
      ring == "centre" ~ 0,
      ring == "inner" ~ 0.64 * sin(2 * pi * (ring_index - 1) / pmax(n_inner, 1) + pi / 6),
      TRUE ~ 1.04 * sin(2 * pi * (ring_index - 1) / pmax(n_outer, 1) + pi / 18)
    ),
    centre_x = centre_x - 0.03,
    module_radius = rescale(sqrt(node_n), to = c(0.100, 0.230)),
    module_radius = if_else(layout_rank == 1, 0.37, module_radius),
    module_class = case_when(core_n > 0 ~ "Core-bearing module", recurrent_n > 0 ~ "Recurrent module", TRUE ~ "Network context"),
    module_code = case_when(
      module_key == "Other" ~ "Other",
      str_detect(coalesce(module_label, module_key), "H_NET_") ~ str_extract(coalesce(module_label, module_key), "H_NET_[0-9]+"),
      str_detect(coalesce(module_label, ""), "^H[0-9]+") ~ str_extract(module_label, "^H[0-9]+"),
      TRUE ~ as.character(module_key)
    ),
    module_display = program_label(paste(module_key, module_label, module_biological_label))
  )

main_nodes <- main_nodes %>%
  left_join(module_summary %>% select(module_key, centre_x, centre_y, module_radius, module_class), by = "module_key") %>%
  group_by(module_key) %>%
  arrange(desc(scaffold_hub), desc(strict_member), desc(relaxed_member), desc(strength), .by_group = TRUE) %>%
  mutate(
    within_rank = row_number(), within_n = n(),
    hub_rank = cumsum(scaffold_hub), hub_n = sum(scaffold_hub),
    base_angle = within_rank * pi * (3 - sqrt(5)),
    angle = if_else(scaffold_hub & hub_n > 1, 2 * pi * (hub_rank - 1) / hub_n - pi / 2, base_angle),
    base_radial = module_radius * 0.95 * sqrt((within_rank - 0.35) / within_n),
    radial = if_else(scaffold_hub & hub_n > 1, pmax(base_radial, module_radius * 0.22), base_radial),
    x_plot = centre_x + radial * cos(angle),
    y_plot = centre_y + radial * sin(angle)
  ) %>%
  ungroup()

circle_angles <- seq(0, 2 * pi, length.out = 100)
module_halos <- bind_rows(lapply(seq_len(nrow(module_summary)), function(i) {
  z <- module_summary[i, ]
  tibble(
    module_key = z$module_key, module_class = z$module_class, angle = circle_angles,
    x = z$centre_x + z$module_radius * cos(circle_angles),
    y = z$centre_y + z$module_radius * sin(circle_angles)
  )
}))

main_edges <- edges %>%
  filter(source %in% main_ids, target %in% main_ids) %>%
  transmute(
    source, target, weight,
    x_from = main_nodes$x_plot[match(source, main_nodes$node_id)],
    y_from = main_nodes$y_plot[match(source, main_nodes$node_id)],
    x_to = main_nodes$x_plot[match(target, main_nodes$node_id)],
    y_to = main_nodes$y_plot[match(target, main_nodes$node_id)]
  )
module_lookup <- main_nodes %>% distinct(node_id, module_key)
module_edges <- edges %>%
  filter(source %in% main_ids, target %in% main_ids) %>%
  left_join(module_lookup %>% rename(source = node_id, source_module = module_key), by = "source") %>%
  left_join(module_lookup %>% rename(target = node_id, target_module = module_key), by = "target") %>%
  filter(source_module != target_module) %>%
  mutate(a = pmin(source_module, target_module), b = pmax(source_module, target_module)) %>%
  count(a, b, wt = weight, name = "edge_weight") %>%
  slice_max(edge_weight, n = 12, with_ties = FALSE) %>%
  left_join(module_summary %>% select(a = module_key, x_from = centre_x, y_from = centre_y), by = "a") %>%
  left_join(module_summary %>% select(b = module_key, x_to = centre_x, y_to = centre_y), by = "b")
core_ids_main <- main_nodes$node_id[main_nodes$strict_member]
core_edges <- main_edges %>%
  filter(source %in% core_ids_main, target %in% core_ids_main)
hub_edges <- main_edges %>%
  filter(source %in% scaffold_hubs | target %in% scaffold_hubs) %>%
  mutate(hub_id = if_else(source %in% scaffold_hubs, source, target)) %>%
  group_by(hub_id) %>%
  slice_max(weight, n = 5, with_ties = FALSE) %>%
  ungroup()
hub_labels <- main_nodes %>% filter(scaffold_hub)
if (
  nrow(hub_labels) != length(scaffold_hubs) ||
  !setequal(hub_labels$node_id, scaffold_hubs) ||
  any(is.na(hub_labels$display_label) | !nzchar(hub_labels$display_label)) ||
  anyDuplicated(hub_labels$node_id)
) {
  stop("Every highlighted host scaffold hub must map to exactly one non-empty label.")
}
module_labels <- module_summary %>%
  filter(core_n > 0 | recurrent_n > 0, module_display != "Other") %>%
  arrange(desc(core_n), desc(recurrent_n), desc(node_n)) %>%
  slice_head(n = 6) %>%
  mutate(
    module_display = recode(module_display, `RNA processing` = "RNA\nprocessing"),
    label_x = centre_x,
    label_y = centre_y + module_radius + 0.055
  )
summary_label <- paste0(
  comma(nrow(main_nodes)), " connected proteins  |  ", comma(nrow(main_edges)), " edges\n",
  comma(sum(main_nodes$relaxed_member)), " recurrent  |  ", comma(sum(main_nodes$strict_member)), " connected core"
)

panel <- ggplot() +
  geom_segment(data = module_edges, aes(x_from, y_from, xend = x_to, yend = y_to, linewidth = edge_weight), colour = "#DCE3E7", alpha = 0.62, lineend = "round") +
  geom_polygon(data = module_halos, aes(x, y, group = module_key, fill = module_class), colour = "#DCE3E7", linewidth = 0.24, alpha = 0.50) +
  geom_segment(data = core_edges, aes(x_from, y_from, xend = x_to, yend = y_to), colour = blue_dark, alpha = 0.32, linewidth = 0.30, lineend = "round") +
  geom_segment(data = hub_edges, aes(x_from, y_from, xend = x_to, yend = y_to), colour = red_dark, alpha = 0.40, linewidth = 0.36, lineend = "round") +
  geom_point(data = filter(main_nodes, network_role == "Network context"), aes(x_plot, y_plot), colour = "#D0DEE5", size = 0.14, alpha = 0.72) +
  geom_point(data = filter(main_nodes, network_role == "Recurrent scaffold"), aes(x_plot, y_plot, size = pmax(strength, 0.1)), shape = 21, fill = blue_pale, colour = blue_dark, stroke = 0.24, alpha = 0.90) +
  geom_point(data = filter(main_nodes, network_role == "Core"), aes(x_plot, y_plot, size = pmax(strength, 0.1)), shape = 21, fill = red, colour = red_dark, stroke = 0.48) +
  geom_point(data = filter(main_nodes, scaffold_hub), aes(x_plot, y_plot), shape = 21, fill = alpha(red_pale, 0.72), colour = alpha(red_dark, 0.64), size = 6.8, stroke = 0.65) +
  geom_point(data = filter(main_nodes, scaffold_hub), aes(x_plot, y_plot, fill = strict_member), shape = 21, colour = ink, size = 1.75, stroke = 0.48) +
  geom_text(
    data = module_labels,
    aes(label_x, label_y, label = module_display),
    colour = ink, size = 1.68, lineheight = 0.86, fontface = "bold", hjust = 0.5, vjust = 0
  ) +
  geom_label_repel(data = hub_labels, aes(x_plot, y_plot, label = display_label), colour = red_dark, fill = alpha("white", 0.92), linewidth = 0.18, label.r = unit(1.2, "mm"), size = 2.20, fontface = "bold", seed = 20260722, max.overlaps = Inf, max.time = 5, max.iter = 20000, min.segment.length = 0, box.padding = 0.32, point.padding = 0.28, force = 1.20, force_pull = 1.80, segment.color = red_dark, segment.size = 0.24) +
  scale_fill_manual(values = c(`Core-bearing module` = red_pale, `Recurrent module` = blue_pale, `Network context` = "#F8F9FA", `FALSE` = blue, `TRUE` = red), guide = "none") +
  scale_linewidth_continuous(range = c(0.45, 1.45), guide = "none") +
  scale_size_continuous(range = c(0.65, 3.5), guide = "none") +
  annotate("text", x = 0, y = -1.40, label = summary_label, hjust = 0.5, vjust = 1, size = 2.05, lineheight = 0.98, colour = ink) +
  coord_equal(xlim = c(-1.36, 1.36), ylim = c(-1.52, 1.36), clip = "off") +
  theme_network +
  theme(plot.margin = margin(4, 4, 2, 4))

save_panel(panel, "Figure_1_panel_c_host_network", 220, 190)
cat(
  "Rendered Figure_1_panel_c_host_network with ", nrow(hub_labels),
  " labeled scaffold hubs: ",
  paste(hub_labels %>% arrange(desc(strength), desc(degree)) %>% pull(display_label), collapse = ", "),
  ".\n", sep = ""
)
