#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(grid)
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
  normalizePath(file.path("01_scripts", "03_panels", "02_figure_2", "03_panel_c_pathogen_network.R"), winslash = "/", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
p <- function(...) file.path(root, ...)

data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
if (nzchar(data_override) && !dir.exists(data_override)) stop("FIGURE_DATA_DIR does not exist: ", data_override)
data_candidates <- c(data_override, p("04_plot_data"))
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
data_dir <- data_candidates[[1]]
read_fig <- function(file) read_tsv(
  file.path(data_dir, "02_figure_2", file),
  show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN")
)

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "02_figure_2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(x, stem, width_mm, height_mm) {
  if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(x))
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), x,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
  ggsave(file.path(output_dir, paste0(stem, ".png")), x,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_png, bg = "white", limitsize = FALSE
  )
}

ink <- "#30363B"
blue_dark <- "#B5C8D1"
blue <- "#D2DEE4"
blue_pale <- "#F1F5F7"
red_dark <- "#CDB0AD"
red <- "#E4D0CD"
red_pale <- "#FAF1F0"
pin1_id <- "5874.Q4UG71"
empty_xy <- tibble(node_id = character(), x_plot = numeric(), y_plot = numeric())

theme_network <- theme_void(base_size = 8, base_family = "sans") +
  theme(
    legend.position = "none",
    plot.margin = margin(6, 8, 8, 8),
    plot.background = element_rect(fill = "white", colour = NA)
  )

as_flag <- function(x) {
  if (is.null(x)) return(FALSE)
  if (is.logical(x)) replace(x, is.na(x), FALSE) else str_to_lower(coalesce(as.character(x), "")) %in% c("true", "t", "1", "yes", "y")
}

sanitize_label <- function(x) {
  z <- as.character(coalesce(x, ""))
  z <- str_replace_all(z, "TA_15795\\?", "TaPin1")
  z <- str_replace_all(z, "TA_15795", "TaPin1")
  z <- str_replace_all(z, "\uFF08", "(")
  z <- str_replace_all(z, "\uFF09", ")")
  z <- str_replace_all(z, "\u2014", "-")
  z <- str_replace_all(z, "\u2013", "-")
  z <- str_replace_all(z, "p104=TA08425[^;]*", "p104")
  z <- str_replace_all(z, "TA08425", "p104")
  z <- str_replace_all(z, "\\(正文确认\\)", "")
  z <- str_replace_all(z, "\\(文中作 TaPin1\\)", "")
  z <- str_replace_all(z, "\\(无 TA 号记录于台账\\)", "")
  z <- str_replace_all(z, "\\(2015 论文[^)]*\\)", "")
  z <- str_replace_all(z, "\\(基因号未给出\\)", "")
  z <- str_replace_all(z, "\\(aa66-413[^)]*\\)", "")
  z <- str_replace_all(z, "\\(XP_952499\\.1\\)", "")
  z <- str_replace_all(z, "\\(XM_950452\\.1\\)", "")
  z <- str_replace_all(z, "（[^）]*）", "")
  str_squish(z)
}

pretty_pathogen_label <- function(display_label, preferred_name = NA_character_, annotation = NA_character_, node_id = NA_character_) {
  dl <- sanitize_label(display_label)
  pn <- sanitize_label(preferred_name)
  ann <- coalesce(as.character(annotation), "")
  nid <- coalesce(as.character(node_id), "")
  dplyr::case_when(
    nid == pin1_id | str_detect(paste(dl, pn), regex("(^|\\b)(Ta)?Pin1(\\b|$)", ignore_case = TRUE)) ~ "TaPin1",
    nid == "5874.Q4U9M9" | str_detect(paste(dl, pn), regex("(^|\\b)p104(\\b|$)", ignore_case = TRUE)) ~ "p104",
    nid == "5874.Q8WPH8" | str_detect(dl, regex("^TaSP$", ignore_case = TRUE)) ~ "TaSP",
    nid == "5874.P25781" | str_detect(dl, regex("^TaCP$", ignore_case = TRUE)) ~ "TaCP",
    str_detect(dl, regex("^TaCyp1$", ignore_case = TRUE)) ~ "TaCyp1",
    str_detect(dl, regex("^TaPHB1$", ignore_case = TRUE)) ~ "TaPHB1",
    str_detect(dl, regex("^MCM7$", ignore_case = TRUE)) | str_detect(ann, "licensing factor MCM7") ~ "MCM7",
    str_detect(dl, regex("^MT-CYB$", ignore_case = TRUE)) | str_detect(ann, regex("^Cytochrome b\\b", ignore_case = TRUE)) ~ "MT-CYB",
    str_detect(ann, regex("Signal recognition particle 54", ignore_case = TRUE)) ~ "SRP54",
    str_detect(ann, regex("Protein transport protein SEC23", ignore_case = TRUE)) ~ "SEC23",
    str_detect(ann, regex("AP-2 complex subunit alpha", ignore_case = TRUE)) ~ "AP2A",
    str_detect(ann, regex("Coatomer subunit beta", ignore_case = TRUE)) ~ "COPB",
    TRUE ~ dl
  )
}

is_ta_locus <- function(x) str_detect(coalesce(x, ""), "^TA[0-9]{4,}$")
is_ribosomal_text <- function(label, annotation) {
  str_detect(
    paste(coalesce(label, ""), coalesce(annotation, "")),
    regex("ribosom|\\bRPL\\d|\\bRPS\\d|\\b40S\\b|\\b60S\\b", ignore_case = TRUE)
  )
}

nodes <- read_fig("05_panel_g_pathogen_network_nodes.tsv")
edges <- read_fig("06_panel_g_pathogen_network_edges.tsv")
if (!"display_label" %in% names(nodes)) {
  label_src <- intersect(c("preferred_name", "canonical_gene", "gene_symbol"), names(nodes))
  nodes$display_label <- if (length(label_src)) dplyr::coalesce(!!!nodes[label_src], nodes$node_id) else nodes$node_id
}
if (!"preferred_name" %in% names(nodes)) nodes$preferred_name <- NA_character_
if (!"annotation" %in% names(nodes)) nodes$annotation <- NA_character_
if (!"weight" %in% names(edges)) edges$weight <- 1
if (!"network_role" %in% names(nodes)) nodes$network_role <- NA_character_
if (!"local_core" %in% names(nodes)) nodes$local_core <- FALSE
if ("scaffold_tier" %in% names(nodes)) {
  nodes <- nodes %>% mutate(
    network_role = case_when(
      network_role %in% c("Core", "Accessible scaffold", "Network context") ~ network_role,
      scaffold_tier %in% c("strict", "recurrent_core", "Core") ~ "Core",
      scaffold_tier %in% c("relaxed", "accessible_scaffold", "Accessible scaffold") ~ "Accessible scaffold",
      TRUE ~ coalesce(na_if(as.character(network_role), ""), "Network context")
    )
  )
}
nodes <- nodes %>%
  mutate(
    degree = as.numeric(degree),
    strength = as.numeric(strength),
    local_core = as_flag(local_core),
    is_pin1 = node_id == pin1_id | str_detect(coalesce(display_label, ""), regex("(^|\\b)(Ta)?Pin1(\\b|$)", ignore_case = TRUE)),
    display_label = pretty_pathogen_label(display_label, preferred_name, annotation, node_id),
    network_role = case_when(
      is_pin1 ~ "Network context",
      local_core ~ "Core",
      network_role %in% c("Core", "strict", "recurrent_core") ~ "Core",
      network_role %in% c("Accessible scaffold", "relaxed", "accessible_scaffold") ~ "Accessible scaffold",
      TRUE ~ coalesce(na_if(as.character(network_role), ""), "Network context")
    ),
    local_core = network_role == "Core",
    is_ribosomal = is_ribosomal_text(display_label, annotation),
    is_named = !is_ta_locus(display_label),
    is_force_effector = display_label %in% c("p104", "TaSP", "TaCP")
  )
if (any(nodes$network_role[nodes$is_pin1] == "Core") || any(nodes$local_core[nodes$is_pin1])) {
  stop("Pin1 (5874.Q4UG71) is never local-core and must not be drawn as Core.")
}
if (sum(nodes$network_role == "Core") != 19L) {
  stop("Fig2 Core must be the 19 local-core nodes; found ", sum(nodes$network_role == "Core"), ".")
}

graph_all <- graph_from_data_frame(
  edges %>% select(source, target, weight),
  directed = FALSE,
  vertices = nodes %>% transmute(name = node_id)
)
component_info <- components(graph_all)
main_component <- which.max(component_info$csize)
main_ids <- names(component_info$membership)[component_info$membership == main_component]

nodes <- nodes %>%
  mutate(
    component = as.integer(component_info$membership[node_id]),
    component_size = component_info$csize[component],
    in_main_component = component == main_component
  )

# Force-directed cloud of the main component. Quantile scaling keeps a few long
# tails from collapsing the biologically useful centre.
graph_main <- induced_subgraph(graph_all, vids = main_ids)
set.seed(20260713)
main_xy <- layout_with_fr(
  graph_main,
  weights = sqrt(pmax(E(graph_main)$weight, 0.05)),
  niter = 5000,
  grid = "nogrid"
)
main_xy <- sweep(main_xy, 2, apply(main_xy, 2, median), FUN = "-")
main_radius <- sqrt(rowSums(main_xy^2))
main_xy <- main_xy * (main_radius^0.58 / pmax(main_radius, 1e-9))
# Uniform radial scale so the force cloud stays circular, not a filled rectangle.
r <- sqrt(rowSums(main_xy^2))
r_cap <- unname(quantile(r, 0.985, na.rm = TRUE, names = FALSE))
main_xy <- main_xy * (1.08 / pmax(r_cap, 1e-9))
r <- sqrt(rowSums(main_xy^2))
too_far <- r > 1.12
if (any(too_far)) {
  main_xy[too_far, ] <- main_xy[too_far, , drop = FALSE] * (1.12 / r[too_far])
}
main_layout <- tibble(
  node_id = V(graph_main)$name,
  x_plot = main_xy[, 1],
  y_plot = main_xy[, 2]
)

peripheral_singletons <- nodes %>%
  filter(!in_main_component, component_size == 1, network_role != "Core") %>%
  arrange(node_id)
set.seed(20260714)
cloud_spec <- tibble::tribble(
  ~cloud, ~cx,  ~cy,   ~sx,  ~sy,  ~prob,
  1L,     -1.04,  0.55, 0.105, 0.145, 0.19,
  2L,     -1.07, -0.37, 0.095, 0.155, 0.20,
  3L,     -0.08, -0.92, 0.220, 0.060, 0.20,
  4L,      1.04, -0.38, 0.095, 0.150, 0.20,
  5L,      1.02,  0.55, 0.105, 0.145, 0.21
)
if (nrow(peripheral_singletons)) {
  cloud_id <- sample(cloud_spec$cloud, nrow(peripheral_singletons), replace = TRUE, prob = cloud_spec$prob)
  cloud_row <- match(cloud_id, cloud_spec$cloud)
  peripheral_layout <- peripheral_singletons %>%
    transmute(
      node_id,
      x_plot = cloud_spec$cx[cloud_row] + rnorm(n(), sd = cloud_spec$sx[cloud_row]),
      y_plot = cloud_spec$cy[cloud_row] + rnorm(n(), sd = cloud_spec$sy[cloud_row])
    )
} else {
  peripheral_layout <- empty_xy
}

minor_nodes <- nodes %>% filter(!in_main_component, component_size > 1)
minor_components <- sort(unique(minor_nodes$component))
minor_centres <- tibble(
  x = c(-1.00, -1.04, -0.86, 0.84, 1.02, 1.00, -0.02),
  y = c( 0.56, -0.34, -0.68, -0.68, -0.34, 0.56, -0.90)
)
if (length(minor_components)) {
  minor_layout <- bind_rows(lapply(seq_along(minor_components), function(i) {
    component_nodes <- minor_nodes %>%
      filter(component == minor_components[[i]]) %>%
      arrange(node_id)
    local_graph <- induced_subgraph(graph_all, vids = component_nodes$node_id)
    local_xy <- if (vcount(local_graph) == 2) {
      matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
    } else {
      layout_in_circle(local_graph)
    }
    local_xy <- local_xy[match(component_nodes$node_id, V(local_graph)$name), , drop = FALSE]
    centre_x <- minor_centres$x[[1 + (i - 1) %% nrow(minor_centres)]]
    centre_y <- minor_centres$y[[1 + (i - 1) %% nrow(minor_centres)]]
    tibble(
      node_id = component_nodes$node_id,
      x_plot = centre_x + 0.052 * local_xy[, 1],
      y_plot = centre_y + 0.052 * local_xy[, 2]
    )
  }))
} else {
  minor_layout <- empty_xy
}

isolated_core <- nodes %>%
  filter(!in_main_component, component_size == 1, network_role == "Core") %>%
  arrange(display_label)
if (nrow(isolated_core)) {
  core_theta <- seq(-1.23, -0.42, length.out = nrow(isolated_core))
  isolated_layout <- isolated_core %>%
    transmute(
      node_id,
      x_plot = 1.13 * cos(core_theta) + 0.04,
      y_plot = 0.94 * sin(core_theta) - 0.02
    )
} else {
  isolated_layout <- empty_xy
}

plot_nodes <- nodes %>%
  left_join(bind_rows(main_layout, peripheral_layout, minor_layout, isolated_layout), by = "node_id")
if (nrow(plot_nodes) != nrow(nodes) || anyNA(plot_nodes$x_plot) || anyNA(plot_nodes$y_plot)) {
  stop("The pathogen network layout must contain every protein exactly once.")
}

# Named effectors first; unnamed TA locus / ribosomal-like hubs are not allowed
# to occupy the label set when those effectors are present.
hub_rank <- plot_nodes %>%
  filter(degree > 0, !is_ribosomal) %>%
  mutate(
    hub_priority = case_when(
      is_force_effector ~ 1L,
      is_pin1 ~ 2L,
      network_role %in% c("Core", "Accessible scaffold") & is_named ~ 3L,
      network_role %in% c("Core", "Accessible scaffold") & !is_ta_locus(display_label) ~ 4L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(hub_priority)) %>%
  arrange(hub_priority, desc(network_role == "Core"), desc(strength), desc(degree), node_id)
n_locked <- sum(hub_rank$hub_priority <= 3L)
hub_ids <- hub_rank %>%
  slice_head(n = max(8L, min(12L, max(n_locked, 10L)))) %>%
  pull(node_id)
if (length(hub_ids) > 12L) hub_ids <- hub_ids[seq_len(12L)]
missing_force <- setdiff(
  plot_nodes$node_id[plot_nodes$is_force_effector],
  hub_ids
)
if (length(missing_force)) hub_ids <- unique(c(missing_force, hub_ids))
if (any(plot_nodes$is_pin1) && !any(hub_ids %in% plot_nodes$node_id[plot_nodes$is_pin1])) {
  hub_ids <- unique(c(plot_nodes$node_id[plot_nodes$is_pin1], hub_ids))
}
if (length(hub_ids) > 12L) {
  keep <- unique(c(
    plot_nodes$node_id[plot_nodes$is_force_effector | plot_nodes$is_pin1],
    hub_ids
  ))
  hub_ids <- keep[seq_len(min(length(keep), 12L))]
}

plot_nodes <- plot_nodes %>%
  mutate(
    is_hub = node_id %in% hub_ids,
    node_size = case_when(
      network_role == "Core" ~ rescale(log1p(strength), to = c(1.15, 2.55)),
      network_role == "Accessible scaffold" ~ rescale(log1p(strength), to = c(0.72, 1.55)),
      TRUE ~ 0.28
    )
  )
if (any(plot_nodes$is_hub & plot_nodes$is_pin1 & plot_nodes$network_role == "Core")) {
  stop("Pin1 may be labeled as accessible/context but must not be labeled as core.")
}

# Keep FR coordinates; only separate labeled hubs that would otherwise stack.
hub_idx <- which(plot_nodes$is_hub)
if (length(hub_idx) > 1L) {
  xy <- cbind(plot_nodes$x_plot[hub_idx], plot_nodes$y_plot[hub_idx])
  min_d <- 0.18
  for (iter in seq_len(120)) {
    moved <- FALSE
    for (i in seq_len(nrow(xy) - 1L)) {
      for (j in (i + 1L):nrow(xy)) {
        dx <- xy[i, 1] - xy[j, 1]
        dy <- xy[i, 2] - xy[j, 2]
        dist <- sqrt(dx * dx + dy * dy)
        if (dist < min_d) {
          push <- (min_d - max(dist, 1e-6)) / 2
          ux <- dx / max(dist, 1e-6)
          uy <- dy / max(dist, 1e-6)
          xy[i, ] <- xy[i, ] + c(push * ux, push * uy)
          xy[j, ] <- xy[j, ] - c(push * ux, push * uy)
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  xy[, 1] <- pmin(pmax(xy[, 1], -1.05), 1.05)
  xy[, 2] <- pmin(pmax(xy[, 2], -0.95), 1.00)
  plot_nodes$x_plot[hub_idx] <- xy[, 1]
  plot_nodes$y_plot[hub_idx] <- xy[, 2]
}

edge_xy <- edges %>%
  transmute(
    source, target, weight,
    x_from = plot_nodes$x_plot[match(source, plot_nodes$node_id)],
    y_from = plot_nodes$y_plot[match(source, plot_nodes$node_id)],
    x_to = plot_nodes$x_plot[match(target, plot_nodes$node_id)],
    y_to = plot_nodes$y_plot[match(target, plot_nodes$node_id)]
  )
n_full_edges <- nrow(edge_xy)

# STRING is too dense to draw (86,106 edges). Keep FR node positions and draw
# only labeled-hub spokes plus the few core-core edges.
core_ids <- plot_nodes$node_id[plot_nodes$network_role == "Core"]
core_edges <- edge_xy %>%
  filter(source %in% core_ids, target %in% core_ids)
if (nrow(core_edges) > 40L) core_edges <- core_edges[0, , drop = FALSE]

hub_neighbour_rows <- bind_rows(
  edges %>% filter(source %in% hub_ids) %>% transmute(hub_id = source, neighbour_id = target, weight),
  edges %>% filter(target %in% hub_ids) %>% transmute(hub_id = target, neighbour_id = source, weight)
) %>%
  filter(hub_id != neighbour_id) %>%
  left_join(
    plot_nodes %>% transmute(
      neighbour_id = node_id,
      neighbour_role = network_role,
      neighbour_strength = strength
    ),
    by = "neighbour_id"
  ) %>%
  mutate(
    role_rank = case_when(
      neighbour_role == "Core" ~ 2L,
      neighbour_role == "Accessible scaffold" ~ 1L,
      TRUE ~ 0L
    )
  )
hub_spokes <- hub_neighbour_rows %>%
  group_by(hub_id) %>%
  arrange(desc(role_rank), desc(neighbour_strength), desc(weight), neighbour_id, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    source = pmin(hub_id, neighbour_id),
    target = pmax(hub_id, neighbour_id)
  ) %>%
  distinct(source, target, .keep_all = TRUE) %>%
  transmute(
    source, target, weight, hub_id,
    x_from = plot_nodes$x_plot[match(source, plot_nodes$node_id)],
    y_from = plot_nodes$y_plot[match(source, plot_nodes$node_id)],
    x_to = plot_nodes$x_plot[match(target, plot_nodes$node_id)],
    y_to = plot_nodes$y_plot[match(target, plot_nodes$node_id)]
  )

context_nodes_main <- plot_nodes %>% filter(network_role == "Network context", in_main_component, !is_hub)
context_nodes_peripheral <- plot_nodes %>% filter(network_role == "Network context", !in_main_component, !is_hub)
scaffold_nodes <- plot_nodes %>% filter(network_role == "Accessible scaffold")
core_nodes <- plot_nodes %>% filter(network_role == "Core")
hub_nodes <- plot_nodes %>% filter(is_hub)
hub_centre <- hub_nodes %>% summarise(x = mean(x_plot), y = mean(y_plot))
hub_labels <- hub_nodes %>%
  mutate(natural_angle = atan2(y_plot - hub_centre$y, x_plot - hub_centre$x)) %>%
  arrange(natural_angle) %>%
  mutate(
    label_angle = seq(-pi + 0.26, pi + 0.26, length.out = n() + 1)[seq_len(n())],
    label_x = hub_centre$x + 0.47 * cos(label_angle),
    label_y = hub_centre$y + 0.39 * sin(label_angle)
  )
if (any(is.na(hub_nodes$display_label) | !nzchar(hub_nodes$display_label)) || anyDuplicated(hub_nodes$node_id)) {
  stop("Every highlighted pathogen hub must map to exactly one non-empty label.")
}
if (any(str_detect(hub_nodes$display_label, "TA_15795\\?|（"))) {
  stop("Hub labels still contain TA_15795? or Chinese parentheticals.")
}

hub_halo <- hub_neighbour_rows %>%
  distinct(hub_id, neighbour_id) %>%
  count(hub_id, name = "n_one_hop") %>%
  right_join(hub_nodes %>% select(hub_id = node_id, x_plot, y_plot), by = "hub_id") %>%
  mutate(
    n_one_hop = coalesce(n_one_hop, 0L),
    halo_radius = rescale(log1p(n_one_hop), to = c(0.048, 0.070))
  ) %>%
  tidyr::crossing(theta = seq(0, 2 * pi, length.out = 80)) %>%
  mutate(
    halo_x = x_plot + halo_radius * cos(theta),
    halo_y = y_plot + halo_radius * sin(theta)
  )

n_scaffold <- sum(plot_nodes$network_role %in% c("Accessible scaffold", "Core"))
n_core <- sum(plot_nodes$network_role == "Core")
summary_label <- paste0(
  comma(nrow(plot_nodes)), " proteins  |  ", comma(n_full_edges), " edges  |  ",
  comma(n_scaffold), " scaffold  |  ", comma(n_core), " core"
)

panel <- ggplot() +
  geom_polygon(
    data = hub_halo,
    aes(halo_x, halo_y, group = hub_id),
    fill = red_pale, colour = red, alpha = 0.18, linewidth = 0.14
  ) +
  geom_segment(
    data = core_edges,
    aes(x_from, y_from, xend = x_to, yend = y_to),
    colour = red, alpha = 0.42, linewidth = 0.28, lineend = "round"
  ) +
  geom_segment(
    data = hub_spokes,
    aes(x_from, y_from, xend = x_to, yend = y_to),
    colour = red_dark, alpha = 0.48, linewidth = 0.34, lineend = "round"
  ) +
  geom_point(
    data = context_nodes_peripheral, aes(x_plot, y_plot),
    colour = blue_dark, size = 0.30, alpha = 0.58
  ) +
  geom_point(
    data = context_nodes_main, aes(x_plot, y_plot),
    colour = blue_dark, size = 0.48, alpha = 0.90
  ) +
  geom_point(
    data = scaffold_nodes, aes(x_plot, y_plot, size = node_size),
    shape = 21, fill = red, colour = red_dark, stroke = 0.36, alpha = 0.95
  ) +
  geom_point(
    data = core_nodes, aes(x_plot, y_plot, size = node_size),
    shape = 21, fill = red, colour = red_dark, stroke = 0.52, alpha = 0.96
  ) +
  geom_point(
    data = hub_nodes, aes(x_plot, y_plot),
    shape = 21, fill = alpha(red_pale, 0.72), colour = alpha(red_dark, 0.70),
    size = 5.4, stroke = 0.65
  ) +
  geom_point(
    data = hub_nodes, aes(x_plot, y_plot),
    shape = 21, fill = "white", colour = red_dark, size = 3.25, stroke = 0.55
  ) +
  geom_point(
    data = hub_nodes, aes(x_plot, y_plot, fill = network_role == "Core"),
    shape = 21, colour = ink, size = 1.52, stroke = 0.46
  ) +
  geom_segment(
    data = hub_labels,
    aes(x_plot, y_plot, xend = label_x, yend = label_y),
    colour = red_dark, alpha = 0.60, linewidth = 0.24, lineend = "round"
  ) +
  geom_label(
    data = hub_labels,
    aes(label_x, label_y, label = display_label),
    colour = red_dark, fill = alpha("white", 0.94), linewidth = 0.18,
    label.r = unit(1.1, "mm"), label.padding = unit(1.3, "mm"),
    size = 2.18, fontface = "bold"
  ) +
  scale_size_identity() +
  scale_fill_manual(values = c(`FALSE` = blue_dark, `TRUE` = red), guide = "none") +
  annotate(
    "label", x = 0, y = -1.22, label = summary_label,
    hjust = 0.5, vjust = 0, size = 2.30, linewidth = 0,
    fill = alpha("white", 0.90), colour = ink
  ) +
  coord_equal(xlim = c(-1.27, 1.27), ylim = c(-1.32, 1.10), clip = "off") +
  theme_network

save_panel(panel, "Figure_2_panel_c_pathogen_network", 220, 190)

hub_report <- hub_nodes %>%
  arrange(desc(is_force_effector), desc(network_role == "Core"), desc(strength)) %>%
  transmute(txt = paste0(display_label, " (", network_role, ")")) %>%
  pull(txt)
pin1_row <- plot_nodes %>% filter(is_pin1)
cat(
  "Rendered Figure_2_panel_c_pathogen_network\n",
  "labeled hubs (", length(hub_ids), "): ", paste(hub_report, collapse = ", "), "\n",
  "full graph: ", nrow(plot_nodes), " proteins / ", n_full_edges, " edges / ",
  n_scaffold, " scaffold / ", n_core, " core\n",
  "drawn edges: ", nrow(hub_spokes), " hub spokes + ", nrow(core_edges),
  " core-core (not the full ", n_full_edges, ")\n",
  "Pin1 ", pin1_id, ": network_role=", if (nrow(pin1_row)) pin1_row$network_role[[1]] else "absent",
  "; local_core=", if (nrow(pin1_row)) pin1_row$local_core[[1]] else NA,
  "; labeled=", pin1_id %in% hub_ids, "\n",
  sep = ""
)
