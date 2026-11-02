#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(patchwork)
  library(ragg)
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
  normalizePath(file.path("01_scripts", "03_panels", "05_figure_5", "04_panel_d_program_mechanism_atlas.R"),
                winslash = "/", mustWork = TRUE)
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
has_input <- vapply(data_candidates, function(x) {
  file.exists(file.path(x, "05_figure_5", "08_panel_d_core_mechanisms.tsv")) ||
    file.exists(file.path(x, "05_figure_5", "07_panel_d_program_coupling.tsv"))
}, logical(1))
data_dir <- data_candidates[has_input][1]
if (is.na(data_dir)) stop("Run 01_scripts/02_build_figure_data.R first.")
read_fig <- function(file) {
  fp <- file.path(data_dir, "05_figure_5", file)
  if (!file.exists(fp)) return(tibble())
  read_tsv(fp, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
}

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "05_figure_5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(plot, stem, width_mm = 205, height_mm = 118) {
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width_mm, height = height_mm,
         units = "mm", dpi = 600, device = ragg::agg_png, bg = "white", limitsize = FALSE)
}

mechanism_palette <- mechanism_cols
mechanism_order <- names(mechanism_palette)

axes <- bind_rows(
  read_fig("08_panel_d_core_mechanisms.tsv"),
  read_fig("panel_d_retained_axes.tsv")
)
axis_alt <- file.path(data_dir, "04_figure_4", "07_panel_d_axis_descriptor.tsv")
if (file.exists(axis_alt)) {
  axes <- bind_rows(axes, read_tsv(axis_alt, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN")))
}
if (!nrow(axes)) stop("Missing core/retained mechanism table.")
for (col in c(
  "axis_id", "biological_axis_id", "interaction_name", "pathogen_gene", "host_gene",
  "parasite_factor_name", "host_target_name", "final_tier", "mechanism_class",
  "axis_family", "aggregate_robustness_score", "host_module", "pathogen_module"
)) {
  if (!col %in% names(axes)) axes[[col]] <- NA
}
axes <- axes %>%
  mutate(
    axis_id_key = axis_key(coalesce(biological_axis_id, axis_id), interaction_name, pathogen_gene, host_gene),
    axis_tier = axis_tier(axis_id_key, final_tier),
    pathogen_factor = endpoint_display(sanitize_label(coalesce(pathogen_gene, parasite_factor_name))),
    host_target = endpoint_display(sanitize_label(coalesce(host_gene, host_target_name))),
    mechanism = coalesce(
      mechanism_from_family(axis_family),
      unname(c(
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
        `TaCP--CRBN` = "Protease / cullin"
      )[axis_id_key])
    ),
    interaction = recode(
      axis_id_key,
      `p104--MAPRE1` = "p104-MAPRE1",
      `TaPin1--PKM` = "TaPin1-PKM",
      `TA20980--CD2AP` = "TA20980-CD2AP",
      `TaSP--tubulin` = "TaSP-tubulin",
      `TaPin1--FBXW7` = "TaPin1-FBXW7",
      `TaPHB1--RUVBL1` = "TaPHB1-RUVBL1",
      `TA11965--TFG` = "TA11965-TFG",
      `TaCyp1--MED21` = "TaCyp1-MED21",
      `TA05560--RBM39` = "TA05560-RBM39",
      `TaCP--PPP4C` = "TaCP-PPP4C",
      `TaCP--CRBN` = "TaCP-CRBN",
      .default = sanitize_label(coalesce(interaction_name, axis_id_key))
    )
  ) %>%
  filter(axis_tier %in% c("Core", "Extended"), !axis_id_key %in% context_axis_ids) %>%
  distinct(axis_id_key, .keep_all = TRUE) %>%
  mutate(mechanism = factor(mechanism, levels = mechanism_order)) %>%
  arrange(mechanism, desc(axis_tier == "Core"), desc(num_col(aggregate_robustness_score)), axis_id_key)
if (any(is.na(axes$mechanism))) stop("Unmapped mechanism class in retained axes.")
if (any(str_detect(paste(axes$interaction, axes$pathogen_factor, axes$host_target), "TA_15795"))) {
  stop("Forbidden label TA_15795? reached the mechanism atlas.")
}
locked_ids <- c(core_axis_ids, extended_axis_ids)
if (nrow(axes) != 11L || !setequal(axes$axis_id_key, locked_ids)) {
  stop(sprintf(
    "Expected 11 retained axes, got %d: %s",
    nrow(axes), paste(axes$axis_id_key, collapse = ", ")
  ))
}

pathogen_program_atlas <- function(z) {
  lab <- pretty_module(z, wrap = FALSE)
  if_else(
    lab %in% c("Folding / processing", "Translation / folding"),
    "Folding /\nprocessing",
    "Surface /\nexport"
  )
}
coupling <- axes %>%
  mutate(
    host_program = str_replace(
      pretty_module(host_module, wrap = TRUE),
      "proteo-\nstasis", "proteostasis"
    ),
    pathogen_program = pathogen_program_atlas(pathogen_module)
  ) %>%
  group_by(host_program, pathogen_program) %>%
  summarise(
    axis_n = n(),
    core_axis_n = sum(axis_tier == "Core"),
    .groups = "drop"
  ) %>%
  mutate(core_fraction = if_else(axis_n > 0, core_axis_n / axis_n, 0))
if (sum(coupling$axis_n) != 11L) {
  stop("Heatmap must cover the 11 retained axes, got ", sum(coupling$axis_n))
}
host_order <- coupling %>%
  group_by(host_program) %>%
  summarise(core_n = sum(core_axis_n), axis_n = sum(axis_n), .groups = "drop") %>%
  arrange(core_n, axis_n, host_program) %>%
  pull(host_program)
pathogen_order <- coupling %>%
  group_by(pathogen_program) %>%
  summarise(core_n = sum(core_axis_n), axis_n = sum(axis_n), .groups = "drop") %>%
  arrange(core_n, axis_n, pathogen_program) %>%
  pull(pathogen_program)
coupling <- coupling %>% mutate(
  host_program = factor(host_program, levels = host_order),
  pathogen_program = factor(pathogen_program, levels = pathogen_order)
)
heat_max <- max(coupling$axis_n)
if (heat_max > 3L) stop("Heatmap size scale is locked to n=1..3, got max ", heat_max)
heat_breaks <- 1:3
left_plot <- ggplot(coupling, aes(pathogen_program, host_program)) +
  geom_point(aes(size = axis_n, fill = core_fraction), shape = 21, colour = "white", stroke = 0.72) +
  scale_size_continuous(
    range = c(5.2, 7.4), breaks = heat_breaks, limits = c(1, 3),
    labels = as.character(heat_breaks), name = "Axes per pair",
    guide = guide_legend(
      order = 1, nrow = 1, byrow = TRUE, title.position = "top", title.hjust = 0.5,
      override.aes = list(
        size = scales::rescale(heat_breaks, to = c(2.2, 4.4)),
        shape = 21, fill = NA, colour = "#927D96", stroke = 0.60
      )
    )
  ) +
  scale_fill_gradientn(
    colours = c("#F6F2F6", "#E0D6DE", "#A895AA"),
    limits = c(0, 1), breaks = c(0, 0.5, 1), name = "Core fraction",
    guide = guide_colorbar(
      order = 2, direction = "horizontal", nbin = 60, title.position = "top", title.hjust = 0.5,
      barwidth = grid::unit(13, "mm"), barheight = grid::unit(1.5, "mm")
    )
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Pathogen program", y = "Host program") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 7.0, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.34),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.22, colour = ink),
    axis.text.x = element_text(colour = "#AE7772", size = 4.35, lineheight = 0.80, margin = margin(t = 1.8)),
    axis.text.y = element_text(colour = "#8099AA", size = 5.05, lineheight = 0.84, margin = margin(r = 1.8)),
    axis.title = element_text(colour = ink, size = 5.7),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.justification = "center",
    legend.title = element_text(colour = ink, size = 5.10, margin = margin(b = 1.2, r = 3.0)),
    legend.text = element_text(colour = ink, size = 4.80),
    legend.key = element_blank(),
    legend.key.width = grid::unit(3.2, "mm"),
    legend.key.height = grid::unit(2.4, "mm"),
    legend.spacing.x = grid::unit(3.6, "mm"),
    legend.spacing.y = grid::unit(1.2, "mm"),
    legend.box.spacing = grid::unit(2.2, "mm"),
    legend.margin = margin(t = 2.8, r = 4, b = 1, l = 4),
    plot.margin = margin(4, 2, 2, 4),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

constellation_slots <- bind_rows(
  tibble(
    mechanism = "Microtubule docking", slot_id = 1:2,
    cx = c(0.140, 0.360), cy = c(0.600, 0.600), angle = c(9, -9)
  ),
  tibble(
    mechanism = "Metabolic / proteostasis", slot_id = 1:3,
    cx = c(0.650, 0.860, 0.755), cy = c(0.710, 0.710, 0.510), angle = c(9, -9, 7)
  ),
  tibble(
    mechanism = "Trafficking / signaling", slot_id = 1:4,
    cx = c(0.140, 0.360, 0.140, 0.360), cy = c(0.205, 0.205, 0.085, 0.085),
    angle = c(9, -9, -9, 9)
  ),
  tibble(
    mechanism = "Protease / cullin", slot_id = 1:2,
    cx = c(0.650, 0.860), cy = c(0.150, 0.150), angle = c(9, -9)
  )
)

in_node_label <- function(z) {
  z <- as.character(z)
  if_else(str_detect(z, "^TA[0-9]{5}$"), str_replace(z, "^(TA)([0-9]+)$", "\\1\n\\2"), z)
}
fit_endpoint_size <- function(z) {
  lab <- as.character(z)
  n <- nchar(str_replace(lab, "\n", ""))
  dplyr::case_when(
    str_detect(lab, "\n") ~ 1.12,
    n >= 6L ~ 1.22,
    TRUE ~ 1.34
  )
}

axes_plot <- axes %>%
  group_by(mechanism) %>%
  mutate(slot_id = row_number()) %>%
  ungroup() %>%
  mutate(mechanism = as.character(mechanism)) %>%
  left_join(constellation_slots, by = c("mechanism", "slot_id")) %>%
  mutate(
    half_length = 0.070,
    angle_rad = angle * pi / 180,
    x_host = cx - half_length * cos(angle_rad),
    x_pathogen = cx + half_length * cos(angle_rad),
    y_host = cy - half_length * sin(angle_rad) * 0.83,
    y_pathogen = cy + half_length * sin(angle_rad) * 0.83,
    label_x = cx - 0.010 * sin(angle_rad),
    label_y = cy + 0.050 * cos(angle_rad),
    mechanism_colour = if_else(
      axis_tier == "Core",
      unname(mechanism_palette[mechanism]),
      unname(mechanism_pale[mechanism])
    ),
    line_width = if_else(axis_tier == "Core", 1.20, 0.70),
    host_dot = in_node_label(host_target),
    pathogen_dot = in_node_label(pathogen_factor),
    host_text_size = fit_endpoint_size(host_dot),
    pathogen_text_size = fit_endpoint_size(pathogen_dot),
    host_node_size = if_else(str_detect(host_dot, "\n") | nchar(host_target) >= 6L, 6.6, 5.4),
    pathogen_node_size = if_else(str_detect(pathogen_dot, "\n") | nchar(pathogen_factor) >= 6L, 6.6, 5.4)
  )
if (any(!is.finite(unlist(select(axes_plot, cx, cy, x_host, x_pathogen, y_host, y_pathogen))))) {
  stop("Incomplete constellation coordinates for retained axes.")
}

n_by_mech <- axes %>% count(mechanism, name = "n")
mechanism_titles <- tibble(
  mechanism = mechanism_order,
  x = c(0.250, 0.750, 0.250, 0.750),
  y = c(0.850, 0.850, 0.340, 0.340),
  label = paste0(
    gsub(" / ", " /\n", mechanism_order, fixed = TRUE), " | n=",
    coalesce(n_by_mech$n[match(mechanism_order, as.character(n_by_mech$mechanism))], 0L)
  ),
  colour = ink
)

right_plot <- ggplot(axes_plot) +
  geom_text(
    data = mechanism_titles,
    aes(x = x, y = y, label = label, colour = colour),
    hjust = 0.5, vjust = 1, family = "sans", size = 2.02,
    fontface = "bold", lineheight = 0.92, show.legend = FALSE
  ) +
  geom_segment(
    aes(x = x_host, xend = x_pathogen, y = y_host, yend = y_pathogen, colour = mechanism_colour, linewidth = line_width),
    lineend = "round", show.legend = FALSE
  ) +
  geom_point(
    aes(x = x_host, y = y_host, size = host_node_size), shape = 21, stroke = 0.50,
    fill = host_fill, colour = "#8099AA", show.legend = FALSE
  ) +
  geom_point(
    aes(x = x_pathogen, y = y_pathogen, size = pathogen_node_size), shape = 21, stroke = 0.50,
    fill = pathogen_fill, colour = "#AE7772", show.legend = FALSE
  ) +
  geom_text(
    aes(x = x_host, y = y_host, label = host_dot, size = host_text_size),
    family = "sans", colour = "#8099AA", fontface = "bold", lineheight = 0.76,
    hjust = 0.5, vjust = 0.5, show.legend = FALSE
  ) +
  geom_text(
    aes(x = x_pathogen, y = y_pathogen, label = pathogen_dot, size = pathogen_text_size),
    family = "sans", colour = "#AE7772", fontface = "bold", lineheight = 0.76,
    hjust = 0.5, vjust = 0.5, show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, y = label_y, label = interaction, angle = angle),
    family = "sans", colour = ink, size = 1.55, fontface = "bold",
    lineheight = 0.80, vjust = 0.0, hjust = 0.5, show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_size_identity() +
  scale_linewidth_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.02, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 7, base_family = "sans") +
  theme(
    plot.margin = margin(2, 4, 4, 2),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

left_compact <- wrap_plots(
  plot_spacer(), left_plot, plot_spacer(),
  ncol = 1, heights = c(0.16, 0.83, 0.01)
)
combined <- wrap_plots(
  left_compact, plot_spacer(), right_plot,
  nrow = 1, widths = c(0.64, 0.018, 1.0), guides = "keep"
) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", unset = ""), "1")) {
  save_panel(combined, "Figure_5_panel_d_program_mechanism_atlas")
}
cat(sprintf(
  "Rendered Figure_5_panel_d_program_mechanism_atlas: left %d program pairs/%d axes; right %d retained axes in %d families (core %d, extended %d).\n",
  nrow(coupling), sum(coupling$axis_n), nrow(axes), length(mechanism_order),
  sum(axes$axis_tier == "Core"), sum(axes$axis_tier == "Extended")
))
