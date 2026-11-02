#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(grid)
})

source_files <- vapply(
  sys.frames(),
  function(x) if (is.null(x$ofile)) "" else as.character(x$ofile[[1]]),
  character(1)
)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (any(nzchar(source_files))) {
  normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/")
} else if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
} else {
  normalizePath(
    "01_scripts/03_panels/03_figure_3/03_panel_c_propagated_nodes.R",
    winslash = "/"
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

figure_data_override <- Sys.getenv("FIGURE_DATA_DIR", "")
data_candidates <- c(
  figure_data_override,
  p("04_plot_data")
)
data_dir <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)][1]
if (is.na(data_dir)) stop("Run stage 09 before rendering Figure 3c.")

input_file <- file.path(data_dir, "03_figure_3", "03_panel_b_rwr_landscape.tsv")
if (!file.exists(input_file)) stop("Missing Figure 3c source data: ", input_file)

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", "")
output_dir <- if (nzchar(output_override)) {
  output_override
} else {
  p("06_panel_output", "03_figure_3")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

nodes <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA", "NaN")
)
required <- c(
  "node_id", "layer", "display_label", "module_label", "degree",
  "strict_member", "relaxed_member", "direct_evidence",
  "opposite_signal", "recruited"
)
missing_cols <- setdiff(required, names(nodes))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

locked_pref <- c(
  "p104", "TaPin1", "TaCP", "TaSP", "MAPRE1", "PKM",
  "CD2AP", "TA20980", "MED21", "FBXW7"
)

is_leftover_program <- function(pretty, raw) {
  pretty <- coalesce(as.character(pretty), "")
  raw <- coalesce(as.character(raw), "")
  str_detect(pretty, "^(NET|FUN)[ _-]*[0-9A-Za-z]*$") |
    str_detect(pretty, "^(H_|P_|host_|parasite_|pathogen_)") |
    pretty %in% c("Unassigned", "unassigned") |
    (str_detect(raw, "^[HP]_NET_[0-9]+$") & str_detect(pretty, "^NET "))
}

is_labelable <- function(x) {
  z <- str_squish(coalesce(as.character(x), ""))
  nzchar(z) &
    !str_detect(z, regex("^ENS", ignore_case = TRUE)) &
    !str_detect(z, "^LOC[0-9]") &
    !str_detect(z, "[\u4e00-\u9fff]") &
    (z %in% locked_pref | !str_detect(z, "^TA[_-]?[0-9]+\\??$"))
}

nodes <- nodes %>%
  filter(flag_col(recruited), layer %in% c("host", "pathogen")) %>%
  mutate(
    layer = factor(layer, levels = c("host", "pathogen")),
    strict_member = flag_col(strict_member),
    relaxed_member = flag_col(relaxed_member),
    display_label = endpoint_display(sanitize_label(coalesce(display_label, node_id))),
    functional_module = pretty_module(coalesce(module_label, "Unassigned")),
    functional_module = if_else(
      is_leftover_program(functional_module, module_label),
      "Other programs",
      functional_module
    ),
    prior_scaffold = strict_member | relaxed_member,
    direct_evidence = pmin(pmax(num_col(direct_evidence, 0), 0), 1),
    opposite_signal = pmin(pmax(num_col(opposite_signal, 0), 0), 1),
    degree = num_col(degree, default = NA_real_)
  )

if (!nrow(nodes)) stop("No recruited nodes were found in panel_b_rwr_landscape.tsv")

finite_degree <- nodes$degree[is.finite(nodes$degree) & nodes$degree > 0]
sane_degree <- finite_degree[finite_degree <= 2000]
obs_degree_max <- if (length(sane_degree)) max(sane_degree) else 500
if (obs_degree_max <= 500 * 1.05) {
  degree_max <- 500
  degree_breaks <- c(1, 10, 100, 500)
} else if (obs_degree_max <= 1000) {
  degree_max <- 1000
  degree_breaks <- c(1, 10, 100, 1000)
} else {
  degree_max <- 1000
  degree_breaks <- c(1, 10, 100, 1000)
}

nodes <- nodes %>%
  mutate(
    degree_plot = pmin(pmax(coalesce(degree, 1), 1), degree_max),
    degree_log10 = log10(degree_plot)
  )

wrap_c_y <- function(prefix, module, n) {
  pre <- paste0(prefix, " · ")
  nlab <- paste0("(n = ", n, ")")
  head <- paste0(pre, module)
  if (grepl(" / ", module, fixed = TRUE)) {
    parts <- strsplit(module, " / ", fixed = TRUE)[[1]]
    return(paste0(pre, parts[[1]], " /\n", parts[[2]], " ", nlab))
  }
  if (grepl(" ", module, fixed = TRUE)) {
    return(paste0(pre, sub(" ", "\n", module, fixed = TRUE), " ", nlab))
  }
  if (nchar(head) > 20 && grepl("-", module, fixed = TRUE)) {
    return(paste0(pre, sub("-", "-\n", module, fixed = TRUE), " ", nlab))
  }
  paste0(head, "\n", nlab)
}

module_groups <- nodes %>%
  group_by(layer, functional_module) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(is_other = functional_module == "Other programs") %>%
  arrange(layer, is_other, desc(n), functional_module) %>%
  mutate(
    group_order = row_number(),
    y_group = rev(seq_len(n())),
    y_group = if_else(layer == "pathogen", y_group - 0.22, as.numeric(y_group)),
    band = rep(c("#F8F9FA", "#FFFFFF"), length.out = n()),
    prefix = if_else(layer == "host", "Host", "Pathogen"),
    module_display = vapply(
      seq_along(functional_module),
      function(i) wrap_c_y(prefix[[i]], functional_module[[i]], n[[i]]),
      character(1)
    )
  ) %>%
  select(-is_other, -prefix)

nodes <- nodes %>%
  left_join(
    module_groups %>% select(layer, functional_module, y_group),
    by = c("layer", "functional_module")
  ) %>%
  group_by(layer, functional_module) %>%
  arrange(node_id, .by_group = TRUE) %>%
  mutate(
    jitter_phase = (row_number() * ((sqrt(5) - 1) / 2)) %% 1,
    y_plot = y_group + (jitter_phase - 0.5) * 0.56
  ) %>%
  ungroup()

locked_nodes <- nodes %>% filter(display_label %in% locked_pref)
row_keys <- nodes %>%
  filter(is_labelable(display_label)) %>%
  group_by(layer, functional_module) %>%
  arrange(
    desc(display_label %in% locked_pref),
    desc(opposite_signal),
    desc(direct_evidence),
    desc(degree_plot),
    node_id,
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()
key_nodes <- bind_rows(locked_nodes, row_keys) %>%
  distinct(node_id, .keep_all = TRUE)

host_colour <- blue_dark
pathogen_colour <- red_dark
separator <- line
row_limits <- range(module_groups$y_group) + c(-0.48, 0.48)
layer_separator <- mean(c(
  min(module_groups$y_group[module_groups$layer == "host"]),
  max(module_groups$y_group[module_groups$layer == "pathogen"])
))
row_breaks <- module_groups$y_group
row_labels <- module_groups$module_display
direct_obs <- max(nodes$direct_evidence, na.rm = TRUE)
direct_hi <- if (is.finite(direct_obs) && direct_obs > 0.42) {
  max(0.52, ceiling(direct_obs * 20) / 20)
} else {
  0.42
}
direct_breaks <- if (direct_hi <= 0.42) c(0, 0.2, 0.4) else c(0, 0.25, 0.5)

theme_metric <- theme_classic(base_size = 7.0, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = frame, linewidth = 0.34),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.05, colour = ink, lineheight = 0.86, hjust = 1),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.30, colour = ink),
    axis.text.x = element_text(size = 6.2, colour = ink),
    axis.title.x = element_text(size = 6.7, colour = ink, margin = margin(t = 4)),
    legend.title = element_blank(),
    legend.text = element_text(size = 6.0),
    legend.key = element_blank(),
    legend.key.width = unit(4.4, "mm"),
    legend.key.height = unit(3.2, "mm"),
    legend.spacing.x = unit(3.2, "mm"),
    legend.box.spacing = unit(4.0, "mm"),
    plot.margin = margin(4, 5, 4, 6),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

add_module_bands <- function(plot) {
  plot +
    geom_rect(
      data = module_groups,
      aes(xmin = -Inf, xmax = Inf, ymin = y_group - 0.38, ymax = y_group + 0.38),
      inherit.aes = FALSE,
      fill = module_groups$band,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_hline(
      yintercept = layer_separator,
      colour = separator,
      linewidth = 0.45
    )
}

layer_scale <- scale_colour_manual(
  values = c(host = host_colour, pathogen = pathogen_colour),
  labels = c(host = "Host", pathogen = "Pathogen"),
  name = NULL
)

add_node_layers <- function(plot, label_keys = FALSE) {
  plot +
    geom_point(
      data = nodes,
      aes(colour = layer),
      size = 0.94,
      alpha = 0.82,
      stroke = 0,
      show.legend = TRUE
    ) +
    geom_point(
      data = filter(nodes, prior_scaffold),
      aes(fill = layer, shape = "Prior accessible scaffold"),
      colour = ink,
      size = 1.72,
      alpha = 0.96,
      stroke = 0.64,
      show.legend = TRUE
    ) +
    {
      if (label_keys) {
        geom_label_repel(
          data = key_nodes,
          aes(label = display_label, colour = layer),
          direction = "y",
          hjust = 1,
          nudge_x = -0.0015,
          size = 1.48,
          fontface = "bold",
          fill = alpha("white", 0.88),
          label.size = 0,
          label.padding = unit(0.55, "mm"),
          box.padding = 0.14,
          point.padding = 0.14,
          force = 0.55,
          force_pull = 1.6,
          min.segment.length = 0,
          segment.colour = separator,
          segment.size = 0.20,
          seed = 23,
          max.overlaps = Inf,
          show.legend = FALSE
        )
      }
    } +
    layer_scale +
    scale_fill_manual(values = c(host = host_colour, pathogen = pathogen_colour), guide = "none") +
    scale_shape_manual(values = c(`Prior accessible scaffold` = 21), name = NULL) +
    guides(
      colour = guide_legend(
        order = 1, nrow = 1,
        override.aes = list(
          shape = 21,
          fill = c(host_colour, pathogen_colour),
          colour = c(host_colour, pathogen_colour),
          size = 2.8,
          stroke = 0.45,
          alpha = 1
        )
      ),
      shape = guide_legend(
        order = 2, nrow = 1,
        override.aes = list(
          shape = 21, size = 3.2, fill = "white",
          colour = ink, stroke = 1.05, alpha = 1
        )
      )
    )
}

direct_panel <- add_module_bands(
  ggplot(nodes, aes(direct_evidence, y_plot))
)
direct_panel <- add_node_layers(direct_panel) +
  scale_x_continuous(
    limits = c(0, direct_hi), breaks = direct_breaks, labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = row_limits, breaks = row_breaks, labels = row_labels,
    expand = expansion(mult = 0)
  ) +
  labs(x = "Direct evidence\npercentile") +
  theme_metric

opposite_panel <- add_module_bands(
  ggplot(nodes, aes(opposite_signal, y_plot))
)
opposite_panel <- add_node_layers(opposite_panel, label_keys = TRUE) +
  scale_x_continuous(
    limits = c(0.948, 1.002), breaks = c(0.95, 0.975, 1),
    labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(limits = row_limits, expand = expansion(mult = 0)) +
  labs(x = "Opposite-side RWR\npercentile") +
  theme_metric +
  theme(axis.text.y = element_blank())

degree_panel <- add_module_bands(
  ggplot(nodes, aes(degree_log10, y_plot))
)
degree_panel <- add_node_layers(degree_panel) +
  scale_x_continuous(
    limits = c(0, log10(degree_max)),
    breaks = log10(degree_breaks),
    labels = comma(degree_breaks),
    expand = expansion(mult = c(0.01, 0.02)),
    oob = scales::oob_squish
  ) +
  scale_y_continuous(limits = row_limits, expand = expansion(mult = 0)) +
  labs(x = "Network degree\n(log scale)") +
  theme_metric +
  theme(axis.text.y = element_blank())

plot <- direct_panel + opposite_panel + degree_panel +
  plot_layout(
    widths = c(1.10, 1.05, 0.95),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(t = 1.5, r = 3.5, b = 0, l = 3.5)
  )

stem <- file.path(output_dir, "Figure_3_panel_c_propagated_nodes")
if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(
  paste0(stem, ".pdf"), plot,
  width = 190, height = 122, units = "mm",
  device = cairo_pdf, bg = "white", limitsize = FALSE
)
ggsave(
  paste0(stem, ".png"), plot,
  width = 190, height = 122, units = "mm", dpi = 600,
  device = ragg::agg_png, bg = "white", limitsize = FALSE
)
}

cat(
  "Rendered Figure 3c with all ", comma(nrow(nodes)), " recruited nodes: ",
  comma(sum(nodes$layer == "host")), " host and ",
  comma(sum(nodes$layer == "pathogen")), " pathogen; ",
  comma(sum(nodes$prior_scaffold)), " retain prior accessible-scaffold membership.\n",
  "Rows: ", paste(module_groups$module_display, collapse = " | "), ".\n",
  "degree_max=", degree_max, " breaks=", paste(degree_breaks, collapse = ","),
  " obs_degree_max=", obs_degree_max, ".\n",
  "Labeled: ", paste(sort(unique(key_nodes$display_label)), collapse = ", "), ".\n",
  sep = ""
)
