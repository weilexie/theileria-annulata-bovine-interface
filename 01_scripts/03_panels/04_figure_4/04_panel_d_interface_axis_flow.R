#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
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
    "01_scripts/03_panels/04_figure_4/04_panel_d_interface_axis_flow.R",
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

muted <- "#6F767D"
# Gold 4d grey wash. White alluvium edges keep the 11 ribbons from merging.
interface_flow <- "#DDD7DE"

axes <- read_fig("07_panel_d_axis_descriptor.tsv")
for (col in c(
  "biological_axis_id", "axis_id", "interaction_name", "pathogen_gene", "host_gene",
  "parasite_factor_name", "host_target_name", "final_tier", "pathogen_program",
  "host_program", "proximity_pathogen_module", "proximity_host_module",
  "context_only", "retained"
)) {
  if (!col %in% names(axes)) axes[[col]] <- NA
}
# Gold two-line program labels: wrap at " / " or "-broad".
# 1-axis names that already fit (Transcription, Chromatin, Metabolism) stay one line.
wrap_program <- function(z) {
  z <- pretty_module(z, wrap = FALSE)
  z <- case_when(
    z == "Ubiquitin / proteostasis" ~ "Ubiquitin /\nproteostasis",
    z == "Cellular process-broad" ~ "Cellular\nprocess-broad",
    z == "Microtubule / cytoskeleton" ~ "Microtubule /\ncytoskeleton",
    z == "Signaling-broad" ~ "Signaling-\nbroad",
    z == "Translation / folding" ~ "Translation /\nfolding",
    z == "Surface family" ~ "Surface\nfamily",
    z == "RNA processing" ~ "RNA\nprocessing",
    z == "Ubiquitin / cullin" ~ "Ubiquitin /\ncullin",
    TRUE ~ z
  )
  needs <- !str_detect(z, "\n")
  z[needs] <- str_replace(z[needs], " / ", " /\n")
  still <- !str_detect(z, "\n")
  z[still] <- str_wrap(z[still], width = 11, whitespace_only = TRUE)
  still <- !str_detect(z, "\n")
  z[still] <- str_replace(z[still], "-broad", "-\nbroad")
  z
}

axes <- axes %>%
  mutate(
    axis_id_key = axis_key(coalesce(biological_axis_id, axis_id), interaction_name, pathogen_gene, host_gene),
    final_tier = axis_tier(axis_id_key, final_tier),
    context_only = flag_col(context_only),
    retained_flag = if (all(is.na(retained))) TRUE else flag_col(retained)
  ) %>%
  filter(
    final_tier %in% c("Core", "Extended"),
    !context_only,
    retained_flag,
    axis_id_key != "p104--CLASP1",
    !str_detect(axis_id_key, "CLASP1")
  ) %>%
  mutate(
    pathogen_gene = endpoint_display(sanitize_label(coalesce(pathogen_gene, parasite_factor_name, "unresolved"))),
    host_gene = endpoint_display(sanitize_label(coalesce(host_gene, host_target_name, "unresolved"))),
    pathogen_program = wrap_program(coalesce(pathogen_program, proximity_pathogen_module)),
    host_program = wrap_program(coalesce(host_program, proximity_host_module)),
    axis_weight = 1
  )

retained_axis_ids <- c(core_axis_ids, extended_axis_ids)
stopifnot(nrow(axes) == 11L)
stopifnot(setequal(axes$axis_id_key, retained_axis_ids))
stopifnot(!any(str_detect(paste(axes$pathogen_program, axes$host_program), "H_NET|P_NET|PIase|CLASP1")))
stopifnot(!any(str_detect(paste(axes$pathogen_gene, axes$host_gene), "TA_15795|CLASP1|正文|文中")))

# Silk-ribbon layout: bundle on the left, barycentric host programs so the
# middle fans instead of crossing the page, targets stay inside their host
# program so right-side curves stay local.
pathogen_program_order <- axes %>%
  count(pathogen_program, name = "n") %>%
  arrange(desc(n), pathogen_program) %>%
  pull(pathogen_program)
pathogen_pos <- setNames(seq_along(pathogen_program_order), pathogen_program_order)

ordered_outer <- function(data, gene_field, program_field, program_order) {
  primary <- data %>%
    count(.data[[gene_field]], .data[[program_field]], name = "within_n") %>%
    mutate(program_order = match(.data[[program_field]], program_order)) %>%
    group_by(.data[[gene_field]]) %>%
    arrange(desc(within_n), program_order, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(gene = .data[[gene_field]], program = .data[[program_field]])

  data %>%
    count(.data[[gene_field]], name = "axis_n") %>%
    transmute(gene = .data[[gene_field]], axis_n) %>%
    left_join(primary, by = "gene") %>%
    mutate(
      program_order = match(program, program_order),
      program_order = coalesce(program_order, length(program_order) + 1L)
    ) %>%
    arrange(program_order, desc(axis_n), gene) %>%
    pull(gene)
}

factor_order <- ordered_outer(
  axes, "pathogen_gene", "pathogen_program", pathogen_program_order
)
factor_idx <- setNames(seq_along(factor_order), factor_order)
host_program_order <- axes %>%
  group_by(host_program) %>%
  summarise(
    n = dplyr::n(),
    bary = mean(pathogen_pos[as.character(pathogen_program)]),
    .groups = "drop"
  ) %>%
  arrange(desc(n), bary, host_program) %>%
  pull(host_program)
target_order <- axes %>%
  mutate(factor_order = unname(factor_idx[as.character(pathogen_gene)])) %>%
  arrange(factor_order, host_gene) %>%
  distinct(host_gene) %>%
  pull(host_gene)
pathogen_programs <- pathogen_program_order[
  pathogen_program_order %in% axes$pathogen_program
]
host_programs <- host_program_order[host_program_order %in% axes$host_program]
axes <- axes %>%
  mutate(
    .p_ord = match(pathogen_program, pathogen_program_order),
    .h_ord = match(host_program, host_program_order),
    .f_ord = unname(factor_idx[as.character(pathogen_gene)])
  ) %>%
  arrange(.p_ord, .h_ord, .f_ord, pathogen_gene, host_gene) %>%
  select(-.p_ord, -.h_ord, -.f_ord)

# Tight, even pastels so one side is not a dark slab and the other paper-white.
# Ribbons are a single warm lilac that sits between terracotta and blue.
pathogen_fill_steps <- c("#E4D0CC", "#EBD9D6", "#F0E4E1")
host_fill_steps <- c("#C9D6DE", "#D3DEE5", "#DCE6EB")
pastel_family <- function(names, counts, steps) {
  n <- length(names)
  if (!n) return(setNames(character(), character()))
  ranked <- names[order(counts[names], decreasing = TRUE, na.last = TRUE)]
  pal <- steps[pmin(seq_len(n), length(steps))]
  setNames(pal, ranked)[names]
}
pathogen_n <- as.numeric(table(axes$pathogen_program)[pathogen_programs])
host_n <- as.numeric(table(axes$host_program)[host_programs])
names(pathogen_n) <- pathogen_programs
names(host_n) <- host_programs
pathogen_program_colours <- pastel_family(
  pathogen_programs, pathogen_n, pathogen_fill_steps
)
host_program_colours <- pastel_family(
  host_programs, host_n, host_fill_steps
)
ribbon <- "#D4C8D2"

dominant_program <- function(data, gene_field, program_field, program_order) {
  data %>%
    count(.data[[gene_field]], .data[[program_field]], name = "axis_n") %>%
    mutate(program_order = match(.data[[program_field]], program_order)) %>%
    group_by(.data[[gene_field]]) %>%
    arrange(desc(axis_n), program_order, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    transmute(gene = .data[[gene_field]], program = .data[[program_field]]) %>%
    { setNames(.$program, .$gene) }
}

factor_program <- dominant_program(
  axes, "pathogen_gene", "pathogen_program", pathogen_program_order
)
target_program <- dominant_program(
  axes, "host_gene", "host_program", host_program_order
)

axis_levels <- c(
  "Pathogen factors",
  "Pathogen programs",
  "Host programs",
  "Host targets"
)
stratum_levels <- c(
  paste0("factor::", rev(factor_order)),
  paste0("pathogen_program::", rev(pathogen_programs)),
  paste0("host_program::", rev(host_programs)),
  paste0("target::", rev(target_order))
)

lodes <- bind_rows(
  axes %>% transmute(
    axis_id = axis_id_key, axis = axis_levels[1], role = "factor",
    stratum_label = pathogen_gene, pathogen_program, host_program, axis_weight
  ),
  axes %>% transmute(
    axis_id = axis_id_key, axis = axis_levels[2], role = "pathogen_program",
    stratum_label = pathogen_program, pathogen_program, host_program, axis_weight
  ),
  axes %>% transmute(
    axis_id = axis_id_key, axis = axis_levels[3], role = "host_program",
    stratum_label = host_program, pathogen_program, host_program, axis_weight
  ),
  axes %>% transmute(
    axis_id = axis_id_key, axis = axis_levels[4], role = "target",
    stratum_label = host_gene, pathogen_program, host_program, axis_weight
  )
) %>%
  mutate(
    axis = factor(axis, levels = axis_levels),
    stratum = factor(paste0(role, "::", stratum_label), levels = stratum_levels),
    flow_fill = ribbon,
    stratum_fill = case_when(
      role == "factor" ~ unname(
        pathogen_program_colours[factor_program[stratum_label]]
      ),
      role == "pathogen_program" ~ unname(
        pathogen_program_colours[stratum_label]
      ),
      role == "host_program" ~ unname(host_program_colours[stratum_label]),
      TRUE ~ unname(host_program_colours[target_program[stratum_label]])
    )
  )

stopifnot(
  nrow(lodes) == 4L * nrow(axes),
  n_distinct(lodes$axis_id) == nrow(axes),
  all(!is.na(lodes$stratum)),
  all(!is.na(lodes$flow_fill)),
  all(!is.na(lodes$stratum_fill))
)

plot <- ggplot(
  lodes,
  aes(x = axis, stratum = stratum, alluvium = axis_id, y = axis_weight)
) +
  ggalluvial::geom_alluvium(
    aes(fill = flow_fill),
    width = 0.22,
    knot.pos = 0.48,
    alpha = 0.38,
    colour = "white",
    linewidth = 0.16,
    lode.guidance = "frontback",
    aes.bind = "flows"
  ) +
  ggalluvial::geom_stratum(
    aes(fill = stratum_fill),
    width = 0.26,
    colour = "white",
    linewidth = 0.32,
    alpha = 0.96
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(sub("^[^:]+::", "", stratum))),
    size = 1.95,
    lineheight = 0.82,
    colour = ink
  ) +
  scale_x_discrete(
    limits = axis_levels,
    position = "top",
    expand = c(0.055, 0.055)
  ) +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0.015, 0.015))) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = NULL,
    caption = paste(
      "Eleven retained biological axes connect parasite factors to host targets",
      "through paired biological programs; block height is axis count."
    )
  ) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 6.7, face = "bold", colour = ink, margin = margin(b = 5)
    ),
    axis.text.y = element_blank(),
    plot.caption = element_text(
      size = 5.7, colour = muted, hjust = 0.5, margin = margin(t = 5)
    ),
    plot.margin = margin(5, 8, 4, 8)
  )

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) {
ggsave(
  file.path(output_dir, "Figure_4_panel_d_interface_axis_flow.pdf"),
  plot,
  width = 235, height = 158, units = "mm",
  device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure_4_panel_d_interface_axis_flow.png"),
  plot,
  width = 235, height = 158, units = "mm", dpi = 600,
  device = ragg::agg_png, bg = "white"
)
}

cat("Rendered Figure_4_panel_d_interface_axis_flow as a four-axis alluvial map.\n")
