#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(circlize)
  library(dplyr)
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
  normalizePath(file.path("01_scripts", "03_panels", "05_figure_5", "02_panel_b_module_topology_chord.R"), winslash = "/", mustWork = TRUE)
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
data_dir <- data_candidates[file.exists(file.path(data_candidates, "05_figure_5", "03_panel_b_module_nodes.tsv"))][1]
if (is.na(data_dir)) stop("Run 01_scripts/02_build_figure_data.R first.")
read_fig <- function(file) read_tsv(file.path(data_dir, "05_figure_5", file), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) output_override else p("06_panel_output", "05_figure_5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#30363B"
host <- "#D0DEE5"; host_dark <- "#8099AA"; host_link <- alpha("#AFC3CE", 0.25)
pathogen <- "#E2C2BF"; pathogen_dark <- "#AE7772"; pathogen_link <- alpha("#D5A7A3", 0.24)
interface <- "#927D96"; interface_pale <- "#E0D6DE"; interface_link <- alpha("#B6A4B9", 0.55)

# Named programs only; unnamed NET IDs collapse to Other / unresolved (Figure 4c honesty).
is_raw_net <- function(z) {
  z <- str_squish(coalesce(as.character(z), ""))
  !nzchar(z) |
    z %in% c("Unassigned", "unassigned") |
    str_detect(z, "^[HP]_NET(_|\\s*)?[0-9]+$") |
    str_detect(z, "^[HP] NET [0-9]+$") |
    str_detect(z, "^NET[ _]?[0-9]+$")
}

pretty_display <- function(id, label) {
  id <- coalesce(as.character(id), "")
  label <- coalesce(as.character(label), "")
  from_id <- pretty_module(id, wrap = FALSE)
  from_label <- pretty_module(label, wrap = FALSE)
  leftover <- function(z) {
    z <- str_squish(coalesce(z, ""))
    !nzchar(z) | z %in% c("Unassigned", id, label) | is_raw_net(z)
  }
  dplyr::case_when(
    !leftover(from_id) ~ from_id,
    !leftover(from_label) ~ from_label,
    TRUE ~ "Other / unresolved"
  )
}

sector_slug <- function(layer, display_base) {
  paste(layer, str_replace_all(str_to_lower(display_base), "[^a-z0-9]+", "_"), sep = "::")
}

wrap_sector_label <- function(z) {
  z <- str_replace(z, " / ", " /\n")
  z <- str_replace(z, "-broad", "-\nbroad")
  z <- str_replace(z, "RNA processing", "RNA\nprocessing")
  z
}

table10_path <- p("05_extra_tables", "01_module_display_names.tsv")
table10 <- if (file.exists(table10_path)) {
  read_tsv(table10_path, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
} else {
  tibble(module_id = character(), layer = character(), display_name = character())
}

named_order <- function(layer_name) {
  ids <- table10 %>% filter(.data$layer == layer_name)
  labs <- unique(c(ids$display_name, ids$module_id))
  labs <- pretty_module(labs, wrap = FALSE)
  unique(labs[!is_raw_net(labs) & labs != "Unassigned" & nzchar(str_squish(labs))])
}

host_named <- named_order("host")
pathogen_named <- named_order("pathogen")

modules_raw <- read_fig("03_panel_b_module_nodes.tsv") %>%
  mutate(
    display_base = pretty_display(module_label, module_label),
    module_key_new = sector_slug(layer, display_base)
  )

modules <- modules_raw %>%
  mutate(
    order_key = case_when(
      display_base == "Other / unresolved" & layer == "host" ~ length(host_named) + 100L,
      display_base == "Other / unresolved" & layer == "pathogen" ~ length(pathogen_named) + 100L,
      layer == "host" ~ match(display_base, host_named),
      TRUE ~ match(display_base, pathogen_named)
    ),
    order_key = if_else(is.na(order_key), 50L, as.integer(order_key))
  ) %>%
  group_by(layer, display_base, module_key_new) %>%
  summarise(
    scaffold_nodes = sum(scaffold_nodes),
    core_endpoint_n = sum(core_endpoint_n),
    hub_n = sum(hub_n),
    order_key = min(order_key),
    .groups = "drop"
  ) %>%
  mutate(
    module_key = module_key_new,
    sector_span = sqrt(scaffold_nodes) + 2,
    sector_fill = if_else(layer == "host", host, pathogen),
    sector_dark = if_else(layer == "host", host_dark, pathogen_dark)
  )

dup_names <- unique(modules$display_base[duplicated(modules$display_base) | duplicated(modules$display_base, fromLast = TRUE)])

modules <- modules %>%
  mutate(
    module_text = case_when(
      display_base %in% dup_names & layer == "host" ~ paste("Host", display_base),
      display_base %in% dup_names & layer == "pathogen" ~ paste("Pathogen", display_base),
      TRUE ~ display_base
    ),
    module_text = wrap_sector_label(module_text)
  ) %>%
  arrange(factor(layer, levels = c("host", "pathogen")), order_key)

key_map <- modules_raw %>%
  distinct(module_key, module_key_new)

module_edges <- read_fig("04_panel_b_module_edges.tsv") %>%
  left_join(key_map %>% transmute(from_module = module_key, from_new = module_key_new), by = "from_module") %>%
  left_join(key_map %>% transmute(to_module = module_key, to_new = module_key_new), by = "to_module") %>%
  mutate(
    from_module = coalesce(from_new, from_module),
    to_module = coalesce(to_new, to_module)
  ) %>%
  filter(from_module != to_module) %>%
  mutate(
    a = if_else(from_module < to_module, from_module, to_module),
    b = if_else(from_module < to_module, to_module, from_module)
  ) %>%
  group_by(from_module = a, to_module = b) %>%
  summarise(
    edge_n = sum(edge_n),
    weight_sum = sum(weight_sum),
    cross_axis_n = sum(cross_axis_n),
    .groups = "drop"
  ) %>%
  left_join(modules %>% select(from_module = module_key, layer_from = layer), by = "from_module") %>%
  left_join(modules %>% select(to_module = module_key, layer_to = layer), by = "to_module") %>%
  mutate(
    edge_class = case_when(
      layer_from != layer_to ~ "Interface axes",
      layer_from == "host" ~ "Host scaffold",
      TRUE ~ "Pathogen scaffold"
    ),
    link_col = case_when(
      edge_class == "Interface axes" ~ interface_link,
      edge_class == "Host scaffold" ~ host_link,
      TRUE ~ pathogen_link
    ),
    link_border = if_else(edge_class == "Interface axes", alpha(interface, 0.55), NA_character_)
  )

other_n <- modules %>%
  filter(display_base == "Other / unresolved") %>%
  count(layer, name = "n")

stopifnot(nrow(modules) > 0)
stopifnot(n_distinct(modules$module_key) == nrow(modules))
stopifnot(sum(modules$scaffold_nodes) == sum(modules_raw$scaffold_nodes))
stopifnot(all(other_n$n <= 1L))
stopifnot(all(module_edges$from_module %in% modules$module_key), all(module_edges$to_module %in% modules$module_key))
stopifnot(!any(str_detect(modules$module_text, "NET[ _]?[0-9]")))
stopifnot(!any(is.na(module_edges$layer_from)), !any(is.na(module_edges$layer_to)))

# Allocate every ribbon within a sector. Ribbon width is proportional to retained edge count,
# while sector span independently reflects scaffold-node burden.
incidence <- bind_rows(
  module_edges %>% transmute(edge_id = row_number(), module = from_module, edge_n),
  module_edges %>% transmute(edge_id = row_number(), module = to_module, edge_n)
) %>%
  group_by(module) %>%
  arrange(desc(edge_n), edge_id, .by_group = TRUE) %>%
  mutate(total = sum(edge_n), frac_start = lag(cumsum(edge_n), default = 0) / total, frac_end = cumsum(edge_n) / total) %>%
  ungroup() %>%
  left_join(modules %>% select(module = module_key, sector_span), by = "module") %>%
  mutate(pos_start = frac_start * sector_span, pos_end = frac_end * sector_span)

edge_pos <- module_edges %>% mutate(edge_id = row_number()) %>%
  left_join(incidence %>% select(edge_id, from_module = module, from_start = pos_start, from_end = pos_end), by = c("edge_id", "from_module")) %>%
  left_join(incidence %>% select(edge_id, to_module = module, to_start = pos_start, to_end = pos_end), by = c("edge_id", "to_module"))

stopifnot(!any(is.na(edge_pos$from_start)), !any(is.na(edge_pos$to_start)))

draw_panel <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit({circos.clear(); par(old_par)}, add = TRUE)
  par(mar = c(4.75, 1.6, 1.0, 1.6), xpd = NA, family = "Arial")
  circos.clear()

  sector_order <- modules$module_key
  n_host <- sum(modules$layer == "host")
  gap_after <- rep(1.4, length(sector_order))
  gap_after[c(n_host, length(sector_order))] <- 8
  circos.par(
    start.degree = 94,
    clock.wise = TRUE,
    gap.after = gap_after,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.004, 0.004),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-1.22, 1.22),
    canvas.ylim = c(-1.18, 1.18)
  )
  circos.initialize(
    factors = sector_order,
    xlim = cbind(rep(0, nrow(modules)), modules$sector_span)
  )

  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.13, bg.border = NA, panel.fun = function(x, y) {})

  for (cls in c("Host scaffold", "Pathogen scaffold", "Interface axes")) {
    idx <- which(edge_pos$edge_class == cls)
    if (!length(idx)) next
    for (i in idx) {
      e <- edge_pos[i, ]
      circos.link(
        e$from_module, c(e$from_start, e$from_end),
        e$to_module, c(e$to_end, e$to_start),
        col = e$link_col, border = e$link_border,
        lwd = if (cls == "Interface axes") 0.45 else 0.2,
        h.ratio = if (cls == "Interface axes") 0.72 else 0.56
      )
    }
  }

  # Redraw the sector track after ribbons so module burden and endpoint rims remain legible.
  for (j in seq_len(nrow(modules))) {
    m <- modules[j, ]
    key <- m$module_key
    xlim <- c(0, m$sector_span)
    xc <- mean(xlim)
    circos.rect(xlim[1], 0.02, xlim[2], 0.77, col = m$sector_fill, border = alpha(m$sector_dark, 0.72), lwd = 0.55, sector.index = key, track.index = 1)
    if (m$core_endpoint_n > 0) {
      circos.rect(xlim[1], 0.78, xlim[2], 1.00, col = interface_pale, border = interface, lwd = 1.0, sector.index = key, track.index = 1)
      circos.lines(xlim, c(0.995, 0.995), col = interface, lwd = 2.2, sector.index = key, track.index = 1)
    }
    circos.text(xc, 0.38, as.character(m$scaffold_nodes), cex = 0.47, col = ink, font = 2, facing = "bending.inside", niceFacing = TRUE, sector.index = key, track.index = 1)
    circos.text(
      xc, 1.10, m$module_text,
      cex = 0.45, col = m$sector_dark, font = 1,
      facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
      sector.index = key, track.index = 1
    )
  }

  # Fixed-column legend.  Base legend() sizes columns from text width, which
  # causes the two rows to drift and leaves the group label detached.  Drawing
  # the glyphs directly keeps all four columns aligned and the whole block centered.
  usr <- par("usr")
  x0 <- usr[1]; x1 <- usr[2]; y0 <- usr[3]; y1 <- usr[4]
  W <- diff(c(x0, x1)); H <- diff(c(y0, y1))
  # Four compact, fixed columns: group label + three keyed items.
  # The block is deliberately narrower than the old text-width legend.
  col_x <- x0 + W * c(0.165, 0.350, 0.545, 0.735)
  y_module <- y0 - H * 0.028
  y_links <- y0 - H * 0.088
  legend_cex <- 0.50
  txt <- function(x, y, label, font = 1) {
    text(x, y, label, adj = c(0, 0.5), cex = legend_cex, col = ink, font = font, xpd = NA)
  }
    glyph_text <- function(x, y, label, kind, col) {
      if (kind == "square") {
        # Compact sector swatches: the main chord already carries the area
        # encoding, so the legend should only provide a small color cue.
        rect(x, y - H * 0.0075, x + W * 0.0082, y + H * 0.0075,
          col = col, border = alpha(col, 0.9), lwd = 0.7, xpd = NA)
        txt(x + W * 0.023, y, label)
      } else if (kind == "ring") {
        # A short curved band mirrors the continuous outer rim drawn around
        # core-endpoint sectors; a dot would imply node-level encoding.
        theta <- seq(-pi / 2, pi / 2, length.out = 32)
        rx <- W * 0.0056
        ry <- H * 0.015
        xc <- x + W * 0.0068
        lines(xc + rx * cos(theta), y + ry * sin(theta),
          col = interface, lwd = 2.1, lend = 1, xpd = NA)
        txt(x + W * 0.023, y, label)
      } else {
      segments(x, y, x + W * 0.033, y, col = col, lwd = if (label == "Interface axes") 3.8 else 2.8, lend = 1, xpd = NA)
      txt(x + W * 0.045, y, label)
    }
  }
  txt(col_x[1], y_module, "Module sectors", font = 2)
  glyph_text(col_x[2], y_module, "Host", "square", host)
  glyph_text(col_x[3], y_module, "Pathogen", "square", pathogen)
  glyph_text(col_x[4], y_module, "Core-endpoint ring", "ring", interface)
  txt(col_x[1], y_links, "Links", font = 2)
  glyph_text(col_x[2], y_links, "Host scaffold", "line", host_dark)
  glyph_text(col_x[3], y_links, "Pathogen scaffold", "line", pathogen_dark)
  glyph_text(col_x[4], y_links, "Interface axes", "line", interface)
}

save_panel <- function(stem, width_mm = 160, height_mm = 162, dpi = 600) {
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  grDevices::cairo_pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial", bg = "white")
  draw_panel()
  dev.off()
  ragg::agg_png(png_path, width = width_mm, height = height_mm, units = "mm", res = dpi, background = "white")
  draw_panel()
  dev.off()
}

if (!identical(Sys.getenv("FIGURE_RETURN_ONLY", unset = ""), "1")) {
  save_panel("Figure_5_panel_b_module_topology_chord")
}
cat("Rendered Figure_5_panel_b_module_topology_chord as a ", nrow(modules), "-module chord diagram.\n", sep = "")
