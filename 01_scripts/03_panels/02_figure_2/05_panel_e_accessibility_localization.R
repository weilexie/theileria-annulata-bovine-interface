#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(scales)
  library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- vapply(
  sys.frames(),
  function(x) {
    z <- x$ofile
    if (is.null(z) || !length(z)) "" else as.character(z[[1]])
  },
  character(1)
)
script_path <- if (any(nzchar(source_files))) {
  normalizePath(tail(source_files[nzchar(source_files)], 1), winslash = "/", mustWork = TRUE)
} else if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path("01_scripts", "03_panels", "02_figure_2", "05_panel_e_accessibility_localization.R"),
    winslash = "/", mustWork = TRUE
  )
}
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
p <- function(...) file.path(root, ...)

data_override <- Sys.getenv("FIGURE_DATA_DIR", unset = "")
data_candidates <- c(data_override, p("04_plot_data"))
data_candidates <- data_candidates[nzchar(data_candidates) & dir.exists(data_candidates)]
if (!length(data_candidates)) stop("Run 01_scripts/02_build_figure_data.R first.")
input_file <- file.path(data_candidates[[1]], "02_figure_2", "03_panel_e_evidence_geometry.tsv")
if (!file.exists(input_file)) stop("Missing panel plot data: ", input_file)

output_override <- Sys.getenv("FIGURE_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) {
  output_override
} else {
  p("06_panel_output", "02_figure_2")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(plot, stem, width_mm = 183, height_mm = 102) {
  if (identical(Sys.getenv("FIGURE_RETURN_ONLY", ""), "1")) return(invisible(plot))
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width_mm, height = height_mm, units = "mm",
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot,
    width = width_mm, height = height_mm, units = "mm", dpi = 600,
    device = ragg::agg_png, bg = "white", limitsize = FALSE
  )
}

wilson_ci <- function(success, total, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- success / total
  den <- 1 + z^2 / total
  centre <- (phat + z^2 / (2 * total)) / den
  half <- z * sqrt(phat * (1 - phat) / total + z^2 / (4 * total^2)) / den
  c(low = max(0, centre - half), high = min(1, centre + half))
}

theileria_feature_levels <- c(
  "accessible_annot",
  "surface_export_trafficking",
  "membrane_bioenergetics",
  "cytoskeletal_engagement",
  "genome_expression_processing",
  "protein_folding_processing",
  "secreted_or_exported",
  "membrane_or_surface",
  "microneme_rhoptry"
)
theileria_feature_labels <- c(
  accessible_annot = "Accessible\nannotation",
  surface_export_trafficking = "Surface / export\ntrafficking",
  membrane_bioenergetics = "Membrane\nbioenergetics",
  cytoskeletal_engagement = "Cytoskeletal\nengagement",
  genome_expression_processing = "Genome-expression\nprocessing",
  protein_folding_processing = "Protein folding /\nprocessing",
  secreted_or_exported = "Secreted or\nexported",
  membrane_or_surface = "Membrane or\nsurface",
  microneme_rhoptry = "Microneme /\nrhoptry"
)

nodes <- read_tsv(
  input_file, show_col_types = FALSE, progress = FALSE,
  na = c("", "NA", "NaN")
)
if (!"scaffold_tier" %in% names(nodes) && "network_role" %in% names(nodes)) {
  nodes$scaffold_tier <- case_when(
    nodes$network_role %in% c("Core", "strict", "recurrent_core") ~ "strict",
    nodes$network_role %in% c("Accessible scaffold", "relaxed", "accessible_scaffold") ~ "relaxed",
    TRUE ~ NA_character_
  )
}
if (!"scaffold_tier" %in% names(nodes)) stop("Missing required column: scaffold_tier")
if (!"localization_class" %in% names(nodes)) nodes$localization_class <- NA_character_
if (!"interface_function" %in% names(nodes)) {
  nodes$interface_function <- if ("accessible_annot" %in% names(nodes)) {
    str_to_lower(coalesce(as.character(nodes$accessible_annot), "")) %in% c("true", "t", "1", "yes", "y")
  } else {
    FALSE
  }
}

nodes <- nodes %>%
  filter(scaffold_tier %in% c("strict", "relaxed", "recurrent_core", "accessible_scaffold", "Core", "Accessible scaffold")) %>%
  mutate(
    scaffold_tier = case_when(
      scaffold_tier %in% c("strict", "recurrent_core", "Core") ~ "strict",
      TRUE ~ "relaxed"
    ),
    tier = factor(
      scaffold_tier,
      levels = c("strict", "relaxed"),
      labels = c("Strict core", "Accessible scaffold")
    ),
    interface_function = as.logical(interface_function)
  )
nodes$interface_function[is.na(nodes$interface_function)] <- FALSE

tier_n <- table(nodes$tier)
core_n <- as.integer(tier_n[["Strict core"]])
scaffold_n <- as.integer(tier_n[["Accessible scaffold"]])
if (!is.finite(core_n) || !is.finite(scaffold_n) || min(core_n, scaffold_n) < 1) {
  stop("Need both strict-core and accessible-scaffold proteins for panel E.")
}
if (!identical(c(core_n, scaffold_n), c(19L, 16L))) {
  warning(
    "Expected 19 strict-core and 16 accessible-scaffold proteins; observed ",
    paste(names(tier_n), as.integer(tier_n), collapse = ", ")
  )
}
if ("node_id" %in% names(nodes)) {
  pin1 <- grepl("Q4UG71", nodes$node_id)
  if (any(pin1 & nodes$tier == "Strict core")) {
    stop("Pin1 is never local-core and must not appear in the strict core.")
  }
}

observed_classes <- unique(na.omit(as.character(nodes$localization_class)))
feature_levels <- intersect(theileria_feature_levels, observed_classes)
feature_levels <- unique(c(feature_levels, setdiff(observed_classes, theileria_feature_levels)))
feature_levels <- feature_levels[nzchar(feature_levels) & !is.na(feature_levels)]
feature_levels <- setdiff(feature_levels, c(
  "interface_function", "cell_wall_or_anchored",
  "experimental_surface_possible_moonlighter",
  "secreted_or_extracellular", "experimentally_surface_or_released"
))
if (!length(feature_levels)) stop("No Theileria localization classes in panel E data.")
feature_labels <- ifelse(
  feature_levels %in% names(theileria_feature_labels),
  unname(theileria_feature_labels[feature_levels]),
  str_replace_all(feature_levels, "_", " ")
)
names(feature_labels) <- feature_levels

feature_member <- function(dat, feature) {
  if (identical(feature, "interface_function")) {
    !is.na(dat$interface_function) & dat$interface_function
  } else {
    !is.na(dat$localization_class) & dat$localization_class == feature
  }
}

prop_rows <- lapply(feature_levels, function(feature) {
  member <- feature_member(nodes, feature)
  bind_rows(lapply(levels(nodes$tier), function(group) {
    in_group <- nodes$tier == group
    count <- sum(member & in_group)
    total <- sum(in_group)
    ci <- wilson_ci(count, total)
    tibble(
      feature = feature,
      tier = factor(group, levels = levels(nodes$tier)),
      count = count,
      total = total,
      proportion = count / total,
      prop_low = ci[["low"]],
      prop_high = ci[["high"]]
    )
  }))
}) %>% bind_rows()

enrichment <- lapply(feature_levels, function(feature) {
  member <- feature_member(nodes, feature)
  core_yes <- sum(member & nodes$tier == "Strict core")
  core_no <- sum(!member & nodes$tier == "Strict core")
  scaffold_yes <- sum(member & nodes$tier == "Accessible scaffold")
  scaffold_no <- sum(!member & nodes$tier == "Accessible scaffold")
  ft <- fisher.test(matrix(
    c(core_yes, core_no, scaffold_yes, scaffold_no),
    nrow = 2, byrow = TRUE
  ))
  # Zero cells make Fisher OR infinite; Haldane-Anscombe keeps a finite
  # gold-style point+CI. P-values stay on the original 2x2.
  if (min(core_yes, core_no, scaffold_yes, scaffold_no) == 0) {
    aa <- core_yes + 0.5
    bb <- core_no + 0.5
    cc <- scaffold_yes + 0.5
    dd <- scaffold_no + 0.5
    odds_ratio <- (aa * dd) / (bb * cc)
    se <- sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    z <- qnorm(0.975)
    or_low <- exp(log(odds_ratio) - z * se)
    or_high <- exp(log(odds_ratio) + z * se)
  } else {
    odds_ratio <- unname(ft$estimate)
    or_low <- unname(ft$conf.int[[1]])
    or_high <- unname(ft$conf.int[[2]])
  }
  tibble(
    feature = feature,
    core_yes = core_yes,
    scaffold_yes = scaffold_yes,
    odds_ratio = odds_ratio,
    or_low = or_low,
    or_high = or_high,
    p_value = ft$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    log2_or = log2(odds_ratio),
    log2_or_low = log2(or_low),
    log2_or_high = log2(or_high)
  )

row_key <- tibble(
  feature = feature_levels,
  row = rev(seq_along(feature_levels)),
  label = unname(feature_labels[feature_levels])
)
prop_rows <- left_join(prop_rows, row_key, by = "feature") %>%
  mutate(
    y_plot = row + if_else(tier == "Strict core", 0.13, -0.13)
  )
enrichment <- left_join(enrichment, row_key, by = "feature") %>%
  mutate(
    stat_label = paste0("FDR ", pvalue(fdr, accuracy = 0.001)),
    significant = fdr < 0.05
  )
fdr_x <- max(enrichment$log2_or_high, enrichment$log2_or, na.rm = TRUE) + 1.15

ink <- "#30363B"
border <- "#DCE3E7"
blue_dark <- "#B5C8D1"
blue <- "#D2DEE4"
red_dark <- "#CDB0AD"
red <- "#E4D0CD"
core_fill <- red
core_dark <- red_dark
scaffold_fill <- blue
scaffold_dark <- blue_dark

theme_panel <- theme_classic(base_size = 8.2, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(fill = NA, colour = border, linewidth = 0.38),
    axis.ticks = element_line(linewidth = 0.32, colour = ink),
    axis.text = element_text(colour = ink),
    axis.title = element_text(colour = ink),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.key.height = unit(2.8, "mm"),
    legend.key.width = unit(4.2, "mm"),
    legend.text = element_text(size = 6.7),
    legend.spacing.x = unit(1.4, "mm"),
    legend.box.spacing = unit(2.0, "mm"),
    legend.margin = margin(t = 1.5),
    plot.margin = margin(5, 7, 6, 7),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

panel_prop <- ggplot(prop_rows, aes(x = proportion, y = y_plot, colour = tier)) +
  geom_errorbar(
    aes(xmin = prop_low, xmax = prop_high),
    orientation = "y", width = 0, linewidth = 0.72
  ) +
  geom_point(
    aes(fill = tier), shape = 21, size = 3.25, stroke = 0.62,
    show.legend = TRUE
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1.08), expand = expansion(mult = c(0.01, 0))
  ) +
  scale_y_continuous(
    breaks = row_key$row, labels = row_key$label,
    limits = c(0.45, length(feature_levels) + 0.55), expand = expansion(mult = 0)
  ) +
  scale_colour_manual(
    values = c("Strict core" = core_dark, "Accessible scaffold" = scaffold_dark),
    labels = c(
      sprintf("Strict core (n = %d)", core_n),
      sprintf("Accessible scaffold (n = %d)", scaffold_n)
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Strict core" = core_fill, "Accessible scaffold" = scaffold_fill),
    labels = c(
      sprintf("Strict core (n = %d)", core_n),
      sprintf("Accessible scaffold (n = %d)", scaffold_n)
    ),
    name = NULL
  ) +
  guides(
    colour = guide_legend(override.aes = list(fill = NA, shape = 21, size = 3.0, stroke = 0.75)),
    fill = "none"
  ) +
  labs(x = "Proteins in tier (Wilson 95% CI)", y = NULL) +
  theme_panel +
  theme(
    axis.text.y = element_text(size = 7.1, hjust = 1, lineheight = 0.92),
    axis.text.x = element_text(size = 6.9),
    axis.title.x = element_text(size = 7.3, margin = margin(t = 5)),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 4, 6, 7)
  )

panel_or <- ggplot(enrichment, aes(x = log2_or, y = row)) +
  geom_vline(xintercept = 0, colour = "#6F767D", linewidth = 0.48, linetype = "22") +
  geom_errorbar(
    aes(xmin = log2_or_low, xmax = log2_or_high, colour = significant),
    orientation = "y", width = 0, linewidth = 0.78
  ) +
  geom_point(
    aes(fill = significant, colour = significant),
    shape = 21, size = 3.35, stroke = 0.68
  ) +
  geom_text(
    aes(x = fdr_x, label = stat_label, colour = significant),
    hjust = 0, size = 2.25, fontface = ifelse(enrichment$significant, "bold", "plain")
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(6),
    limits = c(
      min(enrichment$log2_or_low, enrichment$log2_or, na.rm = TRUE) - 0.45,
      fdr_x + 7.4
    ),
    expand = expansion(mult = c(0.01, 0)),
    oob = squish
  ) +
  scale_y_continuous(
    breaks = row_key$row, labels = NULL,
    limits = c(0.45, length(feature_levels) + 0.55), expand = expansion(mult = 0)
  ) +
  scale_colour_manual(values = c(`FALSE` = "#6F767D", `TRUE` = core_dark), guide = "none") +
  scale_fill_manual(values = c(`FALSE` = "#F4F6F7", `TRUE` = core_fill), guide = "none") +
  labs(
    x = expression("Enrichment in strict core, log"[2] * " odds ratio (95% CI)"),
    y = NULL
  ) +
  theme_panel +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6.9),
    axis.title.x = element_text(size = 7.3, margin = margin(t = 5)),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 7, 6, 2)
  )

panel <- panel_prop + panel_or +
  plot_layout(widths = c(1.28, 1.08), guides = "collect") &
  theme(legend.position = "bottom")

save_panel(panel, "Figure_2_panel_e_accessibility_localization")

cat("Rendered Figure_2_panel_e_accessibility_localization.\n")
report <- enrichment %>%
  left_join(
    prop_rows %>%
      filter(tier == "Strict core") %>%
      transmute(feature, core_count = count, core_total = total, core_prop = proportion),
    by = "feature"
  ) %>%
  left_join(
    prop_rows %>%
      filter(tier == "Accessible scaffold") %>%
      transmute(feature, scaffold_count = count, scaffold_total = total, scaffold_prop = proportion),
    by = "feature"
  ) %>%
  transmute(
    feature = label,
    strict_core = sprintf("%d/%d (%.1f%%)", core_count, core_total, 100 * core_prop),
    accessible_scaffold = sprintf("%d/%d (%.1f%%)", scaffold_count, scaffold_total, 100 * scaffold_prop),
    log2_or = round(log2_or, 3),
    fdr = signif(fdr, 3)
  )
print(report)
