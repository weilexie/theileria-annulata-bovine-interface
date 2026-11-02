# Shared lollipop grammar for Figure 1a / 2a.

ink <- "#30363B"
line <- "#E4E9EC"
blue_dark <- "#B5C8D1"
red_dark <- "#CDB0AD"
red <- "#E4D0CD"
phospho_dark <- "#AE7772"

theme_studies_a <- theme_classic(base_size = 8, base_family = "sans") + theme(
  panel.grid = element_blank(),
  axis.line = element_blank(),
  panel.border = element_rect(fill = NA, colour = "#DCE3E7", linewidth = 0.35),
  axis.ticks = element_line(linewidth = 0.3, colour = ink),
  axis.text = element_text(colour = ink),
  axis.title = element_text(colour = ink),
  legend.position = "bottom",
  legend.justification = "center",
  legend.box.just = "center",
  legend.direction = "horizontal",
  legend.box = "horizontal",
  legend.key = element_blank(),
  legend.key.height = unit(2.8, "mm"),
  legend.key.width = unit(4.2, "mm"),
  legend.spacing.x = unit(1.8, "mm"),
  legend.box.spacing = unit(1.6, "mm"),
  legend.margin = margin(t = 1.5),
  plot.margin = margin(6, 10, 5, 8),
  panel.background = element_rect(fill = "white", colour = NA),
  plot.background = element_rect(fill = "white", colour = NA)
)

studies_a_palettes <- function(levels) {
  colour <- c(
    `Independent proteome` = red_dark,
    `Accessible proteome` = red_dark,
    Phosphoproteome = phospho_dark,
    `Interface-proximal` = blue_dark,
    `Proximity proteome` = blue_dark
  )
  fill <- c(
    `Independent proteome` = red,
    `Accessible proteome` = red,
    Phosphoproteome = red_dark,
    `Interface-proximal` = "white",
    `Proximity proteome` = "white"
  )
  list(
    colour = colour[levels],
    fill = fill[levels],
    shape = setNames(rep(21, length(levels)), levels)
  )
}

draw_studies_lollipop <- function(studies, y, group, x_lab, y_labels) {
  pal <- studies_a_palettes(levels(studies[[group]]))
  ggplot(studies, aes(effect_rows, .data[[y]], colour = .data[[group]], fill = .data[[group]], shape = .data[[group]])) +
    geom_segment(aes(x = 1, xend = effect_rows, yend = .data[[y]]), colour = line, linewidth = 0.5) +
    geom_point(size = 2.7, stroke = 0.7) +
    geom_text(aes(label = comma(effect_rows)), hjust = -0.35, size = 2.3, colour = ink) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale()), expand = expansion(mult = c(0.02, 0.32))) +
    scale_colour_manual(values = pal$colour, breaks = names(pal$colour)) +
    scale_fill_manual(values = pal$fill, breaks = names(pal$fill)) +
    scale_shape_manual(values = pal$shape, breaks = names(pal$shape)) +
    scale_y_discrete(labels = y_labels) +
    labs(x = x_lab, y = NULL, colour = NULL, fill = NULL, shape = NULL) +
    theme_studies_a +
    theme(axis.text.y = element_text(size = 6.4, lineheight = 0.92), legend.text = element_text(size = 6.0))
}
