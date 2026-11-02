pretty_gene <- function(x) {
  x <- coalesce(as.character(x), "")
  x <- str_replace(x, "^ENSBTAP[0-9]+$", "unmapped")
  x <- str_replace(x, "P68402$", "")
  x
}

axis_genes <- c("MAPRE1", "PKM", "CD2AP")
independent_families <- c("pxd010713_ev", "pxd056906_ev", "dual_phospho")
proximal_families <- c("schizont_copur", "mc4404", "bioid")
host_family_levels <- c(independent_families, proximal_families)
host_family_labels <- c(
  pxd010713_ev = "TBL20 EV",
  pxd056906_ev = "Serum EV",
  dual_phospho = "Dual phospho",
  schizont_copur = "Schizont copur",
  mc4404 = "eIF5A pulldown",
  bioid = "CLASP1 BioID"
)
host_family_role <- c(
  pxd010713_ev = "Independent",
  pxd056906_ev = "Independent",
  dual_phospho = "Independent",
  schizont_copur = "Proximal",
  mc4404 = "Proximal",
  bioid = "Proximal"
)

has_source <- function(col, key) {
  str_detect(paste0(";", coalesce(as.character(col), ""), ";"), paste0(";", key, ";"))
}

display_core_catalog <- function(read_fig, p) {
  statistics <- read_fig("04_panel_d_host_statistics.tsv")
  if (!"gene_label" %in% names(statistics) && "display_label" %in% names(statistics)) {
    statistics$gene_label <- statistics$display_label
  }
  nodes <- read_fig("06_panel_g_host_network_nodes.tsv")
  rec_candidates <- unique(c(
    p("02_processed_inputs", "02_evidence", "01_host_recurrence.tsv"),
    file.path(getwd(), "02_processed_inputs", "02_evidence", "01_host_recurrence.tsv")
  ))
  rec_path <- rec_candidates[file.exists(rec_candidates)][1]
  rec <- if (length(rec_path) && !is.na(rec_path) && file.exists(rec_path)) {
    read_tsv(rec_path, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
  } else {
    NULL
  }

  out <- statistics %>%
    mutate(
      gene_label = pretty_gene(gene_label),
      independent_n = as.integer(coalesce(independent_evidence_count, 0))
    ) %>%
    left_join(
      nodes %>% transmute(node_id, degree = as.numeric(degree), net_sources = as.character(sources)),
      by = "node_id"
    )

  if (!is.null(rec) && all(c("string_id", "sources") %in% names(rec))) {
    rec_small <- rec %>% transmute(
      node_id = string_id,
      rec_sources = as.character(sources),
      rec_proximal_n = as.integer(coalesce(interface_proximal_evidence_count, 0)),
      rec_proximal_sources = as.character(interface_proximal_sources)
    )
    out <- out %>% left_join(rec_small, by = "node_id")
  } else {
    out$rec_sources <- NA_character_
    out$rec_proximal_n <- NA_integer_
    out$rec_proximal_sources <- NA_character_
  }

  out %>%
    mutate(
      all_sources = coalesce(na_if(rec_sources, ""), na_if(net_sources, ""), as.character(independent_sources), ""),
      proximal_n = coalesce(
        rec_proximal_n,
        pmax(as.integer(coalesce(evidence_count, 0)) - independent_n, 0L)
      ),
      degree = coalesce(as.numeric(degree), 0),
      axis_end = gene_label %in% axis_genes,
      independent_score = pmin(independent_n, 3) / 3,
      proximal_score = pmin(proximal_n, 3) / 3,
      degree_score = percent_rank(degree),
      axis_score = as.numeric(axis_end)
    ) %>%
    arrange(desc(axis_end), desc(independent_n), desc(proximal_n), desc(degree), gene_label) %>%
    mutate(gene_label = factor(gene_label, levels = rev(unique(gene_label))))
}

find_proj_file <- function(p, ...) {
  rel <- file.path(...)
  cands <- unique(c(do.call(p, list(...)), file.path(getwd(), rel)))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[[1]] else NA_character_
}

geo_study_levels <- c("GSE250088_Cheeseman", "GSE250088_Villares", "GSE119138")
f1d_study_levels <- c(host_family_levels, geo_study_levels)

read_gse250088_dea <- function(path) {
  read_tsv(
    path,
    col_names = c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
    skip = 1,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(
      gene_key = str_to_upper(str_trim(as.character(gene))),
      log2FoldChange = as.numeric(log2FoldChange),
      best_fdr = as.numeric(padj)
    ) %>%
    filter(nzchar(gene_key), is.finite(log2FoldChange), is.finite(best_fdr)) %>%
    group_by(gene_key) %>%
    slice_min(best_fdr, n = 1, with_ties = FALSE) %>%
    ungroup()
}

read_gse119138_dea <- function(path) {
  raw <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  names(raw) <- str_replace_all(str_replace_all(names(raw), "[^A-Za-z0-9]+", "_"), "^_+|_+$", "")
  gene_col <- intersect(c("Gene_Name", "GeneName", "gene_name"), names(raw))
  fdr_col <- intersect(c("Q_value", "Qvalue", "padj", "FDR"), names(raw))
  lfc_col <- intersect(c("log2_fold_change", "log2_fold_change_", "log2FoldChange"), names(raw))
  dir_col <- intersect(c("Regulation", "regulation"), names(raw))
  if (!length(gene_col) || !length(fdr_col)) {
    return(tibble(gene_key = character(), log2FoldChange = numeric(), best_fdr = numeric()))
  }
  lfc <- if (length(lfc_col)) {
    parsed <- suppressWarnings(as.numeric(str_replace(as.character(raw[[lfc_col[[1]]]]), regex("inf", ignore_case = TRUE), "Inf")))
    parsed
  } else {
    rep(NA_real_, nrow(raw))
  }
  if (length(dir_col)) {
    signed <- case_when(
      str_to_upper(coalesce(as.character(raw[[dir_col[[1]]]]), "")) == "UP" ~ 1,
      str_to_upper(coalesce(as.character(raw[[dir_col[[1]]]]), "")) == "DOWN" ~ -1,
      TRUE ~ NA_real_
    )
    lfc <- coalesce(lfc, signed)
  }
  tibble(
    gene_key = str_to_upper(str_trim(as.character(raw[[gene_col[[1]]]]))),
    log2FoldChange = lfc,
    best_fdr = suppressWarnings(as.numeric(raw[[fdr_col[[1]]]]))
  ) %>%
    filter(nzchar(gene_key), is.finite(log2FoldChange), is.finite(best_fdr)) %>%
    group_by(gene_key) %>%
    slice_min(best_fdr, n = 1, with_ties = FALSE) %>%
    ungroup()
}

dea_to_votes <- function(dea, study_family, gene_keys) {
  if (is.null(dea) || !nrow(dea)) {
    return(tibble(gene_key = character(), study_family = character(), family_direction = character(), best_fdr = numeric()))
  }
  dea %>%
    filter(gene_key %in% gene_keys, is.finite(log2FoldChange)) %>%
    transmute(
      gene_key,
      study_family = study_family,
      family_direction = if_else(log2FoldChange >= 0, "up", "down"),
      best_fdr
    )
}

display_core_recurrence_gold <- function(core, p) {
  gene_keys <- str_to_upper(as.character(core$gene_label))
  proteome <- bind_rows(lapply(host_family_levels, function(fam) {
    src <- if (fam %in% independent_families) core$independent_sources else core$all_sources
    tibble(
      gene_label = core$gene_label,
      study_family = fam,
      family_direction = if_else(has_source(src, fam), "up", NA_character_),
      best_fdr = NA_real_
    )
  }))
  geo_specs <- list(
    list(id = "GSE250088_Cheeseman", kind = "gse250088", parts = c("01_data", "01_raw", "geo_transcriptomes", "GSE250088", "GSE250088_dea_Cheeseman_BL3_TBL3.tsv.gz")),
    list(id = "GSE250088_Villares", kind = "gse250088", parts = c("01_data", "01_raw", "geo_transcriptomes", "GSE250088", "GSE250088_dea_Villares_BL3_TBL3.tsv.gz")),
    list(id = "GSE119138", kind = "gse119138", parts = c("01_data", "01_raw", "geo_transcriptomes", "GSE119138", "GSE119138_Bos_taurus_processed_data_file.txt.gz"))
  )
  geo <- bind_rows(lapply(geo_specs, function(spec) {
    path <- do.call(find_proj_file, c(list(p), as.list(spec$parts)))
    if (is.na(path)) return(NULL)
    dea <- if (spec$kind == "gse250088") read_gse250088_dea(path) else read_gse119138_dea(path)
    dea_to_votes(dea, spec$id, gene_keys)
  }))
  if (nrow(geo)) {
    geo <- geo %>%
      mutate(gene_label = core$gene_label[match(gene_key, gene_keys)]) %>%
      filter(!is.na(gene_label)) %>%
      select(gene_label, study_family, family_direction, best_fdr)
  } else {
    geo <- proteome[0, ]
  }
  bind_rows(proteome, geo) %>%
    complete(gene_label = core$gene_label, study_family = f1d_study_levels, fill = list(family_direction = NA_character_, best_fdr = NA_real_)) %>%
    mutate(
      gene_label = factor(gene_label, levels = levels(core$gene_label)),
      study_family = factor(study_family, f1d_study_levels),
      vote = case_when(family_direction == "up" ~ "Up", family_direction == "down" ~ "Down", TRUE ~ "No vote"),
      vote = factor(vote, c("Up", "Down", "No vote")),
      sig = case_when(
        vote == "No vote" ~ 0,
        is.na(best_fdr) ~ 2.2,
        best_fdr <= 0 ~ 6,
        TRUE ~ pmin(-log10(best_fdr), 6)
      )
    )
}

display_core_statistics_long <- function(core) {
  bind_rows(
    core %>% transmute(gene_label, metric = "Independent", value = independent_score),
    core %>% transmute(gene_label, metric = "Proximal", value = proximal_score),
    core %>% transmute(gene_label, metric = "Degree", value = degree_score),
    core %>% transmute(gene_label, metric = "Axis", value = axis_score)
  ) %>%
    mutate(metric = factor(metric, c("Independent", "Proximal", "Degree", "Axis")))
}
