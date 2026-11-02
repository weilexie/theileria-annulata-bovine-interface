#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(igraph)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path("01_scripts", "02_build_figure_data.R"), winslash = "/", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
p <- function(...) file.path(root, ...)
read_tab <- function(...) {
  read_tsv(p(...), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"),
           col_types = cols(.default = col_character()))
}
read_net <- function(...) {
  read_tsv(p(...), col_names = c("source", "target"), show_col_types = FALSE, progress = FALSE,
           col_types = cols(.default = col_character()))
}
flag <- function(x) {
  if (is.logical(x)) replace(x, is.na(x), FALSE)
  else str_to_lower(coalesce(as.character(x), "")) %in% c("true", "t", "1", "yes", "y")
}
num <- function(x, default = 0) {
  z <- suppressWarnings(as.numeric(x))
  z[!is.finite(z)] <- default
  z
}
chr <- function(x) replace(as.character(x), is.na(x), "")
clean_plot_text <- function(x) {
  z <- coalesce(as.character(x), "")
  z <- str_replace_all(z, "TA_15795\\?", "TaPin1")
  z <- str_replace_all(z, "（[^）]*）", "")
  z <- str_replace_all(z, "\\([^)]*正文[^)]*\\)", "")
  z <- str_replace_all(z, "\\s+", " ")
  str_trim(z)
}
split_ids <- function(x) {
  z <- chr(x)
  z <- z[nzchar(z)]
  if (!length(z)) return(character())
  unique(unlist(str_split(z, ";|\\|"), use.names = FALSE))
}
ensure_cols <- function(x, cols) {
  for (nm in cols) if (!nm %in% names(x)) x[[nm]] <- NA
  extras <- setdiff(names(x), cols)
  x[, c(intersect(cols, names(x)), extras), drop = FALSE]
}
first_known <- function(x, fallback = NA_character_) {
  known <- unique(na.omit(as.character(x)))
  known <- known[nzchar(known)]
  if (length(known)) known[[1]] else fallback
}

layout_network_by_module <- function(nodes, edges, seed = 1L) {
  z <- nodes %>%
    mutate(
      layout_module = coalesce(as.character(module_id), "Unassigned"),
      degree = num(degree),
      strength = num(strength)
    )
  modules <- z %>% count(layout_module, name = "module_n") %>%
    arrange(desc(module_n), layout_module) %>%
    mutate(
      module_index = row_number(),
      grid_n = ceiling(sqrt(n())),
      module_x = ((module_index - 1) %% grid_n) * 3.2,
      module_y = -floor((module_index - 1) / grid_n) * 3.2
    )
  set.seed(seed)
  z <- z %>% left_join(modules, by = "layout_module") %>%
    group_by(layout_module) %>%
    arrange(desc(degree), desc(strength), node_id, .by_group = TRUE) %>%
    mutate(
      within_rank = row_number(),
      within_n = n(),
      angle = within_rank * pi * (3 - sqrt(5)),
      radius = 1.25 * sqrt((within_rank - 1) / pmax(within_n - 1, 1)),
      x = module_x + radius * cos(angle),
      y = module_y + radius * sin(angle),
      hub = within_rank <= pmax(1L, ceiling(0.015 * within_n))
    ) %>% ungroup()
  e <- edges %>%
    left_join(z %>% select(source = node_id, x_from = x, y_from = y), by = "source") %>%
    left_join(z %>% select(target = node_id, x_to = x, y_to = y), by = "target") %>%
    filter(is.finite(x_from), is.finite(y_from), is.finite(x_to), is.finite(y_to))
  list(nodes = z, edges = e)
}

degree_from_edges <- function(edges) {
  bind_rows(
    edges %>% transmute(node_id = source),
    edges %>% transmute(node_id = target)
  ) %>% count(node_id, name = "degree")
}

family_display_map <- c(
  microtubule_docking = "Microtubule docking",
  metabolic_proteostasis = "Metabolic / proteostasis",
  trafficking_signaling = "Trafficking / signaling",
  protease_cullin = "Protease / cullin"
)
locked_axes <- tribble(
  ~biological_axis_id, ~parasite_label, ~host_label, ~family_display, ~context_only,
  "p104--MAPRE1", "p104", "MAPRE1", "Microtubule docking", FALSE,
  "TaSP--tubulin", "TaSP", "tubulin", "Microtubule docking", FALSE,
  "p104--CLASP1", "p104", "CLASP1", "Microtubule docking (context only)", TRUE,
  "TaPin1--PKM", "TaPin1", "PKM", "Metabolic / proteostasis", FALSE,
  "TaPin1--FBXW7", "TaPin1", "FBXW7", "Metabolic / proteostasis", FALSE,
  "TaPHB1--RUVBL1", "TaPHB1", "RUVBL1", "Metabolic / proteostasis", FALSE,
  "TA20980--CD2AP", "TA20980", "CD2AP", "Trafficking / signaling", FALSE,
  "TA11965--TFG", "TA11965", "TFG", "Trafficking / signaling", FALSE,
  "TaCyp1--MED21", "TaCyp1", "MED21", "Trafficking / signaling", FALSE,
  "TA05560--RBM39", "TA05560", "RBM39", "Trafficking / signaling", FALSE,
  "TaCP--PPP4C", "TaCP", "PPP4C", "Protease / cullin", FALSE,
  "TaCP--CRBN", "TaCP", "CRBN", "Protease / cullin", FALSE
)

gold_fig1_a <- c("study_family", "dataset_n", "effect_rows", "omics_type", "biological_context",
                 "integration_role", "analyzable", "side", "primary_family", "omics_group")
gold_fig1_c <- c("node_id", "study_family", "family_direction", "family_effect",
                 "family_standardized_effect", "best_fdr", "gene_label", "plot_order",
                 "module_label", "shadow_high_confidence_member")
gold_fig1_d <- c("node_id", "gene_label", "plot_order", "module_label", "shadow_high_confidence_member",
                 "support_fraction", "agreement", "agreement_low", "agreement_high", "posterior",
                 "prevalence", "prevalence_low", "prevalence_high", "grid_frequency",
                 "cross_omics_support", "median_standardized_effect")
gold_fig1_e <- c("module_id", "scaffold_n", "core_n", "module_label", "module_biological_label",
                 "module_theme_signature", "module_size", "measured_n", "context_n",
                 "primary_evidence_n", "replicated_primary_n", "support_only_n",
                 "significant_host_n", "median_abs_effect", "top_nodes", "total_strength",
                 "total_seed_weight", "study_family_n")
gold_fig1_g_nodes <- c(
  "node_id", "gene_symbol", "display_label", "protein_name", "layer", "node_type",
  "degree", "strength", "betweenness", "eigenvector", "module_id", "module_label",
  "module_biological_label", "module_theme_signature", "module_status",
  "strict_member", "relaxed_member", "network_role", "layout_module", "module_n",
  "module_index", "grid_n", "module_x", "module_y", "within_rank", "within_n",
  "angle", "radius", "x", "y", "hub"
)
gold_fig1_g_edges <- c(
  "source", "target", "combined_score", "string_source_ids", "string_target_ids",
  "edge_scope", "edge_type", "weight", "confidence", "evidence_source",
  "species_basis", "network_basis", "x_from", "y_from", "x_to", "y_to"
)
gold_fig2_a <- c("study_family", "dataset_n", "effect_rows", "omics_type", "biological_context",
                 "integration_role", "analyzable", "side", "dataset_n_evidence", "evidence_rows",
                 "voting_family", "omics_group")
gold_fig2_e <- c("module_id", "scaffold_n", "core_n", "module_label", "module_biological_label",
                 "dominant_function_class", "module_function_signature", "module_size",
                 "measured_n", "context_n", "string_mapped_n", "primary_evidence_n",
                 "support_evidence_n", "median_abs_surface_effect", "top_nodes", "total_strength")
gold_fig2_g_nodes <- c(
  "node_id", "display_label", "canonical_gene", "layer", "node_type", "degree", "strength",
  "betweenness", "eigenvector", "module_id", "module_label", "module_biological_label",
  "module_function_signature", "module_status", "scaffold_tier", "network_role",
  "layout_module", "module_n", "module_index", "grid_n", "module_x", "module_y",
  "within_rank", "within_n", "angle", "radius", "x", "y", "hub"
)
gold_fig2_g_edges <- c(
  "source", "target", "edge_scope", "edge_type", "weight", "confidence", "evidence_source",
  "species_basis", "network_basis", "source_rank", "target_rank", "mutual_topk",
  "x_from", "y_from", "x_to", "y_to"
)
gold_fig3_nodes <- c("module_id", "module_label", "module_biological_label",
                     "module_theme_signature", "scaffold_n", "core_n", "side")
gold_fig3_bridges <- c("host_network_module", "pathogen_network_module", "bridge_weight", "axis_n")
gold_fig3_rwr <- c("node_id", "layer", "display_label", "module_label", "degree",
                   "strict_member", "relaxed_member", "direct_evidence", "opposite_signal", "recruited")
gold_fig3_prox <- c(
  "host_module", "pathogen_module", "host_size", "pathogen_size", "host_reachable_n",
  "pathogen_reachable_n", "host_reachable_fraction", "pathogen_reachable_fraction",
  "closest_distance", "host_to_path_distance", "path_to_host_distance", "null_mean",
  "null_sd", "null_q025", "null_q975", "proximity_z", "empirical_p", "module_family",
  "proximity_fdr", "control_type", "proximity_support"
)
gold_fig4_union_nodes <- c("node_id", "selected_runs", "selection_frequency", "layer",
                           "display_label", "module_id", "module_label", "strict_member",
                           "relaxed_member", "side_seed_percentile")
gold_fig4_burden <- c("layer", "module_id", "module_label", "selected_nodes",
                      "mean_selection_frequency", "endpoint_n", "hub_n", "total_strength")
gold_fig4_persist <- c("node_id", "selected_runs", "selection_frequency", "layer", "display_label",
                       "module_id", "module_label", "strict_member", "relaxed_member",
                       "side_seed_percentile", "network_degree", "network_strength",
                       "forest_degree", "weighted_degree", "endpoint", "hub_score", "hub")
gold_fig4_axis <- c(
  "axis_id", "interaction_name", "pathogen_gene", "host_gene", "pathogen_nodes", "host_nodes",
  "interaction_type", "evidence_strength", "bridge_evidence_tier", "main_algorithm_bridge_gate",
  "evidence_source", "evidence_reference", "host_species_evidence", "species_adjusted_confidence",
  "component_rows", "host_measured_gate", "pathogen_measured_gate", "host_recurrence_gate",
  "pathogen_support_gate", "localization_gate", "strain_gate", "species_gate", "direct_physical_edge",
  "pathogen_cross_gain_percentile", "host_cross_gain_percentile", "rwr_endpoint_min_percentile",
  "rwr_support", "pathogen_pcsf_frequency", "host_pcsf_frequency", "pcsf_endpoint_min_frequency",
  "pcsf_edge_frequency", "pcsf_support", "pcsf_successful_runs", "pcsf_total_runs",
  "pcsf_success_fraction", "proximity_support", "proximity_host_module", "proximity_pathogen_module",
  "proximity_module_family", "proximity_distance", "proximity_z", "proximity_p", "proximity_fdr",
  "pcsf_component_path_frequency", "pcsf_all_component_paths_frequency",
  "pcsf_mean_connected_component_fraction", "pcsf_component_support", "strict_evidence_gate",
  "relaxed_evidence_gate", "algorithm_count", "consensus_rule_pass", "consensus_tier",
  "evidence_class", "known_anchor", "discovery_semantics", "primary_bridge_axis_score",
  "all_curated_bridge_axis_score", "primary_bridge_axis_rank", "all_curated_bridge_axis_rank",
  "all_curated_score_shift", "all_curated_rank_shift", "parameter_run_n", "rwr_top10_frequency",
  "rwr_support_frequency", "rwr_median_rank", "rwr_rank_iqr", "rwr_min_parameter_score",
  "rwr_max_parameter_score", "threshold_run_n", "threshold_consensus_frequency",
  "threshold_strict_selection_frequency", "threshold_rwr_support_frequency",
  "threshold_proximity_support_frequency", "threshold_pcsf_support_frequency",
  "proximity_effect", "pcsf_continuous_strength", "rwr_strength_percentile",
  "proximity_strength_percentile", "pcsf_strength_percentile", "continuous_algorithm_score",
  "stability_score", "aggregate_robustness_score", "aggregate_robustness_rank", "robustness_tier",
  "algorithm_combination", "final_tier", "host_network_module", "pathogen_network_module",
  "host_program", "pathogen_program", "interface_score", "confidence", "persistence",
  "axis_label", "rank"
)
gold_fig5_nodes <- c(
  "node_id", "x0", "y0", "layer", "node_type", "display_label", "canonical_gene", "protein_name",
  "biological_theme", "module_id", "module_label", "module_signature", "seed_weight",
  "effect_value", "effect_direction", "p_value", "fdr", "mapped_to_reference_network",
  "degree", "strength", "source_dataset", "study_families", "n_study_families",
  "omics_modalities", "evidence_role", "source_entity_id", "endpoint_status",
  "species_provenance", "network_origin", "side_seed_raw", "side_seed_percentile",
  "strict_member", "relaxed_member", "joint_rwr_percentile", "cross_layer_gain",
  "cross_layer_gain_percentile", "rwr_top1pct_parameter_frequency",
  "rwr_median_parameter_percentile", "pcsf_selected_runs", "pcsf_node_frequency",
  "core_axis_endpoint", "rwr_parameter_stable", "pcsf_parameter_stable", "connector_node",
  "scaffold_role", "x", "y", "label", "hub"
)
gold_fig5_edges <- c(
  "source", "target", "pair_a", "pair_b", "edge_scope", "edge_type", "weight", "confidence",
  "evidence_source", "species_basis", "network_basis", "interaction_name", "interaction_type",
  "evidence_reference", "host_species_evidence", "species_adjusted_confidence", "axis_id",
  "pcsf_edge_selected_runs", "pcsf_edge_frequency", "edge_semantics", "known_cross_layer_anchor",
  "algorithm_derived_cross_layer_edge", "mechanism_class", "x_from", "y_from", "x_to", "y_to"
)
gold_fig5_mod_nodes <- c("module_key", "layer", "module_label", "scaffold_nodes", "core_endpoint_n",
                         "hub_n", "degree_sum", "strength_sum", "odds_ratio", "ci_low", "ci_high", "fdr")
gold_fig5_mod_edges <- c("from_module", "to_module", "edge_n", "weight_sum", "cross_axis_n")
gold_fig5_enrich <- c("layer", "module_label", "selected_n", "selected_total", "background_n",
                      "background_total", "odds_ratio", "ci_low", "ci_high", "p_value", "fdr")
gold_fig5_mech <- c("axis_id", "interaction_name", "pathogen_gene", "host_gene", "final_tier",
                    "aggregate_robustness_score", "host_module", "pathogen_module",
                    "direct_physical_edge", "interaction_type", "edge_type",
                    "species_adjusted_confidence", "mechanism_class")
gold_fig5_coup <- c("host_module", "pathogen_module", "axis_n", "core_axis_n", "mean_score")

figure_data_final <- p("04_plot_data")
stage <- tempfile(pattern = "theileria_figure_data_")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
write_fig <- function(x, fig, file) {
  out <- file.path(stage, fig)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write_tsv(x, file.path(out, file), na = "")
}
host_layout_mode <- "full_700"
contract_notes <- character()

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
host <- read_tab("02_processed_inputs", "02_evidence", "01_host_recurrence.tsv") %>%
  mutate(
    node_id = string_id,
    evidence_count = num(evidence_count),
    independent_evidence_count = num(independent_evidence_count),
    interface_proximal_evidence_count = num(interface_proximal_evidence_count),
    in_700_network = flag(in_700_network)
  )
para <- read_tab("02_processed_inputs", "02_evidence", "02_parasite_accessible.tsv") %>%
  mutate(
    node_id = string_id,
    evidence_count = num(evidence_count),
    accessible_annot = flag(accessible_annot),
    in_para_network = flag(in_para_network),
    preferred_name = clean_plot_text(preferred_name)
  )
annot <- read_tab("02_processed_inputs", "02_evidence", "03_string_annotation_cache.tsv") %>%
  mutate(
    preferred_name = if_else(preferred_name %in% c("9913", "5874"), NA_character_, preferred_name),
    protein_name = annotation
  )
inventory <- read_tab("02_processed_inputs", "04_inventory", "01_dataset_inventory.tsv")
net_modules <- read_tab("03_algorithm_tables", "24_network_module_membership.tsv") %>%
  mutate(module_size = num(module_size))
fun_modules <- read_tab("03_algorithm_tables", "14_functional_module_membership.tsv")
module_names <- read_tab("05_extra_tables", "01_module_display_names.tsv") %>%
  mutate(module_size = num(module_size), top_fdr = num(top_fdr))
enrichment <- read_tab("03_algorithm_tables", "23_network_module_enrichment.tsv") %>%
  mutate(module_size = num(module_size), rank = num(rank), number_of_genes = num(number_of_genes), fdr = num(fdr))
axis_cons <- read_tab("03_algorithm_tables", "06_biological_axis_consensus.tsv")
axis_map <- read_tab("03_algorithm_tables", "07_biological_axis_map.tsv")
axis_records <- read_tab("03_algorithm_tables", "18_interface_axis_algorithm_support.tsv")
rwr <- read_tab("03_algorithm_tables", "36_rwr_cross_layer_gain.tsv") %>%
  mutate(
    degree = num(degree),
    evidence_n = num(evidence_n),
    accessible_annot = flag(accessible_annot),
    strict_member_all_source = flag(strict_member),
    relaxed_member_all_source = flag(relaxed_member),
    side_seed_raw = num(side_seed_raw),
    side_seed_percentile = num(side_seed_percentile, default = NA),
    joint_rwr = num(joint_rwr, default = NA),
    opposite_side_rwr = num(opposite_side_rwr, default = NA),
    joint_rwr_percentile = num(joint_rwr_percentile, default = NA),
    opposite_side_rwr_percentile = num(opposite_side_rwr_percentile, default = NA),
    same_side_diffusion_percentile = num(same_side_diffusion_percentile, default = NA),
    cross_layer_gain_percentile = num(cross_layer_gain_percentile, default = NA)
  )
proximity <- read_tab("03_algorithm_tables", "22_module_proximity_all_pairs.tsv") %>%
  mutate(
    host_size = num(host_size), pathogen_size = num(pathogen_size),
    host_reachable_n = num(host_reachable_n), pathogen_reachable_n = num(pathogen_reachable_n),
    closest_distance = num(closest_distance, default = NA),
    host_to_path_distance = num(host_to_path_distance, default = NA),
    path_to_host_distance = num(path_to_host_distance, default = NA),
    null_mean = num(null_mean, default = NA), null_sd = num(null_sd, default = NA),
    proximity_z = num(proximity_z, default = NA), empirical_p = num(empirical_p, default = NA),
    proximity_fdr = num(proximity_fdr, default = NA), proximity_support = flag(proximity_support)
  )
program_groups <- read_tab("03_algorithm_tables", "21_module_program_groups.tsv")
pcsf_runs <- read_tab("03_algorithm_tables", "35_pcsf_parameter_stability.tsv") %>%
  mutate(w = num(w), b = num(b), mu = num(mu, default = NA), success = flag(success),
         node_n = num(node_n), edge_n = num(edge_n))
pcsf_nodes <- read_tab("03_algorithm_tables", "34_pcsf_node_selection_frequency.tsv") %>%
  mutate(
    selected_runs = num(selected_runs), selection_frequency = num(selection_frequency),
    host_evidence_n = num(host_evidence_n), pathogen_evidence_n = num(pathogen_evidence_n),
    evidence_n = num(evidence_n), accessible_annot = flag(accessible_annot),
    strict_member = flag(strict_member), relaxed_member = flag(relaxed_member),
    seed_raw = num(seed_raw), side_seed_percentile = num(side_seed_percentile, default = NA),
    degree = num(degree)
  )
pcsf_edges <- read_tab("03_algorithm_tables", "32_pcsf_edge_selection_frequency.tsv") %>%
  mutate(selected_runs = num(selected_runs), selection_frequency = num(selection_frequency))
pcsf_node_runs <- read_tab("03_algorithm_tables", "33_pcsf_node_runs.tsv")
local_nodes <- read_tab("03_algorithm_tables", "28_parasite_local_program_nodes.tsv")
axis_coupling <- read_tab("03_algorithm_tables", "25_parasite_axis_program_coupling.tsv") %>%
  mutate(axis_n = num(axis_n), strict_axis_n = num(strict_axis_n),
         mean_rwr_endpoint_min = num(mean_rwr_endpoint_min, default = NA),
         mean_threshold_consensus_frequency = num(mean_threshold_consensus_frequency, default = NA))
host_indep_list <- read_tab("03_algorithm_tables", "16_host_independent_scaffold.tsv")
host_map <- read_tab("02_processed_inputs", "03_id_map", "01_host_string_id_map.tsv")

host_edges_raw <- read_net("02_processed_inputs", "01_networks", "01_host_string_network_700.tsv") %>%
  filter(source != target) %>% distinct(source, target)
para_edges_raw <- read_net("02_processed_inputs", "01_networks", "02_parasite_string_network_400.tsv") %>%
  filter(source != target) %>% distinct(source, target)

# ---------------------------------------------------------------------------
# Display names
# ---------------------------------------------------------------------------
enrich_top <- enrichment %>%
  filter(num(rank) == 1) %>%
  group_by(module_id) %>%
  slice_min(fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(module_id, top_description = description, top_term = term, top_category = category, top_fdr = fdr)

module_ann <- net_modules %>%
  distinct(module_id, layer, module_size) %>%
  left_join(module_names %>% select(module_id, display_name, naming_note), by = "module_id") %>%
  left_join(enrich_top, by = "module_id") %>%
  mutate(
    module_label = coalesce(display_name, module_id),
    module_biological_label = coalesce(display_name, top_description, module_id),
    module_theme_signature = coalesce(top_term, module_id)
  )

fun_ann <- fun_modules %>%
  count(module_id, layer, name = "module_size") %>%
  mutate(
    module_label = module_id,
    module_biological_label = str_replace(module_id, "^[HP]_FUN_", ""),
    module_theme_signature = module_id
  )

module_lookup <- bind_rows(
  module_ann %>% select(module_id, layer, module_size, module_label, module_biological_label, module_theme_signature),
  fun_ann %>% select(module_id, layer, module_size, module_label, module_biological_label, module_theme_signature)
) %>% distinct(module_id, .keep_all = TRUE)

axis_cons <- axis_cons %>%
  mutate(
    retained = flag(retained),
    manuscript_demoted = flag(manuscript_demoted),
    host_independent_recurrence_gate = flag(host_independent_recurrence_gate),
    main_algorithm_bridge = flag(main_algorithm_bridge),
    direct_physical_edge = flag(direct_physical_edge),
    rwr_support = flag(rwr_support),
    pcsf_support = flag(pcsf_support),
    proximity_support = flag(proximity_support),
    consensus_rule_pass = flag(consensus_rule_pass),
    independent_evidence_count = num(independent_evidence_count),
    species_adjusted_confidence = num(species_adjusted_confidence),
    rwr_endpoint_min_percentile = num(rwr_endpoint_min_percentile, default = NA),
    pathogen_cross_gain_percentile = num(pathogen_cross_gain_percentile, default = NA),
    host_cross_gain_percentile = num(host_cross_gain_percentile, default = NA),
    pathogen_pcsf_frequency = num(pathogen_pcsf_frequency),
    host_pcsf_frequency = num(host_pcsf_frequency),
    pcsf_edge_frequency = num(pcsf_edge_frequency),
    pcsf_endpoint_min_frequency = num(pcsf_endpoint_min_frequency),
    proximity_distance = num(proximity_distance, default = NA),
    proximity_z = num(proximity_z, default = NA),
    proximity_p = num(proximity_p, default = NA),
    proximity_fdr = num(proximity_fdr, default = NA),
    algorithm_count = num(algorithm_count),
    threshold_consensus_frequency = num(threshold_consensus_frequency, default = NA)
  ) %>%
  left_join(locked_axes, by = "biological_axis_id") %>%
  mutate(
    parasite_label = coalesce(parasite_label, clean_plot_text(parasite_factor_name)),
    host_label = coalesce(host_label, clean_plot_text(host_target_name)),
    family_display = coalesce(family_display, unname(family_display_map[axis_family]), axis_family),
    interaction_name_clean = biological_axis_id,
    final_tier = recode(consensus_tier, strict = "Core", relaxed = "Extended", context_only = "Context", .default = consensus_tier)
  )

axis_records <- axis_records %>%
  mutate(
    rwr_support = flag(rwr_support), pcsf_support = flag(pcsf_support),
    proximity_support = flag(proximity_support),
    direct_physical_edge = flag(direct_physical_edge),
    main_algorithm_bridge = flag(main_algorithm_bridge),
    species_adjusted_confidence = num(species_adjusted_confidence)
  )

component_nodes <- axis_map %>%
  left_join(
    axis_records %>% select(axis_id, host_node, pathogen_node),
    by = "axis_id"
  ) %>%
  left_join(locked_axes, by = "biological_axis_id") %>%
  mutate(
    parasite_label = coalesce(parasite_label, clean_plot_text(parasite_factor_name)),
    host_label = case_when(
      biological_axis_id == "TaSP--tubulin" ~ clean_plot_text(host_target_name),
      TRUE ~ coalesce(host_label, clean_plot_text(host_target_name))
    )
  )

locked_node_labels <- bind_rows(
  component_nodes %>% transmute(node_id = pathogen_node, display_label = parasite_label),
  component_nodes %>% transmute(node_id = host_node, display_label = host_label)
) %>%
  filter(!is.na(node_id), nzchar(node_id), !is.na(display_label), nzchar(display_label)) %>%
  distinct(node_id, .keep_all = TRUE)

uniprot_fasta <- p("01_data", "01_raw", "PXD070746", "UP_Bos_taurus_tx9913_37504_231027.fasta")
fasta_gn <- tibble(accession = character(), gene_symbol = character())
if (file.exists(uniprot_fasta)) {
  hdr <- grep("^>", readLines(uniprot_fasta, warn = FALSE), value = TRUE)
  fasta_gn <- tibble(header = hdr) %>%
    transmute(
      accession = str_match(header, "^>..\\|([^|]+)\\|")[, 2],
      gene_symbol = str_match(header, " GN=([^ ]+)")[, 2]
    ) %>%
    filter(!is.na(accession), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
    distinct(accession, .keep_all = TRUE)
}
is_uniprot_acc <- function(x) {
  str_detect(x, "^A0A") |
    str_detect(x, "^[OPQ][0-9][A-Z0-9]{3}[0-9]$") |
    str_detect(x, "^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$")
}
host_symbol_from_map <- host_map %>%
  filter(flag(mapped), !is.na(string_id), nzchar(string_id)) %>%
  mutate(
    input_id = chr(input_id),
    acc = str_replace(input_id, "_BOVIN$", ""),
    acc = str_replace(acc, "^sp\\|", ""),
    acc = str_replace(acc, "^tr\\|", ""),
    acc = str_replace(acc, "\\|.*$", ""),
    gene_like = str_detect(input_id, "^[A-Za-z][A-Za-z0-9-]{1,10}$") &
      !str_detect(input_id, "^ENS|^NP_|^XP_|^YP_|^WP_|^NM_|^XM_") &
      !is_uniprot_acc(input_id)
  )
host_symbols <- bind_rows(
  host_symbol_from_map %>%
    inner_join(fasta_gn, by = c("acc" = "accession")) %>%
    transmute(node_id = string_id, gene_symbol, priority = 1L),
  host_symbol_from_map %>%
    filter(gene_like) %>%
    transmute(node_id = string_id, gene_symbol = input_id, priority = 2L)
) %>%
  filter(!is.na(gene_symbol), nzchar(gene_symbol), !is_uniprot_acc(gene_symbol)) %>%
  arrange(priority) %>%
  group_by(node_id) %>%
  summarise(gene_symbol = first_known(gene_symbol), .groups = "drop")

para_names <- bind_rows(
  para %>% transmute(node_id, gene_symbol = preferred_name),
  local_nodes %>% transmute(node_id, gene_symbol = clean_plot_text(preferred_name)),
  pcsf_nodes %>% filter(layer == "pathogen") %>% transmute(node_id, gene_symbol = clean_plot_text(pathogen_name))
) %>%
  filter(!is.na(gene_symbol), nzchar(gene_symbol), !gene_symbol %in% c("9913", "5874")) %>%
  group_by(node_id) %>%
  summarise(gene_symbol = first_known(gene_symbol), .groups = "drop")

name_lookup <- bind_rows(
  locked_node_labels %>% transmute(node_id, gene_symbol = display_label, priority = 1L),
  host_symbols %>% mutate(priority = 2L),
  para_names %>% mutate(priority = 3L)
) %>%
  filter(!is.na(node_id), nzchar(node_id), !is.na(gene_symbol), nzchar(gene_symbol)) %>%
  mutate(gene_symbol = recode(gene_symbol, PIN1 = "TaPin1", TACP = "TaCP", sp = "TaSP", TA08425 = "p104", .default = gene_symbol)) %>%
  filter(!is_uniprot_acc(gene_symbol)) %>%
  arrange(priority) %>%
  distinct(node_id, .keep_all = TRUE) %>%
  select(node_id, gene_symbol)

resolve_name <- function(ids, fallback = ids) {
  tibble(node_id = ids, fallback = fallback) %>%
    left_join(name_lookup, by = "node_id") %>%
    mutate(label = coalesce(gene_symbol, str_remove(fallback, "^(9913|5874)\\."))) %>%
    pull(label)
}

# ---------------------------------------------------------------------------
# Membership flags
# ---------------------------------------------------------------------------
independent_families <- c("pxd010713_ev", "pxd056906_ev", "dual_phospho")
proximal_families <- c("schizont_copur", "mc4404", "bioid")
parasite_families <- c("stage_proteome", "surface_phospho", "schizont_ref",
                       "turboid_export", "drug_target", "surface_bioid")

strict_axis <- axis_cons %>% filter(consensus_tier == "strict")
retained_axis <- axis_cons %>% filter(flag(retained), consensus_tier %in% c("strict", "relaxed"))
strict_host_ends <- unique(strict_axis$host_node)
strict_path_ends <- unique(strict_axis$pathogen_node)
retained_host_ends <- unique(c(retained_axis$host_node, component_nodes$host_node[component_nodes$biological_axis_id %in% retained_axis$biological_axis_id]))
retained_path_ends <- unique(c(retained_axis$pathogen_node, component_nodes$pathogen_node[component_nodes$biological_axis_id %in% retained_axis$biological_axis_id]))
local_core_ids <- unique(local_nodes$node_id)
pin1_id <- "5874.Q4UG71"
local_core_ids <- setdiff(local_core_ids, pin1_id)

host <- host %>%
  mutate(
    independent_scaffold = independent_evidence_count >= 2 & in_700_network,
    display_core = independent_evidence_count >= 3 | node_id %in% strict_host_ends,
    relaxed_member = independent_scaffold,
    strict_member = display_core,
    gene_symbol = resolve_name(node_id),
    display_label = gene_symbol
  )
stopifnot(sum(host$independent_scaffold) == 968L)
stopifnot(!any(str_detect(coalesce(host$display_label, ""), "TA_15795\\?|（")))

para <- para %>%
  mutate(
    accessible_in_network = accessible_annot & in_para_network,
    recurrent_core = layer == "recurrent_core",
    local_core = node_id %in% local_core_ids,
    gene_symbol = resolve_name(node_id, preferred_name),
    display_label = recode(gene_symbol, PIN1 = "TaPin1", TACP = "TaCP", sp = "TaSP", TA08425 = "p104", .default = gene_symbol)
  )
para$local_core[para$node_id == pin1_id] <- FALSE
stopifnot(sum(para$local_core) == 19L)
stopifnot(!isTRUE(para$local_core[match(pin1_id, para$node_id)]))

host_mod <- net_modules %>% filter(layer == "host") %>% select(node_id, module_id, module_size)
para_mod <- net_modules %>% filter(layer == "pathogen") %>% select(node_id, module_id, module_size)
annot_prot <- annot %>% select(node_id, protein_name)

# ---------------------------------------------------------------------------
# Figure 1
# ---------------------------------------------------------------------------
count_host_source <- function(key, col) {
  sum(str_detect(paste0(";", chr(host[[col]]), ";"), paste0(";", key, ";")))
}
host_study_meta <- tribble(
  ~study_family, ~omics_type, ~biological_context, ~integration_role, ~primary_family, ~omics_group,
  "pxd010713_ev", "ev_proteome", "TBL20 vs BL20 extracellular vesicles", "independent_host_state", TRUE, "Proteome",
  "pxd056906_ev", "ev_proteome", "Infected cattle serum extracellular vesicles", "independent_host_state", TRUE, "Proteome",
  "dual_phospho", "dual_phosphoproteome", "BL3 vs TA dual phosphoproteome protein coverage", "independent_host_state", TRUE, "Proteome",
  "schizont_copur", "proteome", "Schizont Bos-space co-purification", "interface_proximal", FALSE, "Proteome",
  "mc4404", "proteome", "eIF5A / MC4404 pulldown host layer", "interface_proximal", FALSE, "Proteome",
  "bioid", "proximity_proteome", "Schizont-surface CLASP1 BioID host neighbours", "interface_proximal", FALSE, "Proteome"
) %>%
  mutate(
    dataset_n = 1L,
    effect_rows = vapply(study_family, function(k) {
      col <- if (k %in% independent_families) "independent_sources" else "sources"
      as.integer(count_host_source(k, col))
    }, integer(1)),
    analyzable = TRUE,
    side = "host"
  )
write_fig(ensure_cols(host_study_meta, gold_fig1_a), "01_figure_1", "01_panel_a_host_studies.tsv")

n_mapped <- nrow(host)
n_any_vote <- sum(host$evidence_count >= 1)
n_indep_vote <- sum(host$independent_evidence_count >= 1)
n_scaffold <- sum(host$independent_scaffold)
n_display_core <- sum(host$display_core)
write_fig(tibble(
  stage = c(
    "Mapped proteins",
    "Independent vote",
    "Independent scaffold",
    "Display core"
  ),
  value = c(n_mapped, n_indep_vote, n_scaffold, n_display_core)
), "01_figure_1", "02_panel_b_host_attrition.tsv")

host_core <- host %>%
  filter(display_core) %>%
  left_join(host_mod, by = "node_id") %>%
  left_join(module_lookup %>% select(module_id, module_label), by = "module_id") %>%
  arrange(desc(independent_evidence_count), desc(node_id %in% strict_host_ends), display_label) %>%
  mutate(
    plot_order = row_number(),
    gene_label = display_label,
    shadow_high_confidence_member = independent_evidence_count >= 3
  )
host_vote_matrix <- lapply(independent_families, function(fam) {
  present <- str_detect(paste0(";", chr(host_core$independent_sources), ";"), paste0(";", fam, ";"))
  tibble(
    node_id = host_core$node_id,
    study_family = fam,
    family_direction = if_else(present, "present", NA_character_),
    family_effect = if_else(present, 1, NA_real_),
    family_standardized_effect = NA_real_,
    best_fdr = NA_real_
  )
}) %>% bind_rows() %>%
  left_join(host_core %>% select(node_id, gene_label, plot_order, module_label, shadow_high_confidence_member), by = "node_id")
write_fig(ensure_cols(host_vote_matrix, gold_fig1_c), "01_figure_1", "03_panel_c_host_recurrence.tsv")

host_stats <- host_core %>%
  transmute(
    node_id, gene_label, plot_order, module_label, shadow_high_confidence_member,
    support_fraction = independent_evidence_count / 3,
    agreement = NA_real_, agreement_low = NA_real_, agreement_high = NA_real_,
    posterior = NA_real_,
    prevalence = independent_evidence_count / 3,
    prevalence_low = NA_real_, prevalence_high = NA_real_,
    grid_frequency = NA_real_,
    cross_omics_support = str_detect(chr(independent_sources), "dual_phospho") &
      str_detect(chr(independent_sources), "pxd0"),
    median_standardized_effect = NA_real_,
    independent_evidence_count, evidence_count, independent_sources
  )
write_fig(ensure_cols(host_stats, gold_fig1_d), "01_figure_1", "04_panel_d_host_statistics.tsv")

host_with_mod <- host %>%
  left_join(host_mod %>% select(node_id, module_id, host_module_size = module_size), by = "node_id") %>%
  mutate(module_id = coalesce(module_id, "Unassigned")) %>%
  left_join(module_lookup %>% select(module_id, module_label, module_biological_label, module_theme_signature),
            by = "module_id")
host_module_burden <- host_with_mod %>%
  filter(independent_scaffold) %>%
  group_by(module_id) %>%
  summarise(
    scaffold_n = n(),
    core_n = sum(display_core),
    module_label = first_known(module_label, first(module_id)),
    module_biological_label = first_known(module_biological_label, first(module_id)),
    module_theme_signature = first_known(module_theme_signature, first(module_id)),
    module_size = max(num(host_module_size), na.rm = TRUE),
    measured_n = n(),
    context_n = 0L,
    primary_evidence_n = n(),
    replicated_primary_n = sum(independent_evidence_count >= 2),
    support_only_n = 0L,
    significant_host_n = sum(display_core),
    median_abs_effect = NA_real_,
    top_nodes = paste(head(display_label[order(-independent_evidence_count, display_label)], 8), collapse = ";"),
    total_strength = NA_real_,
    total_seed_weight = sum(independent_evidence_count),
    study_family_n = n_distinct(unlist(str_split(independent_sources, ";"))),
    .groups = "drop"
  ) %>%
  mutate(module_size = if_else(is.finite(module_size), module_size, as.numeric(scaffold_n))) %>%
  arrange(desc(core_n), desc(scaffold_n), module_id)
write_fig(ensure_cols(host_module_burden, gold_fig1_e), "01_figure_1", "05_panel_e_host_functions.tsv")

host_degree <- degree_from_edges(host_edges_raw)
host_graph_nodes <- sort(unique(c(host_edges_raw$source, host_edges_raw$target)))
host_node_tbl <- tibble(node_id = host_graph_nodes) %>%
  left_join(host %>% select(node_id, gene_symbol, display_label, independent_scaffold, display_core,
                            independent_evidence_count, evidence_count, in_700_network,
                            independent_sources, sources), by = "node_id") %>%
  left_join(host_mod, by = "node_id") %>%
  left_join(module_lookup %>% select(module_id, module_label, module_biological_label, module_theme_signature),
            by = "module_id") %>%
  left_join(host_degree, by = "node_id") %>%
  left_join(annot_prot, by = "node_id") %>%
  mutate(
    gene_symbol = coalesce(gene_symbol, resolve_name(node_id)),
    display_label = coalesce(display_label, gene_symbol),
    protein_name = protein_name,
    layer = "host",
    node_type = "measured",
    degree = coalesce(degree, 0L),
    strength = degree,
    betweenness = NA_real_,
    eigenvector = NA_real_,
    module_status = if_else(!is.na(module_id) & nzchar(module_id), "network_module", "unassigned"),
    strict_member = coalesce(display_core, FALSE),
    relaxed_member = coalesce(independent_scaffold, FALSE),
    network_role = case_when(
      strict_member ~ "Core",
      relaxed_member ~ "Recurrent scaffold",
      TRUE ~ "Network context"
    )
  )

host_edge_tbl <- host_edges_raw %>%
  mutate(
    combined_score = NA_real_,
    string_source_ids = source,
    string_target_ids = target,
    edge_scope = "intra_host",
    edge_type = "string_700",
    weight = 1,
    confidence = NA_real_,
    evidence_source = "STRING",
    species_basis = "Bos_taurus",
    network_basis = "host_string_network_700"
  )

n_scaffold_endpoints <- sum(host_node_tbl$relaxed_member | host_node_tbl$strict_member)
if (nrow(host_edge_tbl) > 380000L) {
  host_layout_mode <- "scaffold_induced"
  keep_ids <- host_node_tbl$node_id[host_node_tbl$relaxed_member | host_node_tbl$strict_member]
  scaf_edges <- host_edge_tbl %>% filter(source %in% keep_ids | target %in% keep_ids)
  ctx_edges <- host_edge_tbl %>% filter(!source %in% keep_ids, !target %in% keep_ids)
  set.seed(20260903)
  ctx_keep <- if (nrow(ctx_edges) > 25000L) slice_sample(ctx_edges, n = 25000L) else ctx_edges
  host_edge_tbl <- bind_rows(scaf_edges, ctx_keep)
  contract_notes <- c(contract_notes, "Host panel_g edges restricted to scaffold-incident edges plus 25k context edges.")
}

host_full_layout <- layout_network_by_module(host_node_tbl, host_edge_tbl, seed = 20260903L)
write_fig(ensure_cols(host_full_layout$nodes, gold_fig1_g_nodes), "01_figure_1", "06_panel_g_host_network_nodes.tsv")
write_fig(ensure_cols(host_full_layout$edges, gold_fig1_g_edges), "01_figure_1", "07_panel_g_host_network_edges.tsv")

# ---------------------------------------------------------------------------
# Figure 2
# ---------------------------------------------------------------------------
count_para_source <- function(key) {
  sum(str_detect(paste0(";", chr(para$sources), ";"), paste0(";", key, ";")))
}
para_study_meta <- tribble(
  ~study_family, ~omics_type, ~biological_context, ~integration_role, ~omics_group,
  "stage_proteome", "proteome", "Schizont vs merozoite stage proteome", "parasite_accessible_vote", "Proteome",
  "surface_phospho", "phosphoproteome", "Schizont surface phosphoproteome", "parasite_accessible_vote", "Proteome",
  "schizont_ref", "proteome", "Witschi 2013 schizont proteome reference", "parasite_accessible_vote", "Proteome",
  "turboid_export", "proximity_proteome", "TurboID secreted / exportome", "parasite_accessible_vote", "Proteome",
  "drug_target", "proteome", "MC4404 / eIF5A parasite proteins", "parasite_accessible_vote", "Proteome",
  "surface_bioid", "proximity_proteome", "CLASP1 BioID parasite neighbours", "parasite_accessible_vote", "Proteome"
) %>%
  mutate(
    dataset_n = 1L,
    evidence_rows = vapply(study_family, count_para_source, integer(1)),
    dataset_n_evidence = 1L,
    effect_rows = evidence_rows,
    analyzable = TRUE,
    side = "pathogen",
    voting_family = TRUE
  )
write_fig(ensure_cols(para_study_meta, gold_fig2_a), "02_figure_2", "01_panel_a_pathogen_studies.tsv")

write_fig(tibble(
  stage = c(
    "Parasite evidence proteins",
    "In parasite network",
    "Accessible annot",
    "Accessible in network",
    "Local core"
  ),
  value = c(
    nrow(para),
    sum(para$in_para_network),
    sum(para$accessible_annot),
    sum(para$accessible_in_network),
    sum(para$local_core)
  )
), "02_figure_2", "02_panel_b_pathogen_attrition.tsv")

para_degree <- degree_from_edges(para_edges_raw)
para_graph_nodes <- sort(unique(c(para_edges_raw$source, para_edges_raw$target)))
para_cent <- tibble(node_id = para_graph_nodes, betweenness = NA_real_, eigenvector = NA_real_)
if (length(para_graph_nodes) >= 2L && nrow(para_edges_raw) >= 1L) {
  g_para <- graph_from_data_frame(
    para_edges_raw %>% select(source, target),
    directed = FALSE,
    vertices = para_graph_nodes
  )
  para_cent$betweenness <- as.numeric(betweenness(g_para, normalized = TRUE))
  para_cent$eigenvector <- as.numeric(eigen_centrality(g_para, directed = FALSE)$vector)
}

para_node_tbl <- tibble(node_id = para_graph_nodes) %>%
  left_join(para %>% select(node_id, gene_symbol, display_label, preferred_name, annotation,
                            evidence_count, sources, accessible_annot, in_para_network,
                            layer, accessible_in_network, recurrent_core, local_core),
            by = "node_id") %>%
  left_join(para_mod, by = "node_id") %>%
  left_join(module_lookup %>% select(module_id, module_label, module_biological_label, module_theme_signature),
            by = "module_id") %>%
  left_join(local_nodes %>% select(node_id, parasite_program), by = "node_id") %>%
  left_join(para_degree, by = "node_id") %>%
  left_join(para_cent, by = "node_id") %>%
  mutate(
    display_label = coalesce(display_label, resolve_name(node_id)),
    canonical_gene = display_label,
    layer = "pathogen",
    node_type = "measured",
    degree = coalesce(degree, 0L),
    strength = degree,
    betweenness = coalesce(betweenness, 0),
    eigenvector = coalesce(eigenvector, 0),
    module_status = if_else(!is.na(module_id) & nzchar(module_id), "network_module", "unassigned"),
    module_function_signature = coalesce(parasite_program, module_theme_signature),
    local_core = coalesce(local_core, FALSE) & node_id != pin1_id,
    accessible_in_network = coalesce(accessible_in_network, FALSE),
    recurrent_core = coalesce(recurrent_core, FALSE),
    scaffold_tier = case_when(
      local_core ~ "strict",
      accessible_in_network ~ "relaxed",
      TRUE ~ "context"
    ),
    network_role = case_when(
      local_core ~ "Core",
      accessible_in_network ~ "Accessible scaffold",
      TRUE ~ "Network context"
    ),
    strict_member = local_core,
    relaxed_member = accessible_in_network
  )

para_edge_tbl <- para_edges_raw %>%
  mutate(
    edge_scope = "intra_pathogen",
    edge_type = "string",
    weight = 1,
    confidence = NA_real_,
    evidence_source = "STRING",
    species_basis = "Theileria_annulata",
    network_basis = "parasite_string_network_v1",
    source_rank = NA_real_,
    target_rank = NA_real_,
    mutual_topk = NA
  )
para_full_layout <- layout_network_by_module(para_node_tbl, para_edge_tbl, seed = 20260904L)
write_fig(ensure_cols(para_full_layout$nodes, gold_fig2_g_nodes), "02_figure_2", "05_panel_g_pathogen_network_nodes.tsv")
write_fig(ensure_cols(para_full_layout$edges, gold_fig2_g_edges), "02_figure_2", "06_panel_g_pathogen_network_edges.tsv")

para_functions <- para_node_tbl %>%
  filter(network_role %in% c("Core", "Accessible scaffold")) %>%
  mutate(module_id = coalesce(module_id, "Unassigned")) %>%
  group_by(module_id) %>%
  summarise(
    scaffold_n = n(),
    core_n = sum(local_core),
    module_label = first_known(module_label, first(module_id)),
    module_biological_label = first_known(module_biological_label, first(module_id)),
    dominant_function_class = first_known(parasite_program, "unassigned"),
    module_function_signature = first_known(module_function_signature, first(module_id)),
    module_size = max(num(module_size), na.rm = TRUE),
    measured_n = n(),
    context_n = 0L,
    string_mapped_n = n(),
    primary_evidence_n = n(),
    support_evidence_n = sum(recurrent_core),
    median_abs_surface_effect = NA_real_,
    top_nodes = paste(head(display_label[order(-degree, display_label)], 8), collapse = ";"),
    total_strength = sum(strength),
    .groups = "drop"
  ) %>%
  mutate(module_size = if_else(is.finite(module_size), module_size, as.numeric(scaffold_n))) %>%
  arrange(desc(core_n), desc(scaffold_n))
write_fig(ensure_cols(para_functions, gold_fig2_e), "02_figure_2", "04_panel_e_pathogen_functions.tsv")

write_fig(
  para_node_tbl %>%
    filter(scaffold_tier %in% c("strict", "relaxed")) %>%
    transmute(
      node_id,
      scaffold_tier,
      localization_class = coalesce(parasite_program, if_else(accessible_annot, "accessible_annot", "other")),
      interface_function = local_core,
      local_core,
      recurrent_core,
      accessible_annot,
      parasite_program
    ),
  "figure_2", "03_panel_e_evidence_geometry.tsv"
)

# ---------------------------------------------------------------------------
# Figure 3
# ---------------------------------------------------------------------------
supported_prox <- proximity %>% filter(proximity_support)
axis_mod_pairs <- retained_axis %>%
  transmute(
    host_network_module = proximity_host_module,
    pathogen_network_module = proximity_pathogen_module,
    biological_axis_id,
    species_adjusted_confidence
  ) %>%
  filter(!is.na(host_network_module), !is.na(pathogen_network_module),
         nzchar(host_network_module), nzchar(pathogen_network_module))

bridge_from_prox <- supported_prox %>%
  transmute(
    host_network_module = host_module,
    pathogen_network_module = pathogen_module,
    proximity_z, proximity_fdr, module_family
  ) %>%
  left_join(
    axis_mod_pairs %>%
      group_by(host_network_module, pathogen_network_module) %>%
      summarise(axis_n = n(), conf_sum = sum(species_adjusted_confidence, na.rm = TRUE), .groups = "drop"),
    by = c("host_network_module", "pathogen_network_module")
  ) %>%
  mutate(
    axis_n = coalesce(axis_n, 0L),
    bridge_weight = pmax(-num(proximity_z), 0) + coalesce(conf_sum, 0)
  )

axis_only <- axis_mod_pairs %>%
  anti_join(bridge_from_prox, by = c("host_network_module", "pathogen_network_module")) %>%
  group_by(host_network_module, pathogen_network_module) %>%
  summarise(axis_n = n(), bridge_weight = sum(species_adjusted_confidence, na.rm = TRUE), .groups = "drop")

module_bridges <- bind_rows(
  bridge_from_prox %>% select(host_network_module, pathogen_network_module, bridge_weight, axis_n),
  axis_only
) %>%
  group_by(host_network_module, pathogen_network_module) %>%
  summarise(bridge_weight = max(bridge_weight, na.rm = TRUE), axis_n = max(axis_n, na.rm = TRUE), .groups = "drop")
write_fig(ensure_cols(module_bridges, gold_fig3_bridges), "03_figure_3", "02_panel_a_module_bridges.tsv")

host_scaf_ids <- host$node_id[host$independent_scaffold]
host_core_ids <- host$node_id[host$display_core]
para_scaf_ids <- para$node_id[para$in_para_network]
para_core_ids2 <- para$node_id[para$local_core]

count_in_module <- function(mod_id, ids, layer_name) {
  if (str_detect(mod_id, "_FUN_")) {
    sum(fun_modules$module_id == mod_id & fun_modules$node_id %in% ids)
  } else {
    sum(net_modules$module_id == mod_id & net_modules$node_id %in% ids)
  }
}

fig3_module_ids <- unique(c(module_bridges$host_network_module, module_bridges$pathogen_network_module))
module_nodes_fig3 <- tibble(module_id = fig3_module_ids) %>%
  left_join(module_lookup, by = "module_id") %>%
  mutate(
    side = case_when(
      str_starts(module_id, "H_") | layer == "host" ~ "host",
      str_starts(module_id, "P_") | layer == "pathogen" ~ "pathogen",
      TRUE ~ if_else(str_starts(module_id, "H"), "host", "pathogen")
    ),
    module_label = coalesce(module_label, module_id),
    module_biological_label = coalesce(module_biological_label, module_id),
    module_theme_signature = coalesce(module_theme_signature, module_id),
    scaffold_n = vapply(seq_len(n()), function(i) {
      if (side[i] == "host") count_in_module(module_id[i], host_scaf_ids, "host")
      else count_in_module(module_id[i], para_scaf_ids, "pathogen")
    }, integer(1)),
    core_n = vapply(seq_len(n()), function(i) {
      if (side[i] == "host") count_in_module(module_id[i], host_core_ids, "host")
      else count_in_module(module_id[i], para_core_ids2, "pathogen")
    }, integer(1))
  ) %>%
  mutate(scaffold_n = pmax(scaffold_n, 1L))
write_fig(ensure_cols(module_nodes_fig3, gold_fig3_nodes), "03_figure_3", "01_panel_a_module_nodes.tsv")

rwr_landscape <- rwr %>%
  left_join(name_lookup, by = "node_id") %>%
  left_join(net_modules %>% select(node_id, module_id), by = "node_id") %>%
  left_join(module_lookup %>% select(module_id, module_label), by = "module_id") %>%
  left_join(host %>% select(node_id, host_strict = strict_member, host_relaxed = relaxed_member), by = "node_id") %>%
  left_join(para %>% select(node_id, para_strict = local_core, para_relaxed = accessible_in_network), by = "node_id") %>%
  transmute(
    node_id, layer,
    display_label = coalesce(gene_symbol, resolve_name(node_id)),
    module_label,
    degree,
    strict_member = case_when(
      layer == "host" ~ coalesce(host_strict, FALSE),
      TRUE ~ coalesce(para_strict, FALSE)
    ),
    relaxed_member = case_when(
      layer == "host" ~ coalesce(host_relaxed, FALSE),
      TRUE ~ coalesce(para_relaxed, FALSE)
    ),
    direct_evidence = side_seed_percentile,
    opposite_signal = cross_layer_gain_percentile,
    recruited = is.finite(side_seed_percentile) & is.finite(cross_layer_gain_percentile) &
      side_seed_percentile < 0.5 & cross_layer_gain_percentile >= 0.95,
    opposite_side_rwr_percentile, same_side_diffusion_percentile
  ) %>%
  filter(is.finite(direct_evidence), is.finite(opposite_signal))
write_fig(ensure_cols(rwr_landscape, gold_fig3_rwr), "03_figure_3", "03_panel_b_rwr_landscape.tsv")

proximity_out <- proximity %>%
  left_join(
    program_groups %>%
      select(host_module, pathogen_module, program_group_id, program_group, display_label,
             host_display_name, pathogen_display_name),
    by = c("host_module", "pathogen_module")
  ) %>%
  mutate(
    host_reachable_fraction = host_reachable_n / pmax(host_size, 1),
    pathogen_reachable_fraction = pathogen_reachable_n / pmax(pathogen_size, 1),
    null_q025 = NA_real_,
    null_q975 = NA_real_,
    control_type = NA_character_
  )
write_fig(ensure_cols(proximity_out, gold_fig3_prox), "03_figure_3", "04_panel_e_module_proximity.tsv")

# ---------------------------------------------------------------------------
# Figure 4
# ---------------------------------------------------------------------------
write_fig(pcsf_runs, "04_figure_4", "01_panel_a_pcsf_parameters.tsv")

pcsf_nodes_ann <- pcsf_nodes %>%
  left_join(name_lookup, by = "node_id") %>%
  left_join(net_modules %>% select(node_id, module_id), by = "node_id") %>%
  left_join(module_lookup %>% select(module_id, module_label), by = "module_id") %>%
  left_join(host %>% select(node_id, host_strict = strict_member, host_relaxed = relaxed_member), by = "node_id") %>%
  left_join(para %>% select(node_id, para_strict = local_core, para_relaxed = accessible_in_network), by = "node_id") %>%
  mutate(
    display_label = coalesce(gene_symbol, resolve_name(node_id), clean_plot_text(pathogen_name), node_id),
    module_label = coalesce(module_label, module_id),
    strict_member = case_when(
      layer == "host" ~ coalesce(host_strict, strict_member),
      TRUE ~ coalesce(para_strict, strict_member)
    ),
    relaxed_member = case_when(
      layer == "host" ~ coalesce(host_relaxed, relaxed_member),
      TRUE ~ coalesce(para_relaxed, relaxed_member)
    )
  )
write_fig(ensure_cols(pcsf_nodes_ann, gold_fig4_union_nodes), "04_figure_4", "02_panel_a_pcsf_union_nodes.tsv")
write_fig(pcsf_edges, "04_figure_4", "03_panel_a_pcsf_union_edges.tsv")

stopifnot(sum(pcsf_nodes$selection_frequency >= 0.5) == 104L)
stopifnot(sum(pcsf_edges$selection_frequency >= 0.5) == 95L)
stable_pcsf_nodes <- pcsf_nodes_ann %>%
  filter(selection_frequency >= 0.5) %>%
  arrange(layer, desc(selection_frequency))
stable_pcsf_edges <- pcsf_edges %>%
  filter(selection_frequency >= 0.5,
         source %in% stable_pcsf_nodes$node_id,
         target %in% stable_pcsf_nodes$node_id)

endpoint_nodes <- unique(c(retained_axis$host_node, retained_axis$pathogen_node))
if (nrow(stable_pcsf_nodes) >= 2 && nrow(stable_pcsf_edges)) {
  g <- graph_from_data_frame(
    stable_pcsf_edges %>% select(source, target, selection_frequency),
    directed = FALSE,
    vertices = stable_pcsf_nodes %>% transmute(name = node_id)
  )
  set.seed(20260905)
  xy <- layout_with_fr(g, weights = pmax(E(g)$selection_frequency, 0.01), niter = 1200, grid = "nogrid")
  pcsf_xy <- tibble(node_id = V(g)$name, x0 = xy[, 1], y = xy[, 2]) %>%
    left_join(stable_pcsf_nodes, by = "node_id") %>%
    mutate(x = if_else(layer == "host", -0.5 - abs(x0), 0.5 + abs(x0)))
  pcsf_edge_xy <- stable_pcsf_edges %>%
    left_join(pcsf_xy %>% select(source = node_id, x_from = x, y_from = y), by = "source") %>%
    left_join(pcsf_xy %>% select(target = node_id, x_to = x, y_to = y), by = "target")
} else {
  pcsf_xy <- tibble()
  pcsf_edge_xy <- tibble()
}
write_fig(pcsf_xy, "04_figure_4", "04_panel_b_pcsf_nodes.tsv")
write_fig(pcsf_edge_xy, "04_figure_4", "05_panel_b_pcsf_edges.tsv")

stable_incidence <- bind_rows(
  stable_pcsf_edges %>% transmute(node_id = source, edge_weight = selection_frequency),
  stable_pcsf_edges %>% transmute(node_id = target, edge_weight = selection_frequency)
) %>%
  group_by(node_id) %>%
  summarise(forest_degree = n(), weighted_degree = sum(edge_weight), .groups = "drop")
pcsf_node_annotated <- stable_pcsf_nodes %>%
  left_join(stable_incidence, by = "node_id") %>%
  mutate(
    forest_degree = coalesce(forest_degree, 0L),
    weighted_degree = coalesce(weighted_degree, 0),
    network_degree = coalesce(degree, 0),
    network_strength = coalesce(degree, 0),
    endpoint = node_id %in% endpoint_nodes
  ) %>%
  group_by(layer) %>%
  mutate(
    hub_score = 0.50 * percent_rank(forest_degree) +
      0.35 * percent_rank(weighted_degree) +
      0.15 * percent_rank(selection_frequency),
    hub = min_rank(desc(hub_score)) <= 4
  ) %>%
  ungroup()

pcsf_module_burden <- pcsf_node_annotated %>%
  mutate(module_id = coalesce(module_id, "Unassigned"),
         module_label = coalesce(module_label, module_id)) %>%
  group_by(layer, module_id, module_label) %>%
  summarise(
    selected_nodes = n(),
    mean_selection_frequency = mean(selection_frequency, na.rm = TRUE),
    endpoint_n = sum(endpoint, na.rm = TRUE),
    hub_n = sum(hub, na.rm = TRUE),
    total_strength = sum(weighted_degree, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(layer, desc(selected_nodes), desc(mean_selection_frequency), module_label)
write_fig(ensure_cols(pcsf_module_burden, gold_fig4_burden), "04_figure_4", "06_panel_c_pcsf_module_burden.tsv")

axis_descriptor <- axis_cons %>%
  mutate(
    axis_id = biological_axis_id,
    interaction_name = biological_axis_id,
    pathogen_gene = parasite_label,
    host_gene = host_label,
    pathogen_nodes = pathogen_node,
    host_nodes = host_node,
    interaction_type = bridge_evidence_tier,
    evidence_strength = physical_claim,
    main_algorithm_bridge_gate = main_algorithm_bridge,
    pathogen_support_gate = pathogen_recurrence_gate,
    algorithm_combination = NA_character_,
    host_network_module = proximity_host_module,
    pathogen_network_module = proximity_pathogen_module,
    host_program = coalesce(
      module_lookup$module_label[match(proximity_host_module, module_lookup$module_id)],
      proximity_host_module,
      "Host program unresolved"
    ),
    pathogen_program = coalesce(
      module_lookup$module_label[match(proximity_pathogen_module, module_lookup$module_id)],
      proximity_pathogen_module,
      family_display,
      "Pathogen program unresolved"
    ),
    interface_score = rwr_endpoint_min_percentile,
    confidence = species_adjusted_confidence,
    persistence = threshold_consensus_frequency,
    aggregate_robustness_score = rwr_endpoint_min_percentile,
    stability_score = threshold_consensus_frequency,
    axis_label = paste0(parasite_label, " -> ", host_label),
    pcsf_component_path_frequency = num(component_path_frequency),
    pcsf_successful_runs = 18,
    pcsf_total_runs = 18,
    pcsf_success_fraction = 1
  ) %>%
  rowwise() %>%
  mutate(
    algorithm_combination = {
      z <- c("RWR", "Proximity", "PCSF")[c(isTRUE(rwr_support), isTRUE(proximity_support), isTRUE(pcsf_support))]
      if (length(z)) paste(z, collapse = "+") else "No algorithm gate"
    }
  ) %>%
  ungroup() %>%
  arrange(desc(num(interface_score)), biological_axis_id) %>%
  mutate(rank = row_number())
stopifnot(!any(str_detect(coalesce(axis_descriptor$axis_label, ""), "TA_15795\\?|（")))
write_fig(ensure_cols(axis_descriptor, gold_fig4_axis), "04_figure_4", "07_panel_d_axis_descriptor.tsv")
write_fig(ensure_cols(
  pcsf_node_annotated %>% arrange(layer, desc(selection_frequency), desc(forest_degree), display_label),
  gold_fig4_persist
), "04_figure_4", "08_panel_f_pcsf_node_persistence.tsv")

# ---------------------------------------------------------------------------
# Figure 5
# ---------------------------------------------------------------------------
rwr_join <- rwr %>%
  select(node_id, joint_rwr_percentile, cross_layer_gain_percentile, side_seed_raw, side_seed_percentile,
         same_side_diffusion_percentile)
layer_of <- function(id) if_else(str_starts(id, "9913."), "host", "pathogen")

final_node_ids <- unique(c(stable_pcsf_nodes$node_id, retained_host_ends, retained_path_ends))
final_node_ids <- final_node_ids[!is.na(final_node_ids) & nzchar(final_node_ids)]

core_end_ids <- unique(c(strict_host_ends, strict_path_ends))

intra_edges <- stable_pcsf_edges %>%
  mutate(source_layer = layer_of(source), target_layer = layer_of(target)) %>%
  filter(source_layer == target_layer) %>%
  transmute(
    source, target,
    pair_a = pmin(source, target), pair_b = pmax(source, target),
    edge_scope = paste0("intra_", source_layer),
    edge_type = "pcsf_persistent",
    weight = selection_frequency,
    confidence = selection_frequency,
    evidence_source = "PCSF",
    species_basis = if_else(source_layer == "host", "Bos_taurus", "Theileria_annulata"),
    network_basis = "pcsf_persistent",
    interaction_name = NA_character_,
    interaction_type = "pcsf_intra",
    evidence_reference = NA_character_,
    host_species_evidence = NA_character_,
    species_adjusted_confidence = NA_real_,
    axis_id = NA_character_,
    pcsf_edge_selected_runs = selected_runs,
    pcsf_edge_frequency = selection_frequency,
    edge_semantics = "pcsf_persistent_intra",
    known_cross_layer_anchor = FALSE,
    algorithm_derived_cross_layer_edge = FALSE,
    mechanism_class = NA_character_
  )

cross_edges <- retained_axis %>%
  transmute(
    source = pathogen_node,
    target = host_node,
    pair_a = pmin(pathogen_node, host_node),
    pair_b = pmax(pathogen_node, host_node),
    edge_scope = "cross",
    edge_type = "biological_axis",
    weight = if_else(consensus_tier == "strict", 3, 1),
    confidence = species_adjusted_confidence,
    evidence_source = bridge_evidence_tier,
    species_basis = "cross_layer",
    network_basis = "biological_axis_consensus",
    interaction_name = biological_axis_id,
    interaction_type = bridge_evidence_tier,
    evidence_reference = claim_boundary,
    host_species_evidence = "Bos_taurus",
    species_adjusted_confidence,
    axis_id = biological_axis_id,
    pcsf_edge_selected_runs = num(selected_runs),
    pcsf_edge_frequency,
    edge_semantics = if_else(consensus_tier == "strict", "core_axis", "extended_axis"),
    known_cross_layer_anchor = TRUE,
    algorithm_derived_cross_layer_edge = FALSE,
    mechanism_class = family_display
  )
final_edges <- bind_rows(intra_edges, cross_edges)

final_nodes <- tibble(node_id = final_node_ids) %>%
  mutate(layer = layer_of(node_id)) %>%
  left_join(pcsf_nodes_ann %>% select(node_id, selected_runs, selection_frequency, pcsf_degree = degree,
                                      pcsf_strict = strict_member, pcsf_relaxed = relaxed_member,
                                      pcsf_seed_percentile = side_seed_percentile,
                                      display_label_pcsf = display_label,
                                      module_id_pcsf = module_id, module_label_pcsf = module_label),
            by = "node_id") %>%
  left_join(net_modules %>% select(node_id, module_id, net_module_size = module_size), by = "node_id") %>%
  left_join(module_lookup %>% select(module_id, module_label, module_biological_label, module_theme_signature),
            by = "module_id") %>%
  left_join(rwr_join, by = "node_id") %>%
  left_join(host %>% select(node_id, host_strict = strict_member, host_relaxed = relaxed_member,
                            host_evidence = evidence_count, independent_evidence_count, independent_sources),
            by = "node_id") %>%
  left_join(para %>% select(node_id, para_strict = local_core, para_relaxed = accessible_in_network,
                            para_evidence = evidence_count, sources, annotation),
            by = "node_id") %>%
  left_join(annot_prot, by = "node_id") %>%
  left_join(name_lookup, by = "node_id") %>%
  left_join(bind_rows(host_degree, para_degree) %>% distinct(node_id, .keep_all = TRUE), by = "node_id") %>%
  mutate(
    display_label = coalesce(gene_symbol, display_label_pcsf, resolve_name(node_id)),
    canonical_gene = display_label,
    protein_name = coalesce(protein_name, annotation),
    module_id = coalesce(module_id, module_id_pcsf, "Unassigned"),
    module_label = coalesce(module_label, module_label_pcsf, module_id),
    module_signature = module_theme_signature,
    node_type = "measured",
    biological_theme = if_else(layer == "host", "independent_host_state", "parasite_accessible_state"),
    mapped_to_reference_network = TRUE,
    degree = coalesce(pcsf_degree, degree, 0),
    strength = degree,
    strict_member = coalesce(host_strict, para_strict, pcsf_strict, FALSE),
    relaxed_member = coalesce(host_relaxed, para_relaxed, pcsf_relaxed, FALSE),
    pcsf_selected_runs = coalesce(selected_runs, 0),
    pcsf_node_frequency = coalesce(selection_frequency, 0),
    core_axis_endpoint = node_id %in% core_end_ids,
    pcsf_parameter_stable = coalesce(selection_frequency, 0) >= 0.5,
    rwr_parameter_stable = FALSE,
    connector_node = !core_axis_endpoint & coalesce(selection_frequency, 0) >= 0.5,
    scaffold_role = case_when(
      core_axis_endpoint ~ "core_axis_endpoint",
      pcsf_parameter_stable ~ "pcsf_persistent",
      TRUE ~ "axis_endpoint"
    ),
    n_study_families = if_else(layer == "host", independent_evidence_count, para_evidence),
    study_families = if_else(layer == "host", independent_sources, sources),
    seed_weight = coalesce(side_seed_percentile, pcsf_seed_percentile)
  )

if (nrow(final_nodes) >= 2 && nrow(final_edges)) {
  g5 <- graph_from_data_frame(
    final_edges %>% select(source, target, weight),
    directed = FALSE,
    vertices = final_nodes %>% transmute(name = node_id)
  )
  set.seed(20260906)
  initial <- cbind(
    ifelse(final_nodes$layer[match(V(g5)$name, final_nodes$node_id)] == "host", -1, 1),
    runif(vcount(g5), -1, 1)
  )
  xy <- layout_with_fr(g5, coords = initial, weights = pmax(E(g5)$weight, 0.01), niter = 1800, grid = "nogrid")
  final_xy <- tibble(node_id = V(g5)$name, x0 = xy[, 1], y0 = xy[, 2]) %>%
    left_join(final_nodes, by = "node_id") %>%
    mutate(
      x = if_else(layer == "host", -0.6 - abs(x0), 0.6 + abs(x0)),
      y = y0,
      label = display_label,
      hub = percent_rank(desc(num(degree))) <= 0.05
    )
  final_edge_xy <- final_edges %>%
    left_join(final_xy %>% select(source = node_id, x_from = x, y_from = y), by = "source") %>%
    left_join(final_xy %>% select(target = node_id, x_to = x, y_to = y), by = "target")
} else {
  final_xy <- tibble()
  final_edge_xy <- tibble()
}
stopifnot(sum(final_xy$core_axis_endpoint) == length(core_end_ids))
stopifnot(sum(final_edge_xy$known_cross_layer_anchor %in% TRUE) == nrow(retained_axis))
stopifnot(!any(str_detect(coalesce(final_xy$display_label, ""), "TA_15795\\?|（")))

module_enrichment <- function(layer_name, selected, background) {
  bg <- background %>% filter(layer == layer_name, !is.na(module_label), nzchar(module_label))
  sel <- selected %>% filter(layer == layer_name, !is.na(module_label), nzchar(module_label))
  modules <- sort(unique(bg$module_label))
  if (!length(modules) || !nrow(sel)) {
    return(tibble(
      layer = character(), module_label = character(), selected_n = integer(),
      selected_total = integer(), background_n = integer(), background_total = integer(),
      odds_ratio = numeric(), ci_low = numeric(), ci_high = numeric(), p_value = numeric(), fdr = numeric()
    ))
  }
  bind_rows(lapply(modules, function(module_name) {
    selected_n <- sum(sel$module_label == module_name)
    background_n <- sum(bg$module_label == module_name)
    other_selected <- nrow(sel) - selected_n
    other_background <- nrow(bg) - nrow(sel) - (background_n - selected_n)
    other_background <- max(other_background, 0)
    contingency <- matrix(
      c(selected_n, other_selected, max(background_n - selected_n, 0), other_background),
      nrow = 2, byrow = TRUE
    )
    if (any(!is.finite(contingency)) || any(contingency < 0)) {
      fit <- list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = NA_real_)
    } else {
      fit <- tryCatch(fisher.test(contingency), error = function(e) {
        list(estimate = NA_real_, conf.int = c(NA_real_, NA_real_), p.value = NA_real_)
      })
    }
    tibble(
      layer = layer_name,
      module_label = module_name,
      selected_n = as.integer(selected_n),
      selected_total = nrow(sel),
      background_n = as.integer(background_n),
      background_total = nrow(bg),
      odds_ratio = unname(fit$estimate),
      ci_low = unname(fit$conf.int[[1]]),
      ci_high = unname(fit$conf.int[[2]]),
      p_value = fit$p.value
    )
  })) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    arrange(fdr, desc(odds_ratio), module_label)
}

background_nodes <- bind_rows(
  host_node_tbl %>% transmute(node_id, layer, module_label = coalesce(module_label, "Unassigned")),
  para_node_tbl %>% transmute(node_id, layer, module_label = coalesce(module_label, "Unassigned"))
)
host_enrichment <- module_enrichment("host", final_xy, background_nodes)
pathogen_enrichment <- module_enrichment("pathogen", final_xy, background_nodes)

module_nodes5 <- final_xy %>%
  mutate(
    module_label = coalesce(module_label, "Unassigned"),
    module_key = paste(layer, module_label, sep = "::")
  ) %>%
  group_by(module_key, layer, module_label) %>%
  summarise(
    scaffold_nodes = n(),
    core_endpoint_n = sum(flag(core_axis_endpoint)),
    hub_n = sum(flag(hub)),
    degree_sum = sum(num(degree)),
    strength_sum = sum(num(strength)),
    .groups = "drop"
  ) %>%
  left_join(
    bind_rows(host_enrichment, pathogen_enrichment) %>%
      select(layer, module_label, odds_ratio, ci_low, ci_high, fdr),
    by = c("layer", "module_label")
  ) %>%
  arrange(layer, desc(scaffold_nodes), module_label)

module_lookup5 <- final_xy %>%
  transmute(node_id, module_key = paste(layer, coalesce(module_label, "Unassigned"), sep = "::"))
module_edges5 <- final_edges %>%
  left_join(module_lookup5 %>% rename(source = node_id, source_module = module_key), by = "source") %>%
  left_join(module_lookup5 %>% rename(target = node_id, target_module = module_key), by = "target") %>%
  filter(!is.na(source_module), !is.na(target_module), source_module != target_module) %>%
  mutate(from_module = pmin(source_module, target_module), to_module = pmax(source_module, target_module)) %>%
  group_by(from_module, to_module) %>%
  summarise(
    edge_n = n(),
    weight_sum = sum(pmax(num(weight), 0.01)),
    cross_axis_n = sum(flag(known_cross_layer_anchor)),
    .groups = "drop"
  ) %>%
  arrange(desc(cross_axis_n), desc(edge_n), from_module, to_module)

core_mechanisms <- axis_descriptor %>%
  filter(final_tier %in% c("Core", "Extended")) %>%
  transmute(
    axis_id, interaction_name, pathogen_gene, host_gene, final_tier,
    aggregate_robustness_score, host_module = host_program, pathogen_module = pathogen_program,
    direct_physical_edge, interaction_type, edge_type = bridge_evidence_tier,
    species_adjusted_confidence, mechanism_class = family_display
  ) %>%
  arrange(mechanism_class, desc(num(aggregate_robustness_score)), interaction_name)
stopifnot(nrow(core_mechanisms) == 11L)
stopifnot(sum(core_mechanisms$final_tier == "Core") == 3L)
stopifnot(all(!is.na(core_mechanisms$mechanism_class)))
stopifnot(!any(str_detect(paste(coalesce(core_mechanisms$pathogen_gene, ""), coalesce(core_mechanisms$host_gene, ""), coalesce(core_mechanisms$interaction_name, "")), "TA_15795\\?|（")))
stopifnot(!any(core_mechanisms$axis_id == "p104--CLASP1"))

program_coupling <- axis_coupling %>%
  transmute(
    host_module = host_program,
    pathogen_module = parasite_program,
    axis_n,
    core_axis_n = strict_axis_n,
    mean_score = mean_rwr_endpoint_min
  ) %>%
  arrange(desc(core_axis_n), desc(axis_n), desc(mean_score), host_module, pathogen_module)

write_fig(ensure_cols(final_xy, gold_fig5_nodes), "05_figure_5", "01_panel_a_final_nodes.tsv")
write_fig(ensure_cols(final_edge_xy, gold_fig5_edges), "05_figure_5", "02_panel_a_final_edges.tsv")
write_fig(ensure_cols(module_nodes5, gold_fig5_mod_nodes), "05_figure_5", "03_panel_b_module_nodes.tsv")
write_fig(ensure_cols(module_edges5, gold_fig5_mod_edges), "05_figure_5", "04_panel_b_module_edges.tsv")
write_fig(ensure_cols(host_enrichment, gold_fig5_enrich), "05_figure_5", "05_panel_c_host_module_enrichment.tsv")
write_fig(ensure_cols(pathogen_enrichment, gold_fig5_enrich), "05_figure_5", "06_panel_c_pathogen_module_enrichment.tsv")
write_fig(ensure_cols(program_coupling, gold_fig5_coup), "05_figure_5", "07_panel_d_program_coupling.tsv")
write_fig(ensure_cols(core_mechanisms, gold_fig5_mech), "05_figure_5", "08_panel_d_core_mechanisms.tsv")

expected <- c(
  "01_figure_1/01_panel_a_host_studies.tsv", "01_figure_1/02_panel_b_host_attrition.tsv",
  "01_figure_1/03_panel_c_host_recurrence.tsv", "01_figure_1/04_panel_d_host_statistics.tsv",
  "01_figure_1/05_panel_e_host_functions.tsv",
  "01_figure_1/06_panel_g_host_network_nodes.tsv", "01_figure_1/07_panel_g_host_network_edges.tsv",
  "02_figure_2/01_panel_a_pathogen_studies.tsv", "02_figure_2/02_panel_b_pathogen_attrition.tsv",
  "02_figure_2/03_panel_e_evidence_geometry.tsv", "02_figure_2/04_panel_e_pathogen_functions.tsv",
  "02_figure_2/05_panel_g_pathogen_network_nodes.tsv", "02_figure_2/06_panel_g_pathogen_network_edges.tsv",
  "03_figure_3/01_panel_a_module_nodes.tsv", "03_figure_3/02_panel_a_module_bridges.tsv",
  "03_figure_3/03_panel_b_rwr_landscape.tsv", "03_figure_3/04_panel_e_module_proximity.tsv",
  "04_figure_4/01_panel_a_pcsf_parameters.tsv", "04_figure_4/02_panel_a_pcsf_union_nodes.tsv",
  "04_figure_4/03_panel_a_pcsf_union_edges.tsv", "04_figure_4/04_panel_b_pcsf_nodes.tsv",
  "04_figure_4/05_panel_b_pcsf_edges.tsv", "04_figure_4/06_panel_c_pcsf_module_burden.tsv",
  "04_figure_4/07_panel_d_axis_descriptor.tsv", "04_figure_4/08_panel_f_pcsf_node_persistence.tsv",
  "05_figure_5/01_panel_a_final_nodes.tsv", "05_figure_5/02_panel_a_final_edges.tsv",
  "05_figure_5/03_panel_b_module_nodes.tsv", "05_figure_5/04_panel_b_module_edges.tsv",
  "05_figure_5/05_panel_c_host_module_enrichment.tsv", "05_figure_5/06_panel_c_pathogen_module_enrichment.tsv",
  "05_figure_5/panel_d_program_coupling.tsv", "05_figure_5/panel_d_core_mechanisms.tsv"
)
info <- file.info(file.path(stage, expected))
if (any(!file.exists(file.path(stage, expected))) || any(info$size == 0)) {
  stop("Five-figure data staging validation failed.")
}

dir.create(figure_data_final, recursive = TRUE, showWarnings = FALSE)
for (fig in paste0("figure_", 1:5)) {
  dest <- file.path(figure_data_final, fig)
  src <- file.path(stage, fig)
  if (dir.exists(dest)) unlink(dest, recursive = TRUE, force = TRUE)
  renamed <- file.rename(src, dest)
  if (!isTRUE(renamed) || !dir.exists(dest)) {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(list.files(src, full.names = TRUE), dest, overwrite = TRUE)
    if (any(!ok)) stop("Could not publish ", fig)
  }
}

contract <- c(
  "# 09 Figure-data contract",
  "",
  "Gold-shaped panel TSVs live under `02_results/02_figure_data/figure_1` … `figure_5`.",
  "Parent-level `fig01_*` copies from `08_export_publication_data.py` are left in place.",
  "Algorithm tables under `02_results/04_intermediate/01_algorithms/` are not modified.",
  "",
  "## Funnel",
  "",
  paste0("- Fig1 stages: mapped host proteins ", n_mapped,
         "; any-source vote ", n_any_vote,
         "; independent vote ", n_indep_vote,
         "; independent recurrent in STRING700 ", n_scaffold,
         "; display core ", n_display_core, "."),
  "- Display core = `independent_evidence_count >= 3` OR host endpoint of a manuscript-strict axis (p104--MAPRE1, TaPin1--PKM, TA20980--CD2AP).",
  "- Independent scaffold (Fig1 highlight / `relaxed_member`) = `independent_evidence_count >= 2` AND `in_700_network`.",
  "- All-source `evidence_count` remains the RWR seed and is not presented as Fig1 core.",
  paste0("- Fig2 stages: parasite evidence proteins ", nrow(para),
         "; in parasite network ", sum(para$in_para_network),
         "; accessible annot ", sum(para$accessible_annot),
         "; accessible in network ", sum(para$accessible_in_network),
         "; local core ", sum(para$local_core), "."),
  "- Fig2 `network_role` Core = 19 local-core nodes. Pin1 (`5874.Q4UG71`) is never local-core.",
  "- Fig2 Accessible scaffold = `accessible_annot` proteins in the parasite network that are not local-core.",
  "",
  "## Networks",
  "",
  paste0("- Host panel_g layout mode: `", host_layout_mode, "`."),
  paste0("- Host STRING 700 file actually contains ", length(host_graph_nodes), " nodes / ",
         nrow(host_full_layout$edges), " edges (network_summary 8835/350k is the broader host_intra universe, not this 700 file)."),
  paste0("- Parasite STRING nodes/edges: ", length(para_graph_nodes), " / ", nrow(para_full_layout$edges), "."),
  "- Layout function is gold `layout_network_by_module` with seeds 20260903 (host), 20260904 (pathogen), 20260905 (PCSF), 20260906 (atlas).",
  "",
  "## Atlas",
  "",
  paste0("- PCSF persistent (>=0.5): ", nrow(stable_pcsf_nodes), " nodes / ", nrow(stable_pcsf_edges), " edges."),
  paste0("- Fig5 nodes: ", nrow(final_xy), " (PCSF persistent union 11 retained axis endpoints)."),
  paste0("- Fig5 cross-layer axes: ", sum(flag(final_edge_xy$known_cross_layer_anchor)),
         " retained; 3 Core thick (`weight=3`); `core_axis_endpoint` on strict ends only."),
  "- `p104--CLASP1` is context-only and is not in Fig5 atlas / core-mechanism tables.",
  "- Plot labels use locked axis names. `TA_15795?` and Chinese parentheticals are stripped.",
  "",
  if (length(contract_notes)) paste(c("## Notes", "", paste("-", contract_notes)), collapse = "\n") else NULL
)
writeLines(contract[!is.na(contract)], p("00_docs", "09_FIGURE_DATA_CONTRACT.md"))

manifest <- list.files(figure_data_final, pattern = "^panel_.*\\.tsv$", recursive = TRUE)
manifest <- sort(manifest[str_detect(manifest, "^figure_[1-5]/")])
cat("Published panel plot data at:", figure_data_final, "\n")
cat(sprintf("%-55s %8s\n", "file", "rows"))
for (rel in manifest) {
  path <- file.path(figure_data_final, rel)
  n <- length(readLines(path, warn = FALSE)) - 1L
  cat(sprintf("%-55s %8d\n", rel, n))
}
cat("Host layout mode:", host_layout_mode, "\n")
cat("Fig1 display core / independent scaffold:", n_display_core, "/", n_scaffold, "\n")
cat("Fig2 local core / accessible in network:", sum(para$local_core), "/", sum(para$accessible_in_network), "\n")
cat("Fig5 nodes / cross axes / core endpoints:", nrow(final_xy), "/",
    sum(flag(final_edge_xy$known_cross_layer_anchor)), "/", sum(final_xy$core_axis_endpoint), "\n")
