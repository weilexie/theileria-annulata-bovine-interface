# S. aureus-aligned host-parasite interface analysis.
# The separation of roles is deliberate:
#   continuous side evidence -> bidirectional RWR
#   independently defined network modules -> cross-species proximity
#   heterogeneous graph -> PCSF persistence
#   curated interaction axes -> evidence gate plus algorithm consensus

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(Matrix)
  library(matrixStats)
  library(PCSF)
  library(readr)
  library(stringr)
  library(tidyr)
})

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
p <- function(...) file.path(root, ...)
read_tab <- function(...) read_tsv(p(...), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
flag <- function(x) str_to_lower(coalesce(as.character(x), "")) %in% c("yes", "true", "t", "1", "y")
rank01 <- function(x) {
  z <- rep(0, length(x)); ok <- is.finite(x) & x > 0
  if (sum(ok) == 1L) z[ok] <- 1
  if (sum(ok) > 1L) z[ok] <- rank(x[ok], ties.method = "average") / sum(ok)
  z
}
collapse_values <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")
safe_z <- function(observed, null) if (sd(null) > 0) (observed - mean(null)) / sd(null) else NA_real_

out_dir <- p("03_algorithm_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_network <- function(path, layer) {
  read_tsv(path, col_names = c("source", "target"), show_col_types = FALSE, progress = FALSE) %>%
    filter(source != target) %>% distinct(source, target) %>%
    mutate(weight = 1, edge_scope = paste0("within_", layer))
}

host_edges <- read_network(p("02_processed_inputs", "01_networks", "01_host_string_network_700.tsv"), "host")
path_edges <- read_network(p("02_processed_inputs", "01_networks", "02_parasite_string_network_400.tsv"), "pathogen")
bridge_raw <- read_tab("02_processed_inputs", "01_networks", "03_bridge_edges_string_mapped.tsv")
host_ev <- read_tab("02_processed_inputs", "02_evidence", "01_host_recurrence.tsv")
path_ev <- read_tab("02_processed_inputs", "02_evidence", "02_parasite_accessible.tsv")
string_annotations <- read_tab("02_processed_inputs", "02_evidence", "03_string_annotation_cache.tsv")

tier_confidence <- c(
  direct_validated = 1.00, direct_literature = 0.85, cocomplex = 0.80,
  y2h_hit_unconfirmed = 0.60, paralog_inferred_cocomplex = 0.45
)
main_tiers <- c("direct_validated", "direct_literature", "cocomplex")

cross_all <- bridge_raw %>%
  filter(flag(both_mapped), !is.na(parasite_string_id), !is.na(host_string_id)) %>%
  transmute(
    axis_id = edge_id,
    interaction_name = paste(coalesce(parasite_gene_id_orig, parasite_string_id),
                             coalesce(host_gene_symbol_orig, host_string_id), sep = "--"),
    pathogen_node = parasite_string_id, host_node = host_string_id,
    bridge_evidence_tier = tier,
    species_adjusted_confidence = unname(tier_confidence[tier]),
    main_algorithm_bridge = tier %in% main_tiers,
    direct_physical_edge = tier == "direct_validated",
    claim_boundary
  ) %>%
  filter(!is.na(species_adjusted_confidence)) %>%
  mutate(
    host_endpoint_in_network = host_node %in% c(host_edges$source, host_edges$target),
    pathogen_endpoint_in_network = pathogen_node %in% c(path_edges$source, path_edges$target),
    double_measured = host_endpoint_in_network & pathogen_endpoint_in_network,
    main_algorithm_bridge = main_algorithm_bridge & double_measured
  )
cross <- cross_all %>% filter(main_algorithm_bridge)
if (!nrow(cross)) stop("No main bridge axes remain after evidence gating.", call. = FALSE)

cross_edges_all <- cross_all %>%
  transmute(source = pathogen_node, target = host_node,
            weight = species_adjusted_confidence, edge_scope = "cross",
            axis_id, main_algorithm_bridge)
cross_edges <- cross_edges_all %>% filter(main_algorithm_bridge)
edges <- bind_rows(host_edges, path_edges, cross_edges)
edges_all <- bind_rows(host_edges, path_edges, cross_edges_all)

node_ids <- sort(unique(c(edges$source, edges$target)))
node_layer <- if_else(str_starts(node_ids, "9913."), "host", "pathogen")
nodes <- tibble(node_id = node_ids, layer = node_layer) %>%
  left_join(host_ev %>% transmute(node_id = string_id, host_evidence_n = as.numeric(evidence_count)), by = "node_id") %>%
  left_join(path_ev %>% transmute(node_id = string_id, pathogen_name = preferred_name,
                                  pathogen_evidence_n = as.numeric(evidence_count),
                                  accessible_annot = flag(accessible_annot)), by = "node_id") %>%
  mutate(
    host_evidence_n = coalesce(host_evidence_n, 0),
    pathogen_evidence_n = coalesce(pathogen_evidence_n, 0),
    accessible_annot = coalesce(accessible_annot, FALSE),
    evidence_n = if_else(layer == "host", host_evidence_n, pathogen_evidence_n),
    strict_member = if_else(layer == "host", host_evidence_n >= 2, pathogen_evidence_n >= 2),
    relaxed_member = evidence_n >= 1,
    seed_raw = if_else(
      layer == "host",
      (host_evidence_n / max(host_evidence_n))^2,
      (pathogen_evidence_n / max(pathogen_evidence_n))^2
    )
  ) %>%
  group_by(layer) %>% mutate(side_seed_percentile = rank01(seed_raw)) %>% ungroup()

g_degree <- graph_from_data_frame(edges %>% select(source, target), directed = FALSE,
                                  vertices = nodes %>% transmute(name = node_id))
nodes$degree <- degree(g_degree, v = nodes$node_id)
node_index <- setNames(seq_along(node_ids), node_ids)

build_transition <- function(delta, edge_frame = edges) {
  directed <- bind_rows(
    edge_frame %>% transmute(from = source, to = target, weight = pmax(as.numeric(weight), 1e-6), is_cross = edge_scope == "cross"),
    edge_frame %>% transmute(from = target, to = source, weight = pmax(as.numeric(weight), 1e-6), is_cross = edge_scope == "cross")
  ) %>%
    filter(from %in% node_ids, to %in% node_ids) %>%
    group_by(from) %>%
    mutate(
      has_cross = any(is_cross), has_within = any(!is_cross),
      type_sum = ave(weight, is_cross, FUN = sum),
      allocation = case_when(
        has_cross & has_within & is_cross ~ delta,
        has_cross & has_within & !is_cross ~ 1 - delta,
        TRUE ~ 1
      ),
      probability = allocation * weight / type_sum
    ) %>% ungroup()
  M <- sparseMatrix(i = unname(node_index[directed$to]), j = unname(node_index[directed$from]),
                    x = directed$probability, dims = c(length(node_ids), length(node_ids)))
  isolated <- which(Matrix::colSums(M) == 0)
  if (length(isolated)) M <- M + sparseMatrix(i = isolated, j = isolated, x = 1, dims = dim(M))
  M
}

run_rwr <- function(M, p0, restart, tol = 1e-11, max_iter = 1000L) {
  if (sum(p0) <= 0) stop("RWR seed vector is empty.", call. = FALSE)
  p0 <- p0 / sum(p0); state <- p0
  for (iter in seq_len(max_iter)) {
    next_state <- as.numeric((1 - restart) * (M %*% state) + restart * p0)
    if (sum(abs(next_state - state)) < tol) return(list(score = next_state, iterations = iter))
    state <- next_state
  }
  list(score = state, iterations = max_iter)
}

seed_host <- ifelse(nodes$layer == "host", nodes$side_seed_percentile^2, 0)
seed_path <- ifelse(nodes$layer == "pathogen", nodes$side_seed_percentile^2, 0)
seed_joint <- 0.5 * seed_host / sum(seed_host) + 0.5 * seed_path / sum(seed_path)
restart_grid <- c(0.2, 0.3, 0.5, 0.7)
delta_grid <- c(0.25, 0.5, 0.75)
transition_cache <- setNames(lapply(delta_grid, build_transition), as.character(delta_grid))
rwr_grid <- list(); k <- 1L
for (delta in delta_grid) {
  M <- transition_cache[[as.character(delta)]]
  for (restart in restart_grid) {
    joint_fit <- run_rwr(M, seed_joint, restart)
    host_fit <- run_rwr(M, seed_host, restart)
    path_fit <- run_rwr(M, seed_path, restart)
    rwr_grid[[k]] <- tibble(
      node_id = node_ids, restart, delta, joint_rwr = joint_fit$score,
      from_host_rwr = host_fit$score, from_pathogen_rwr = path_fit$score,
      iterations = joint_fit$iterations
    )
    k <- k + 1L
  }
}
rwr_all <- bind_rows(rwr_grid) %>% left_join(nodes, by = "node_id") %>%
  group_by(restart, delta, layer) %>%
  mutate(
    opposite_side_rwr = if_else(layer == "host", from_pathogen_rwr, from_host_rwr),
    opposite_side_rwr_percentile = rank01(opposite_side_rwr),
    opposite_side_rwr_rank = min_rank(desc(opposite_side_rwr))
  ) %>% ungroup()

primary_restart <- 0.3; primary_delta <- 0.5
M_primary <- transition_cache[[as.character(primary_delta)]]
fit_joint <- run_rwr(M_primary, seed_joint, primary_restart)
fit_host <- run_rwr(M_primary, seed_host, primary_restart)
fit_path <- run_rwr(M_primary, seed_path, primary_restart)
M_uncoupled <- build_transition(0)
fit_host_uncoupled <- run_rwr(M_uncoupled, seed_host, primary_restart)
fit_path_uncoupled <- run_rwr(M_uncoupled, seed_path, primary_restart)

rwr_primary <- nodes %>% transmute(
  node_id, layer, degree, strict_member, relaxed_member, evidence_n, accessible_annot,
  side_seed_raw = seed_raw, side_seed_percentile,
  joint_rwr = fit_joint$score, from_host_rwr = fit_host$score, from_pathogen_rwr = fit_path$score,
  uncoupled_same_side_rwr = if_else(layer == "host", fit_host_uncoupled$score, fit_path_uncoupled$score),
  opposite_side_rwr = if_else(layer == "host", fit_path$score, fit_host$score)
) %>% group_by(layer) %>% mutate(
  joint_rwr_percentile = rank01(joint_rwr),
  opposite_side_rwr_percentile = rank01(opposite_side_rwr),
  same_side_diffusion_percentile = rank01(uncoupled_same_side_rwr),
  cross_layer_gain_percentile = opposite_side_rwr_percentile
) %>% ungroup()

recruitment <- rwr_primary %>%
  filter(side_seed_percentile < 0.50, opposite_side_rwr_percentile >= 0.95) %>%
  arrange(layer, desc(opposite_side_rwr_percentile))

primary_ranks <- rwr_primary %>% group_by(layer) %>%
  mutate(primary_rank = min_rank(desc(opposite_side_rwr))) %>% ungroup() %>%
  select(node_id, layer, primary_percentile = opposite_side_rwr_percentile, primary_rank)
rwr_rank_stability <- rwr_all %>% group_by(restart, delta) %>% group_modify(~{
  x <- .x %>% left_join(primary_ranks, by = c("node_id", "layer"))
  z <- lapply(c("host", "pathogen"), function(s) {
    y <- x %>% filter(layer == s)
    a <- y$node_id[y$opposite_side_rwr_rank <= 200]; b <- y$node_id[y$primary_rank <= 200]
    tibble(layer = s,
           spearman = cor(y$opposite_side_rwr_percentile, y$primary_percentile, method = "spearman"),
           top200_jaccard = length(intersect(a, b)) / length(union(a, b)))
  }) %>% bind_rows()
  tibble(spearman_to_primary = min(z$spearman), top200_jaccard = min(z$top200_jaccard),
         host_spearman = z$spearman[z$layer == "host"], pathogen_spearman = z$spearman[z$layer == "pathogen"])
}) %>% ungroup()

M_all <- build_transition(primary_delta, edges_all)
all_host <- run_rwr(M_all, seed_host, primary_restart)$score
all_path <- run_rwr(M_all, seed_path, primary_restart)$score
rwr_bridge_tier_sensitivity <- nodes %>% transmute(
  node_id, layer,
  primary_score = rwr_primary$opposite_side_rwr,
  all_curated_score = if_else(layer == "host", all_path, all_host)
) %>% group_by(layer) %>% mutate(
  primary_percentile = rank01(primary_score), all_curated_percentile = rank01(all_curated_score),
  percentile_shift = all_curated_percentile - primary_percentile
) %>% ungroup()

write_tsv(rwr_all, file.path(out_dir, "38_rwr_node_scores_all_parameters.tsv"), na = "")
write_tsv(rwr_primary, file.path(out_dir, "36_rwr_cross_layer_gain.tsv"), na = "")
write_tsv(recruitment, file.path(out_dir, "37_rwr_node_recruitment.tsv"), na = "")
write_tsv(rwr_rank_stability, file.path(out_dir, "40_rwr_rank_stability.tsv"), na = "")
write_tsv(rwr_bridge_tier_sensitivity, file.path(out_dir, "39_rwr_primary_vs_all_curated_bridges.tsv"), na = "")

# Degree- and measurement-matched convergence. Measurement availability is used
# as a stratum; evidence count itself is not, which avoids a deterministic null.
conv_nodes <- rwr_primary %>% group_by(layer) %>%
  mutate(degree_bin = ntile(rank(degree, ties.method = "average"), 8L)) %>% ungroup() %>%
  mutate(coverage_class = if_else(evidence_n > 0, "measured", "not_measured"),
         localization_group = if_else(layer == "pathogen" & accessible_annot, "accessible", "other"),
         match_stratum = paste(degree_bin, coverage_class, localization_group, sep = "|"))

sample_matched_mean <- function(frame, selected, n_perm, seed) {
  selected_frame <- frame[selected, , drop = FALSE]
  counts <- table(selected_frame$match_stratum)
  pools <- split(seq_len(nrow(frame)), frame$match_stratum)
  stopifnot(all(vapply(names(counts), function(s) length(pools[[s]]) >= counts[[s]], logical(1))))
  set.seed(seed); observed <- mean(selected_frame$opposite_side_rwr_percentile)
  null <- vapply(seq_len(n_perm), function(i) {
    idx <- unlist(lapply(names(counts), function(s) sample(pools[[s]], counts[[s]], replace = FALSE)), use.names = FALSE)
    mean(frame$opposite_side_rwr_percentile[idx])
  }, numeric(1))
  null_sd <- sd(null)
  list(summary = tibble(observed, null_mean = mean(null), null_sd = sd(null),
                        null_q025 = quantile(null, .025), null_q975 = quantile(null, .975),
                        z = safe_z(observed, null),
                        empirical_p = if_else(null_sd < 1e-12, 1,
                          (1 + sum(null >= observed - 1e-12)) / (n_perm + 1))),
       null = null)
}

n_conv_perm <- as.integer(Sys.getenv("N_CONV_PERM", "10000"))
conv_summary_rows <- list(); conv_null_rows <- list(); kk <- 1L
for (s in c("host", "pathogen")) {
  frame <- conv_nodes %>% filter(layer == s)
  for (tier_name in c("strict", "relaxed")) {
    selected <- if (tier_name == "strict") frame$strict_member else frame$relaxed_member
    ans <- sample_matched_mean(frame, selected, n_conv_perm, 20260902L + kk)
    conv_summary_rows[[kk]] <- ans$summary %>% mutate(layer = s, tier = tier_name,
                                                      selected_n = sum(selected), universe_n = nrow(frame))
    conv_null_rows[[kk]] <- tibble(layer = s, tier = tier_name, iteration = seq_along(ans$null), null_statistic = ans$null)
    kk <- kk + 1L
  }
}
conv_summary <- bind_rows(conv_summary_rows)
conv_null <- bind_rows(conv_null_rows)
host_strict <- conv_null %>% filter(layer == "host", tier == "strict") %>% pull(null_statistic)
path_strict <- conv_null %>% filter(layer == "pathogen", tier == "strict") %>% pull(null_statistic)
joint_observed <- mean(conv_summary$observed[conv_summary$tier == "strict"])
joint_null <- (host_strict + path_strict) / 2
conv_summary <- bind_rows(conv_summary, tibble(
  observed = joint_observed, null_mean = mean(joint_null), null_sd = sd(joint_null),
  null_q025 = quantile(joint_null, .025), null_q975 = quantile(joint_null, .975),
  z = safe_z(joint_observed, joint_null), empirical_p = (1 + sum(joint_null >= joint_observed)) / (length(joint_null) + 1),
  layer = "joint_two_sided", tier = "strict",
  selected_n = sum(conv_summary$selected_n[conv_summary$tier == "strict"]),
  universe_n = sum(conv_summary$universe_n[conv_summary$tier == "strict"])
))
conv_null <- bind_rows(conv_null, tibble(layer = "joint_two_sided", tier = "strict",
                                         iteration = seq_along(joint_null), null_statistic = joint_null))
write_tsv(conv_summary, file.path(out_dir, "41_two_sided_convergence.tsv"), na = "")
write_tsv(conv_null, file.path(out_dir, "42_two_sided_convergence_null.tsv"), na = "")

# Axis-level products use uncoupled same-side diffusion. They test whether the
# curated pairing aligns two independently recovered state fields.
axis_state <- cross %>%
  left_join(rwr_primary %>% select(host_node = node_id,
                                   host_diffusion = same_side_diffusion_percentile,
                                   host_direct = side_seed_percentile), by = "host_node") %>%
  left_join(rwr_primary %>% select(pathogen_node = node_id,
                                   pathogen_diffusion = same_side_diffusion_percentile,
                                   pathogen_direct = side_seed_percentile), by = "pathogen_node") %>%
  mutate(weighted_diffusion = species_adjusted_confidence * host_diffusion * pathogen_diffusion,
         weighted_direct = species_adjusted_confidence * host_direct * pathogen_direct)
write_tsv(axis_state %>% select(axis_id, interaction_name, bridge_evidence_tier,
                                species_adjusted_confidence, host_node, pathogen_node,
                                weighted_diffusion, weighted_direct),
          file.path(out_dir, "02_axis_collapsed_weighted_convergence.tsv"), na = "")

node_null_frame <- conv_nodes %>% select(node_id, layer, match_stratum,
                                         same_side_diffusion_percentile, side_seed_percentile)
make_endpoint_sampler <- function(endpoint_ids, layer_name) {
  frame <- node_null_frame[node_null_frame$layer == layer_name, , drop = FALSE]
  endpoint_ids <- unique(endpoint_ids)
  endpoint_strata <- frame$match_stratum[match(endpoint_ids, frame$node_id)]
  endpoint_groups <- split(endpoint_ids, endpoint_strata)
  pools <- split(frame$node_id, frame$match_stratum)
  stopifnot(all(vapply(names(endpoint_groups), function(s) length(pools[[s]]) >= length(endpoint_groups[[s]]), logical(1))))
  function() {
    mapped <- unlist(lapply(names(endpoint_groups), function(s) {
      setNames(sample(pools[[s]], length(endpoint_groups[[s]]), replace = FALSE), endpoint_groups[[s]])
    }), use.names = TRUE)
    mapped[endpoint_ids]
  }
}
weighted_axis_stat <- function(path_ids, host_ids, metric) {
  values <- setNames(node_null_frame[[metric]], node_null_frame$node_id)
  product <- unname(values[path_ids] * values[host_ids])
  weighted.mean(product, axis_state$species_adjusted_confidence)
}

n_label_perm <- as.integer(Sys.getenv("N_LABEL_PERM", "10000"))
observed_diffusion <- weighted_axis_stat(axis_state$pathogen_node, axis_state$host_node,
                                         "same_side_diffusion_percentile")
observed_direct <- weighted_axis_stat(axis_state$pathogen_node, axis_state$host_node,
                                      "side_seed_percentile")
set.seed(20260903)
label_diffusion <- numeric(n_label_perm); label_direct <- numeric(n_label_perm)
sample_host_endpoints <- make_endpoint_sampler(unique(axis_state$host_node), "host")
sample_path_endpoints <- make_endpoint_sampler(unique(axis_state$pathogen_node), "pathogen")
for (i in seq_len(n_label_perm)) {
  hm <- sample_host_endpoints()
  pm <- sample_path_endpoints()
  h <- unname(hm[axis_state$host_node]); q <- unname(pm[axis_state$pathogen_node])
  label_diffusion[[i]] <- weighted_axis_stat(q, h, "same_side_diffusion_percentile")
  label_direct[[i]] <- weighted_axis_stat(q, h, "side_seed_percentile")
}
label_summary <- bind_rows(
  tibble(null_model = "matched_label", statistic = "weighted_diffusion",
         observed = observed_diffusion, null_mean = mean(label_diffusion), null_sd = sd(label_diffusion),
         standardized_effect = safe_z(observed_diffusion, label_diffusion),
         empirical_p = (1 + sum(label_diffusion >= observed_diffusion)) / (n_label_perm + 1), n_perm = n_label_perm),
  tibble(null_model = "matched_label", statistic = "weighted_direct",
         observed = observed_direct, null_mean = mean(label_direct), null_sd = sd(label_direct),
         standardized_effect = safe_z(observed_direct, label_direct),
         empirical_p = (1 + sum(label_direct >= observed_direct)) / (n_label_perm + 1), n_perm = n_label_perm)
)

rewire_once <- function(path_ids, host_ids, n_swap = 200L) {
  q <- path_ids; h <- host_ids; accepted <- 0L
  for (attempt in seq_len(n_swap * 20L)) {
    ij <- sample.int(length(q), 2L); i <- ij[[1]]; j <- ij[[2]]
    if (q[[i]] == q[[j]] || h[[i]] == h[[j]]) next
    if (any(q == q[[i]] & h == h[[j]]) || any(q == q[[j]] & h == h[[i]])) next
    tmp <- h[[i]]; h[[i]] <- h[[j]]; h[[j]] <- tmp; accepted <- accepted + 1L
    if (accepted >= n_swap) break
  }
  list(path = q, host = h, accepted = accepted)
}
n_rewire_perm <- as.integer(Sys.getenv("N_REWIRE_PERM", "5000"))
set.seed(20260904)
rewire_diffusion <- numeric(n_rewire_perm); rewire_direct <- numeric(n_rewire_perm); accepted_swaps <- integer(n_rewire_perm)
for (i in seq_len(n_rewire_perm)) {
  rw <- rewire_once(axis_state$pathogen_node, axis_state$host_node)
  rewire_diffusion[[i]] <- weighted_axis_stat(rw$path, rw$host, "same_side_diffusion_percentile")
  rewire_direct[[i]] <- weighted_axis_stat(rw$path, rw$host, "side_seed_percentile")
  accepted_swaps[[i]] <- rw$accepted
}
rewire_summary <- bind_rows(
  tibble(null_model = "degree_preserving_rewire", statistic = "weighted_diffusion",
         observed = observed_diffusion, null_mean = mean(rewire_diffusion), null_sd = sd(rewire_diffusion),
         standardized_effect = safe_z(observed_diffusion, rewire_diffusion),
         empirical_p = (1 + sum(rewire_diffusion >= observed_diffusion)) / (n_rewire_perm + 1), n_perm = n_rewire_perm),
  tibble(null_model = "degree_preserving_rewire", statistic = "weighted_direct",
         observed = observed_direct, null_mean = mean(rewire_direct), null_sd = sd(rewire_direct),
         standardized_effect = safe_z(observed_direct, rewire_direct),
         empirical_p = (1 + sum(rewire_direct >= observed_direct)) / (n_rewire_perm + 1), n_perm = n_rewire_perm)
)
write_tsv(bind_rows(label_summary, rewire_summary), file.path(out_dir, "11_bridge_weighted_convergence_summary.tsv"), na = "")
write_tsv(tibble(iteration = seq_len(n_rewire_perm), weighted_diffusion = rewire_diffusion,
                 weighted_direct = rewire_direct, accepted_swaps),
          file.path(out_dir, "10_bridge_rewiring_null.tsv"), na = "")

# Re-derive modules independently on each selected within-species graph. The old
# 63/30 module universe came from a mixed score-400 graph and over-fragmented the
# proximity family relative to the S. aureus analysis contract.
host_graph <- graph_from_data_frame(host_edges %>% select(source, target), directed = FALSE)
path_graph <- graph_from_data_frame(path_edges %>% select(source, target), directed = FALSE)
derive_modules <- function(graph, layer_name, prefix, min_size = 5L) {
  set.seed(if (layer_name == "host") 20260905 else 20260906)
  fit <- cluster_leiden(graph, objective_function = "modularity", resolution = 1, n_iterations = -1)
  raw <- tibble(node_id = V(graph)$name, raw_module = as.integer(membership(fit))) %>%
    add_count(raw_module, name = "module_size") %>% filter(module_size >= min_size)
  module_order <- raw %>% distinct(raw_module, module_size) %>%
    arrange(desc(module_size), raw_module) %>% mutate(module_id = paste0(prefix, row_number()))
  raw %>% left_join(module_order %>% select(raw_module, module_id), by = "raw_module") %>%
    transmute(node_id, layer = layer_name, module_id, module_size)
}
module_membership <- bind_rows(
  derive_modules(host_graph, "host", "H_NET_"),
  derive_modules(path_graph, "pathogen", "P_NET_")
)
host_sets <- split(module_membership$node_id[module_membership$layer == "host"],
                   module_membership$module_id[module_membership$layer == "host"])
path_sets <- split(module_membership$node_id[module_membership$layer == "pathogen"],
                   module_membership$module_id[module_membership$layer == "pathogen"])
write_tsv(module_membership, file.path(out_dir, "24_network_module_membership.tsv"), na = "")

# Predefined functional programs use cached STRING names/annotations and are
# independent of the observed proximity P values.
annotated_nodes <- nodes %>% left_join(string_annotations, by = "node_id") %>%
  mutate(text = str_to_lower(paste(coalesce(as.character(preferred_name), ""),
                                   coalesce(as.character(annotation), ""))))
host_patterns <- c(
  proteasome = "proteasome",
  ribosome_translation = "ribosomal|translation|elongation factor",
  spliceosome = "splicing|spliceosome|small nuclear ribonucleoprotein|srsf|rbm",
  endocytosis_phagosome = "endocyt|phagos|clathrin|dynamin|cd2ap",
  glycolysis_metabolism = "glycolys|pyruvate|glyceraldehyde|glucose|metabol",
  cytoskeleton_microtubule = "tubulin|actin|kinesin|microtubule|clasp|mapre",
  vesicle_COPII = "sec24|sec23|copii|vesicle|rab|gtpase activ",
  chaperone_folding = "heat shock|chaperon|hsp|prohibitin",
  mediator_transcription = "mediator|polymerase|polr2|transcription factor",
  ubiquitin_cullin = "ubiquitin|cullin|e3 ligase|nedd8|cereblon|fbxw"
)
path_patterns <- c(
  surface_family = "surface|tpr|svsp|tams1|gp34|tasp|hypervariable|subtelomer",
  secretion_exported = "secreted|exported|rhoptry|microneme|dense granule",
  kinase_phosphatase = "kinase|phosphatase",
  translation_folding = "ribosomal|elongation|chaperon|heat shock|prohibitin|cyclophilin|isomerase",
  metabolism_ATP = "atp synth|dehydrogenase|oxidoreduct|metabol|proton pump",
  DNA_RNA_binding = "helicase|dna-binding|rna-binding|nucleic acid|histone",
  membrane_transport = "transporter|atp-binding cassette|membrane protein|v-type",
  PIase_PPI = "isomerase|ppiase"
)
make_functional_sets <- function(frame, patterns, prefix) {
  out <- lapply(names(patterns), function(nm) frame$node_id[str_detect(frame$text, patterns[[nm]])])
  names(out) <- paste0(prefix, names(patterns)); out[lengths(out) >= 5]
}
host_functional <- make_functional_sets(annotated_nodes %>% filter(layer == "host"), host_patterns, "H_FUN_")
path_functional <- make_functional_sets(annotated_nodes %>% filter(layer == "pathogen"), path_patterns, "P_FUN_")
functional_membership <- bind_rows(c(
  lapply(names(host_functional), function(nm) tibble(node_id = host_functional[[nm]], layer = "host", module_id = nm)),
  lapply(names(path_functional), function(nm) tibble(node_id = path_functional[[nm]], layer = "pathogen", module_id = nm))
))
all_module_membership <- bind_rows(module_membership %>% select(node_id, layer, module_id), functional_membership) %>% distinct()
host_sets <- c(host_sets, host_functional); path_sets <- c(path_sets, path_functional)
write_tsv(functional_membership, file.path(out_dir, "14_functional_module_membership.tsv"), na = "")
bridge_host <- cross$host_node; bridge_path <- cross$pathogen_node
Dh <- distances(host_graph, v = bridge_host, to = V(host_graph), weights = NA)
Dp <- distances(path_graph, v = bridge_path, to = V(path_graph), weights = NA)
degree_strata <- function(graph, reachable) {
  tibble(node_id = V(graph)$name, degree = degree(graph)) %>%
    mutate(degree_bin = ntile(rank(degree, ties.method = "average"), 10L)) %>% filter(node_id %in% reachable)
}
host_reachable <- colnames(Dh)[apply(is.finite(Dh), 2, any)]
path_reachable <- colnames(Dp)[apply(is.finite(Dp), 2, any)]
host_universe <- degree_strata(host_graph, host_reachable)
path_universe <- degree_strata(path_graph, path_reachable)
host_col <- setNames(seq_len(ncol(Dh)), colnames(Dh)); path_col <- setNames(seq_len(ncol(Dp)), colnames(Dp))
host_pools <- split(unname(host_col[host_universe$node_id]), host_universe$degree_bin)
path_pools <- split(unname(path_col[path_universe$node_id]), path_universe$degree_bin)

module_distance <- function(A, B) {
  if (!length(A) || !length(B)) return(c(distance = NA, host_to_path = NA, path_to_host = NA))
  p_to_B <- matrixStats::rowMins(Dp[, B, drop = FALSE])
  h_to_A <- matrixStats::rowMins(Dh[, A, drop = FALSE])
  dA <- matrixStats::colMins(Dh[, A, drop = FALSE] + p_to_B) + 1
  dB <- matrixStats::colMins(Dp[, B, drop = FALSE] + h_to_A) + 1
  c(distance = mean(c(mean(dA), mean(dB))), host_to_path = mean(dA), path_to_host = mean(dB))
}
make_sampler <- function(ids, universe, pools) {
  bins <- universe$degree_bin[match(ids, universe$node_id)]; counts <- table(bins)
  function() unlist(Map(function(b, n) sample(pools[[b]], as.integer(n), replace = FALSE),
                        names(counts), as.integer(counts)), use.names = FALSE)
}

n_prox_perm <- as.integer(Sys.getenv("N_PROX_PERM", "1000"))
pair_grid <- crossing(host_module = names(host_sets), pathogen_module = names(path_sets))
proximity_rows <- vector("list", nrow(pair_grid)); set.seed(20260902)
for (i in seq_len(nrow(pair_grid))) {
  hm <- pair_grid$host_module[[i]]; pm <- pair_grid$pathogen_module[[i]]
  A0 <- unique(host_sets[[hm]]); B0 <- unique(path_sets[[pm]])
  A <- intersect(A0, host_universe$node_id); B <- intersect(B0, path_universe$node_id)
  if (length(A) < 5 || length(B) < 5) next
  Aidx <- unname(host_col[A]); Bidx <- unname(path_col[B])
  sample_host <- make_sampler(A, host_universe, host_pools)
  sample_path <- make_sampler(B, path_universe, path_pools)
  obs <- module_distance(Aidx, Bidx)
  null <- vapply(seq_len(n_prox_perm), function(j) module_distance(sample_host(), sample_path())[["distance"]], numeric(1))
  proximity_rows[[i]] <- tibble(
    host_module = hm, pathogen_module = pm, host_size = length(A0), pathogen_size = length(B0),
    host_reachable_n = length(A), pathogen_reachable_n = length(B),
    closest_distance = obs[["distance"]], host_to_path_distance = obs[["host_to_path"]],
    path_to_host_distance = obs[["path_to_host"]], null_mean = mean(null), null_sd = sd(null),
    proximity_z = safe_z(obs[["distance"]], null),
    empirical_p = (1 + sum(null <= obs[["distance"]])) / (n_prox_perm + 1)
  )
}
proximity <- bind_rows(proximity_rows) %>%
  mutate(module_family = if_else(str_detect(host_module, "_NET_") & str_detect(pathogen_module, "_NET_"),
                                 "network_modules", "functional_modules")) %>%
  group_by(module_family) %>% mutate(proximity_fdr = p.adjust(empirical_p, method = "BH")) %>% ungroup() %>% mutate(
  proximity_support = !is.na(proximity_z) & proximity_z < 0 & proximity_fdr <= 0.05
)
write_tsv(proximity, file.path(out_dir, "22_module_proximity_all_pairs.tsv"), na = "")

# A real PCSF sweep on the heterogeneous graph, including main cross edges.
g_pcsf <- graph_from_data_frame(edges %>% transmute(source, target, cost = 1 / pmax(weight, 0.05)),
                                directed = FALSE, vertices = nodes %>% transmute(name = node_id))
E(g_pcsf)$weight <- E(g_pcsf)$cost
terminal_frame <- nodes %>% filter(seed_raw > 0, evidence_n >= 4) %>%
  mutate(prize_raw = side_seed_percentile + 0.75 * strict_member + 0.25 * relaxed_member)
terminals <- setNames(terminal_frame$prize_raw, terminal_frame$node_id)
pcsf_grid <- crossing(w = c(2, 4), b = c(1, 2, 4), mu = c(1e-4, 5e-4, 1e-3))
pcsf_nodes <- list(); pcsf_edges <- list(); pcsf_status <- list()
for (i in seq_len(nrow(pcsf_grid))) {
  pars <- pcsf_grid[i, ]; run_id <- sprintf("w%s_b%s_mu%s", pars$w, pars$b, pars$mu)
  fit <- tryCatch(PCSF::PCSF(g_pcsf, terminals, w = pars$w, b = pars$b, mu = pars$mu), error = function(e) e)
  if (inherits(fit, "error")) {
    pcsf_status[[i]] <- tibble(run_id, w = pars$w, b = pars$b, mu = pars$mu,
                               success = FALSE, node_n = 0L, edge_n = 0L, message = conditionMessage(fit))
  } else {
    en <- as_data_frame(fit, what = "edges")
    pcsf_nodes[[i]] <- tibble(run_id, node_id = as.character(V(fit)$name))
    pcsf_edges[[i]] <- en %>% transmute(run_id, source = pmin(as.character(from), as.character(to)),
                                        target = pmax(as.character(from), as.character(to)))
    pcsf_status[[i]] <- tibble(run_id, w = pars$w, b = pars$b, mu = pars$mu,
                               success = TRUE, node_n = vcount(fit), edge_n = ecount(fit), message = "")
  }
}
pcsf_status <- bind_rows(pcsf_status); successful_n <- sum(pcsf_status$success)
if (!successful_n) stop("All PCSF parameter runs failed.", call. = FALSE)
pcsf_node_runs <- bind_rows(pcsf_nodes); pcsf_edge_runs <- bind_rows(pcsf_edges)
pcsf_node_frequency <- pcsf_node_runs %>% dplyr::count(node_id, name = "selected_runs") %>%
  mutate(selection_frequency = selected_runs / successful_n) %>% left_join(nodes, by = "node_id")
pcsf_edge_frequency <- pcsf_edge_runs %>% dplyr::count(source, target, name = "selected_runs") %>%
  mutate(selection_frequency = selected_runs / successful_n)

pcsf_axis_runs <- lapply(pcsf_status$run_id[pcsf_status$success], function(run_id) {
  vn <- pcsf_node_runs %>% filter(.data$run_id == !!run_id) %>% pull(node_id)
  en <- pcsf_edge_runs %>% filter(.data$run_id == !!run_id)
  rg <- graph_from_data_frame(en %>% select(source, target), directed = FALSE, vertices = tibble(name = vn))
  membership <- components(rg)$membership
  cross %>% transmute(
    axis_id, run_id,
    connected = pathogen_node %in% names(membership) & host_node %in% names(membership) &
      membership[pathogen_node] == membership[host_node]
  )
}) %>% bind_rows()
pcsf_axis_connection <- pcsf_axis_runs %>% group_by(axis_id) %>%
  summarise(component_path_selected_runs = sum(connected),
            component_path_frequency = mean(connected), .groups = "drop")
write_tsv(pcsf_status, file.path(out_dir, "35_pcsf_parameter_stability.tsv"), na = "")
write_tsv(pcsf_node_runs, file.path(out_dir, "33_pcsf_node_runs.tsv"), na = "")
write_tsv(pcsf_edge_runs, file.path(out_dir, "31_pcsf_edge_runs.tsv"), na = "")
write_tsv(pcsf_node_frequency, file.path(out_dir, "34_pcsf_node_selection_frequency.tsv"), na = "")
write_tsv(pcsf_edge_frequency, file.path(out_dir, "32_pcsf_edge_selection_frequency.tsv"), na = "")
write_tsv(pcsf_axis_connection, file.path(out_dir, "30_pcsf_axis_connection_frequency.tsv"), na = "")

# Map independent algorithm outputs back to the curated axes.
axis_modules <- cross %>% select(axis_id, host_node, pathogen_node) %>%
  left_join(all_module_membership %>% filter(layer == "host") %>% select(host_node = node_id, host_module = module_id),
            by = "host_node", relationship = "many-to-many") %>%
  left_join(all_module_membership %>% filter(layer == "pathogen") %>% select(pathogen_node = node_id, pathogen_module = module_id),
            by = "pathogen_node", relationship = "many-to-many") %>%
  left_join(proximity, by = c("host_module", "pathogen_module"))
best_proximity <- axis_modules %>% mutate(valid = coalesce(proximity_support, FALSE), rank_value = if_else(valid, proximity_fdr, Inf)) %>%
  arrange(axis_id, rank_value, empirical_p) %>% group_by(axis_id) %>% slice_head(n = 1) %>% ungroup() %>%
  transmute(axis_id, proximity_support = valid, proximity_host_module = host_module,
            proximity_pathogen_module = pathogen_module, proximity_distance = closest_distance,
            proximity_z, proximity_p = empirical_p, proximity_fdr)

node_support <- rwr_primary %>% select(node_id, side_seed_raw, side_seed_percentile, opposite_side_rwr_percentile)
pcsf_node_support <- pcsf_node_frequency %>% select(node_id, selection_frequency)
pcsf_edge_support <- pcsf_edge_frequency %>% rename(pcsf_edge_frequency = selection_frequency)
pcsf_success_fraction <- successful_n / nrow(pcsf_grid)

axis_rows <- cross %>%
  left_join(node_support %>% rename(pathogen_node = node_id, pathogen_seed_raw = side_seed_raw,
                                    pathogen_seed_percentile = side_seed_percentile,
                                    pathogen_cross_gain_percentile = opposite_side_rwr_percentile), by = "pathogen_node") %>%
  left_join(node_support %>% rename(host_node = node_id, host_seed_raw = side_seed_raw,
                                    host_seed_percentile = side_seed_percentile,
                                    host_cross_gain_percentile = opposite_side_rwr_percentile), by = "host_node") %>%
  left_join(pcsf_node_support %>% rename(pathogen_node = node_id, pathogen_pcsf_frequency = selection_frequency), by = "pathogen_node") %>%
  left_join(pcsf_node_support %>% rename(host_node = node_id, host_pcsf_frequency = selection_frequency), by = "host_node") %>%
  mutate(pair_a = pmin(pathogen_node, host_node), pair_b = pmax(pathogen_node, host_node)) %>%
  left_join(pcsf_edge_support, by = c("pair_a" = "source", "pair_b" = "target")) %>%
  left_join(host_ev %>% select(host_node = string_id, host_evidence_n = evidence_count), by = "host_node") %>%
  left_join(path_ev %>% select(pathogen_node = string_id, pathogen_evidence_n = evidence_count), by = "pathogen_node")

axis_support <- axis_rows %>% mutate(
  across(c(pathogen_pcsf_frequency, host_pcsf_frequency, pcsf_edge_frequency), ~coalesce(.x, 0)),
  host_measured_gate = coalesce(host_evidence_n, 0) >= 1,
  pathogen_measured_gate = coalesce(pathogen_evidence_n, 0) >= 1,
  host_recurrence_gate = coalesce(host_evidence_n, 0) >= 2,
  pathogen_recurrence_gate = coalesce(pathogen_evidence_n, 0) >= 2,
  rwr_endpoint_min_percentile = pmin(pathogen_cross_gain_percentile, host_cross_gain_percentile),
  rwr_support = pathogen_cross_gain_percentile >= 0.95 & host_cross_gain_percentile >= 0.95,
  pcsf_endpoint_min_frequency = pmin(pathogen_pcsf_frequency, host_pcsf_frequency),
  pcsf_support = pcsf_success_fraction >= 0.8 & (pcsf_endpoint_min_frequency >= 0.5 | pcsf_edge_frequency >= 0.5),
  strict_evidence_gate = main_algorithm_bridge & host_recurrence_gate & pathogen_recurrence_gate,
  relaxed_evidence_gate = main_algorithm_bridge & host_measured_gate & pathogen_measured_gate
) %>% left_join(best_proximity, by = "axis_id") %>% left_join(pcsf_axis_connection, by = "axis_id") %>%
  mutate(
    proximity_support = coalesce(proximity_support, FALSE),
    component_path_frequency = coalesce(component_path_frequency, 0),
    algorithm_count = as.integer(rwr_support) + as.integer(proximity_support) + as.integer(pcsf_support),
    consensus_rule_pass = algorithm_count >= 2 | (algorithm_count >= 1 & direct_physical_edge),
    consensus_tier = case_when(
      strict_evidence_gate & consensus_rule_pass ~ "strict",
      relaxed_evidence_gate & consensus_rule_pass ~ "relaxed",
      TRUE ~ "context_only"
    )
  ) %>% arrange(factor(consensus_tier, levels = c("strict", "relaxed", "context_only")),
                desc(algorithm_count), desc(rwr_endpoint_min_percentile))

# Joint threshold grid follows the same axis-level statistic as the primary call.
threshold_grid <- crossing(rwr_threshold = c(0.90, 0.95, 0.975),
                           proximity_fdr_threshold = c(0.01, 0.05, 0.10),
                           pcsf_threshold = c(0.30, 0.50, 0.70),
                           confidence_threshold = c(0.70, 0.75, 0.80, 0.90)) %>%
  mutate(threshold_id = row_number())
threshold_rows <- crossing(axis_id = axis_support$axis_id, threshold_id = threshold_grid$threshold_id) %>%
  left_join(threshold_grid, by = "threshold_id") %>% left_join(axis_support, by = "axis_id") %>%
  mutate(
    threshold_rwr_support = rwr_endpoint_min_percentile >= rwr_threshold,
    threshold_proximity_support = !is.na(proximity_fdr) & proximity_fdr <= proximity_fdr_threshold & proximity_z < 0,
    threshold_pcsf_support = pcsf_success_fraction >= 0.8 &
      (pcsf_endpoint_min_frequency >= pcsf_threshold | pcsf_edge_frequency >= pcsf_threshold),
    threshold_algorithm_count = as.integer(threshold_rwr_support) + as.integer(threshold_proximity_support) + as.integer(threshold_pcsf_support),
    threshold_consensus_pass = threshold_algorithm_count >= 2 | (threshold_algorithm_count >= 1 & direct_physical_edge),
    threshold_strict_evidence_gate = strict_evidence_gate & species_adjusted_confidence >= confidence_threshold,
    threshold_strict_selected = threshold_strict_evidence_gate & threshold_consensus_pass
  )
robustness <- threshold_rows %>% group_by(axis_id) %>% summarise(
  threshold_run_n = n(), threshold_consensus_frequency = mean(threshold_consensus_pass),
  threshold_strict_selection_frequency = mean(threshold_strict_selected),
  threshold_rwr_support_frequency = mean(threshold_rwr_support),
  threshold_proximity_support_frequency = mean(threshold_proximity_support),
  threshold_pcsf_support_frequency = mean(threshold_pcsf_support), .groups = "drop"
) %>% left_join(
  rwr_all %>% filter(node_id %in% c(cross$host_node, cross$pathogen_node)) %>%
    select(node_id, restart, delta, opposite_side_rwr_percentile) %>%
    inner_join(cross %>% select(axis_id, host_node, pathogen_node), by = c("node_id" = "host_node")) %>%
    rename(host_score = opposite_side_rwr_percentile) %>%
    inner_join(rwr_all %>% select(node_id, restart, delta, opposite_side_rwr_percentile),
               by = c("pathogen_node" = "node_id", "restart", "delta"), suffix = c("", "_path")) %>%
    group_by(axis_id) %>% summarise(rwr_support_frequency = mean(pmin(host_score, opposite_side_rwr_percentile) >= 0.95), .groups = "drop"),
  by = "axis_id"
)
axis_support <- axis_support %>% left_join(robustness, by = "axis_id") %>% mutate(
  robustness_tier = case_when(
    consensus_tier == "strict" & algorithm_count >= 2 & threshold_strict_selection_frequency >= 0.70 & rwr_support_frequency >= 0.50 ~ "robust_strict",
    consensus_tier == "strict" ~ "supported_strict",
    consensus_tier == "relaxed" & threshold_consensus_frequency >= 0.70 & rwr_support_frequency >= 0.50 ~ "stable_relaxed",
    consensus_tier == "relaxed" ~ "supported_relaxed",
    TRUE ~ "context_only"
  )
)

write_tsv(cross_all, file.path(out_dir, "13_curated_bridge_axis_contract.tsv"), na = "")
write_tsv(axis_support, file.path(out_dir, "18_interface_axis_algorithm_support.tsv"), na = "")
write_tsv(axis_support %>% filter(consensus_tier %in% c("strict", "relaxed")),
          file.path(out_dir, "12_consensus_interface_axes.tsv"), na = "")
write_tsv(axis_support %>% filter(consensus_tier == "strict"), file.path(out_dir, "19_interface_core.tsv"), na = "")
write_tsv(threshold_rows, file.path(out_dir, "05_axis_threshold_selection_grid.tsv"), na = "")
write_tsv(robustness, file.path(out_dir, "04_axis_robustness_summary.tsv"), na = "")

primary_check <- rwr_all %>% filter(restart == primary_restart, delta == primary_delta) %>%
  select(node_id, grid = opposite_side_rwr_percentile) %>%
  inner_join(rwr_primary %>% select(node_id, primary = opposite_side_rwr_percentile), by = "node_id")
stopifnot(
  n_distinct(cross$axis_id) == nrow(cross),
  all(cross$bridge_evidence_tier %in% main_tiers),
  max(abs(primary_check$grid - primary_check$primary)) < 1e-12,
  all(between(axis_support$rwr_endpoint_min_percentile, 0, 1)),
  all(between(axis_support$component_path_frequency, 0, 1)),
  nrow(threshold_rows) == nrow(axis_support) * nrow(threshold_grid)
)

message("S. aureus-aligned interface analysis complete")
message("  main / all mapped curated axes: ", nrow(cross), " / ", nrow(cross_all))
message("  heterogeneous graph: ", nrow(nodes), " nodes; ", nrow(edges), " edges")
message("  strict joint convergence empirical P: ", signif(conv_summary$empirical_p[conv_summary$layer == "joint_two_sided"], 4))
message("  module pairs / FDR-supported: ", nrow(proximity), " / ", sum(proximity$proximity_support))
message("  PCSF successful runs: ", successful_n, "/", nrow(pcsf_grid))
message("  strict / relaxed / context axes: ", paste(table(factor(axis_support$consensus_tier,
        levels = c("strict", "relaxed", "context_only"))), collapse = " / "))
