# Revision rerun: drop 1.25 seed bonus. PCSF terminals = evidence_n >= 4.
# Writes 02_results/04_intermediate/03_revision_rerun/. Does not overwrite 01_algorithms.
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(Matrix)
  library(PCSF)
  library(readr)
  library(stringr)
  library(tidyr)
})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1]])), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}
out <- file.path(root, "03_algorithm_tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_tab <- function(path) read_tsv(path, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
flag <- function(x) str_to_lower(coalesce(as.character(x), "")) %in% c("yes", "true", "t", "1", "y")
rank01 <- function(x) {
  z <- rep(0, length(x)); ok <- is.finite(x) & x > 0
  if (sum(ok) == 1L) z[ok] <- 1
  if (sum(ok) > 1L) z[ok] <- rank(x[ok], ties.method = "average") / sum(ok)
  z
}
read_network <- function(path, layer) {
  read_tsv(path, col_names = c("source", "target"), show_col_types = FALSE, progress = FALSE) %>%
    filter(source != target) %>% distinct(source, target) %>%
    mutate(weight = 1, edge_scope = paste0("within_", layer))
}
cat(format(Sys.time()), "loading\n"); flush.console()
host_edges <- read_network(file.path(root, "02_processed_inputs/01_networks/01_host_string_network_700.tsv"), "host")
path_edges <- read_network(file.path(root, "02_processed_inputs/01_networks/02_parasite_string_network_400.tsv"), "pathogen")
bridge_raw <- read_tab(file.path(root, "02_processed_inputs/01_networks/03_bridge_edges_string_mapped.tsv"))
host_ev <- read_tab(file.path(root, "02_processed_inputs/02_evidence/01_host_recurrence.tsv"))
path_ev <- read_tab(file.path(root, "02_processed_inputs/02_evidence/02_parasite_accessible.tsv"))
stored <- read_tab(file.path(root, "03_algorithm_tables/36_rwr_cross_layer_gain.tsv"))
axes <- read_tab(file.path(root, "03_algorithm_tables/06_biological_axis_consensus.tsv")) %>%
  filter(flag(retained))
tier_confidence <- c(direct_validated = 1, direct_literature = 0.85, cocomplex = 0.80,
                     y2h_hit_unconfirmed = 0.60, paralog_inferred_cocomplex = 0.45)
main_tiers <- c("direct_validated", "direct_literature", "cocomplex")
cross <- bridge_raw %>%
  filter(flag(both_mapped), !is.na(parasite_string_id), !is.na(host_string_id)) %>%
  transmute(
    axis_id = edge_id,
    pathogen_node = parasite_string_id, host_node = host_string_id,
    tier = tier,
    species_adjusted_confidence = unname(tier_confidence[tier]),
    main_algorithm_bridge = tier %in% main_tiers
  ) %>%
  filter(!is.na(species_adjusted_confidence)) %>%
  mutate(
    host_endpoint_in_network = host_node %in% c(host_edges$source, host_edges$target),
    pathogen_endpoint_in_network = pathogen_node %in% c(path_edges$source, path_edges$target),
    main_algorithm_bridge = main_algorithm_bridge & host_endpoint_in_network & pathogen_endpoint_in_network
  ) %>% filter(main_algorithm_bridge)
cross_edges <- cross %>% transmute(source = pathogen_node, target = host_node,
                                   weight = species_adjusted_confidence, edge_scope = "cross")
edges <- bind_rows(host_edges, path_edges, cross_edges)
node_ids <- sort(unique(c(edges$source, edges$target)))
node_index <- setNames(seq_along(node_ids), node_ids)
deg <- bind_rows(
  edges %>% transmute(node_id = source),
  edges %>% transmute(node_id = target)
) %>% count(node_id, name = "degree")
base <- tibble(node_id = node_ids, layer = if_else(str_starts(node_ids, "9913."), "host", "pathogen")) %>%
  left_join(host_ev %>% transmute(node_id = string_id, host_evidence_n = as.numeric(evidence_count)), by = "node_id") %>%
  left_join(path_ev %>% transmute(node_id = string_id, pathogen_evidence_n = as.numeric(evidence_count),
                                  accessible_annot = flag(accessible_annot)), by = "node_id") %>%
  left_join(deg, by = "node_id") %>%
  mutate(
    host_evidence_n = coalesce(host_evidence_n, 0),
    pathogen_evidence_n = coalesce(pathogen_evidence_n, 0),
    accessible_annot = coalesce(accessible_annot, FALSE),
    evidence_n = if_else(layer == "host", host_evidence_n, pathogen_evidence_n),
    strict_member = if_else(layer == "host", host_evidence_n >= 2, pathogen_evidence_n >= 2),
    relaxed_member = evidence_n >= 1,
    degree = coalesce(degree, 0)
  )
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
  p0 <- p0 / sum(p0); state <- p0
  for (iter in seq_len(max_iter)) {
    next_state <- as.numeric((1 - restart) * (M %*% state) + restart * p0)
    if (sum(abs(next_state - state)) < tol) return(list(score = next_state, iterations = iter))
    state <- next_state
  }
  list(score = state, iterations = max_iter)
}
nodes <- base %>% mutate(
  seed_raw = if_else(
    layer == "host",
    (host_evidence_n / max(host_evidence_n))^2,
    (pathogen_evidence_n / max(pathogen_evidence_n))^2
  )
) %>% group_by(layer) %>% mutate(side_seed_percentile = rank01(seed_raw)) %>% ungroup()
seed_path <- ifelse(nodes$layer == "pathogen", nodes$side_seed_percentile^2, 0)
seed_host <- ifelse(nodes$layer == "host", nodes$side_seed_percentile^2, 0)
cat(format(Sys.time()), "rwr\n"); flush.console()
M <- build_transition(0.5)
M0 <- build_transition(0)
fit_path <- run_rwr(M, seed_path, 0.3)
fit_host <- run_rwr(M, seed_host, 0.3)
fit_path_u <- run_rwr(M0, seed_path, 0.3)
fit_host_u <- run_rwr(M0, seed_host, 0.3)
rwr <- nodes %>% transmute(
  node_id, layer, degree, strict_member, relaxed_member, evidence_n, accessible_annot,
  side_seed_raw = seed_raw, side_seed_percentile,
  from_host_rwr = fit_host$score, from_pathogen_rwr = fit_path$score,
  uncoupled_same_side_rwr = if_else(layer == "host", fit_host_u$score, fit_path_u$score),
  opposite_side_rwr = if_else(layer == "host", fit_path$score, fit_host$score)
) %>% group_by(layer) %>% mutate(
  opposite_side_rwr_percentile = rank01(opposite_side_rwr),
  same_side_diffusion_percentile = rank01(uncoupled_same_side_rwr)
) %>% ungroup()
write_tsv(rwr, file.path(out, "rerun_rwr_nodes.tsv"))
cmp <- stored %>% select(node_id, stored_opp = opposite_side_rwr_percentile) %>%
  left_join(rwr %>% select(node_id, new_opp = opposite_side_rwr_percentile), by = "node_id")
axis_rwr <- axes %>%
  transmute(biological_axis_id, consensus_tier, host_node, pathogen_node,
            stored_host = host_cross_gain_percentile, stored_path = pathogen_cross_gain_percentile) %>%
  left_join(rwr %>% select(host_node = node_id, host_pct = opposite_side_rwr_percentile), by = "host_node") %>%
  left_join(rwr %>% select(pathogen_node = node_id, path_pct = opposite_side_rwr_percentile), by = "pathogen_node") %>%
  mutate(min_pct = pmin(host_pct, path_pct), rwr_support = min_pct >= 0.95)
write_tsv(axis_rwr, file.path(out, "rerun_axis_rwr.tsv"))
rec_n <- rwr %>% filter(side_seed_percentile < 0.50, opposite_side_rwr_percentile >= 0.95) %>% nrow()
cat(format(Sys.time()), "RWR vs stored max abs", max(abs(cmp$new_opp - cmp$stored_opp), na.rm = TRUE), "\n")
cat("axes losing 0.95:", sum(axis_rwr$rwr_support == FALSE), "recruitment", rec_n, "\n")
print(axis_rwr %>% filter(consensus_tier == "strict") %>%
        select(biological_axis_id, stored_host, host_pct, stored_path, path_pct, min_pct))
flush.console()

terminal_frame <- rwr %>% filter(side_seed_raw > 0, evidence_n >= 4) %>%
  mutate(prize_raw = side_seed_percentile + 0.75 * strict_member + 0.25 * relaxed_member)
cat(format(Sys.time()), "terminals", nrow(terminal_frame),
    "host", sum(terminal_frame$layer == "host"),
    "pathogen", sum(terminal_frame$layer == "pathogen"), "\n"); flush.console()
write_tsv(terminal_frame, file.path(out, "rerun_terminals_ev4.tsv"))
terminals <- setNames(terminal_frame$prize_raw, terminal_frame$node_id)
cat(format(Sys.time()), "pcsf graph\n"); flush.console()
g_pcsf <- graph_from_data_frame(edges %>% transmute(source, target, cost = 1 / pmax(weight, 0.05)),
                                directed = FALSE, vertices = nodes %>% transmute(name = node_id))
E(g_pcsf)$weight <- E(g_pcsf)$cost
pcsf_grid <- crossing(w = c(2, 4), b = c(1, 2, 4), mu = c(1e-4, 5e-4, 1e-3))
pcsf_nodes <- list(); pcsf_edges <- list(); pcsf_status <- list()
for (i in seq_len(nrow(pcsf_grid))) {
  pars <- pcsf_grid[i, ]
  run_id <- sprintf("w%s_b%s_mu%s", pars$w, pars$b, pars$mu)
  cat(format(Sys.time()), "pcsf", run_id, "\n"); flush.console()
  fit <- tryCatch(PCSF::PCSF(g_pcsf, terminals, w = pars$w, b = pars$b, mu = pars$mu), error = function(e) e)
  if (inherits(fit, "error")) {
    pcsf_status[[i]] <- tibble(run_id, w = pars$w, b = pars$b, mu = pars$mu,
                               success = FALSE, node_n = 0L, edge_n = 0L, message = conditionMessage(fit))
  } else {
    en <- as_data_frame(fit, what = "edges")
    pcsf_nodes[[i]] <- tibble(run_id, node_id = as.character(V(fit)$name))
    pcsf_edges[[i]] <- en %>% transmute(run_id, src = pmin(as.character(from), as.character(to)),
                                        tgt = pmax(as.character(from), as.character(to)))
    pcsf_status[[i]] <- tibble(run_id, w = pars$w, b = pars$b, mu = pars$mu,
                               success = TRUE, node_n = vcount(fit), edge_n = ecount(fit), message = "")
  }
  write_tsv(bind_rows(pcsf_status), file.path(out, "rerun_pcsf_status.tsv"))
}
pcsf_status <- bind_rows(pcsf_status)
successful_n <- sum(pcsf_status$success)
pcsf_node_runs <- bind_rows(pcsf_nodes)
pcsf_edge_runs <- bind_rows(pcsf_edges)
pcsf_node_frequency <- pcsf_node_runs %>% count(node_id, name = "selected_runs") %>%
  mutate(selection_frequency = selected_runs / successful_n)
pcsf_edge_frequency <- pcsf_edge_runs %>% count(src, tgt, name = "selected_runs") %>%
  mutate(selection_frequency = selected_runs / successful_n)
write_tsv(pcsf_status, file.path(out, "rerun_pcsf_status.tsv"))
write_tsv(pcsf_node_frequency, file.path(out, "rerun_pcsf_node_frequency.tsv"))
write_tsv(pcsf_edge_frequency, file.path(out, "rerun_pcsf_edge_frequency.tsv"))
core <- axes %>% filter(consensus_tier == "strict") %>%
  mutate(src = pmin(pathogen_node, host_node), tgt = pmax(pathogen_node, host_node))
edge_hits <- core %>% left_join(pcsf_edge_frequency %>% select(src, tgt, freq = selection_frequency),
                               by = c("src", "tgt")) %>%
  mutate(freq = coalesce(freq, 0))
backbone_n <- sum(pcsf_node_frequency$selection_frequency >= 0.5)
backbone_e <- sum(pcsf_edge_frequency$selection_frequency >= 0.5)
every_n <- sum(pcsf_node_frequency$selection_frequency >= 1 - 1e-12)
summary <- c(
  paste("time", format(Sys.time())),
  paste("rwr_vs_stored_max_abs", max(abs(cmp$new_opp - cmp$stored_opp), na.rm = TRUE)),
  paste("strict_axes_below_0.95", sum(axis_rwr$consensus_tier == "strict" & !axis_rwr$rwr_support)),
  paste("retained_axes_below_0.95", sum(!axis_rwr$rwr_support)),
  paste("recruitment_n", rec_n),
  paste("terminals_ev4", nrow(terminal_frame)),
  paste("pcsf_success", successful_n, "of", nrow(pcsf_status)),
  paste("backbone_freq_ge_0.5_nodes", backbone_n, "edges", backbone_e),
  paste("nodes_in_every_run", every_n),
  paste("unique_node_n_values", paste(sort(unique(pcsf_status$node_n[pcsf_status$success])), collapse = ",")),
  paste("core_edge_freq", paste(edge_hits$biological_axis_id, round(edge_hits$freq, 3), collapse = "; "))
)
writeLines(summary, file.path(out, "rerun_no125_term4_summary.txt"))
print(pcsf_status)
print(edge_hits %>% select(biological_axis_id, freq))
cat(paste(summary, collapse = "\n"), "\n")
cat(format(Sys.time()), "done R\n")
