# Patch revision 02_results to no-1.25 primary RWR, 12-grid Spearman, two-sided P.
# Writes only under revision/02_results. Does not touch original/.
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(Matrix)
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
alg <- file.path(root, "03_algorithm_tables")
tmp <- file.path(root, "03_algorithm_tables")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
read_tab <- function(path) read_tsv(path, show_col_types = FALSE, progress = FALSE, na = c("", "NA", "NaN"))
flag <- function(x) str_to_lower(coalesce(as.character(x), "")) %in% c("yes", "true", "t", "1", "y")
rank01 <- function(x) {
  z <- rep(0, length(x)); ok <- is.finite(x) & x > 0
  if (sum(ok) == 1L) z[ok] <- 1
  if (sum(ok) > 1L) z[ok] <- rank(x[ok], ties.method = "average") / sum(ok)
  z
}
safe_z <- function(observed, null) if (sd(null) > 0) (observed - mean(null)) / sd(null) else NA_real_
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
tier_confidence <- c(direct_validated = 1, direct_literature = 0.85, cocomplex = 0.80,
                     y2h_hit_unconfirmed = 0.60, paralog_inferred_cocomplex = 0.45)
main_tiers <- c("direct_validated", "direct_literature", "cocomplex")
cross <- bridge_raw %>%
  filter(flag(both_mapped), !is.na(parasite_string_id), !is.na(host_string_id)) %>%
  transmute(
    pathogen_node = parasite_string_id, host_node = host_string_id, tier = tier,
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
deg <- bind_rows(edges %>% transmute(node_id = source), edges %>% transmute(node_id = target)) %>%
  count(node_id, name = "degree")
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
seed_joint <- 0.5 * seed_host / sum(seed_host) + 0.5 * seed_path / sum(seed_path)
stopifnot(max(nodes$seed_raw[nodes$layer == "pathogen" & nodes$accessible_annot]) <= 1 + 1e-12)

cat(format(Sys.time()), "primary RWR\n"); flush.console()
M <- build_transition(0.5)
M0 <- build_transition(0)
fit_path <- run_rwr(M, seed_path, 0.3)
fit_host <- run_rwr(M, seed_host, 0.3)
fit_joint <- run_rwr(M, seed_joint, 0.3)
fit_path_u <- run_rwr(M0, seed_path, 0.3)
fit_host_u <- run_rwr(M0, seed_host, 0.3)
rwr_primary <- nodes %>% transmute(
  node_id, layer, degree, strict_member, relaxed_member, evidence_n, accessible_annot,
  side_seed_raw = seed_raw, side_seed_percentile,
  joint_rwr = fit_joint$score, from_host_rwr = fit_host$score, from_pathogen_rwr = fit_path$score,
  uncoupled_same_side_rwr = if_else(layer == "host", fit_host_u$score, fit_path_u$score),
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
cat("recruitment", nrow(recruitment),
    "host", sum(recruitment$layer == "host"),
    "pathogen", sum(recruitment$layer == "pathogen"), "\n"); flush.console()

cat(format(Sys.time()), "12-grid\n"); flush.console()
restart_grid <- c(0.2, 0.3, 0.5, 0.7)
delta_grid <- c(0.25, 0.5, 0.75)
primary_ranks <- rwr_primary %>% group_by(layer) %>%
  mutate(primary_rank = min_rank(desc(opposite_side_rwr))) %>% ungroup() %>%
  select(node_id, layer, primary_percentile = opposite_side_rwr_percentile, primary_rank)
stab_rows <- list(); k <- 1L
for (delta in delta_grid) {
  Md <- build_transition(delta)
  for (restart in restart_grid) {
    host_fit <- run_rwr(Md, seed_host, restart)
    path_fit <- run_rwr(Md, seed_path, restart)
    tmp <- nodes %>% transmute(
      node_id, layer,
      opposite_side_rwr = if_else(layer == "host", path_fit$score, host_fit$score)
    ) %>% group_by(layer) %>%
      mutate(opposite_side_rwr_percentile = rank01(opposite_side_rwr),
             opposite_side_rwr_rank = min_rank(desc(opposite_side_rwr))) %>% ungroup() %>%
      left_join(primary_ranks, by = c("node_id", "layer"))
    z <- lapply(c("host", "pathogen"), function(s) {
      y <- tmp %>% filter(layer == s)
      a <- y$node_id[y$opposite_side_rwr_rank <= 200]
      b <- y$node_id[y$primary_rank <= 200]
      tibble(layer = s,
             spearman = cor(y$opposite_side_rwr_percentile, y$primary_percentile, method = "spearman"),
             top200_jaccard = length(intersect(a, b)) / length(union(a, b)))
    }) %>% bind_rows()
    stab_rows[[k]] <- tibble(
      restart, delta,
      spearman_to_primary = min(z$spearman),
      top200_jaccard = min(z$top200_jaccard),
      host_spearman = z$spearman[z$layer == "host"],
      pathogen_spearman = z$spearman[z$layer == "pathogen"]
    )
    k <- k + 1L
    cat(format(Sys.time()), "grid", restart, delta, "spearman", min(z$spearman), "\n"); flush.console()
  }
}
rwr_rank_stability <- bind_rows(stab_rows)
cat("min Spearman", min(rwr_rank_stability$spearman_to_primary), "\n"); flush.console()

cat(format(Sys.time()), "two-sided 10000\n"); flush.console()
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
  list(summary = tibble(observed, null_mean = mean(null), null_sd = sd(null),
                        null_q025 = unname(quantile(null, .025)), null_q975 = unname(quantile(null, .975)),
                        z = safe_z(observed, null),
                        empirical_p = if_else(sd(null) < 1e-12, 1,
                          (1 + sum(null >= observed - 1e-12)) / (n_perm + 1))),
       null = null)
}
n_conv_perm <- 10000L
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
    cat(format(Sys.time()), s, tier_name, "P", ans$summary$empirical_p, "z", ans$summary$z, "\n"); flush.console()
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
  null_q025 = unname(quantile(joint_null, .025)), null_q975 = unname(quantile(joint_null, .975)),
  z = safe_z(joint_observed, joint_null),
  empirical_p = (1 + sum(joint_null >= joint_observed)) / (length(joint_null) + 1),
  layer = "joint_two_sided", tier = "strict",
  selected_n = sum(conv_summary$selected_n[conv_summary$tier == "strict"]),
  universe_n = sum(conv_summary$universe_n[conv_summary$tier == "strict"])
))
conv_null <- bind_rows(conv_null, tibble(layer = "joint_two_sided", tier = "strict",
                                         iteration = seq_along(joint_null), null_statistic = joint_null))

write_tsv(rwr_primary, file.path(alg, "36_rwr_cross_layer_gain.tsv"), na = "")
write_tsv(recruitment, file.path(alg, "37_rwr_node_recruitment.tsv"), na = "")
write_tsv(rwr_rank_stability, file.path(alg, "40_rwr_rank_stability.tsv"), na = "")
write_tsv(conv_summary, file.path(alg, "41_two_sided_convergence.tsv"), na = "")
write_tsv(conv_null, file.path(alg, "42_two_sided_convergence_null.tsv"), na = "")
write_tsv(rwr_rank_stability, file.path(root, "02_results/01_tables/03_publication_tables/02_supplementary/S5_network/S5B_rwr_stability.tsv"), na = "")
writeLines(c(
  paste("time", format(Sys.time())),
  paste("recruitment", nrow(recruitment)),
  paste("min_spearman", min(rwr_rank_stability$spearman_to_primary)),
  paste("host_strict_p", conv_summary$empirical_p[conv_summary$layer == "host" & conv_summary$tier == "strict"]),
  paste("pathogen_strict_p", conv_summary$empirical_p[conv_summary$layer == "pathogen" & conv_summary$tier == "strict"]),
  paste("joint_strict_p", conv_summary$empirical_p[conv_summary$layer == "joint_two_sided"])
), file.path(tmp, "patch_revision_primary_summary.txt"))
cat(format(Sys.time()), "wrote primary RWR tables\n")
print(conv_summary %>% select(layer, tier, observed, z, empirical_p, selected_n))
print(rwr_rank_stability)
