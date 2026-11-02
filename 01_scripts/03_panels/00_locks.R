# Shared Theileria biology locks for Figure 3-5 panel scripts.
# Visual contract stays with the gold S. aureus panels; labels and counts are Theileria-specific.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

ink <- "#30363B"
frame <- "#DCE3E7"
line <- "#E4E9EC"
blue_dark <- "#B5C8D1"
blue <- "#D2DEE4"
blue_pale <- "#F1F5F7"
red_dark <- "#CDB0AD"
red <- "#E4D0CD"
red_pale <- "#FAF1F0"
host_fill <- blue
host_dark <- blue_dark
host_pale <- blue_pale
pathogen_fill <- red
pathogen_dark <- red_dark
pathogen_pale <- red_pale
interface_fill <- "#D5C8D4"
interface_dark <- "#B7A6B8"
interface_mid <- "#C9B8C8"
interface_pale <- "#F6F2F6"

mechanism_cols <- c(
  "Microtubule docking" = "#C2AFBA",
  "Metabolic / proteostasis" = "#DAB3A6",
  "Trafficking / signaling" = "#A5B8C5",
  "Protease / cullin" = "#A2BDB7"
)
mechanism_pale <- c(
  "Microtubule docking" = "#E4D9E0",
  "Metabolic / proteostasis" = "#EBD6D0",
  "Trafficking / signaling" = "#C9D6DE",
  "Protease / cullin" = "#C7D8D4"
)

core_axis_ids <- c("p104--MAPRE1", "TaPin1--PKM", "TA20980--CD2AP")
extended_axis_ids <- c(
  "TaSP--tubulin", "TaPin1--FBXW7", "TaPHB1--RUVBL1", "TA11965--TFG",
  "TaCyp1--MED21", "TA05560--RBM39", "TaCP--PPP4C", "TaCP--CRBN"
)
context_axis_ids <- c("p104--CLASP1")
core_endpoint_labels <- c("p104", "TaPin1", "CD2AP", "MAPRE1", "PKM", "TA20980")

axis_family_label <- c(
  microtubule_docking = "Microtubule docking",
  metabolic_proteostasis = "Metabolic / proteostasis",
  trafficking_signaling = "Trafficking / signaling",
  protease_cullin = "Protease / cullin"
)

ensure_cols <- function(d, cols, fill = NA) {
  for (col in cols) {
    if (!col %in% names(d)) d[[col]] <- fill
  }
  d
}

flag_col <- function(x) {
  if (is.null(x)) return(logical())
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  tolower(as.character(coalesce(x, ""))) %in% c("true", "t", "1", "yes", "y")
}

num_col <- function(x, default = 0) {
  z <- suppressWarnings(as.numeric(x))
  z[!is.finite(z)] <- default
  z
}

sanitize_label <- function(x) {
  x <- as.character(coalesce(x, ""))
  x <- str_replace_all(x, "TA_15795\\?", "TaPin1")
  x <- str_replace_all(x, "TA_15795", "TaPin1")
  x <- str_replace_all(x, "\uFF08", "(")
  x <- str_replace_all(x, "\uFF09", ")")
  x <- str_replace_all(x, "\u2014", "-")
  x <- str_replace_all(x, "\u2013", "-")
  x <- str_replace_all(x, "p104=TA08425[^;]*", "p104")
  x <- str_replace_all(x, "TA08425", "p104")
  x <- str_replace_all(x, "\\(正文确认\\)", "")
  x <- str_replace_all(x, "\\(文中作 TaPin1\\)", "")
  x <- str_replace_all(x, "\\(无 TA 号记录于台账\\)", "")
  x <- str_replace_all(x, "\\(2015 论文[^)]*\\)", "")
  x <- str_replace_all(x, "\\(基因号未给出\\)", "")
  x <- str_replace_all(x, "\\(aa66-413[^)]*\\)", "")
  x <- str_replace_all(x, "\\(XP_952499\\.1\\)", "")
  x <- str_replace_all(x, "\\(XM_950452\\.1\\)", "")
  x <- str_squish(x)
  x
}

pretty_module <- function(z, wrap = FALSE) {
  raw <- coalesce(as.character(z), "Unassigned")
  mapped <- recode(
    raw,
    H_NET_1 = "Signaling-broad",
    host_signaling_broad = "Signaling-broad",
    H_NET_4 = "RNA processing",
    host_rna_processing = "RNA processing",
    H_NET_5 = "Metabolism",
    host_metabolism = "Metabolism",
    H_NET_7 = "Chromatin",
    host_chromatin_organization = "Chromatin",
    H_NET_8 = "Ubiquitin / proteostasis",
    host_ubiquitin_proteostasis = "Ubiquitin / proteostasis",
    H_NET_12 = "Transcription",
    host_transcription_machinery = "Transcription",
    H_NET_14 = "Microtubule / cytoskeleton",
    host_microtubule_cytoskeleton = "Microtubule / cytoskeleton",
    H_FUN_cytoskeleton_microtubule = "Microtubule / cytoskeleton",
    H_FUN_glycolysis_metabolism = "Glycolysis / metabolism",
    H_FUN_mediator_transcription = "Mediator / transcription",
    H_FUN_ubiquitin_cullin = "Ubiquitin / cullin",
    H_FUN_chaperone_folding = "Chaperone / folding",
    H_FUN_endocytosis_phagosome = "Endocytosis",
    H_FUN_proteasome = "Proteasome",
    H_FUN_ribosome_translation = "Ribosome / translation",
    H_FUN_spliceosome = "Spliceosome",
    H_FUN_vesicle_COPII = "COPII vesicle",
    P_NET_1 = "Cellular process-broad",
    parasite_cellular_process_broad = "Cellular process-broad",
    P_NET_3 = "Cytoplasmic metabolism",
    parasite_cytoplasmic_metabolism = "Cytoplasmic metabolism",
    P_NET_4 = "RNA processing",
    parasite_rna_processing = "RNA processing",
    P_FUN_surface_family = "Surface family",
    parasite_surface_family = "Surface family",
    P_FUN_PIase_PPI = "PIase / PPI",
    parasite_PIase_PPI = "PIase / PPI",
    P_FUN_translation_folding = "Translation / folding",
    parasite_translation_folding = "Translation / folding",
    P_FUN_DNA_RNA_binding = "DNA / RNA binding",
    P_FUN_kinase_phosphatase = "Kinase / phosphatase",
    P_FUN_membrane_transport = "Membrane transport",
    P_FUN_metabolism_ATP = "Metabolism / ATP",
    surface_export_trafficking = "Surface / export",
    protein_folding_processing = "Folding / processing",
    Unassigned = "Unassigned",
    `__unassigned__` = "Unassigned",
    .default = NA_character_
  )
  out <- if_else(
    is.na(mapped),
    str_replace_all(str_remove(raw, "^[HP]_"), "_", " "),
    mapped
  )
  if (isTRUE(wrap)) {
    out <- str_replace(out, " / ", " /\n")
    out <- str_replace(out, "-broad", "-\nbroad")
    out <- str_replace(out, "proteostasis", "proteo-\nstasis")
  }
  out
}

mechanism_from_family <- function(x) {
  z <- str_to_lower(str_replace_all(coalesce(as.character(x), ""), "[^a-z0-9]+", "_"))
  z <- str_replace(z, "^_+|_+$", "")
  unname(axis_family_label[z])
}

axis_key <- function(axis_id, interaction_name, pathogen_gene, host_gene) {
  ids <- sanitize_label(coalesce(axis_id, ""))
  names <- sanitize_label(coalesce(interaction_name, ""))
  pairs <- paste(
    sanitize_label(coalesce(pathogen_gene, "")),
    sanitize_label(coalesce(host_gene, "")),
    sep = "--"
  )
  dplyr::case_when(
    ids %in% core_axis_ids | names %in% core_axis_ids | pairs %in% core_axis_ids ~ ids,
    str_detect(ids, "p104") & str_detect(paste(names, pairs), "MAPRE1") ~ "p104--MAPRE1",
    str_detect(paste(ids, names, pairs), "TaPin1") & str_detect(paste(ids, names, pairs), "PKM") ~ "TaPin1--PKM",
    str_detect(paste(ids, names, pairs), "TA20980") & str_detect(paste(ids, names, pairs), "CD2AP") ~ "TA20980--CD2AP",
    str_detect(paste(ids, names, pairs), "TaSP") & str_detect(paste(ids, names, pairs), "tubulin|TUBA|TUBG") ~ "TaSP--tubulin",
    str_detect(paste(ids, names, pairs), "TaPin1") & str_detect(paste(ids, names, pairs), "FBXW7") ~ "TaPin1--FBXW7",
    str_detect(paste(ids, names, pairs), "TaPHB1") ~ "TaPHB1--RUVBL1",
    str_detect(paste(ids, names, pairs), "TA11965") ~ "TA11965--TFG",
    str_detect(paste(ids, names, pairs), "TaCyp1") ~ "TaCyp1--MED21",
    str_detect(paste(ids, names, pairs), "TA05560") ~ "TA05560--RBM39",
    str_detect(paste(ids, names, pairs), "TaCP") & str_detect(paste(ids, names, pairs), "PPP4C") ~ "TaCP--PPP4C",
    str_detect(paste(ids, names, pairs), "TaCP") & str_detect(paste(ids, names, pairs), "CRBN") ~ "TaCP--CRBN",
    str_detect(paste(ids, names, pairs), "CLASP1") ~ "p104--CLASP1",
    TRUE ~ if_else(nzchar(ids), ids, if_else(nzchar(pairs) & pairs != "--", pairs, names))
  )
}

axis_tier <- function(key, final_tier = NULL) {
  key <- as.character(key)
  dplyr::case_when(
    key %in% core_axis_ids ~ "Core",
    key %in% extended_axis_ids ~ "Extended",
    key %in% context_axis_ids ~ "Context",
    str_to_lower(coalesce(as.character(final_tier), "")) %in% c("core", "strict") ~ "Core",
    str_to_lower(coalesce(as.character(final_tier), "")) %in% c("extended", "relaxed") ~ "Extended",
    TRUE ~ NA_character_
  )
}

endpoint_display <- function(x) {
  z <- sanitize_label(x)
  z <- recode(
    z,
    TA08425 = "p104",
    `p104=TA08425` = "p104",
    TUBA1A = "tubulin",
    TUBG1 = "tubulin",
    `TUBA1A|TUBG1` = "tubulin",
    .default = z
  )
  z
}
