#!/usr/bin/env Rscript
# =============================================================================
# Script 10: Priorización multi-criterio v3 → candidatos de repurposing
# HNSCC Drug Repurposing — Fase 3
# =============================================================================
# Objetivo: Calcular un score compuesto de DOS FACTORES para cada fármaco
#           candidato, anclado en la estructura de la red PPI (módulos + hubs).
#
#   composite = w_target · TargetPriority + w_drug · DrugViability
#
#   TargetPriority (nivel-DIANA; compartido por fármacos de la misma diana)
#     - centrality     : centralidad topológica del hub que engancha el fármaco
#                        media(degree_norm, betweenness_norm_scaled, eigenvector)
#     - pi_directional : pi = sign(logFC)·|logFC|·-log10(FDR), min-max DIRECCIONAL
#                        (premia inhibir dianas UP; dianas down puntúan bajo)
#   DrugViability (nivel-DROGA; rompe empates entre fármacos de la misma diana)
#     - reversal       : reversión de firma L2S2/LINCS (cmap_score; NA → 0)
#     - regulatory     : clase regulatoria A/B/C/D (graduado, no saturado)
#     - breadth        : nº bases de datos de soporte / 4
#
# Cada fármaco se atribuye a su diana MÁS CENTRAL (primary_target = argmax
# centrality sobre sus dianas DE), de la que hereda módulo, tier y TargetPriority.
#   tier ∈ {hub_central (diana ∈ druggable hubs), peripheral_diff (diana en red,
#           no hub), off_network (diana fuera de red)}
#
# Pesos en config/analysis_params.yaml → bloque scoring_v3 + regulatory_class_scores.
# Cambios v3: se eliminan score_pathway (redundante con red/módulo, saturado),
#             score_evidence (muerto), score_clinical crudo (→ clase regulatoria),
#             y el bonus. Resuelve los empates EGFR (drogas difieren en DrugViability).
#
# Input:  results/tables/drug_targets/08_multi_source_candidates.tsv
#         results/tables/de_limma/02_TVsS_significant_with_ids.tsv
#         results/tables/network/09_network_node_metrics.tsv  (degree, betweenness,
#                 eigenvector, module_id, module_name)
#         results/tables/network/09_druggable_hubs.tsv        (set de hubs drogables)
#         config/analysis_params.yaml
#
# Output: results/tables/10_all_candidates_scored.tsv
#         results/tables/10_top20_candidates.tsv
#         results/tables/10_module_hub_candidates.tsv  (tabla-resumen de revisión)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(yaml)
})

# --- Working directory (raíz del proyecto vía scripts/_setup.R) --------------
source(here::here("scripts", "_setup.R"))
setup_project()

# --- Log setup ---------------------------------------------------------------
dir.create("logs",   showWarnings = FALSE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)

log_file <- paste0("logs/10_prioritization_scoring_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
con_log <- file(log_file, open = "wt")
sink(con_log, type = "output")
sink(con_log, type = "message", append = TRUE)

cat("=== 10_prioritization_scoring.R (v3: TargetPriority × DrugViability) ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# --- Parametros --------------------------------------------------------------
params <- yaml::read_yaml("config/analysis_params.yaml")
sc <- params$scoring_v3
if (is.null(sc)) stop("FATAL: falta bloque scoring_v3 en config/analysis_params.yaml")

w_target <- sc$weight_target
w_drug   <- sc$weight_drug
w_cent   <- sc$weight_centrality
w_pi     <- sc$weight_pi
w_rev    <- sc$weight_reversal
w_reg    <- sc$weight_regulatory
w_brd    <- sc$weight_breadth
top_n    <- params$candidates$top_n

reg_scores <- setNames(
  as.numeric(unlist(params$regulatory_class_scores)),
  names(params$regulatory_class_scores)
)

WEIGHT_TOLERANCE <- 0.001
check_sum <- function(x, label) {
  if (abs(sum(x) - 1.0) > WEIGHT_TOLERANCE)
    stop(sprintf("FATAL: pesos %s suman %.4f, deben sumar 1.0", label, sum(x)))
}
check_sum(c(w_target, w_drug),        "composite (target+drug)")
check_sum(c(w_cent, w_pi),            "TargetPriority (centrality+pi)")
check_sum(c(w_rev, w_reg, w_brd),     "DrugViability (reversal+regulatory+breadth)")

cat(sprintf("composite = %.2f·TargetPriority + %.2f·DrugViability\n", w_target, w_drug))
cat(sprintf("  TargetPriority = %.2f·centrality + %.2f·pi_directional\n", w_cent, w_pi))
cat(sprintf("  DrugViability  = %.2f·reversal + %.2f·regulatory + %.2f·breadth\n",
            w_rev, w_reg, w_brd))
cat(sprintf("Top N: %d\n\n", top_n))

# =============================================================================
# CARGAR DATOS
# =============================================================================
cat("--- Cargando datos ---\n")

# Candidatos multi-fuente (script 08)
candidates_raw <- read.delim("results/tables/drug_targets/08_multi_source_candidates.tsv",
                              stringsAsFactors = FALSE)
cat(sprintf("Candidatos multi-fuente (raw): %d\n", nrow(candidates_raw)))

# Deduplicar por ChEMBL ID canonico: si varios nombres apuntan al mismo ChEMBL,
# quedarse con el nombre mas corto (nombre base, sin sal/formulacion).
candidates <- candidates_raw %>%
  filter(!is.na(chembl_id) & chembl_id != "") %>%
  group_by(chembl_id) %>%
  slice_min(nchar(drug_name_norm), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  bind_rows(
    candidates_raw %>%
      filter(is.na(chembl_id) | chembl_id == "") %>%
      distinct(drug_name_norm, .keep_all = TRUE)
  ) %>%
  distinct(drug_name_norm, .keep_all = TRUE)

cat(sprintf("Candidatos tras dedup por ChEMBL ID: %d\n", nrow(candidates)))

# Deduplicar por nombre base (eliminar sufijos de sal/formulacion).
SALT_SUFFIXES <- paste0(
  "(\\s+(hydrochloride|hydrobromide|hydroiodide|",
  "hcl|hbr|hcloride|",
  "calcium|sodium|potassium|magnesium|zinc|",
  "hyclate|monohydrate|dihydrate|trihydrate|anhydrous|",
  "phosphate|sulfate|sulphate|acetate|tartrate|maleate|",
  "fumarate|mesylate|besylate|citrate|bromide|chloride|",
  "lactate|succinate|gluconate|malate|nitrate|oxalate|",
  "tosylate|pamoate|xinafoate|edisylate))+$"
)

candidates <- candidates %>%
  mutate(drug_name_base = str_trim(str_remove(drug_name_norm,
                                               regex(SALT_SUFFIXES, ignore_case = TRUE)))) %>%
  group_by(drug_name_base) %>%
  mutate(n_sources = max(n_sources)) %>%
  slice_min(nchar(drug_name_norm), n = 1, with_ties = FALSE) %>%
  mutate(drug_name_norm = drug_name_base) %>%
  ungroup()

cat(sprintf("Candidatos tras dedup por nombre base (sales): %d\n", nrow(candidates)))

# Tabla DE (script 02) — logFC y p-valor por gen → pi_stat
sig <- read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
                  stringsAsFactors = FALSE)
gene_de <- sig %>%
  select(symbol_org, logFC_TVsS, adj.P.Val_TVsS) %>%
  distinct(symbol_org, .keep_all = TRUE) %>%
  mutate(pi_stat = sign(logFC_TVsS) * abs(logFC_TVsS) *
                   (-log10(adj.P.Val_TVsS + 1e-300)))

# Normalización direccional (min-max sobre pi con signo):
#   1.0 = gen más UP-regulado; 0.0 = gen más DOWN-regulado del dataset.
pi_min   <- min(gene_de$pi_stat, na.rm = TRUE)
pi_max   <- max(gene_de$pi_stat, na.rm = TRUE)
pi_range <- pi_max - pi_min
if (pi_range == 0) stop("FATAL: pi_range = 0 — revisar tabla DE.")
gene_de <- gene_de %>% mutate(pi_dir = (pi_stat - pi_min) / pi_range)

# Métricas de red (script 09) — centralidad + módulo por nodo
net <- read.delim("results/tables/network/09_network_node_metrics.tsv",
                  stringsAsFactors = FALSE)

# Hubs de red (is_hub = top 10% por grado, script 09) — definen el tier hub_central.
is_hub_set <- net$gene_symbol[net$is_hub == "TRUE" | net$is_hub == TRUE]
cat(sprintf("Hubs de red (is_hub): %d\n", length(is_hub_set)))

# Aristas CREÍBLES fármaco→diana (master table, script 08).
# CRÍTICO: DGIdb anota dianas de forma promiscua y sus interaction_score no
# identifican el mecanismo de acción (ej. decitabine→NR5A1 > DNMT1; nivolumab→EGFR
# espurio). En cambio, ChEMBL y OpenTargets son fuentes CURADAS de mecanismo/
# asociación y traen gene_symbol. Por eso una arista fármaco→diana se considera
# CREÍBLE solo si proviene de ChEMBL u OpenTargets. Esto recupera correctamente:
# EGFR→{gefitinib,cetuximab,afatinib…}, PSMA2/PSMB3→carfilzomib, NDUFx→metformina,
# DNMT1→{decitabine,azacitidine}, MAOA→tranylcypromine; y excluye el ruido DGIdb.
# (Limitación: fármacos sin arista curada — ej. bortezomib — quedan sin ancla de red;
#  se revisan manualmente en la curación del shortlist.)
master <- read.delim("results/tables/drug_targets/08_drug_target_master_table.tsv",
                     stringsAsFactors = FALSE)
strip_salts <- function(x) str_trim(str_remove(x, regex(SALT_SUFFIXES, ignore_case = TRUE)))
credible_edges <- master %>%
  filter(source_db %in% c("ChEMBL", "OpenTargets"),
         !is.na(gene_symbol) & gene_symbol != "") %>%
  mutate(drug_base = strip_salts(drug_name_norm)) %>%
  distinct(drug_base, gene_symbol)
cat(sprintf("Aristas creíbles (curadas ChEMBL/OpenTargets): %d en %d fármacos\n",
            nrow(credible_edges), n_distinct(credible_edges$drug_base)))
credible_by_drug <- split(credible_edges$gene_symbol, credible_edges$drug_base)

# Centralidad topológica por nodo: media de degree/betweenness/eigenvector
# (cada uno min-max a 0–1 sobre los nodos de la red).
norm01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}
net <- net %>%
  mutate(
    deg_n  = norm01(degree),
    betw_n = norm01(betweenness_norm),
    eig_n  = norm01(eigenvector),
    centrality = (deg_n + betw_n + eig_n) / 3
  )

# Lookups por gen
cent_by_gene   <- setNames(net$centrality, net$gene_symbol)
mod_id_by_gene <- setNames(net$module_id,   net$gene_symbol)
mod_nm_by_gene <- setNames(net$module_name, net$gene_symbol)
pi_by_gene     <- setNames(gene_de$pi_dir,  gene_de$symbol_org)

# TargetPriority por gen (nivel-diana): w_cent·centralidad + w_pi·pi_directional
tp_by_gene <- setNames(
  w_cent * net$centrality + w_pi * ifelse(net$gene_symbol %in% names(pi_by_gene),
                                          pi_by_gene[net$gene_symbol], 0),
  net$gene_symbol
)

# Genes-ancla excluidos (artefactos revisados en checkpoint: EML4 + sangre/ECM)
excl_anchor <- toupper(unlist(params$exclusions$exclude_anchor_genes))
cat(sprintf("Genes-ancla excluidos: %s\n", paste(excl_anchor, collapse = ", ")))

# Ancla un fármaco a su diana CREÍBLE (curada) más central → hereda módulo/tier/TP.
# Desempate: pi_directional, luego alfabético. Sin diana creíble de red → off_network.
anchor_fn <- function(drug_base) {
  genes <- credible_by_drug[[drug_base]]
  genes <- genes[genes %in% net$gene_symbol]            # solo dianas en la red PPI
  genes <- genes[!toupper(genes) %in% excl_anchor]      # quitar genes-ancla excluidos
  if (length(genes) == 0)
    return(list(primary = NA, cent = 0, pi = 0, tp = 0,
                module_id = NA, module_name = NA, tier = "off_network"))
  cent <- unname(cent_by_gene[genes]); pid <- unname(pi_by_gene[genes])
  pid[is.na(pid)] <- 0
  prim <- genes[order(-cent, -pid, genes)[1]]
  list(
    primary     = prim,
    cent        = unname(cent_by_gene[prim]),
    pi          = ifelse(prim %in% names(pi_by_gene), unname(pi_by_gene[prim]), 0),
    tp          = unname(tp_by_gene[prim]),
    module_id   = unname(mod_id_by_gene[prim]),
    module_name = unname(mod_nm_by_gene[prim]),
    tier        = if (prim %in% is_hub_set) "hub_central" else "peripheral_diff"
  )
}

# =============================================================================
# HELPER: dianas declaradas de un fármaco (para descriptivos)
# =============================================================================
get_target_genes <- function(gene_str) {
  if (is.na(gene_str) || gene_str == "") return(character(0))
  unique(str_split(gene_str, "\\|")[[1]])
}

# DrugViability: reversal L2S2
min_cmap <- min(candidates$cmap_score, na.rm = TRUE)  # el mas negativo
score_reversal_fn <- function(cs) {
  if (is.na(cs) || cs >= 0) return(0)
  cs / min_cmap   # ambos negativos → 0–1 (min_cmap → 1)
}
# DrugViability: clase regulatoria
score_regulatory_fn <- function(cls) {
  v <- reg_scores[as.character(cls)]
  if (length(v) == 0 || is.na(v)) 0 else as.numeric(v)
}

# =============================================================================
# APLICAR SCORING
# =============================================================================
cat("\n--- Calculando scores (TargetPriority × DrugViability) ---\n")

an <- lapply(candidates$drug_name_norm, anchor_fn)
scored <- candidates %>%
  mutate(
    primary_target            = vapply(an, function(x) as.character(x$primary), character(1)),
    primary_target_centrality = vapply(an, `[[`, numeric(1), "cent"),
    primary_pi                = vapply(an, `[[`, numeric(1), "pi"),
    TargetPriority            = vapply(an, `[[`, numeric(1), "tp"),
    module_id                 = vapply(an, function(x) as.character(x$module_id), character(1)),
    module_name               = vapply(an, function(x) as.character(x$module_name), character(1)),
    tier                      = vapply(an, `[[`, character(1), "tier"),
    is_network_hub            = primary_target %in% is_hub_set,
    has_credible_target       = tier != "off_network"
  ) %>%
  rowwise() %>%
  mutate(
    s_reversal   = score_reversal_fn(cmap_score),
    s_regulatory = score_regulatory_fn(drug_class),
    s_breadth    = min(n_sources / 4, 1.0)
  ) %>%
  ungroup() %>%
  mutate(
    DrugViability   = w_rev * s_reversal + w_reg * s_regulatory + w_brd * s_breadth,
    composite_score = w_target * TargetPriority + w_drug * DrugViability,
    final_score     = composite_score   # alias backward-compat (scripts 15/18)
  ) %>%
  arrange(desc(composite_score))

cat(sprintf("Composite range: %.4f - %.4f\n",
            min(scored$composite_score), max(scored$composite_score)))

# =============================================================================
# EXCLUSIONES (config -> exclusions)
# =============================================================================
excluded_drugs <- toupper(unlist(params$exclusions$drugs))
n_excl_name <- 0
if (length(excluded_drugs) > 0) {
  n_excl_name <- sum(toupper(scored$drug_name_norm) %in% excluded_drugs)
  scored <- scored %>% filter(!toupper(drug_name_norm) %in% excluded_drugs)
  cat(sprintf("Exclusiones por nombre: %d fármacos removidos\n", n_excl_name))
}
excl_targets <- toupper(unlist(params$exclusions$exclude_single_target_genes))
n_excl_target <- 0
if (length(excl_targets) > 0) {
  is_single_excl <- toupper(scored$de_genes) %in% excl_targets
  n_excl_target  <- sum(is_single_excl)
  scored <- scored %>% filter(!is_single_excl)
  cat(sprintf("Exclusiones por target unico: %d fármacos removidos\n", n_excl_target))
}
cat(sprintf("Total excluidos: %d | Candidatos restantes: %d\n",
            n_excl_name + n_excl_target, nrow(scored)))

# =============================================================================
# TOP N
# =============================================================================
top20 <- scored %>% arrange(desc(composite_score)) %>% slice_head(n = top_n)

cat(sprintf("\nDistribución de tiers en top %d:\n", top_n))
print(top20 %>% count(tier))
cat(sprintf("Módulos representados en top %d: %d\n",
            top_n, n_distinct(top20$module_name)))

cat(sprintf("\n=== TOP %d CANDIDATOS ===\n", top_n))
cat(sprintf("  %-28s %6s %6s %6s  %-14s  %-16s  %s\n",
            "Drug", "Comp", "TP", "DV", "Tier", "Module", "Primary"))
cat(paste(rep("-", 100), collapse=""), "\n")
for (i in seq_len(nrow(top20))) {
  r <- top20[i, ]
  cat(sprintf("  %2d. %-24s %.3f  %.3f  %.3f  %-14s  %-16s  %s\n",
              i, substr(r$drug_name_norm, 1, 24),
              r$composite_score, r$TargetPriority, r$DrugViability,
              r$tier, substr(ifelse(is.na(r$module_name), "—", r$module_name), 1, 16),
              r$primary_target))
}

# =============================================================================
# TABLA-RESUMEN DE REVISIÓN: módulo → hub (primary_target) → tier → fármacos
# =============================================================================
cat("\n--- Tabla-resumen módulo→hub→droga (para checkpoint) ---\n")
review <- scored %>%
  group_by(module_id, module_name, primary_target, tier) %>%
  summarise(
    n_drugs        = n(),
    max_composite  = max(composite_score),
    centrality     = first(primary_target_centrality),
    primary_pi     = first(primary_pi),
    top_drugs      = paste(head(str_to_title(str_to_lower(drug_name_norm[order(-composite_score)])), 4),
                           collapse = " | "),
    best_class     = drug_class[which.max(composite_score)],
    .groups = "drop"
  ) %>%
  arrange(factor(tier, levels = c("hub_central","peripheral_diff","off_network")),
          desc(max_composite))

write.table(review, "results/tables/10_module_hub_candidates.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Exportado: 10_module_hub_candidates.tsv (%d filas hub×módulo)\n", nrow(review)))

# =============================================================================
# EXPORTAR tablas de scoring
# =============================================================================
cat("\n--- Exportando ---\n")
out_cols <- c("drug_name_norm", "chembl_id", "drug_class", "drug_class_label",
              "n_sources", "sources", "max_phase", "is_approved", "hnscc_indication",
              "cmap_score", "n_de_genes", "n_up_genes", "n_down_genes", "de_genes",
              "primary_target", "primary_target_centrality", "primary_pi",
              "module_id", "module_name", "tier", "is_network_hub", "has_credible_target",
              "s_reversal", "s_regulatory", "s_breadth",
              "TargetPriority", "DrugViability", "composite_score", "final_score")
out_cols_avail <- intersect(out_cols, colnames(scored))

write.table(scored %>% select(all_of(out_cols_avail)),
            "results/tables/10_all_candidates_scored.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top20 %>% select(all_of(out_cols_avail)),
            "results/tables/10_top20_candidates.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Exportado: 10_all_candidates_scored.tsv (%d)\n", nrow(scored)))
cat(sprintf("Exportado: 10_top20_candidates.tsv (top %d)\n", top_n))

# =============================================================================
# RESUMEN
# =============================================================================
cat("\n=== RESUMEN ===\n")
cat(sprintf("  Candidatos evaluados: %d\n", nrow(scored)))
cat(sprintf("  Tiers (todos):\n")); print(scored %>% count(tier))
cat(sprintf("  Composite top: %.4f (%s)\n",
            max(scored$composite_score), scored$drug_name_norm[which.max(scored$composite_score)]))
cat("\nFin:", format(Sys.time()), "\n")

sink(type = "message"); sink(); close(con_log)
cat("Script completado. Log en:", log_file, "\n")
