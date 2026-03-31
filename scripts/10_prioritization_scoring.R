#!/usr/bin/env Rscript
# =============================================================================
# Script 10: Priorización multi-criterio → Top 20 candidatos
# HNSCC Drug Repurposing — Fase 3
# =============================================================================
# Objetivo: Calcular un score compuesto para cada farmaco candidato usando
#           6 dimensiones independientes. Seleccionar top 20 para validacion.
#
# Dimensiones de scoring (pesos en config/analysis_params.yaml):
#   1. score_pi_stat       (0.325): pi-stat = sign(logFC)*|logFC|*(-log10 adj.P)
#   2. score_clinical      (0.200): fase clinica maxima del farmaco
#   3. score_cmap          (0.130): connectivity reversal score L2S2/LINCS (reversal)
#   4. score_pathway       (0.150): proporcion genes diana en vias enriquecidas
#   5. score_network       (0.195): centralidad + diversidad de modulos en red PPI
# NOTA v2: score_evidence eliminado del composite (sesgo publicacion + circularidad
#          con revision manual del grupo). Se mantiene como columna descriptiva.
#
# Input:  results/tables/drug_targets/08_drug_summary_per_drug.tsv
#         results/tables/drug_targets/08_multi_source_candidates.tsv
#         results/tables/de_limma/02_TVsS_significant_with_ids.tsv
#         results/tables/network/09_network_node_metrics.tsv
#         results/tables/pathway_enrichment/03_GO_BP_ORA_simplified.tsv
#         results/tables/pathway_enrichment/03_KEGG_ORA.tsv
#         results/tables/pathway_enrichment/03_Reactome_ORA.tsv
#         config/analysis_params.yaml
#
# Output: results/tables/10_top20_candidates.tsv
#         results/tables/10_all_candidates_scored.tsv
#         results/figures/10_scoring_heatmap.pdf
#         results/figures/10_top20_barplot.pdf
#         results/figures/10_score_vs_targets.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
})

# --- Working directory -------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  proj_dir    <- dirname(dirname(script_path))
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml")))
    setwd(proj_dir)
}
cat("Working dir:", getwd(), "\n")

# --- Log setup ---------------------------------------------------------------
dir.create("logs",   showWarnings = FALSE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

log_file <- paste0("logs/10_prioritization_scoring_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
con_log <- file(log_file, open = "wt")
sink(con_log, type = "output")
sink(con_log, type = "message", append = TRUE)

cat("=== 10_prioritization_scoring.R ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# --- Parametros --------------------------------------------------------------
params <- yaml::read_yaml("config/analysis_params.yaml")

w_pistat  <- params$scoring$weight_pi_stat           # reemplaza logFC + sig
w_clin    <- params$scoring$weight_clinical_phase
w_cmap    <- params$scoring$weight_cmap_connectivity  # reducido a 0.10
w_path    <- params$scoring$weight_pathway_relevance
w_net     <- params$scoring$weight_network_centrality
w_evid    <- params$scoring$weight_evidence %||% 0.0  # excluido del composite (v2: sesgo publicacion)
top_n     <- params$candidates$top_n

phase_scores <- setNames(
  as.numeric(unlist(params$clinical_phase_scores)),
  names(params$clinical_phase_scores)
)

WEIGHT_TOLERANCE <- 0.001  # acceptable deviation from sum-to-1 for scoring weights

cat(sprintf("Pesos: pi_stat=%.3f | clinical=%.3f | cmap=%.3f | pathway=%.3f | network=%.3f (evidence excluido)\n",
            w_pistat, w_clin, w_cmap, w_path, w_net))
weight_sum <- w_pistat + w_clin + w_cmap + w_path + w_net
cat(sprintf("Sum pesos (sin evidence): %.4f\n", weight_sum))
if (abs(weight_sum - 1.0) > WEIGHT_TOLERANCE)
  stop(sprintf("FATAL: pesos suman %.4f, deben sumar exactamente 1.0", weight_sum))
cat(sprintf("Top N candidatos: %d\n\n", top_n))

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
# Se filtra primero para no agrupar entradas sin ChEMBL junto con las que sí lo tienen.
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
# Algunos compuestos tienen ChEMBL IDs distintos para la sal y la base (ej.
# METFORMIN CHEMBL1431 vs METFORMIN HYDROCHLORIDE CHEMBL1703). Se colapsan
# por nombre base, conservando la entrada con mayor n_sources.
# Drug name normalization: common salt/form suffixes to strip for deduplication
# Source: manually curated list of common pharmaceutical salt suffixes
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
  # Usar el nombre base para display; conservar el maximo n_sources del grupo
  mutate(n_sources = max(n_sources)) %>%
  slice_min(nchar(drug_name_norm), n = 1, with_ties = FALSE) %>%
  mutate(drug_name_norm = drug_name_base) %>%
  ungroup()

cat(sprintf("Candidatos tras dedup por nombre base (sales): %d\n", nrow(candidates)))

# Tabla DE (script 02) — logFC y p-valor por gen
sig <- read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
                  stringsAsFactors = FALSE)
gene_de <- sig %>%
  select(symbol_org, logFC_TVsS, adj.P.Val_TVsS) %>%
  distinct(symbol_org, .keep_all = TRUE)

# Metricas de red (script 09) — degree y betweenness
net_metrics <- read.delim("results/tables/network/09_network_node_metrics.tsv",
                          stringsAsFactors = FALSE) %>%
  select(gene_symbol, degree, betweenness_norm, eigenvector)

# Genes en vias enriquecidas (script 03)
# Extraer todos los gene symbols de los resultados ORA
load_pathway_genes <- function(file) {
  # clusterProfiler con readable=TRUE escribe gene SYMBOLS separados por "/"
  if (!file.exists(file)) return(character(0))
  df <- read.delim(file, stringsAsFactors = FALSE)
  if (!"geneID" %in% colnames(df)) return(character(0))
  genes <- unlist(strsplit(df$geneID, "/"))
  unique(trimws(genes))
}

pathway_syms_go    <- load_pathway_genes("results/tables/pathway_enrichment/03_GO_BP_ORA_simplified.tsv")
pathway_syms_kegg  <- load_pathway_genes("results/tables/pathway_enrichment/03_KEGG_ORA.tsv")
pathway_syms_react <- load_pathway_genes("results/tables/pathway_enrichment/03_Reactome_ORA.tsv")
pathway_syms_all   <- unique(c(pathway_syms_go, pathway_syms_kegg, pathway_syms_react))

cat(sprintf("Genes en vias enriquecidas (ORA): %d\n", length(pathway_syms_all)))

# =============================================================================
# HELPER: genes diana de un farmaco
# =============================================================================
get_target_genes <- function(gene_str) {
  if (is.na(gene_str) || gene_str == "") return(character(0))
  unique(str_split(gene_str, "\\|")[[1]])
}

# =============================================================================
# DIMENSION 1 (nueva): score_pi_stat — combina magnitud y significancia
#             pi = sign(logFC) × |logFC| × -log10(adj.P.Val)
#             Normalizado 0-1 sobre el máximo del dataset.
# =============================================================================
cat("\n--- Calculando scores ---\n")

gene_de <- gene_de %>%
  mutate(
    pi_stat = sign(logFC_TVsS) * abs(logFC_TVsS) *
              (-log10(adj.P.Val_TVsS + 1e-300))
  )
max_pi <- max(abs(gene_de$pi_stat), na.rm = TRUE)

score_pistat_fn <- function(gene_str) {
  genes   <- get_target_genes(gene_str)
  pi_vals <- abs(gene_de$pi_stat[gene_de$symbol_org %in% genes])
  if (length(pi_vals) == 0) return(0)
  mean(pi_vals, na.rm = TRUE) / max_pi
}

# =============================================================================
# DIMENSION 3: score_clinical — fase clinica maxima → score 0-1
# =============================================================================
score_clinical_fn <- function(phase) {
  if (is.na(phase)) return(0)
  key <- as.character(round(phase))
  val <- phase_scores[key]
  if (is.na(val)) 0 else as.numeric(val)
}

# =============================================================================
# DIMENSION 4: score_cmap — connectivity score (mas negativo = mejor)
#              min score en dataset -> 1, 0 -> 0, positivo -> 0
# =============================================================================
min_cmap <- min(candidates$cmap_score, na.rm = TRUE)  # el mas negativo

score_cmap_fn <- function(cs) {
  if (is.na(cs) || cs >= 0) return(0)
  # Escala: cs va de min_cmap a 0, queremos que min_cmap = 1
  (cs - 0) / (min_cmap - 0)   # = cs/min_cmap (ambos negativos -> positivo)
}

# =============================================================================
# DIMENSION 5: score_pathway — proporcion genes diana en vias ORA enriquecidas
# =============================================================================
score_pathway_fn <- function(gene_str) {
  genes <- get_target_genes(gene_str)
  if (length(genes) == 0) return(0)
  sum(genes %in% pathway_syms_all) / length(genes)
}

# =============================================================================
# DIMENSION 5 (red): score_network — centralidad + diversidad de módulos
#              usa combinacion degree + betweenness normalizada + module diversity
# =============================================================================
max_degree <- max(net_metrics$degree, na.rm = TRUE)

# Cargar tabla de módulos (generada por script 09)
modules_tbl <- tryCatch(
  read.delim("results/tables/network/09_modules.tsv", stringsAsFactors = FALSE) %>%
    select(gene_symbol, module_id) %>%
    distinct(),
  error = function(e) {
    cat("  WARN: 09_modules.tsv no encontrado — se usa centralidad sin módulos\n")
    NULL
  }
)

score_network_fn <- function(gene_str) {
  genes <- get_target_genes(gene_str)
  nm    <- net_metrics %>% filter(gene_symbol %in% genes)
  if (nrow(nm) == 0) return(0)

  deg_norm  <- mean(nm$degree / max_degree, na.rm = TRUE)
  betw_norm <- mean(nm$betweenness_norm,    na.rm = TRUE)
  centrality_score <- 0.6 * deg_norm + 0.4 * betw_norm

  if (!is.null(modules_tbl)) {
    n_modules_hit <- modules_tbl %>%
      filter(gene_symbol %in% genes) %>%
      pull(module_id) %>%
      n_distinct()
    module_diversity <- min(n_modules_hit / 3, 1.0)
    return(0.7 * centrality_score + 0.3 * module_diversity)
  }
  centrality_score
}

# =============================================================================
# DIMENSION 6 (nueva): score_evidence — evidencia clínica pre-calculada
# =============================================================================
evidence_tbl <- tryCatch(
  read.delim("results/tables/evidence/11_clinical_evidence.tsv",
             stringsAsFactors = FALSE) %>%
    mutate(drug_name_norm_ev = toupper(trimws(drug_name))) %>%
    select(drug_name_norm_ev, evidence_score),
  error = function(e) {
    cat("  INFO: 11_clinical_evidence.tsv no encontrado — score_evidence = 0\n")
    NULL
  }
)

score_evidence_fn <- function(drug_name) {
  if (is.null(evidence_tbl)) return(0)
  hit <- evidence_tbl$evidence_score[
    toupper(trimws(drug_name)) == evidence_tbl$drug_name_norm_ev
  ]
  if (length(hit) == 0 || all(is.na(hit))) return(0)
  max(hit, na.rm = TRUE)
}

# =============================================================================
# APLICAR SCORING A TODOS LOS CANDIDATOS
# =============================================================================
cat("Calculando scores para todos los candidatos...\n")

scored <- candidates %>%
  rowwise() %>%
  mutate(
    s_pi_stat  = score_pistat_fn(de_genes),
    s_clinical = score_clinical_fn(max_phase),
    s_cmap     = score_cmap_fn(cmap_score),
    s_pathway  = score_pathway_fn(de_genes),
    s_network  = score_network_fn(de_genes),
    s_evidence = score_evidence_fn(drug_name_norm)
  ) %>%
  ungroup() %>%
  mutate(
    composite_score = w_pistat * s_pi_stat   +
                      w_clin   * s_clinical  +
                      w_cmap   * s_cmap      +
                      w_path   * s_pathway   +
                      w_net    * s_network,
    # Bonus: candidatos con CMap + gene-target o Clase C multi-fuente
    # Clase A removida del bonus (se reporta como control positivo, no como candidato)
    bonus = case_when(
      !is.na(cmap_score) & n_sources >= 3 ~ 0.03,
      drug_class == "C" & n_sources >= 2  ~ 0.01,
      TRUE                                ~ 0
    ),
    final_score = pmin(composite_score + bonus, 1.0)
  ) %>%
  arrange(desc(final_score))

cat(sprintf("Score range: %.4f - %.4f\n",
            min(scored$final_score), max(scored$final_score)))

# =============================================================================
# APLICAR EXCLUSIONES (config -> exclusions$drugs)
# =============================================================================
# Exclusion por nombre de fármaco
excluded_drugs <- toupper(unlist(params$exclusions$drugs))
n_excl_name <- 0
if (length(excluded_drugs) > 0) {
  n_excl_name <- sum(toupper(scored$drug_name_norm) %in% excluded_drugs)
  scored <- scored %>% filter(!toupper(drug_name_norm) %in% excluded_drugs)
  cat(sprintf("Exclusiones por nombre: %d fármacos removidos\n", n_excl_name))
}

# Exclusion por gen diana unico (si el farmaco solo tiene 1 target y es excluido)
excl_targets <- toupper(unlist(params$exclusions$exclude_single_target_genes))
n_excl_target <- 0
if (length(excl_targets) > 0) {
  is_single_excl <- scored$de_genes %in% excl_targets
  n_excl_target  <- sum(is_single_excl)
  scored <- scored %>% filter(!is_single_excl)
  cat(sprintf("Exclusiones por target unico: %d fármacos removidos\n", n_excl_target))
}

cat(sprintf("Total excluidos: %d | Candidatos restantes: %d\n",
            n_excl_name + n_excl_target, nrow(scored)))

# =============================================================================
# TOP N — sin límite por target (criterio: LOD stability aplicado post-hoc)
# =============================================================================
# El límite de diversidad por target fue removido (v2): era arbitrario y
# ocultaba candidatos LOD-stable legítimos (ej. cluster EGFR).
# La agrupación por mecanismo se hace editorialmente en la discusión.
scored <- scored %>%
  mutate(primary_target = sapply(de_genes, function(g) {
    if (is.na(g) || g == "") return("unknown")
    str_split(g, "\\|")[[1]][1]
  }))

top20 <- scored %>%
  arrange(desc(final_score)) %>%
  slice_head(n = top_n)

cat(sprintf("Targets primarios representados en top %d: %d\n",
            top_n, n_distinct(top20$primary_target)))

cat(sprintf("\n=== TOP %d CANDIDATOS ===\n", top_n))
cat(sprintf("  %-35s %6s %5s %5s %4s  %s\n",
            "Drug", "Score", "Clin", "CMap", "Src", "Class"))
cat(paste(rep("-", 80), collapse=""), "\n")

for (i in seq_len(nrow(top20))) {
  r <- top20[i, ]
  cat(sprintf("  %2d. %-31s %6.4f %5s %5s %4d  %s\n",
              i,
              substr(r$drug_name_norm, 1, 31),
              r$final_score,
              ifelse(is.na(r$max_phase), "?", sprintf("%.0f", r$max_phase)),
              ifelse(is.na(r$cmap_score), "  -", sprintf("%.2f", r$cmap_score)),
              r$n_sources,
              r$drug_class))
}

# Score components breakdown para top 20
cat(sprintf("\n--- Desglose de scores (top %d) ---\n", top_n))
cat(sprintf("  %-30s  %6s  %5s  %5s  %5s  %5s  %5s  %5s\n",
            "Drug", "PI-stat", "Clin", "CMap", "Path", "Net", "Evid", "FINAL"))
for (i in seq_len(nrow(top20))) {
  r <- top20[i, ]
  cat(sprintf("  %-30s  %6.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f  %5.3f\n",
              substr(r$drug_name_norm, 1, 30),
              r$s_pi_stat, r$s_clinical,
              r$s_cmap, r$s_pathway, r$s_network, r$s_evidence,
              r$final_score))
}

# =============================================================================
# FIGURAS
# =============================================================================
cat("\n--- Generando figuras ---\n")

safe_pdf <- function(file, expr, w = 10, h = 7) {
  tryCatch({
    pdf(file, width = w, height = h)
    force(expr)
    dev.off()
    cat(sprintf("  %s\n", file))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("  WARN %s: %s\n", file, e$message))
  })
}

# 1. Barplot top 20 con color por clase
safe_pdf("results/figures/10_top20_barplot.pdf", w = 11, h = 8, {
  p_data <- top20 %>%
    mutate(rank = seq_len(n()),
           label = paste0(rank, ". ", str_to_title(str_to_lower(drug_name_norm))),
           label = factor(label, levels = rev(label)))

  p <- ggplot(p_data, aes(x = label, y = final_score, fill = drug_class_label)) +
    geom_col(width = 0.75, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.3f", final_score)),
              hjust = -0.1, size = 3.2) +
    scale_fill_manual(values = c(
      "A: Approved + HNSCC evidence"   = "#d62728",
      "B: Approved other cancer"       = "#ff7f0e",
      "C: Approved non-oncology"       = "#2ca02c",
      "D: Not approved / Experimental" = "#7f7f7f"
    )) +
    coord_flip() +
    expand_limits(y = 1.1) +
    labs(title = sprintf("Top %d drug repurposing candidates — HNSCC", top_n),
         subtitle = "Multi-criteria composite score (pi-stat + clinical phase + L2S2 + pathways + network modules + clinical evidence)",
         x = NULL, y = "Composite score (0-1)",
         fill = "Classification") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "bottom")
  print(p)
})

# 2. Heatmap de componentes de score (top 20)
safe_pdf("results/figures/10_scoring_heatmap.pdf", w = 12, h = 8, {
  score_cols <- c("s_pi_stat", "s_clinical", "s_cmap", "s_pathway", "s_network", "s_evidence")
  col_labels <- c("PI-stat\n(logFC×sig)", "Clinical\nPhase",
                  "L2S2\nReversal", "Pathway\nRelevance",
                  "Network\n(modules)", "Clinical\nEvidence")

  mat_data <- top20 %>%
    mutate(drug_label = sprintf("%d. %s", seq_len(n()),
                                str_sub(str_to_title(str_to_lower(drug_name_norm)), 1, 28))) %>%
    select(drug_label, all_of(score_cols)) %>%
    pivot_longer(-drug_label, names_to = "component", values_to = "score") %>%
    mutate(component = factor(component, levels = score_cols, labels = col_labels),
           drug_label = factor(drug_label, levels = rev(unique(drug_label))))

  p <- ggplot(mat_data, aes(x = component, y = drug_label, fill = score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", score)), size = 3) +
    scale_fill_gradient2(low = "#f7f7f7", mid = "#74c476", high = "#00441b",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Score\n(0-1)") +
    labs(title = sprintf("Score components — Top %d HNSCC drug candidates", top_n),
         x = "Scoring dimension", y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
          panel.grid = element_blank())
  print(p)
})

# 3. Scatter: score vs n_de_genes, coloreado por clase
safe_pdf("results/figures/10_score_vs_targets.pdf", w = 10, h = 7, {
  p_data <- scored %>%
    mutate(label = ifelse(seq_len(n()) <= top_n,
                          str_to_title(str_to_lower(drug_name_norm)), NA))
  p <- ggplot(p_data, aes(x = n_de_genes, y = final_score,
                           color = drug_class_label)) +
    geom_point(aes(size = n_sources), alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.8,
                             max.overlaps = 20, na.rm = TRUE) +
    geom_hline(yintercept = min(top20$final_score),
               linetype = "dashed", color = "gray40", linewidth = 0.7) +
    annotate("text", x = max(scored$n_de_genes, na.rm=TRUE),
             y = min(top20$final_score) + 0.01,
             label = sprintf("Top %d cutoff", top_n),
             hjust = 1, color = "gray40", size = 3) +
    scale_color_manual(values = c(
      "A: Approved + HNSCC evidence"   = "#d62728",
      "B: Approved other cancer"       = "#ff7f0e",
      "C: Approved non-oncology"       = "#2ca02c",
      "D: Not approved / Experimental" = "#7f7f7f"
    )) +
    scale_size_continuous(range = c(2, 6), name = "N databases") +
    labs(title = "Drug candidates: composite score vs. DE target genes",
         subtitle = "Labeled = top 20 | dashed line = top 20 cutoff",
         x = "Number of DE target genes", y = "Composite score",
         color = "Classification") +
    theme_bw(base_size = 12)
  print(p)
})

# =============================================================================
# EXPORTAR
# =============================================================================
cat("\n--- Exportando ---\n")

# Seleccionar columnas relevantes para output
out_cols <- c("drug_name_norm", "chembl_id", "drug_class", "drug_class_label",
              "n_sources", "sources", "max_phase", "is_approved", "hnscc_indication",
              "cmap_score", "n_de_genes", "n_up_genes", "n_down_genes", "de_genes",
              "primary_target",
              "s_pi_stat", "s_clinical", "s_cmap", "s_pathway", "s_network", "s_evidence",
              "composite_score", "bonus", "final_score")
out_cols_avail <- intersect(out_cols, colnames(scored))

write.table(scored %>% select(all_of(out_cols_avail)),
            "results/tables/10_all_candidates_scored.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(top20 %>% select(all_of(out_cols_avail)),
            "results/tables/10_top20_candidates.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Exportado: 10_all_candidates_scored.tsv (%d candidatos)\n", nrow(scored)))
cat(sprintf("Exportado: 10_top20_candidates.tsv (top %d)\n", top_n))

# =============================================================================
# RESUMEN FINAL
# =============================================================================
cat("\n=== RESUMEN ===\n")
cat(sprintf("  Candidatos evaluados: %d\n", nrow(scored)))
cat(sprintf("  Top %d seleccionados:\n", top_n))
class_dist <- top20 %>% count(drug_class_label)
for (i in seq_len(nrow(class_dist))) {
  cat(sprintf("    %s: %d\n", class_dist$drug_class_label[i], class_dist$n[i]))
}
cat(sprintf("  Score minimo top %d: %.4f\n", top_n, min(top20$final_score)))
cat(sprintf("  Score maximo:        %.4f\n", max(top20$final_score)))
cat(sprintf("\nSiguiente: scripts/11_clinicaltrials_pubmed.py\n"))
cat("Fin:", format(Sys.time()), "\n")

sink(type = "message"); sink(); close(con_log)
cat("Script completado. Log en:", log_file, "\n")
