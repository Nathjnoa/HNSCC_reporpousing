#!/usr/bin/env Rscript
# =============================================================================
# Script 08: Integrar fuentes de drug targets
# HNSCC Drug Repurposing — Fase 3
# =============================================================================
# Objetivo: Unificar 4 fuentes (DGIdb, ChEMBL, Open Targets, L2S2) en una
#           tabla maestra. Clasificar farmacos:
#           A = aprobado + indicacion HNSCC directa
#           B = aprobado + otra indicacion oncologica (reposicionamiento intra-oncologico)
#           C = aprobado + indicacion no-oncologica (reposicionamiento clasico)
#           D = no aprobado (fase 3, experimental)
#
# Input:  results/tables/drug_targets/04_dgidb_raw.tsv
#         results/tables/drug_targets/05_chembl_drugs.tsv
#         results/tables/drug_targets/06_opentargets_gene_drugs.tsv
#         results/tables/drug_targets/07_cmap_top_reversors.tsv
#         results/tables/de_limma/02_TVsS_significant_with_ids.tsv
#         config/analysis_params.yaml
#
# Output: results/tables/drug_targets/08_drug_target_master_table.tsv
#         results/tables/drug_targets/08_drug_summary_per_drug.tsv
#         results/figures/08_drug_sources_upset.pdf
#         results/figures/08_drug_class_barplot.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(yaml)
})

# --- Log setup ---------------------------------------------------------------
dir.create("logs",    showWarnings = FALSE)
dir.create("results/tables/drug_targets", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures",             showWarnings = FALSE, recursive = TRUE)

log_file <- paste0("logs/08_integrate_drug_targets_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
con_log <- file(log_file, open = "wt")
sink(con_log, type = "output")
sink(con_log, type = "message", append = TRUE)

cat("=== 08_integrate_drug_targets.R ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# --- Parametros --------------------------------------------------------------
params     <- yaml::read_yaml("config/analysis_params.yaml")
min_db     <- params$candidates$min_databases  # 2

# --- Cargar tabla DE base (para gene metadata) --------------------------------
cat("--- Cargando tabla DE ---\n")
sig <- read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
                  stringsAsFactors = FALSE)
gene_meta <- sig %>%
  select(symbol_org, uniprot_id, logFC_TVsS, adj.P.Val_TVsS) %>%
  mutate(direction = ifelse(logFC_TVsS > 0, "up", "down")) %>%
  distinct(uniprot_id, .keep_all = TRUE)
cat(sprintf("Proteinas DE: %d\n\n", nrow(gene_meta)))

# =============================================================================
# 1. CARGAR Y NORMALIZAR CADA FUENTE
# =============================================================================
cat("--- 1. Cargando y normalizando fuentes ---\n\n")

# ---------------------------------------------------------------------------
# 1a. DGIdb
# ---------------------------------------------------------------------------
cat("  [DGIdb]\n")
dgidb_raw <- read.delim("results/tables/drug_targets/04_dgidb_raw.tsv",
                        stringsAsFactors = FALSE)

dgidb <- dgidb_raw %>%
  filter(!is.na(drug_name), drug_name != "") %>%
  mutate(
    drug_name_norm = str_to_upper(str_trim(drug_name)),
    # Extraer ChEMBL ID si viene en drug_concept_id o drug_chembl_id
    chembl_id = case_when(
      !is.na(drug_chembl_id) & str_detect(drug_chembl_id, "CHEMBL") ~
        str_extract(str_to_upper(drug_chembl_id), "CHEMBL[0-9]+"),
      !is.na(drug_concept_id) & str_detect(str_to_upper(drug_concept_id), "CHEMBL") ~
        str_extract(str_to_upper(drug_concept_id), "CHEMBL[0-9]+"),
      TRUE ~ NA_character_
    ),
    source_db      = "DGIdb",
    max_phase_db   = NA_real_,   # DGIdb no provee fase clinica directamente
    is_approved_db = NA,
    cmap_score     = NA_real_,
    gene_symbol    = str_to_upper(str_trim(query_gene))
  ) %>%
  select(gene_symbol, drug_name_norm, chembl_id, source_db,
         max_phase_db, is_approved_db, cmap_score,
         interaction_score, direction)

cat(sprintf("    DGIdb: %d filas, %d drugs unicos, %d genes\n",
            nrow(dgidb), n_distinct(dgidb$drug_name_norm), n_distinct(dgidb$gene_symbol)))

# ---------------------------------------------------------------------------
# 1b. ChEMBL
# ---------------------------------------------------------------------------
cat("  [ChEMBL]\n")
chembl_raw <- read.delim("results/tables/drug_targets/05_chembl_drugs.tsv",
                         stringsAsFactors = FALSE)

chembl <- chembl_raw %>%
  filter(!is.na(pref_name), pref_name != "") %>%
  mutate(
    drug_name_norm    = str_to_upper(str_trim(pref_name)),
    chembl_id         = str_to_upper(str_trim(molecule_chembl_id)),
    source_db         = "ChEMBL",
    max_phase_db      = as.numeric(max_phase),
    is_approved_db    = (max_phase >= 4),
    cmap_score        = NA_real_,
    interaction_score = NA_real_,
    gene_symbol       = str_to_upper(str_trim(
                          coalesce(symbol_org, gene_symbol)
                        )),
    direction         = ifelse(logFC_TVsS > 0, "up", "down")
  ) %>%
  select(gene_symbol, drug_name_norm, chembl_id, source_db,
         max_phase_db, is_approved_db, cmap_score,
         interaction_score, direction)

cat(sprintf("    ChEMBL: %d filas, %d drugs, %d genes\n",
            nrow(chembl), n_distinct(chembl$drug_name_norm), n_distinct(chembl$gene_symbol)))

# ---------------------------------------------------------------------------
# Keywords oncologicos para detectar indicaciones cancer (Open Targets)
# Se aplica sobre indication_name (texto libre de OT)
# ---------------------------------------------------------------------------
CANCER_KW <- paste(
  c("cancer", "carcinoma", "tumor", "tumour", "lymphoma",
    "leukemia", "leukaemia", "sarcoma", "melanoma",
    "neoplasm", "glioma", "myeloma", "adenocarcinoma",
    "mesothelioma", "blastoma", "malignant"),
  collapse = "|"
)

# ---------------------------------------------------------------------------
# 1c. Open Targets
# ---------------------------------------------------------------------------
cat("  [Open Targets]\n")
ot_raw <- read.delim("results/tables/drug_targets/06_opentargets_gene_drugs.tsv",
                     stringsAsFactors = FALSE)

ot <- ot_raw %>%
  filter(!is.na(drug_name), drug_name != "") %>%
  mutate(
    drug_name_norm      = str_to_upper(str_trim(drug_name)),
    chembl_id           = str_to_upper(str_trim(drug_id)),
    source_db           = "OpenTargets",
    max_phase_db        = as.numeric(max_phase_global),
    # Python escribe "True"/"False" — R necesita tolower() para parsear
    is_approved_db      = (tolower(as.character(is_approved)) == "true"),
    is_hnscc_indication = (tolower(as.character(is_hnscc_indication)) == "true"),
    # Flag: indicacion oncologica (cualquier cancer, incluye HNSCC)
    has_cancer_indication = str_detect(str_to_lower(indication_name), CANCER_KW),
    has_cancer_indication = ifelse(is.na(has_cancer_indication), FALSE, has_cancer_indication),
    cmap_score          = NA_real_,
    interaction_score   = NA_real_,
    gene_symbol         = str_to_upper(str_trim(symbol_org)),
    direction           = direction
  ) %>%
  select(gene_symbol, drug_name_norm, chembl_id, source_db,
         max_phase_db, is_approved_db, cmap_score,
         interaction_score, direction,
         is_hnscc_indication, has_cancer_indication)

cat(sprintf("    OpenTargets: %d filas, %d drugs, %d genes\n",
            nrow(ot), n_distinct(ot$drug_name_norm), n_distinct(ot$gene_symbol)))

# ---------------------------------------------------------------------------
# 1d. L2S2
# ---------------------------------------------------------------------------
cat("  [L2S2 — LINCS L1000 Signature Search]\n")
cmap_raw <- read.delim("results/tables/drug_targets/07_l2s2_top_reversors.tsv",
                       stringsAsFactors = FALSE)

# L2S2 usa scaled_score (negativo = reversor)
cmap_score_col <- if ("scaled_score" %in% colnames(cmap_raw)) "scaled_score" else
                  if ("raw_score"    %in% colnames(cmap_raw)) "raw_score"    else
                  colnames(cmap_raw)[which(sapply(cmap_raw, is.numeric))[1]]
cat(sprintf("    L2S2 score column detectada: %s\n", cmap_score_col))

# Solo tomamos los reversores (score < 0)
cmap <- cmap_raw %>%
  filter(type == "trt_cp", .data[[cmap_score_col]] < 0) %>%
  mutate(
    drug_name_norm      = str_to_upper(str_trim(pert)),
    chembl_id           = NA_character_,
    gene_symbol         = NA_character_,  # evidencia a nivel de firma, no gen-fármaco
    source_db           = "L2S2",
    max_phase_db        = NA_real_,
    is_approved_db      = TRUE,           # L2S2 filtrò a FDA-approved (filterFda=TRUE)
    cmap_score          = .data[[cmap_score_col]],
    direction           = NA_character_,
    interaction_score   = NA_real_,
    is_hnscc_indication = FALSE,
    has_cancer_indication = FALSE
  ) %>%
  select(gene_symbol, drug_name_norm, chembl_id, source_db,
         max_phase_db, is_approved_db, cmap_score,
         interaction_score, direction,
         is_hnscc_indication, has_cancer_indication)

cat(sprintf("    L2S2: %d compuestos reversores\n", nrow(cmap)))

# =============================================================================
# 2. COMBINAR EN TABLA MAESTRA (long format)
# =============================================================================
cat("\n--- 2. Combinando fuentes ---\n")

# Columnas ausentes en DGIdb y ChEMBL — asignar valores por defecto
dgidb$is_hnscc_indication   <- FALSE
dgidb$has_cancer_indication <- FALSE
chembl$is_hnscc_indication  <- FALSE
chembl$has_cancer_indication <- FALSE

# Asegurar tipos consistentes antes del bind
for (col in c("is_hnscc_indication", "has_cancer_indication")) {
  dgidb[[col]]  <- as.logical(dgidb[[col]])
  chembl[[col]] <- as.logical(chembl[[col]])
  ot[[col]]     <- as.logical(ot[[col]])
  cmap[[col]]   <- as.logical(cmap[[col]])
}

all_cols <- c("gene_symbol", "drug_name_norm", "chembl_id", "source_db",
              "max_phase_db", "is_approved_db", "cmap_score",
              "interaction_score", "direction",
              "is_hnscc_indication", "has_cancer_indication")

master_long <- bind_rows(
  dgidb  %>% select(any_of(all_cols)),
  chembl %>% select(any_of(all_cols)),
  ot     %>% select(any_of(all_cols)),
  cmap   %>% select(any_of(all_cols))
)

cat(sprintf("Filas totales (long): %d\n", nrow(master_long)))

# =============================================================================
# 3. CONSTRUIR TABLA RESUMEN POR FARMACO
# =============================================================================
cat("\n--- 3. Construyendo resumen por farmaco ---\n")

drug_summary <- master_long %>%
  group_by(drug_name_norm) %>%
  summarise(
    # Fuentes que lo detectaron
    sources           = paste(sort(unique(source_db)), collapse = "|"),
    n_sources         = n_distinct(source_db),

    # ChEMBL ID (el primero no-NA)
    chembl_id         = first(na.omit(chembl_id)),

    # Fase clinica maxima
    max_phase         = suppressWarnings(max(max_phase_db, na.rm = TRUE)),
    max_phase         = ifelse(is.infinite(max_phase), NA_real_, max_phase),

    # Aprobado en alguna fuente
    is_approved       = any(is_approved_db == TRUE, na.rm = TRUE),

    # Indicacion HNSCC directa (Open Targets)
    hnscc_indication      = any(is_hnscc_indication == TRUE, na.rm = TRUE),

    # Indicacion oncologica (cualquier cancer, desde OT)
    has_cancer_indication = any(has_cancer_indication == TRUE, na.rm = TRUE),

    # CMap connectivity score (negativo = mejor reversor)
    cmap_score        = suppressWarnings(min(cmap_score, na.rm = TRUE)),
    cmap_score        = ifelse(is.infinite(cmap_score), NA_real_, cmap_score),

    # Genes diana afectados (conteo de genes unicos, no filas)
    n_de_genes        = n_distinct(na.omit(gene_symbol)),
    de_genes          = paste(sort(unique(na.omit(gene_symbol))), collapse = "|"),

    # Direccion: genes unicos up vs down (no filas, para evitar doble conteo)
    n_up_genes        = n_distinct(gene_symbol[direction == "up"   & !is.na(gene_symbol)]),
    n_down_genes      = n_distinct(gene_symbol[direction == "down" & !is.na(gene_symbol)]),

    # Max interaction score (DGIdb)
    max_interaction_score = suppressWarnings(max(interaction_score, na.rm = TRUE)),
    max_interaction_score = ifelse(is.infinite(max_interaction_score),
                                   NA_real_, max_interaction_score),

    .groups = "drop"
  )

cat(sprintf("Farmacos unicos totales: %d\n", nrow(drug_summary)))

# =============================================================================
# 4. CLASIFICACION A / B / C / D
# =============================================================================
cat("\n--- 4. Clasificando farmacos ---\n")

drug_summary <- drug_summary %>%
  mutate(
    drug_class = case_when(
      # A: aprobado + evidencia directa en HNSCC (uso actual o ensayo avanzado)
      is_approved & hnscc_indication                           ~ "A",
      # B: aprobado para otro cancer — reposicionamiento intra-oncologico
      is_approved & has_cancer_indication & !hnscc_indication  ~ "B",
      # C: aprobado para indicacion no-oncologica — reposicionamiento clasico
      is_approved & !has_cancer_indication & !hnscc_indication ~ "C",
      # D: no aprobado (fase 3, experimental, solo CMap)
      TRUE                                                      ~ "D"
    ),
    drug_class_label = case_when(
      drug_class == "A" ~ "A: Approved + HNSCC evidence",
      drug_class == "B" ~ "B: Approved other cancer",
      drug_class == "C" ~ "C: Approved non-oncology",
      drug_class == "D" ~ "D: Not approved / Experimental"
    )
  )

class_counts <- drug_summary %>%
  count(drug_class, drug_class_label) %>%
  arrange(drug_class)

cat("\nDistribucion por clase:\n")
for (i in seq_len(nrow(class_counts))) {
  cat(sprintf("  %s: %d farmacos\n",
              class_counts$drug_class_label[i], class_counts$n[i]))
}

# ---------------------------------------------------------------------------
# 4b. Corrección manual: drogas aprobadas para HNSCC con indicación conocida
#     que pueden no estar correctamente clasificadas por EFO matching.
# ---------------------------------------------------------------------------
HNSCC_APPROVED_OVERRIDE <- c(
  "CETUXIMAB",      # FDA/EMA: SCCHNeck + colorectal; EFO mismatch frecuente
  "PEMBROLIZUMAB",  # FDA: HNSCC R/M (2nd line)
  "NIVOLUMAB"       # FDA: HNSCC R/M
)

n_overridden <- sum(
  toupper(drug_summary$drug_name_norm) %in% HNSCC_APPROVED_OVERRIDE &
    drug_summary$drug_class != "A"
)
if (n_overridden > 0) {
  cat(sprintf("\n  Corrigiendo %d drogas a Clase A (HNSCC override manual):\n",
              n_overridden))
  drug_summary <- drug_summary %>%
    mutate(
      drug_class = ifelse(
        toupper(drug_name_norm) %in% HNSCC_APPROVED_OVERRIDE,
        "A", drug_class
      ),
      drug_class_label = case_when(
        drug_class == "A" ~ "A: Approved + HNSCC evidence",
        drug_class == "B" ~ "B: Approved other cancer",
        drug_class == "C" ~ "C: Approved non-oncology",
        drug_class == "D" ~ "D: Not approved / Experimental"
      ),
      hnscc_indication = ifelse(
        toupper(drug_name_norm) %in% HNSCC_APPROVED_OVERRIDE,
        TRUE, hnscc_indication
      )
    )
  cat(paste(
    drug_summary$drug_name_norm[toupper(drug_summary$drug_name_norm) %in%
                                   HNSCC_APPROVED_OVERRIDE],
    collapse = ", "
  ), "\n")
}

# =============================================================================
# 5. FILTRAR CANDIDATOS CON SOPORTE EN >= min_db FUENTES
# =============================================================================
cat(sprintf("\n--- 5. Filtrando candidatos con soporte >= %d fuentes ---\n", min_db))

candidates <- drug_summary %>%
  filter(n_sources >= min_db | drug_class %in% c("A", "B")) %>%
  arrange(drug_class, desc(n_sources), desc(n_de_genes))

cat(sprintf("Candidatos multi-fuente: %d\n", nrow(candidates)))

cat("\nTop 30 candidatos:\n")
cat(sprintf("  %-35s %-5s %-5s %-10s %-30s\n",
            "Drug", "Class", "nSrc", "MaxPhase", "Sources"))
for (i in seq_len(min(30, nrow(candidates)))) {
  r <- candidates[i, ]
  cat(sprintf("  %-35s %-5s %-5s %-10s %-30s\n",
              r$drug_name_norm, r$drug_class,
              r$n_sources,
              ifelse(is.na(r$max_phase), "?", as.character(r$max_phase)),
              r$sources))
}

# =============================================================================
# 6. TABLA MAESTRA DETALLADA (con metadatos de genes)
# =============================================================================
cat("\n--- 6. Construyendo tabla maestra detallada ---\n")

# Join con metadatos DE por gen
master_detail <- master_long %>%
  left_join(gene_meta, by = c("gene_symbol" = "symbol_org")) %>%
  left_join(drug_summary %>%
              select(drug_name_norm, drug_class, drug_class_label,
                     n_sources, sources, max_phase, is_approved,
                     hnscc_indication, cmap_score, chembl_id),
            by = "drug_name_norm") %>%
  arrange(drug_class, drug_name_norm, gene_symbol)

cat(sprintf("Tabla maestra detallada: %d filas\n", nrow(master_detail)))

# =============================================================================
# 7. FIGURAS
# =============================================================================
cat("\n--- 7. Generando figuras ---\n")

safe_pdf <- function(file, expr, w = 8, h = 6) {
  tryCatch({
    pdf(file, width = w, height = h)
    force(expr)
    dev.off()
    cat(sprintf("    Figura guardada: %s\n", file))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("    WARN figura %s: %s\n", file, e$message))
  })
}

# 7a. Barplot de clases
safe_pdf("results/figures/08_drug_class_barplot.pdf", {
  p <- drug_summary %>%
    count(drug_class_label) %>%
    mutate(drug_class_label = factor(drug_class_label,
                                     levels = rev(sort(unique(drug_class_label))))) %>%
    ggplot(aes(x = drug_class_label, y = n,
               fill = drug_class_label)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = n), hjust = -0.2, size = 4) +
    scale_fill_manual(values = c(
      "A: Approved + HNSCC evidence"   = "#d62728",
      "B: Approved other cancer"       = "#ff7f0e",
      "C: Approved non-oncology"       = "#2ca02c",
      "D: Not approved / Experimental" = "#7f7f7f"
    )) +
    coord_flip() +
    expand_limits(y = max(drug_summary %>% count(drug_class_label) %>% pull(n)) * 1.15) +
    labs(title = "Drug candidates by classification",
         subtitle = "HNSCC drug repurposing — integrated 4 sources",
         x = NULL, y = "Number of drugs") +
    theme_bw(base_size = 13) +
    theme(panel.grid.major.y = element_blank())
  print(p)
})

# 7b. Barplot fuentes de soporte (n_sources)
safe_pdf("results/figures/08_drug_sources_support.pdf", {
  p <- drug_summary %>%
    count(n_sources) %>%
    mutate(n_sources = factor(n_sources)) %>%
    ggplot(aes(x = n_sources, y = n, fill = n_sources)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = n), vjust = -0.3, size = 4) +
    scale_fill_brewer(palette = "Blues") +
    labs(title = "Drug support across databases",
         subtitle = "Number of databases detecting each drug",
         x = "Number of supporting databases", y = "Number of drugs") +
    theme_bw(base_size = 13)
  print(p)
})

# 7c. Top 20 candidatos multi-fuente: lollipop por n_sources + n_de_genes
top_plot <- candidates %>%
  filter(drug_class != "D") %>%
  slice_head(n = 20)

if (nrow(top_plot) > 0) {
  safe_pdf("results/figures/08_top_candidates_lollipop.pdf", w = 10, h = 7, {
    p <- top_plot %>%
      mutate(drug_name_norm = factor(drug_name_norm,
                                     levels = rev(drug_name_norm))) %>%
      ggplot(aes(x = drug_name_norm, y = n_de_genes,
                 color = drug_class_label)) +
      geom_segment(aes(xend = drug_name_norm, y = 0, yend = n_de_genes),
                   linewidth = 1) +
      geom_point(aes(size = n_sources), alpha = 0.9) +
      scale_color_manual(values = c(
        "A: Approved + HNSCC evidence"   = "#d62728",
        "B: Approved other cancer"       = "#ff7f0e",
        "C: Approved non-oncology"       = "#2ca02c",
        "D: Not approved / Experimental" = "#7f7f7f"
      )) +
      scale_size_continuous(range = c(3, 8), breaks = 1:4) +
      coord_flip() +
      labs(title = "Top multi-source drug candidates",
           subtitle = "Size = number of supporting databases",
           x = NULL, y = "Number of DE target genes",
           color = "Classification", size = "N databases") +
      theme_bw(base_size = 12) +
      theme(legend.position = "right",
            panel.grid.major.y = element_blank())
    print(p)
  })
}

# 7d. Scatter: CMap score vs n_de_genes para candidatos con CMap
cmap_cands <- candidates %>%
  filter(!is.na(cmap_score))

if (nrow(cmap_cands) >= 3) {
  safe_pdf("results/figures/08_cmap_vs_targets.pdf", w = 9, h = 6, {
    p <- cmap_cands %>%
      ggplot(aes(x = n_de_genes, y = cmap_score,
                 color = drug_class_label,
                 label = drug_name_norm)) +
      geom_point(size = 2.5, alpha = 0.8) +
      ggrepel::geom_text_repel(size = 2.8, max.overlaps = 15) +
      scale_color_manual(values = c(
        "A: Approved + HNSCC evidence"   = "#d62728",
        "B: Approved other cancer"       = "#ff7f0e",
        "C: Approved non-oncology"       = "#2ca02c",
        "D: Not approved / Experimental" = "#7f7f7f"
      )) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      labs(title = "CMap connectivity score vs. number of DE target genes",
           subtitle = "Negative score = tumor signature reversal; more targets = broader effect",
           x = "DE target genes", y = "CMap scaled score (negative = reversal)",
           color = "Classification") +
      theme_bw(base_size = 12)
    print(p)
  })
}

cat("\nVerificación clasificación Clase A:\n")
drug_summary %>%
  filter(drug_class == "A") %>%
  select(drug_name_norm, drug_class, hnscc_indication) %>%
  print()

# =============================================================================
# 8. EXPORTAR
# =============================================================================
cat("\n--- 8. Exportando ---\n")

write.table(master_detail,
            "results/tables/drug_targets/08_drug_target_master_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(drug_summary %>% arrange(drug_class, desc(n_sources)),
            "results/tables/drug_targets/08_drug_summary_per_drug.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(candidates,
            "results/tables/drug_targets/08_multi_source_candidates.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Exportado: 08_drug_target_master_table.tsv (%d filas)\n",
            nrow(master_detail)))
cat(sprintf("Exportado: 08_drug_summary_per_drug.tsv (%d farmacos)\n",
            nrow(drug_summary)))
cat(sprintf("Exportado: 08_multi_source_candidates.tsv (%d candidatos)\n",
            nrow(candidates)))

# =============================================================================
# RESUMEN FINAL
# =============================================================================
cat("\n=== RESUMEN ===\n")
cat(sprintf("  Fuentes integradas:          4 (DGIdb, ChEMBL, OpenTargets, L2S2)\n"))
cat(sprintf("  Farmacos totales:            %d\n", nrow(drug_summary)))
for (i in seq_len(nrow(class_counts))) {
  cat(sprintf("    Clase %s: %d\n", class_counts$drug_class[i], class_counts$n[i]))
}
cat(sprintf("  Candidatos multi-fuente (>=%d): %d\n", min_db, nrow(candidates)))
cat(sprintf("  Con CMap score:              %d\n", sum(!is.na(drug_summary$cmap_score))))
cat(sprintf("\nSiguiente: scripts/09_string_network.R\n"))
cat("Fin:", format(Sys.time()), "\n")

sink(type = "message"); sink(); close(con_log)
cat("Script completado. Log en:", log_file, "\n")
