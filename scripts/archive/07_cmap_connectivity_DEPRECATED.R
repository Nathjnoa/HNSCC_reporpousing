#!/usr/bin/env Rscript
# =============================================================================
# Script 07: CMap/LINCS Connectivity Analysis
# HNSCC Drug Repurposing — Fase 3
# =============================================================================
# ============================================================
# DEPRECADO: Este script usa CMap2 (signatureSearch/ExperimentHub).
# Fue reemplazado por scripts/07_l2s2_connectivity.py que usa
# L2S2 (LINCS L1000, 248 lineas celulares, API GraphQL, NAR 2025).
# Conservado como referencia. NO ejecutar en el pipeline actual.
# ============================================================
# Objetivo: buscar compuestos que reviertan la firma proteomica tumoral TVsS.
# Usa signatureSearch + CMap2 via ExperimentHub (descarga ~3.5GB si no esta cacheada).
#
# Input:  results/tables/de_limma/02_TVsS_significant_with_ids.tsv
#         results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv
#         config/analysis_params.yaml
# Output: results/tables/drug_targets/07_cmap_results.tsv
#         results/tables/drug_targets/07_cmap_top_reversors.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(signatureSearch)
  library(ExperimentHub)
  library(SummarizedExperiment)
  library(yaml)
  library(dplyr)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1]) && a[1] != "") a[1] else b

# --- Log setup ---------------------------------------------------------------
dir.create("logs", showWarnings = FALSE)
dir.create("results/tables/drug_targets", showWarnings = FALSE, recursive = TRUE)
log_file <- paste0("logs/07_cmap_connectivity_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
con_log  <- file(log_file, open = "wt")
sink(con_log, type = "output")
sink(con_log, type = "message", append = TRUE)

cat("=== 07_cmap_connectivity.R ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# --- Parametros --------------------------------------------------------------
params <- yaml::read_yaml("config/analysis_params.yaml")
cmap_score_threshold <- params$cmap$connectivity_score_threshold  # -90

# --- Cargar firma DE ---------------------------------------------------------
cat("--- Cargando firma DE TVsS ---\n")
sig <- read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
                  stringsAsFactors = FALSE)
all_prot <- read.delim("results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv",
                       stringsAsFactors = FALSE)

# CMap2 usa Entrez IDs (como character). Usamos entrez_id de script 02.
up_genes <- sig %>%
  filter(logFC_TVsS > 0, !is.na(entrez_id)) %>%
  mutate(entrez_char = as.character(as.integer(entrez_id))) %>%
  pull(entrez_char) %>% unique()

down_genes <- sig %>%
  filter(logFC_TVsS < 0, !is.na(entrez_id)) %>%
  mutate(entrez_char = as.character(as.integer(entrez_id))) %>%
  pull(entrez_char) %>% unique()

cat(sprintf("Genes upregulated:   %d (Entrez IDs)\n", length(up_genes)))
cat(sprintf("Genes downregulated: %d (Entrez IDs)\n", length(down_genes)))

# --- Paso 1: Cargar CMap2 reference database ---------------------------------
cat("\n--- Paso 1: Cargando CMap2 via ExperimentHub ---\n")
cat("(Si no esta cacheado descargara ~3.5GB — puede tardar varios minutos)\n")
cat("Inicio descarga/carga:", format(Sys.time()), "\n")

eh <- ExperimentHub()

# Usamos cmap_expr (EH3224) = matrix de expresion CMap2 normalizada
# Alternativa: cmap_rank (EH3225) para metodo de ranking
tryCatch({
  cmap_db_path <- eh[["EH3225"]]  # CMap2 cmap_rank HDF5 — requerido por gess_cmap()
  cat("CMap2 DB cargado:", format(Sys.time()), "\n")
  cat("Path:", cmap_db_path, "\n")
}, error = function(e) {
  cat("ERROR cargando CMap2:", e$message, "\n")
  cat("Intentando metodo alternativo...\n")
})

# Verificar el archivo HDF5
if (exists("cmap_db_path") && file.exists(cmap_db_path)) {
  cat(sprintf("Tamano archivo CMap2: %.1f MB\n",
              file.info(cmap_db_path)$size / 1e6))
} else {
  cat("ADVERTENCIA: No se pudo cargar CMap2 DB\n")
  sink(type = "message"); sink(); close(con_log)
  quit(status = 1)
}

# --- Paso 2: Verificar overlap con landmark genes ----------------------------
cat("\n--- Paso 2: Overlap con genes CMap2 ---\n")

# Obtener los nombres de genes en CMap2 usando rhdf5
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  cat("ERROR: el paquete 'rhdf5' no esta instalado.\n")
  cat("Instalar con: BiocManager::install('rhdf5')\n")
  sink(type = "message"); sink(); close(con_log)
  quit(status = 1)
}
suppressPackageStartupMessages(library(rhdf5))

h5_contents <- tryCatch(
  h5ls(cmap_db_path),
  error = function(e) { cat("Error leyendo HDF5:", e$message, "\n"); NULL }
)

if (!is.null(h5_contents)) {
  cat("Estructura HDF5:\n")
  print(head(h5_contents, 20))

  # Intentar leer row names (genes) del HDF5
  tryCatch({
    gene_ids <- h5read(cmap_db_path, "rownames")
    cat(sprintf("\nGenes en CMap2: %d\n", length(gene_ids)))
    cat("Primeros 10 genes:", paste(head(gene_ids, 10), collapse=", "), "\n")

    overlap_up   <- intersect(toupper(up_genes),   toupper(gene_ids))
    overlap_down <- intersect(toupper(down_genes), toupper(gene_ids))

    cat(sprintf("Overlap up (nuestros / CMap): %d / %d\n",
                length(overlap_up), length(up_genes)))
    cat(sprintf("Overlap down (nuestros / CMap): %d / %d\n",
                length(overlap_down), length(down_genes)))
  }, error = function(e) {
    cat("Error leyendo rownames del HDF5:", e$message, "\n")
  })
}

# --- Paso 3: Construir qSig y correr GESS ------------------------------------
cat("\n--- Paso 3: Construyendo query signature y corriendo GESS ---\n")

# La funcion qSig construye la firma de busqueda
# up   = genes que se sobre-expresan en TUMOR vs NORMAL
# dn   = genes que se sub-expresan en TUMOR vs NORMAL
# Compuesto ideal: deberia BAJAR los genes up y SUBIR los genes down

# Usar top 150 up y top 150 down (mas genes = mas especificidad)
# Rankeados por |logFC|
top_up <- sig %>%
  filter(logFC_TVsS > 0, !is.na(entrez_id)) %>%
  arrange(desc(logFC_TVsS)) %>%
  slice_head(n = 150) %>%
  mutate(entrez_char = as.character(as.integer(entrez_id))) %>%
  pull(entrez_char)

top_down <- sig %>%
  filter(logFC_TVsS < 0, !is.na(entrez_id)) %>%
  arrange(logFC_TVsS) %>%
  slice_head(n = 150) %>%
  mutate(entrez_char = as.character(as.integer(entrez_id))) %>%
  pull(entrez_char)

cat(sprintf("Firma query: %d up + %d down genes\n",
            length(top_up), length(top_down)))

# Construir qSig
qsig <- tryCatch({
  qSig(
    query    = list(upset = top_up, downset = top_down),
    gess_method = "CMAP",
    refdb    = cmap_db_path
  )
}, error = function(e) {
  cat("ERROR construyendo qSig:", e$message, "\n")
  NULL
})

if (is.null(qsig)) {
  cat("No se pudo construir qSig. Revisa el formato del HDF5.\n")
  sink(type = "message"); sink(); close(con_log)
  quit(status = 1)
}

cat("qSig construido. Corriendo gess_cmap()...\n")
cat("Inicio busqueda CMap:", format(Sys.time()), "\n")

cmap_res <- tryCatch({
  gess_cmap(
    qSig    = qsig,
    chunk_size = 5000,
    ref_trts = NULL,
    workers  = 1
  )
}, error = function(e) {
  cat("ERROR en gess_cmap:", e$message, "\n")
  NULL
})

if (is.null(cmap_res)) {
  cat("gess_cmap fallo. Abortando.\n")
  sink(type = "message"); sink(); close(con_log)
  quit(status = 1)
}

cat("gess_cmap completado:", format(Sys.time()), "\n")

# --- Paso 4: Procesar resultados --------------------------------------------
cat("\n--- Paso 4: Procesando resultados ---\n")

df_res <- result(cmap_res)
cat(sprintf("Total compuestos evaluados: %d\n", nrow(df_res)))

# Columnas clave: WTCS (Weighted Tau Connectivity Score) - negativo = reversor
# NCS: normalized connectivity score
cat("Columnas disponibles:", paste(colnames(df_res), collapse = ", "), "\n")

# gess_cmap() genera: pert, cell, type, trend, raw_score, scaled_score, ...
# scaled_score: -100 a +100 (negativo = reversal potente)
# raw_score: score crudo
score_col <- if ("scaled_score" %in% colnames(df_res)) "scaled_score" else
             if ("raw_score"    %in% colnames(df_res)) "raw_score"    else
             colnames(df_res)[which(sapply(df_res, is.numeric))[1]]

cat(sprintf("Score column: %s\n", score_col))
cat(sprintf("Score range: %.3f to %.3f\n",
            min(df_res[[score_col]], na.rm=TRUE),
            max(df_res[[score_col]], na.rm=TRUE)))

score_vals <- as.numeric(df_res[[score_col]])
cat(sprintf("Score 5th percentile: %.3f\n", quantile(score_vals, 0.05, na.rm=TRUE)))
cat(sprintf("Score 1st percentile: %.3f\n", quantile(score_vals, 0.01, na.rm=TRUE)))

# Filtrar: solo small molecules (type = "trt_cp") + bottom 5% de scaled_score
threshold_5pct <- quantile(score_vals, 0.05, na.rm = TRUE)
reversors <- df_res %>%
  filter(type == "trt_cp",
         as.numeric(.data[[score_col]]) <= threshold_5pct) %>%
  arrange(as.numeric(.data[[score_col]]))

cat(sprintf("\nReversores potentes (bottom 5%%): %d compuestos\n", nrow(reversors)))

# Top 20 reversores
top20 <- head(reversors, 20)
cat("\nTop 20 compuestos reversores:\n")
for (i in seq_len(nrow(top20))) {
  r <- top20[i, ]
  pert_name <- if ("pert" %in% colnames(r)) r$pert else as.character(r[1, 1])
  cat(sprintf("  %2d. %-35s score=%.4f\n", i, pert_name, r[[score_col]]))
}

# --- Paso 5: Enriquecimiento funcional (opcional) ----------------------------
cat("\n--- Paso 5: Drug set enrichment analysis (DSEA) ---\n")
dsea_res <- tryCatch({
  dsea_hyperG(
    drugs    = reversors$pert[seq_len(min(100, nrow(reversors)))],
    type     = "GO",
    ont      = "BP",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500
  )
}, error = function(e) {
  cat("DSEA no disponible o error:", e$message, "\n")
  NULL
})

if (!is.null(dsea_res)) {
  df_dsea <- result(dsea_res)
  cat(sprintf("DSEA: %d terminos GO BP significativos\n", nrow(df_dsea)))
  write.table(df_dsea,
              "results/tables/drug_targets/07_cmap_dsea_go_bp.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Exportado: 07_cmap_dsea_go_bp.tsv\n")
} else {
  cat("DSEA sin resultados — tabla no generada\n")
}

# --- Paso 6: Exportar --------------------------------------------------------
cat("\n--- Paso 6: Exportando resultados ---\n")

# Todos los resultados ordenados
df_export <- df_res %>% arrange(.data[[score_col]])
write.table(df_export,
            "results/tables/drug_targets/07_cmap_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Top reversores
write.table(head(df_export, 200),
            "results/tables/drug_targets/07_cmap_top_reversors.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Exportado: 07_cmap_results.tsv (%d compuestos)\n", nrow(df_export)))
cat(sprintf("Exportado: 07_cmap_top_reversors.tsv (top 200)\n"))
if (!is.null(dsea_res)) cat("Exportado: 07_cmap_dsea_go_bp.tsv\n")

# Resumen
cat("\n=== RESUMEN ===\n")
cat(sprintf("  Compuestos evaluados: %d\n", nrow(df_res)))
cat(sprintf("  Reversores (bottom 5%%): %d\n", nrow(reversors)))
cat(sprintf("  Score minimo (mejor reversor): %.3f\n",
            min(df_res[[score_col]], na.rm=TRUE)))
cat(sprintf("  Score maximo: %.3f\n",
            max(df_res[[score_col]], na.rm=TRUE)))

cat(sprintf("\nSiguiente: scripts/08_integrate_drug_targets.R\n"))
cat("Fin:", format(Sys.time()), "\n")

sink(type = "message"); sink(); close(con_log)
cat("Script completado. Log en:", log_file, "\n")
