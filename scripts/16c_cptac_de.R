#!/usr/bin/env Rscript
# =============================================================================
# 16c_cptac_de.R — DE proteómico CPTAC-HNSCC (limma) + diagnóstico de concordancia
# =============================================================================
# Lee la matriz proteómica (log2-ratios TMT) exportada por 16b_cptac_fetch.py y
# calcula DE tumor-vs-normal con limma, el MISMO método que nuestra proteómica de
# descubrimiento (proteoDA/limma), garantizando comparabilidad metodológica.
#
# Diseño: PAREADO (66 pares T/N identificados por patient_id) →
#   duplicateCorrelation(block = patient_id) + lmFit(~condition) + eBayes
# Si por algún motivo no hay pares suficientes, cae a limma de dos grupos.
#
# Outputs:
#   data/intermediate/cptac/16c_cptac_hnscc_de.tsv
#     gene_symbol, logFC_cptac, adjP_cptac, n_tumor, n_normal
#
# Diagnóstico impreso al log (Fase 1 del plan — decide inclusión en Fig6):
#   • Pearson r  logFC_our vs logFC_cptac (genes solapantes)
#   • % concordancia direccional
#   • Tabla 14 dianas-ancla del shortlist con logFC_our, logFC_cptac, FDR, Y/N
#
# Ambiente: omics-R
# Ejecución (desde raíz del proyecto, DESPUÉS de 16b_cptac_fetch.py):
#   Rscript scripts/16c_cptac_de.R
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
  library(yaml)
})

cat("=== 16c_cptac_de.R ===\n")
source(here::here("scripts", "_setup.R")); setup_project()

# =============================================================================
# SECCIÓN 1: PARÁMETROS
# =============================================================================
CPTAC_DIR  <- "data/intermediate/cptac"
OUT_DIR    <- CPTAC_DIR
OUT_FILE   <- file.path(OUT_DIR, "16c_cptac_hnscc_de.tsv")
MIN_PAIRS  <- 10   # mínimo de pares para usar diseño pareado

params     <- yaml::read_yaml(here::here("config", "analysis_params.yaml"))
N_SHORTLIST <- 14

# =============================================================================
# SECCIÓN 2: LEER DATOS CPTAC (output de 16b)
# =============================================================================
cat("\n--- Leyendo datos proteómicos CPTAC ---\n")

prot_path <- file.path(CPTAC_DIR, "cptac_hnscc_proteomics.tsv")
samp_path <- file.path(CPTAC_DIR, "cptac_hnscc_samples.tsv")

if (!file.exists(prot_path) || !file.exists(samp_path))
  stop("Faltan outputs de 16b_cptac_fetch.py en ", CPTAC_DIR,
       "\n  Ejecutar antes: python scripts/16b_cptac_fetch.py")

prot_mat <- read_tsv(prot_path, show_col_types = FALSE) |>
  column_to_rownames("gene_symbol") |>
  as.matrix()

samples <- read_tsv(samp_path, show_col_types = FALSE)

cat(sprintf("  Matriz proteómica: %d genes × %d muestras\n",
            nrow(prot_mat), ncol(prot_mat)))
cat(sprintf("  Tumores:  %d  |  Normales: %d\n",
            sum(samples$condition == "Tumor"),
            sum(samples$condition == "Normal")))

# Alinear orden de columnas de la matriz con el orden de samples
samples    <- samples |> filter(sample_id %in% colnames(prot_mat))
prot_mat   <- prot_mat[, samples$sample_id]

# =============================================================================
# SECCIÓN 3: DISEÑO EXPERIMENTAL Y LIMMA
# =============================================================================
cat("\n--- Diseño limma ---\n")

condition  <- factor(samples$condition, levels = c("Normal", "Tumor"))
patient_id <- samples$patient_id

paired_ids <- intersect(
  patient_id[condition == "Normal"],
  patient_id[condition == "Tumor"]
)
n_pairs <- length(paired_ids)
cat(sprintf("  Pares completos T/N: %d\n", n_pairs))

design <- model.matrix(~ condition)

if (n_pairs >= MIN_PAIRS) {
  cat("  Diseño: PAREADO (duplicateCorrelation + block = patient_id)\n")
  # Primera pasada para estimar la correlación inter-bloque
  fit0    <- lmFit(prot_mat, design)
  corfit  <- duplicateCorrelation(prot_mat, design, block = patient_id)
  cat(sprintf("  Correlación entre-bloque (inter-paciente): %.3f\n",
              corfit$consensus.correlation))
  # Segunda pasada con la correlación estimada
  fit     <- lmFit(prot_mat, design, block = patient_id,
                   correlation = corfit$consensus.correlation)
  de_mode <- "paired (duplicateCorrelation)"
} else {
  cat("  Diseño: DOS GRUPOS (sin emparejamiento, n_pairs < MIN_PAIRS)\n")
  fit     <- lmFit(prot_mat, design)
  de_mode <- "two-group"
}

fit   <- eBayes(fit)
tt    <- topTable(fit, coef = "conditionTumor", number = Inf, sort.by = "none")

# =============================================================================
# SECCIÓN 4: EXPORTAR RESULTADOS DE DE
# =============================================================================
cat("\n--- Exportando DE CPTAC ---\n")

n_tumor  <- sum(samples$condition == "Tumor")
n_normal <- sum(samples$condition == "Normal")

de_cptac <- tibble(
  gene_symbol  = rownames(tt),
  logFC_cptac  = tt$logFC,
  adjP_cptac   = tt$adj.P.Val,
  n_tumor      = n_tumor,
  n_normal     = n_normal
) |>
  filter(!is.na(gene_symbol), gene_symbol != "") |>
  distinct(gene_symbol, .keep_all = TRUE)

write_tsv(de_cptac, OUT_FILE)
cat(sprintf("  Exportado: %s  (%d genes)\n", OUT_FILE, nrow(de_cptac)))

# =============================================================================
# SECCIÓN 5: DIAGNÓSTICO DE CONCORDANCIA (decide Fase 2)
# =============================================================================
cat("\n========================================\n")
cat("DIAGNÓSTICO DE CONCORDANCIA — CPTAC\n")
cat("========================================\n")

# Cargar nuestro proteoma (igual que script 16)
de_our_raw <- read_tsv("results/tables/de_limma/01_TVsS_significant.tsv",
                        show_col_types = FALSE)
id_map     <- read_tsv("data/intermediate/id_mapping/02_uniprot_to_ids.tsv",
                        show_col_types = FALSE)

# Añadir gene_symbol vía id_map solo si no está ya en la tabla (igual que script 16)
if (!"gene_symbol" %in% names(de_our_raw)) {
  de_our_raw <- de_our_raw |>
    left_join(id_map |> dplyr::select(uniprot_id, gene_symbol = symbol_org),
              by = "uniprot_id")
}

de_our <- de_our_raw |>
  filter(!is.na(gene_symbol), gene_symbol != "") |>
  dplyr::select(gene_symbol, logFC_our = logFC_TVsS, adjP_our = adj.P.Val_TVsS) |>
  distinct(gene_symbol, .keep_all = TRUE)

# Concordancia global
conc <- de_our |>
  inner_join(de_cptac |> select(gene_symbol, logFC_cptac, adjP_cptac),
             by = "gene_symbol")

n_overlap   <- nrow(conc)
pearson_r   <- cor(conc$logFC_our, conc$logFC_cptac, method = "pearson")
pearson_p   <- cor.test(conc$logFC_our, conc$logFC_cptac)$p.value
pct_concord <- mean(sign(conc$logFC_our) == sign(conc$logFC_cptac)) * 100

cat(sprintf("  Genes solapantes (nuestro vs CPTAC): %d\n", n_overlap))
cat(sprintf("  Pearson r:                           %.3f  (p = %.2e)\n",
            pearson_r, pearson_p))
cat(sprintf("  Concordancia direccional:            %.1f%%\n", pct_concord))
cat(sprintf("  Método DE CPTAC:                     %s\n", de_mode))

# Dianas del shortlist
cat("\n--- 14 dianas-ancla del shortlist en CPTAC ---\n")

scored <- read_tsv("results/tables/10_all_candidates_scored.tsv",
                   show_col_types = FALSE)

w_target <- params$scoring_v3$weight_target
w_drug   <- params$scoring_v3$weight_drug

shortlist <- scored |>
  filter(tier != "off_network") |>
  group_by(primary_target) |>
  slice_max(composite_score, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(composite_score, n = N_SHORTLIST, with_ties = FALSE)

targets_diag <- shortlist |>
  select(primary_target, drug = drug_name_norm, composite_score) |>
  left_join(de_our |> select(gene_symbol, logFC_our),
            by = c("primary_target" = "gene_symbol")) |>
  left_join(de_cptac |> select(gene_symbol, logFC_cptac, adjP_cptac),
            by = c("primary_target" = "gene_symbol")) |>
  mutate(
    concord = case_when(
      !is.na(logFC_our) & !is.na(logFC_cptac) &
        sign(logFC_our) == sign(logFC_cptac) ~ "Y",
      TRUE ~ "N"
    ),
    sig_cptac = case_when(
      !is.na(adjP_cptac) & adjP_cptac < 0.05 ~ "Y",
      !is.na(adjP_cptac)                       ~ "N (ns)",
      TRUE                                      ~ "N (ND)"
    )
  )

# Imprimir tabla
cat(sprintf("%-12s %-20s %7s %8s %8s %10s %10s\n",
            "Target", "Drug", "Score", "logFC_ours", "logFC_cptac", "FDR_cptac", "Concord"))
cat(strrep("-", 82), "\n")
for (i in seq_len(nrow(targets_diag))) {
  r <- targets_diag[i, ]
  cat(sprintf("%-12s %-20s %7.3f %8.3f %8s %10s %10s\n",
              r$primary_target,
              substr(r$drug, 1, 20),
              r$composite_score,
              ifelse(is.na(r$logFC_our),   NA_real_, r$logFC_our),
              ifelse(is.na(r$logFC_cptac), "ND", sprintf("%.3f", r$logFC_cptac)),
              ifelse(is.na(r$adjP_cptac),  "ND", sprintf("%.3e", r$adjP_cptac)),
              r$concord))
}

n_data     <- sum(!is.na(targets_diag$logFC_cptac))
n_concord  <- sum(targets_diag$concord == "Y", na.rm = TRUE)
n_sig      <- sum(targets_diag$sig_cptac == "Y", na.rm = TRUE)

cat(strrep("-", 82), "\n")
cat(sprintf("Con datos en CPTAC:      %d / %d\n", n_data, N_SHORTLIST))
cat(sprintf("Concordantes prot/prot:  %d / %d  (%.0f%%)\n",
            n_concord, n_data, 100 * n_concord / max(n_data, 1)))
cat(sprintf("Significativas FDR<0.05: %d / %d\n", n_sig, n_data))

cat("\n========================================\n")
cat("RESUMEN 16c_cptac_de.R\n")
cat("========================================\n")
cat(sprintf("Fuente:                CPTAC-HNSCC (umich, TMT proteomics)\n"))
cat(sprintf("Método DE:             %s\n", de_mode))
cat(sprintf("Genes con DE:          %d\n", nrow(de_cptac)))
cat(sprintf("Genes solapantes:      %d\n", n_overlap))
cat(sprintf("Pearson r (global):    %.3f  (p = %.2e)\n", pearson_r, pearson_p))
cat(sprintf("Concordancia dir.:     %.1f%%\n", pct_concord))
cat(sprintf("Dianas shortlist:      %d / 14 con datos en CPTAC\n", n_data))
cat(sprintf("  Concordantes:        %d / %d\n", n_concord, n_data))
cat(sprintf("  FDR<0.05 en CPTAC:  %d / %d\n", n_sig, n_data))
cat("\nOutput:\n")
cat(sprintf("  %s\n", OUT_FILE))
cat("\n>>> DECISIÓN (Fase 2): revisar r y concordancia de dianas con el usuario <<<\n")
cat("\nFin:", format(Sys.time()), "\n")
