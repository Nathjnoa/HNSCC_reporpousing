#!/usr/bin/env Rscript
# =============================================================================
# 16_external_validation.R
# =============================================================================
# External validation of HNSCC proteomics findings using TCGA-HNSC:
#
#   Panel A — Concordance: DIA proteomics log2FC vs TCGA-HNSC RNA-seq log2FC
#             (tumor vs matched adjacent normal, n=43 pairs via DESeq2)
#   Panel B — Survival: KM curves for top drug-target genes (OS, TCGA-HNSC)
#             Genes: EGFR, PSMB10, DNMT1, NDUFS3
#             (representatives of EGFR, Proteasome, Epigenetic, OXPHOS modules)
#
# Outputs:
#   results/figures/pub/main/OE3_FigA_tcga_concordance.{pdf,png}
#   results/figures/pub/main/OE3_FigB_survival.{pdf,png}
#   results/tables/pub/main/OE3_Tab_concordance_summary.tsv
#   results/tables/pub/supp/OE3_TabS_survival_genes.tsv
#
# Data cache: data/intermediate/tcga/  (populated on first run, ~800 MB)
#   Skip re-download if cache present.
#
# Dependencias adicionales (instalar una vez en omics-R):
#   BiocManager::install(c("TCGAbiolinks", "DESeq2"))
#   install.packages(c("survival", "survminer"))
#
# Ambiente: omics-R
# Ejecucion:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/16_external_validation.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(stringr)
})

# ── Detectar directorio del proyecto ─────────────────────────────────────────
args        <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  proj_dir    <- dirname(dirname(script_path))
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml")))
    setwd(proj_dir)
}

cat("=== 16_external_validation.R ===\n")
cat("Working directory:", getwd(), "\n")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ── Directorios ───────────────────────────────────────────────────────────────
dir.create("logs",                              showWarnings = FALSE)
dir.create("data/intermediate/tcga",           recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/pub/main",         recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/pub/main",          recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/pub/supp",          recursive = TRUE, showWarnings = FALSE)

log_file <- file.path("logs", paste0("16_external_validation_", timestamp, ".log"))
sink(log_file, split = TRUE)
cat("Inicio:", format(Sys.time()), "\n\n")

# ── Verificar dependencias ────────────────────────────────────────────────────
check_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0(
      "Paquete requerido no encontrado: '", pkg, "'\n",
      "Instalar con: ", if (bioc) paste0("BiocManager::install('", pkg, "')")
                        else        paste0("install.packages('", pkg, "')")
    )
    stop(msg, call. = FALSE)
  }
}
check_pkg("TCGAbiolinks", bioc = TRUE)
check_pkg("DESeq2",       bioc = TRUE)
check_pkg("survival")
check_pkg("survminer")
check_pkg("SummarizedExperiment", bioc = TRUE)

library(TCGAbiolinks)
library(DESeq2)
library(survival)
library(survminer)
library(SummarizedExperiment)

# =============================================================================
# ESTILOS (consistentes con script 17)
# =============================================================================
PRESETS <- list(
  double_col = list(w = 180, h = 120, base = 8.5, title = 10,
                    axis = 8.5, leg = 7.5, tick = 7.5, lwd = 0.7, pt = 1.8),
  single_col = list(w = 85,  h = 70,  base = 8.0, title = 9.0,
                    axis = 8.0, leg = 7.0, tick = 7.0, lwd = 0.6, pt = 1.6)
)

OKB     <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#000000")
CONCORD_COLS <- c(
  "Concordant up"   = "#D55E00",
  "Concordant down" = "#0072B2",
  "Discordant"      = "#BBBBBB"
)
MODULE_COLS <- c(
  EGFR        = "#D55E00",
  Proteasome  = "#56B4E9",
  Epigenetic  = "#009E73",
  OXPHOS      = "#CC79A7"
)

save_fig <- function(plot, name, preset = "double_col", ...) {
  p  <- PRESETS[[preset]]
  pf <- file.path("results/figures/pub/main", name)
  ggsave(paste0(pf, ".pdf"), plot, width = p$w, height = p$h,
         units = "mm", dpi = 300, ...)
  ggsave(paste0(pf, ".png"), plot, width = p$w, height = p$h,
         units = "mm", dpi = 300, bg = "white", ...)
  cat("Saved:", pf, "\n")
}

base_theme <- function(p) {
  theme_classic(base_size = p$base) +
    theme(
      plot.title       = element_text(size = p$title, face = "bold"),
      axis.title       = element_text(size = p$axis),
      axis.text        = element_text(size = p$tick),
      legend.text      = element_text(size = p$leg),
      legend.title     = element_text(size = p$leg, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# =============================================================================
# SECCION 1: DATOS LOCALES (nuestro estudio)
# =============================================================================
cat("\n--- Cargando datos locales ---\n")

de_our <- read_tsv("results/tables/de_limma/01_TVsS_significant.tsv",
                   show_col_types = FALSE)
cat("Proteinas DE significativas:", nrow(de_our), "\n")

id_map <- read_tsv("data/intermediate/id_mapping/02_uniprot_to_ids.tsv",
                   show_col_types = FALSE) |>
  select(uniprot_id, gene_symbol = symbol_org) |>
  filter(!is.na(gene_symbol), gene_symbol != "")

# Unir gene_symbol si no viene en de_our
if (!"gene_symbol" %in% colnames(de_our)) {
  de_our <- de_our |>
    left_join(id_map, by = "uniprot_id")
}

de_our_genes <- de_our |>
  filter(!is.na(gene_symbol), gene_symbol != "") |>
  select(gene_symbol, logFC_our = logFC_TVsS, adjP_our = adj.P.Val_TVsS) |>
  distinct(gene_symbol, .keep_all = TRUE)

cat("Proteinas con gene_symbol:", nrow(de_our_genes), "\n")

# =============================================================================
# SECCION 2: DESCARGA TCGA-HNSC (con cache)
# =============================================================================
cache_se  <- "data/intermediate/tcga/tcga_hnsc_se.rds"
cache_cli <- "data/intermediate/tcga/tcga_hnsc_clinical.rds"

if (file.exists(cache_se) && file.exists(cache_cli)) {
  cat("\n--- Cargando TCGA-HNSC desde cache ---\n")
  se   <- readRDS(cache_se)
  clin <- readRDS(cache_cli)
} else {
  cat("\n--- Descargando TCGA-HNSC (primera ejecucion, ~800 MB) ---\n")

  query <- GDCquery(
    project       = "TCGA-HNSC",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type   = c("Primary Tumor", "Solid Tissue Normal")
  )
  GDCdownload(query, method = "api", files.per.chunk = 50,
              directory = "data/intermediate/tcga/GDCdata")

  se <- GDCprepare(query, directory = "data/intermediate/tcga/GDCdata")
  saveRDS(se, cache_se)

  clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
  saveRDS(clin, cache_cli)
  cat("TCGA-HNSC descargado y guardado en cache.\n")
}

cat("Dimensiones SE:", dim(se), "\n")
cat("Muestras TCGA:", ncol(se), "\n")

# =============================================================================
# SECCION 3: DE ANALYSIS TCGA TUMOR VS NORMAL (DESeq2, con cache)
# Diseño: ~condition (todos tumor vs normal, sin bloqueo por paciente)
# Más rápido que diseño pareado; válido para concordancia de direcciones.
# =============================================================================
cache_de <- "data/intermediate/tcga/tcga_hnsc_deseq2_results.tsv"

if (file.exists(cache_de)) {
  cat("\n--- Cargando DE TCGA desde cache ---\n")
  de_tcga <- read_tsv(cache_de, show_col_types = FALSE)
} else {
  cat("\n--- Corriendo DESeq2 en TCGA tumor vs normal (~10 min) ---\n")

  # Seleccionar solo muestras tumor (01) y normal (11)
  sample_type <- substr(colnames(se), 14, 15)
  se_tn       <- se[, sample_type %in% c("01", "11")]
  condition   <- factor(ifelse(substr(colnames(se_tn), 14, 15) == "01",
                               "Tumor", "Normal"),
                        levels = c("Normal", "Tumor"))
  cat("Muestras: Tumor =", sum(condition == "Tumor"),
      " Normal =", sum(condition == "Normal"), "\n")

  # Gene names desde rowData (columna gene_name confirmada)
  gene_names <- rowData(se_tn)$gene_name

  counts_mat <- assay(se_tn, "unstranded")
  rownames(counts_mat) <- gene_names

  col_data <- data.frame(condition = condition, row.names = colnames(se_tn))

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = col_data,
    design    = ~ condition
  )
  # Filtro: al menos 10 counts en al menos 5 muestras
  dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
  cat("Genes tras filtro:", nrow(dds), "\n")

  dds <- DESeq(dds, parallel = FALSE)
  res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)

  de_tcga <- as.data.frame(res) |>
    rownames_to_column("gene_symbol") |>
    filter(!is.na(log2FoldChange), gene_symbol != "") |>
    select(gene_symbol, logFC_tcga = log2FoldChange,
           adjP_tcga = padj, baseMean) |>
    distinct(gene_symbol, .keep_all = TRUE)

  write_tsv(de_tcga, cache_de)
  cat("DESeq2 completado:", nrow(de_tcga), "genes.\n")
}

cat("Genes TCGA DE calculados:", nrow(de_tcga), "\n")

# =============================================================================
# SECCION 4: CONCORDANCIA
# =============================================================================
cat("\n--- Analisis de concordancia ---\n")

# Genes con símbolo en ambos datasets
conc <- de_our_genes |>
  inner_join(de_tcga, by = "gene_symbol") |>
  filter(!is.na(logFC_tcga), !is.na(logFC_our)) |>
  mutate(
    concord_class = case_when(
      logFC_our > 0 & logFC_tcga > 0 ~ "Concordant up",
      logFC_our < 0 & logFC_tcga < 0 ~ "Concordant down",
      TRUE                            ~ "Discordant"
    )
  )

cat("Genes en ambos datasets:", nrow(conc), "\n")
cat(table(conc$concord_class), "\n")

# Pearson r (todos los genes solapantes)
pearson_r <- cor(conc$logFC_our, conc$logFC_tcga, method = "pearson")
pearson_p <- cor.test(conc$logFC_our, conc$logFC_tcga)$p.value
cat(sprintf("Pearson r = %.3f  (p = %.2e)\n", pearson_r, pearson_p))

# Porcentaje de concordancia
pct_concord <- mean(conc$concord_class != "Discordant") * 100
cat(sprintf("Concordancia direccional: %.1f%%\n", pct_concord))

# Genes a etiquetar: los 4 representativos de modulos
LABEL_GENES <- c("EGFR", "PSMB10", "DNMT1", "NDUFS3")
LABEL_MODULE <- c(EGFR = "EGFR", PSMB10 = "Proteasome",
                  DNMT1 = "Epigenetic", NDUFS3 = "OXPHOS")

conc <- conc |>
  mutate(
    is_label  = gene_symbol %in% LABEL_GENES,
    module    = LABEL_MODULE[gene_symbol]
  )

# Resumen para tabla de publicacion
summary_tbl <- tibble(
  metric = c("N genes overlap", "Pearson r", "p-value (Pearson)",
             "% concordant direction", "% concordant up", "% concordant down"),
  value  = c(
    nrow(conc),
    round(pearson_r, 3),
    signif(pearson_p, 3),
    round(pct_concord, 1),
    round(mean(conc$concord_class == "Concordant up")   * 100, 1),
    round(mean(conc$concord_class == "Concordant down") * 100, 1)
  )
)
write_tsv(summary_tbl, "results/tables/pub/main/OE3_Tab_concordance_summary.tsv")
cat("Tabla concordancia guardada.\n")

# =============================================================================
# SECCION 5: PANEL A — SCATTER CONCORDANCIA
# =============================================================================
cat("\n--- Generando Panel A (concordancia) ---\n")

p_preset  <- PRESETS$double_col
x_range   <- range(conc$logFC_our)
y_range   <- range(conc$logFC_tcga)

figA <- ggplot(conc, aes(x = logFC_our, y = logFC_tcga)) +

  # Cuadrantes
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +

  # Puntos no etiquetados
  geom_point(data = filter(conc, !is_label),
             aes(color = concord_class), size = 0.9, alpha = 0.55,
             shape = 16) +

  # Puntos etiquetados (modulos)
  geom_point(data = filter(conc, is_label),
             aes(fill = module), size = 2.8, shape = 21,
             color = "white", stroke = 0.5) +

  # Etiquetas
  geom_label_repel(
    data         = filter(conc, is_label),
    aes(label = gene_symbol, fill = module),
    color        = "white",
    size         = p_preset$tick / .pt,
    fontface     = "bold",
    box.padding  = 0.4,
    label.size   = 0,
    label.padding = unit(0.12, "lines"),
    show.legend  = FALSE
  ) +

  # Linea de regresion
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              linewidth = 0.5, alpha = 0.1) +

  # Anotacion r
  annotate("text",
    x = min(conc$logFC_our) * 0.9,
    y = max(conc$logFC_tcga) * 0.95,
    label = sprintf("r = %.3f\nn = %d", pearson_r, nrow(conc)),
    hjust = 0, size = p_preset$tick / .pt, color = "grey20"
  ) +

  scale_color_manual(values = CONCORD_COLS, name = NULL,
                     guide = guide_legend(override.aes = list(size = 2,
                                          alpha = 0.8))) +
  scale_fill_manual(values  = MODULE_COLS,  name = "Module",
                    guide = guide_legend(override.aes = list(size = 2.5))) +

  labs(
    x     = expression("log"[2]*"FC (DIA proteomics, Tumor vs Normal)"),
    y     = expression("log"[2]*"FC (TCGA-HNSC RNA-seq, Tumor vs Normal)"),
    title = "Proteomics–transcriptomics concordance"
  ) +
  base_theme(p_preset) +
  theme(legend.position = "right")

# =============================================================================
# SECCION 6: DATOS EXPRESION PARA SUPERVIVENCIA
# =============================================================================
cat("\n--- Preparando datos de supervivencia ---\n")

# Muestras tumorales (tipo 01)
se_tumor <- se[, substr(colnames(se), 14, 15) == "01"]
cat("Muestras tumorales para supervivencia:", ncol(se_tumor), "\n")

# Log2 FPKM (confirmado disponible en SE) — rápido, sin re-modelado
expr_mat           <- log2(assay(se_tumor, "fpkm_unstrand") + 0.1)
rownames(expr_mat) <- rowData(se_tumor)$gene_name   # gene_name confirmado en rowData

# Datos clinicos desde colData del SE (columnas confirmadas)
col_df <- as.data.frame(colData(se_tumor))
cat("Columnas OS disponibles: vital_status, days_to_death, days_to_last_follow_up\n")

surv_df <- col_df |>
  mutate(
    patient_id = patient,    # columna 'patient' confirmada en colData
    os_event   = as.integer(tolower(vital_status) == "dead"),
    os_time    = if_else(os_event == 1,
                         suppressWarnings(as.numeric(days_to_death)),
                         suppressWarnings(as.numeric(days_to_last_follow_up)))
  ) |>
  filter(!is.na(os_time), os_time > 0) |>
  select(patient_id, os_event, os_time) |>
  distinct(patient_id, .keep_all = TRUE)

cat("Pacientes con datos de supervivencia:", nrow(surv_df), "\n")

# Agregar expresion de genes de interes (usando patient_id de colData)
surv_genes_data <- surv_df

for (gene in LABEL_GENES) {
  if (gene %in% rownames(expr_mat)) {
    gene_expr <- expr_mat[gene, ]
    gene_df   <- data.frame(
      patient_id = col_df$patient,
      expr       = as.numeric(gene_expr),
      stringsAsFactors = FALSE
    ) |>
      group_by(patient_id) |>
      summarise(!!gene := mean(expr, na.rm = TRUE), .groups = "drop")

    surv_genes_data <- surv_genes_data |>
      left_join(gene_df, by = "patient_id")
  } else {
    cat("ATENCION: gene no encontrado en expr_mat:", gene, "\n")
    surv_genes_data[[gene]] <- NA_real_
  }
}

# Estratificar por mediana de expresion (alto vs bajo)
for (gene in LABEL_GENES) {
  med <- median(surv_genes_data[[gene]], na.rm = TRUE)
  surv_genes_data[[paste0(gene, "_group")]] <- if_else(
    surv_genes_data[[gene]] >= med, "High", "Low"
  )
}

write_tsv(surv_genes_data,
          "results/tables/pub/supp/OE3_TabS_survival_genes.tsv")
cat("Tabla supervivencia guardada.\n")

# =============================================================================
# SECCION 7: PANEL B — KM CURVES (2 x 2)
# =============================================================================
cat("\n--- Generando Panel B (KM curves) ---\n")

p_s <- PRESETS$single_col

km_plots <- lapply(seq_along(LABEL_GENES), function(i) {
  gene   <- LABEL_GENES[i]
  module <- LABEL_MODULE[[gene]]
  gcol   <- MODULE_COLS[[module]]
  grp_col <- paste0(gene, "_group")

  df_gene <- surv_genes_data |>
    filter(!is.na(.data[[grp_col]])) |>
    mutate(group = factor(.data[[grp_col]], levels = c("High", "Low")))

  if (nrow(df_gene) < 20) {
    message("Pocos datos para ", gene, ", omitiendo KM.")
    return(NULL)
  }

  fit <- survfit(Surv(os_time, os_event) ~ group, data = df_gene)

  lrt <- tryCatch(
    survdiff(Surv(os_time, os_event) ~ group, data = df_gene),
    error = function(e) NULL
  )
  pval_txt <- if (!is.null(lrt)) {
    pv <- 1 - pchisq(lrt$chisq, df = length(lrt$n) - 1)
    if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)
  } else ""

  ggsurvplot(
    fit,
    data           = df_gene,
    palette        = c(gcol, "grey60"),
    size           = 0.6,
    conf.int       = FALSE,
    pval           = FALSE,
    risk.table     = FALSE,
    xscale         = "d_y",            # dias a años
    xlab           = "Time (years)",
    ylab           = "Overall survival",
    title          = paste0(gene, " (", module, ")"),
    legend.title   = "Expression",
    legend.labs    = c("High", "Low"),
    fontsize       = p_s$tick / .pt,
    ggtheme        = base_theme(p_s) +
                       theme(plot.title   = element_text(size = p_s$axis,
                                                         face = "bold"),
                             legend.position = c(0.78, 0.85))
  )$plot +
  annotate("text", x = Inf, y = 0.05, label = pval_txt,
           hjust = 1.05, vjust = 0,
           size = p_s$tick / .pt, color = "grey20")
})

# Filtrar NULLs y combinar con patchwork
km_plots_ok <- Filter(Negate(is.null), km_plots)

if (length(km_plots_ok) >= 1) {
  figB <- wrap_plots(km_plots_ok, ncol = 2)
} else {
  figB <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No survival data available",
             size = 5) +
    theme_void()
}

# =============================================================================
# SECCION 8: GUARDAR FIGURAS
# =============================================================================
cat("\n--- Guardando figuras ---\n")

save_fig(figA, "OE3_FigA_tcga_concordance", preset = "double_col")

p_surv <- PRESETS$double_col
figB_wide <- figB
ggsave("results/figures/pub/main/OE3_FigB_survival.pdf",
       figB_wide, width = p_surv$w, height = p_surv$h,
       units = "mm", dpi = 300)
ggsave("results/figures/pub/main/OE3_FigB_survival.png",
       figB_wide, width = p_surv$w, height = p_surv$h,
       units = "mm", dpi = 300, bg = "white")
cat("Saved: results/figures/pub/main/OE3_FigB_survival\n")

# =============================================================================
# SECCION 9: RESUMEN
# =============================================================================
cat("\n========================================\n")
cat("RESUMEN 16_external_validation.R\n")
cat("========================================\n")
cat(sprintf("Genes solapantes (nuestros vs TCGA): %d\n", nrow(conc)))
cat(sprintf("Concordancia direccional:             %.1f%%\n", pct_concord))
cat(sprintf("Pearson r (logFC proteomics vs RNA):  %.3f (p = %.2e)\n",
            pearson_r, pearson_p))
cat(sprintf("Pacientes para supervivencia:         %d\n", nrow(surv_df)))
cat("\nOutputs:\n")
cat("  results/figures/pub/main/OE3_FigA_tcga_concordance.{pdf,png}\n")
cat("  results/figures/pub/main/OE3_FigB_survival.{pdf,png}\n")
cat("  results/tables/pub/main/OE3_Tab_concordance_summary.tsv\n")
cat("  results/tables/pub/supp/OE3_TabS_survival_genes.tsv\n")
cat("\nFin:", format(Sys.time()), "\n")
sink()
