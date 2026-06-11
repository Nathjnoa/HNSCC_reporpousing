#!/usr/bin/env Rscript
# =============================================================================
# 19_methylation_tcga.R
# =============================================================================
# Methylation analysis of promoter CpG islands in TCGA-HNSC
# (Illumina HumanMethylation450K, GDC harmonized beta values):
#
#   Panel A (OE4_FigA) — Volcano of DMPs (differentially methylated positions)
#                         in promoter CpG islands, tumor vs adjacent normal.
#                         Highlights known HNSCC suppressors: CDKN2A, CDH1,
#                         DAPK1, RASSF1A.
#   Panel B (OE4_FigB) — DNMT1 overexpression correlates with promoter
#                         hypermethylation and transcriptional silencing of
#                         tumor suppressors (β-value vs RNA, stratified by DNMT1).
#   Panel C (OE4_FigC) — Methylation burden score (mean β top DMPs) predicts
#                         Overall Survival: KM curves, High vs Low (median-split).
#
# Outputs:
#   results/figures/pub/main/OE4_FigA_dmp_volcano.{pdf,png}
#   results/figures/pub/main/OE4_FigB_dnmt1_meth_expr.{pdf,png}
#   results/figures/pub/main/OE4_FigC_survival_burden.{pdf,png}
#   results/tables/pub/main/OE4_Tab1_dmps_promoter.tsv
#   results/tables/pub/main/OE4_Tab2_meth_expr_corr.tsv
#   results/tables/pub/supp/OE4_TabS_survival_burden.tsv
#
# Data cache: data/intermediate/tcga/
#   tcga_hnsc_meth450_se.rds — β-values GDC armonizados (primera descarga ~200 MB)
#   tcga_hnsc_clinical.rds   — reutilizado de script 16 (si existe)
#   tcga_hnsc_se.rds         — SE RNA reutilizado para DNMT1 y expresión
#
# Dependencias adicionales (instalar una vez en omics-R):
#   BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
#
# Ambiente: omics-R
# Ejecución:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/19_methylation_tcga.R
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

cat("=== 19_methylation_tcga.R ===\n")
cat("Working directory:", getwd(), "\n")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ── Directorios ───────────────────────────────────────────────────────────────
dir.create("logs",                          showWarnings = FALSE)
dir.create("data/intermediate/tcga",       recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/pub/main",     recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/pub/main",      recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/pub/supp",      recursive = TRUE, showWarnings = FALSE)

log_file <- file.path("logs", paste0("19_methylation_", timestamp, ".log"))
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
check_pkg("TCGAbiolinks",                                    bioc = TRUE)
check_pkg("limma",                                           bioc = TRUE)
check_pkg("SummarizedExperiment",                            bioc = TRUE)
check_pkg("IlluminaHumanMethylation450kanno.ilmn12.hg19",   bioc = TRUE)
check_pkg("survival")
check_pkg("survminer")

library(TCGAbiolinks)
library(limma)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(survival)
library(survminer)

# =============================================================================
# ESTILOS (consistentes con scripts 16 y 17)
# =============================================================================
PRESETS <- list(
  double_col = list(w = 180, h = 120, base = 8.5, title = 10,
                    axis = 8.5, leg = 7.5, tick = 7.5, lwd = 0.7, pt = 1.8),
  single_col = list(w = 85,  h = 70,  base = 8.0, title = 9.0,
                    axis = 8.0, leg = 7.0, tick = 7.0, lwd = 0.6, pt = 1.6),
  quad_col   = list(w = 180, h = 160, base = 8.0, title = 9.0,
                    axis = 8.0, leg = 7.0, tick = 7.0, lwd = 0.6, pt = 1.6)
)

OKB <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7", "#000000")

DMP_COLS <- c(
  "Hypermethylated" = "#D55E00",
  "Hypomethylated"  = "#0072B2",
  "Not significant" = "#BBBBBB"
)
DNMT1_COLS  <- c("DNMT1 High" = "#E69F00", "DNMT1 Low" = "#56B4E9")
BURDEN_COLS <- c("High burden" = "#D55E00", "Low burden"  = "#56B4E9")

# Supresores a destacar (claim literatura: ⚠️ ver CLINICAL_CRITERIA.md)
SUPPRESSOR_GENES <- c("CDKN2A", "CDH1", "DAPK1", "RASSF1A")

# Umbrales de análisis
DELTA_BETA_THRESH <- 0.1   # Diferencia biológicamente relevante en β
ADJ_P_THRESH      <- 0.05
N_BURDEN_PROBES   <- 100   # Sondas top para construir el burden score

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
# SECCION 1: DESCARGA 450K TCGA-HNSC (con cache)
# =============================================================================
cat("\n--- SECCION 1: Descarga/cache 450K ---\n")

cache_meth <- "data/intermediate/tcga/tcga_hnsc_meth450_se.rds"
cache_cli  <- "data/intermediate/tcga/tcga_hnsc_clinical.rds"
cache_rna  <- "data/intermediate/tcga/tcga_hnsc_se.rds"

if (file.exists(cache_meth)) {
  cat("Cargando 450K desde cache...\n")
  se_meth <- readRDS(cache_meth)
} else {
  cat("Descargando TCGA-HNSC 450K (primera ejecucion, ~200 MB)...\n")

  query_meth <- GDCquery(
    project       = "TCGA-HNSC",
    data.category = "DNA Methylation",
    data.type     = "Methylation Beta Value",
    platform      = "Illumina Human Methylation 450",
    sample.type   = c("Primary Tumor", "Solid Tissue Normal")
  )
  GDCdownload(query_meth, method = "api", files.per.chunk = 50,
              directory = "data/intermediate/tcga/GDCdata")
  se_meth <- GDCprepare(query_meth,
                        directory = "data/intermediate/tcga/GDCdata")
  saveRDS(se_meth, cache_meth)
  cat("450K guardado en cache.\n")
}

# Clinical data (reutiliza cache de script 16 si existe)
if (file.exists(cache_cli)) {
  clin <- readRDS(cache_cli)
} else {
  clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
  saveRDS(clin, cache_cli)
}

cat("Dimensiones SE metilacion:", paste(dim(se_meth), collapse = " x "), "\n")
cat("Muestras totales:", ncol(se_meth), "\n")

# =============================================================================
# SECCION 2: ANOTACION + FILTRO A PROMOTORES EN ISLAS CpG
# =============================================================================
cat("\n--- SECCION 2: Anotacion y filtro de sondas ---\n")

anno <- as.data.frame(
  getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
)
cat("Sondas en anotacion 450K:", nrow(anno), "\n")

# Nombre de columna CpG island (varía entre versiones del paquete)
cpg_col <- grep("Relation_to.*Island", colnames(anno), value = TRUE)[1]
cat("Columna CpG Island detectada:", cpg_col, "\n")

# Definir promotor: TSS1500 o TSS200 en isla CpG, con gen asignado
is_promoter <- grepl("TSS1500|TSS200", anno$UCSC_RefGene_Group) &
               anno[[cpg_col]] == "Island" &
               nchar(as.character(anno$UCSC_RefGene_Name)) > 0

anno_prom <- anno[is_promoter, ]
cat("Sondas promotoras en islas CpG:", nrow(anno_prom), "\n")

# Beta matrix desde el SE y filtro a sondas promotoras
beta_all <- assay(se_meth)
cat("Dimensiones beta completa:", paste(dim(beta_all), collapse = " x "), "\n")

common_probes <- intersect(rownames(beta_all), rownames(anno_prom))
cat("Sondas promotoras presentes en datos:", length(common_probes), "\n")

beta_prom <- beta_all[common_probes, , drop = FALSE]
anno_prom <- anno_prom[common_probes, ]

rm(beta_all); gc()

# Gene symbol por sonda (primer gen si hay múltiples separados por ";")
probe_gene1 <- sapply(
  strsplit(as.character(anno_prom$UCSC_RefGene_Name), ";"),
  function(x) x[1]
)
probe_genes <- data.frame(
  probe_id    = rownames(anno_prom),
  gene_symbol = probe_gene1,
  stringsAsFactors = FALSE
)

# =============================================================================
# SECCION 3: DMPs TUMOR VS NORMAL (limma sobre M-values)
# =============================================================================
cat("\n--- SECCION 3: DMPs tumor vs normal ---\n")

sample_type <- substr(colnames(beta_prom), 14, 15)
is_tumor    <- sample_type == "01"
is_normal   <- sample_type == "11"

beta_tumor  <- beta_prom[, is_tumor,  drop = FALSE]
beta_normal <- beta_prom[, is_normal, drop = FALSE]

cat("Muestras tumorales:", ncol(beta_tumor), "\n")
cat("Muestras normales:",  ncol(beta_normal), "\n")

# Filtrar sondas con >20% NA en cualquier grupo
keep_probes <- rowMeans(is.na(beta_tumor))  <= 0.20 &
               rowMeans(is.na(beta_normal)) <= 0.20

beta_tumor  <- beta_tumor[keep_probes, , drop = FALSE]
beta_normal <- beta_normal[keep_probes, , drop = FALSE]
anno_prom   <- anno_prom[keep_probes, ]
probe_genes <- probe_genes[keep_probes, ]

cat("Sondas tras filtro NA:", nrow(beta_tumor), "\n")

# Beta → M-values (logit), clampeado a [-10, 10] para estabilidad
b2m <- function(b) {
  m <- log2(pmax(b, 1e-6) / pmax(1 - b, 1e-6))
  m[m >  10] <-  10
  m[m < -10] <- -10
  m
}
mval_tumor  <- b2m(beta_tumor)
mval_normal <- b2m(beta_normal)

# limma: modelo lineal tumor vs normal
mval_all <- cbind(mval_tumor, mval_normal)
cond     <- factor(c(rep("Tumor",  ncol(mval_tumor)),
                     rep("Normal", ncol(mval_normal))),
                   levels = c("Normal", "Tumor"))
design   <- model.matrix(~ cond)

cat("Ajustando modelo limma...\n")
fit <- lmFit(mval_all, design)
fit <- eBayes(fit)

# Delta-beta (diferencia de medias de β, Tumor - Normal)
delta_beta_vec <- rowMeans(beta_tumor, na.rm = TRUE) -
                  rowMeans(beta_normal, na.rm = TRUE)
names(delta_beta_vec) <- rownames(beta_tumor)

top_dmp <- topTable(fit, coef = 2, number = nrow(mval_all), sort.by = "none") |>
  rownames_to_column("probe_id") |>
  mutate(delta_beta = delta_beta_vec[probe_id]) |>
  left_join(probe_genes, by = "probe_id") |>
  filter(!is.na(gene_symbol), gene_symbol != "") |>
  select(probe_id, gene_symbol, logFC, P.Value, adj.P.Val, delta_beta) |>
  arrange(adj.P.Val)

n_hyper <- sum(top_dmp$adj.P.Val < ADJ_P_THRESH &
               top_dmp$delta_beta > DELTA_BETA_THRESH, na.rm = TRUE)
n_hypo  <- sum(top_dmp$adj.P.Val < ADJ_P_THRESH &
               top_dmp$delta_beta < -DELTA_BETA_THRESH, na.rm = TRUE)
cat(sprintf("DMPs hipermetilados (adj.P<%.2f, Δβ>%.1f): %d\n",
            ADJ_P_THRESH, DELTA_BETA_THRESH, n_hyper))
cat(sprintf("DMPs hipometilados  (adj.P<%.2f, Δβ<-%.1f): %d\n",
            ADJ_P_THRESH, DELTA_BETA_THRESH, n_hypo))

# Tabla principal: DMPs significativos
dmp_sig <- top_dmp |>
  filter(adj.P.Val < ADJ_P_THRESH, abs(delta_beta) > DELTA_BETA_THRESH) |>
  arrange(desc(delta_beta))

write_tsv(dmp_sig, "results/tables/pub/main/OE4_Tab1_dmps_promoter.tsv")
cat("Tabla DMPs guardada:", nrow(dmp_sig), "sondas.\n")

# ── FigA: Volcán de DMPs ──────────────────────────────────────────────────────
cat("Generando OE4_FigA (volcan DMPs)...\n")

p_dmp <- PRESETS$double_col

top_dmp_plot <- top_dmp |>
  mutate(
    sig_class = case_when(
      adj.P.Val < ADJ_P_THRESH & delta_beta >  DELTA_BETA_THRESH ~ "Hypermethylated",
      adj.P.Val < ADJ_P_THRESH & delta_beta < -DELTA_BETA_THRESH ~ "Hypomethylated",
      TRUE                                                        ~ "Not significant"
    ),
    is_suppressor = gene_symbol %in% SUPPRESSOR_GENES,
    neg_log10_p   = pmin(-log10(pmax(adj.P.Val, 1e-50)), 50)
  )

# Sonda más significativa por supresor (hipermetilada)
label_df <- top_dmp_plot |>
  filter(is_suppressor, delta_beta > 0) |>
  group_by(gene_symbol) |>
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) |>
  ungroup()
cat("Supresores etiquetados en FigA:", paste(label_df$gene_symbol, collapse = ", "), "\n")

figA <- ggplot(top_dmp_plot,
               aes(x = delta_beta, y = neg_log10_p)) +

  geom_hline(yintercept = -log10(ADJ_P_THRESH),
             linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_vline(xintercept =  DELTA_BETA_THRESH,
             linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = -DELTA_BETA_THRESH,
             linetype = "dashed", color = "grey50", linewidth = 0.35) +

  geom_point(data = filter(top_dmp_plot, !is_suppressor),
             aes(color = sig_class), size = 0.5, alpha = 0.35, shape = 16) +

  geom_point(data = filter(top_dmp_plot, is_suppressor),
             aes(fill = sig_class), size = 2.5, shape = 21,
             color = "white", stroke = 0.5) +

  geom_label_repel(
    data          = label_df,
    aes(label = gene_symbol, fill = sig_class),
    color         = "white",
    size          = p_dmp$tick / .pt,
    fontface      = "bold",
    box.padding   = 0.5,
    label.size    = 0,
    label.padding = unit(0.13, "lines"),
    show.legend   = FALSE
  ) +

  scale_color_manual(values = DMP_COLS, name = NULL) +
  scale_fill_manual(values  = DMP_COLS, name = NULL) +

  labs(
    x     = expression(Delta*beta[mean]~"(Tumor - Normal)"),
    y     = expression(-log[10]~"(adj. p-value)"),
    title = "Promoter CpG island methylation: DMPs (tumor vs normal)"
  ) +
  base_theme(p_dmp) +
  theme(legend.position = "top")

save_fig(figA, "OE4_FigA_dmp_volcano", preset = "double_col")

# =============================================================================
# SECCION 4: DNMT1 OVEREXPRESSION vs SILENCIAMIENTO EPIGENETICO
# =============================================================================
cat("\n--- SECCION 4: DNMT1 vs silenciamiento epigenético ---\n")

run_sec4 <- file.exists(cache_rna)
if (!run_sec4) {
  cat("ATENCION: cache RNA no encontrado (", cache_rna, ").\n")
  cat("Ejecutar script 16 primero para generar tcga_hnsc_se.rds.\n")
  cat("Omitiendo seccion 4.\n")
}

if (run_sec4) {
  se_rna       <- readRDS(cache_rna)
  se_rna_tumor <- se_rna[, substr(colnames(se_rna), 14, 15) == "01"]
  expr_mat     <- log2(assay(se_rna_tumor, "fpkm_unstrand") + 0.1)
  rownames(expr_mat) <- rowData(se_rna_tumor)$gene_name
  rm(se_rna); gc()
  cat("Expresion RNA cargada:", nrow(expr_mat), "genes,",
      ncol(expr_mat), "muestras tumorales.\n")

  # IDs de paciente (primeros 12 chars del barcode)
  rna_pts    <- substr(colnames(se_rna_tumor), 1, 12)
  meth_pts   <- substr(colnames(beta_tumor),   1, 12)
  common_pts <- intersect(rna_pts, meth_pts)
  cat("Pacientes comunes meth+RNA:", length(common_pts), "\n")

  rna_idx  <- match(common_pts, rna_pts)
  meth_idx <- match(common_pts, meth_pts)

  if (!"DNMT1" %in% rownames(expr_mat)) {
    cat("ATENCION: DNMT1 no encontrado en expr_mat. Omitiendo seccion 4.\n")
    run_sec4 <- FALSE
  }
}

if (run_sec4) {
  # Estratificar por DNMT1 (median-split en pacientes comunes)
  dnmt1_expr  <- expr_mat["DNMT1", rna_idx]
  dnmt1_med   <- median(dnmt1_expr, na.rm = TRUE)
  dnmt1_group <- factor(
    ifelse(dnmt1_expr >= dnmt1_med, "DNMT1 High", "DNMT1 Low"),
    levels = c("DNMT1 High", "DNMT1 Low")
  )
  cat(sprintf("DNMT1 mediana: %.3f | High: %d | Low: %d\n",
              dnmt1_med, sum(dnmt1_group=="DNMT1 High"), sum(dnmt1_group=="DNMT1 Low")))

  # Media de beta por paciente para sondas promotoras de un gen
  get_promoter_beta <- function(gene) {
    gene_mask   <- grepl(paste0("(^|;)", gene, "(;|$)"),
                         as.character(anno_prom$UCSC_RefGene_Name))
    gene_probes <- rownames(anno_prom)[gene_mask]
    gene_probes <- intersect(gene_probes, rownames(beta_tumor))
    if (length(gene_probes) == 0) return(NULL)
    colMeans(beta_tumor[gene_probes, meth_idx, drop = FALSE], na.rm = TRUE)
  }

  corr_records <- list()
  scatter_data <- list()

  for (gene in SUPPRESSOR_GENES) {
    if (!gene %in% rownames(expr_mat)) {
      cat("Gen no en RNA:", gene, "\n"); next
    }
    beta_gene <- get_promoter_beta(gene)
    if (is.null(beta_gene)) {
      cat("Sin sondas promotoras para:", gene, "\n"); next
    }

    rna_gene <- expr_mat[gene, rna_idx]

    df_gene <- data.frame(
      patient_id  = common_pts,
      beta        = as.numeric(beta_gene),
      expr_rna    = as.numeric(rna_gene),
      dnmt1_group = dnmt1_group,
      gene        = gene,
      stringsAsFactors = FALSE
    ) |> filter(!is.na(beta), !is.na(expr_rna))

    if (nrow(df_gene) < 20) {
      cat("Pocos datos para correlacion:", gene, "\n"); next
    }

    # Spearman global y por estrato
    rho_all <- cor(df_gene$beta, df_gene$expr_rna, method = "spearman")
    p_all   <- suppressWarnings(
      cor.test(df_gene$beta, df_gene$expr_rna, method = "spearman")$p.value)

    strata_corr <- df_gene |>
      group_by(dnmt1_group) |>
      summarise(
        rho = cor(beta, expr_rna, method = "spearman"),
        p   = suppressWarnings(
                cor.test(beta, expr_rna, method = "spearman")$p.value),
        n   = n(),
        .groups = "drop"
      )

    rho_h <- strata_corr |> filter(dnmt1_group == "DNMT1 High") |> pull(rho)
    p_h   <- strata_corr |> filter(dnmt1_group == "DNMT1 High") |> pull(p)
    rho_l <- strata_corr |> filter(dnmt1_group == "DNMT1 Low")  |> pull(rho)
    p_l   <- strata_corr |> filter(dnmt1_group == "DNMT1 Low")  |> pull(p)

    corr_records[[gene]] <- tibble(
      gene       = gene,
      n_patients = nrow(df_gene),
      rho_all    = round(rho_all, 3),
      p_all      = signif(p_all, 3),
      rho_high   = if (length(rho_h) > 0) round(rho_h, 3) else NA_real_,
      p_high     = if (length(p_h)   > 0) signif(p_h,  3) else NA_real_,
      rho_low    = if (length(rho_l) > 0) round(rho_l, 3) else NA_real_,
      p_low      = if (length(p_l)   > 0) signif(p_l,  3) else NA_real_
    )
    scatter_data[[gene]] <- df_gene
    cat(sprintf("  %s: rho_all=%.3f, rho_high=%.3f, rho_low=%.3f\n",
                gene,
                rho_all,
                if (length(rho_h) > 0) rho_h else NA,
                if (length(rho_l) > 0) rho_l else NA))
  }

  if (length(corr_records) > 0) {
    corr_tbl <- bind_rows(corr_records)
    write_tsv(corr_tbl, "results/tables/pub/main/OE4_Tab2_meth_expr_corr.tsv")
    cat("Tabla correlaciones guardada.\n")
  }

  # ── FigB: Scatter facetado por supresor, estratificado por DNMT1 ───────────
  if (length(scatter_data) >= 1) {
    cat("Generando OE4_FigB (DNMT1 correlation)...\n")

    p_q <- PRESETS$quad_col
    scatter_all <- bind_rows(scatter_data) |>
      mutate(gene = factor(gene, levels = SUPPRESSOR_GENES))

    # Posiciones de anotación ρ por panel (top-left de cada faceta)
    rho_annot <- bind_rows(lapply(names(corr_records), function(g) {
      cr <- corr_records[[g]]
      sd <- scatter_data[[g]]
      x0 <- min(sd$beta, na.rm = TRUE)
      y1 <- max(sd$expr_rna, na.rm = TRUE)
      yr <- diff(range(sd$expr_rna, na.rm = TRUE))
      bind_rows(
        tibble(gene = g, dnmt1_group = "DNMT1 High",
               x = x0, y = y1,
               label = sprintf("ρ=%.2f (high)", cr$rho_high)),
        tibble(gene = g, dnmt1_group = "DNMT1 Low",
               x = x0, y = y1 - yr * 0.13,
               label = sprintf("ρ=%.2f (low)",  cr$rho_low))
      )
    })) |>
      mutate(
        gene        = factor(gene, levels = SUPPRESSOR_GENES),
        dnmt1_group = factor(dnmt1_group, levels = c("DNMT1 High", "DNMT1 Low"))
      )

    figB <- ggplot(scatter_all,
                   aes(x = beta, y = expr_rna, color = dnmt1_group)) +
      geom_point(size = 0.7, alpha = 0.5, shape = 16) +
      geom_smooth(aes(fill = dnmt1_group), method = "lm", se = TRUE,
                  linewidth = 0.5, alpha = 0.12, show.legend = FALSE) +
      geom_text(data = rho_annot,
                aes(x = x, y = y, label = label, color = dnmt1_group),
                hjust = 0, size = p_q$tick / .pt,
                fontface = "italic", inherit.aes = FALSE) +
      facet_wrap(~ gene, ncol = 2, scales = "free") +
      scale_color_manual(values = DNMT1_COLS, name = "DNMT1 expression") +
      scale_fill_manual(values  = DNMT1_COLS, name = "DNMT1 expression") +
      labs(
        x     = expression("Promoter mean "*beta*"-value"),
        y     = expression("log"[2]*" FPKM (RNA-seq)"),
        title = "Promoter hypermethylation and transcriptional silencing"
      ) +
      base_theme(p_q) +
      theme(
        legend.position = "top",
        strip.text      = element_text(size = p_q$axis, face = "bold.italic")
      )

    save_fig(figB, "OE4_FigB_dnmt1_meth_expr", preset = "quad_col")
  }
}

# =============================================================================
# SECCION 5: METHYLATION BURDEN SCORE ~ OVERALL SURVIVAL
# =============================================================================
cat("\n--- SECCION 5: Burden score ~ OS ---\n")

# Top DMPs hipermetilados para el burden score
top_hyper <- top_dmp |>
  filter(adj.P.Val < ADJ_P_THRESH, delta_beta > DELTA_BETA_THRESH) |>
  head(N_BURDEN_PROBES)

burden_probes <- intersect(top_hyper$probe_id, rownames(beta_tumor))
cat("Sondas para burden score:", length(burden_probes), "\n")

run_sec5 <- length(burden_probes) >= 10
if (!run_sec5)
  cat("ATENCION: muy pocas sondas para burden score. Omitiendo seccion 5.\n")

if (run_sec5) {
  burden_mat   <- beta_tumor[burden_probes, , drop = FALSE]
  burden_score <- colMeans(burden_mat, na.rm = TRUE)

  burden_df <- data.frame(
    patient_id   = substr(colnames(beta_tumor), 1, 12),
    burden_score = burden_score,
    stringsAsFactors = FALSE
  ) |> distinct(patient_id, .keep_all = TRUE)

  # Supervivencia desde colData del SE metilación (tumor)
  se_meth_tumor <- se_meth[, is_tumor]
  col_meth      <- as.data.frame(colData(se_meth_tumor))

  has_surv_cols <- all(c("patient", "vital_status") %in% colnames(col_meth))

  if (has_surv_cols) {
    surv_df <- col_meth |>
      mutate(
        patient_id = patient,
        os_event   = as.integer(tolower(vital_status) == "dead"),
        os_time    = if_else(os_event == 1,
                             suppressWarnings(as.numeric(days_to_death)),
                             suppressWarnings(as.numeric(days_to_last_follow_up)))
      ) |>
      filter(!is.na(os_time), os_time > 0) |>
      select(patient_id, os_event, os_time) |>
      distinct(patient_id, .keep_all = TRUE)
  } else {
    # Fallback: clinical RDS de script 16
    cat("Columnas de supervivencia no en colData meth SE. Usando clinical RDS...\n")
    pt_col <- if ("submitter_id"        %in% colnames(clin)) "submitter_id"
              else if ("bcr_patient_barcode" %in% colnames(clin)) "bcr_patient_barcode"
              else colnames(clin)[1]
    cat("Columna paciente en clinical:", pt_col, "\n")

    surv_df <- data.frame(
      patient_id = substr(colnames(se_meth_tumor), 1, 12),
      stringsAsFactors = FALSE
    ) |>
      distinct(patient_id) |>
      left_join(
        clin |> select(all_of(c(pt_col, "vital_status",
                                "days_to_death", "days_to_last_follow_up"))),
        by = setNames(pt_col, "patient_id")
      ) |>
      mutate(
        os_event = as.integer(tolower(vital_status) == "dead"),
        os_time  = if_else(os_event == 1,
                           suppressWarnings(as.numeric(days_to_death)),
                           suppressWarnings(as.numeric(days_to_last_follow_up)))
      ) |>
      filter(!is.na(os_time), os_time > 0) |>
      select(patient_id, os_event, os_time) |>
      distinct(patient_id, .keep_all = TRUE)
  }

  cat("Pacientes con datos de supervivencia:", nrow(surv_df), "\n")

  # Unir burden + OS y estratificar por mediana
  surv_burden <- burden_df |>
    inner_join(surv_df, by = "patient_id") |>
    mutate(
      burden_group = factor(
        if_else(burden_score >= median(burden_score, na.rm = TRUE),
                "High burden", "Low burden"),
        levels = c("High burden", "Low burden")
      )
    )

  cat("Pacientes con burden+OS:", nrow(surv_burden), "\n")
  write_tsv(surv_burden, "results/tables/pub/supp/OE4_TabS_survival_burden.tsv")
  cat("Tabla burden/supervivencia guardada.\n")

  # ── FigC: KM curves (patrón idéntico a OE3_FigB de script 16) ─────────────
  if (nrow(surv_burden) >= 20) {
    cat("Generando OE4_FigC (survival burden KM)...\n")

    fit_burden <- survfit(Surv(os_time, os_event) ~ burden_group,
                          data = surv_burden)

    lrt_burden <- tryCatch(
      survdiff(Surv(os_time, os_event) ~ burden_group, data = surv_burden),
      error = function(e) NULL
    )
    pval_txt <- if (!is.null(lrt_burden)) {
      pv <- 1 - pchisq(lrt_burden$chisq, df = length(lrt_burden$n) - 1)
      if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)
    } else ""
    cat("Survival log-rank:", pval_txt, "\n")

    p_s <- PRESETS$single_col

    figC <- ggsurvplot(
      fit_burden,
      data         = surv_burden,
      palette      = BURDEN_COLS,
      size         = 0.6,
      conf.int     = FALSE,
      pval         = FALSE,
      risk.table   = FALSE,
      xscale       = "d_y",
      xlab         = "Time (years)",
      ylab         = "Overall survival",
      title        = "Methylation burden score",
      legend.title = "Burden",
      legend.labs  = c("High", "Low"),
      fontsize     = p_s$tick / .pt,
      ggtheme      = base_theme(p_s) +
                       theme(legend.position = c(0.78, 0.85))
    )$plot +
    annotate("text", x = Inf, y = 0.05, label = pval_txt,
             hjust = 1.05, vjust = 0,
             size = p_s$tick / .pt, color = "grey20")

    save_fig(figC, "OE4_FigC_survival_burden", preset = "single_col")
  } else {
    cat("ATENCION: n <20, no se genera FigC.\n")
  }
}

# =============================================================================
# SECCION 6: RESUMEN
# =============================================================================
cat("\n========================================\n")
cat("RESUMEN 19_methylation_tcga.R\n")
cat("========================================\n")
cat(sprintf("Sondas promotoras analizadas:     %d\n", nrow(beta_tumor)))
cat(sprintf("DMPs hipermetilados significativos: %d\n", n_hyper))
cat(sprintf("DMPs hipometilados significativos:  %d\n", n_hypo))
cat("\nOutputs:\n")
cat("  results/figures/pub/main/OE4_FigA_dmp_volcano.{pdf,png}\n")
if (run_sec4 && length(scatter_data) >= 1)
  cat("  results/figures/pub/main/OE4_FigB_dnmt1_meth_expr.{pdf,png}\n")
if (run_sec5 && nrow(surv_burden) >= 20)
  cat("  results/figures/pub/main/OE4_FigC_survival_burden.{pdf,png}\n")
cat("  results/tables/pub/main/OE4_Tab1_dmps_promoter.tsv\n")
if (exists("corr_tbl"))
  cat("  results/tables/pub/main/OE4_Tab2_meth_expr_corr.tsv\n")
if (run_sec5)
  cat("  results/tables/pub/supp/OE4_TabS_survival_burden.tsv\n")
cat("\nFin:", format(Sys.time()), "\n")
sink()
