#!/usr/bin/env Rscript
# =============================================================================
# 17_pub_figures.R
# =============================================================================
# Figuras de calidad publicación para el proyecto HNSCC Drug Repurposing.
# Aplica tema centralizado (pub-figures double_col / single_col), paleta
# Okabe-Ito (colorblind-safe) y exporta PDF + PNG 300 DPI.
#
# Nota sobre nomenclatura: TVsS = Tumor vs. Sano (tejido adyacente del mismo
# paciente). En inglés: Tumor vs. Adjacent Normal.
#
# Figuras generadas (results/figures/pub/):
#   SECCIÓN A — Proteómica QC/DE
#     A1_volcano          — Volcano plot mejorado (top labels ggrepel)
#     A2_MA_plot          — MA plot (logFC vs intensidad media)
#     A3_PCA              — PCA biplot con líneas de pares
#     A4_heatmap_topDE    — Heatmap top 30 proteínas DE × 20 muestras
#   SECCIÓN B — Enriquecimiento funcional
#     B1_hallmarks_barplot — Hallmarks GSEA horizontal por NES
#   SECCIÓN C — Bases de datos de fármacos
#     C1_drug_sources_bar  — Nº candidatos por fuente
#     C2_drug_phase_dist   — Distribución fases clínicas
#   SECCIÓN D — Red PPI + Scoring
#     D1_degree_vs_logFC   — Scatter grado PPI vs logFC
#     D2_scoring_components_dot — Dot matrix 6 componentes × Top 20
#     D3_top20_lollipop    — Lollipop ranking final mejorado
#   SECCIÓN E — Evidencia final
#     E1_evidence_heatmap  — Heatmap binario de evidencia (ComplexHeatmap)
#   SECCIÓN F — Análisis de sensibilidad
#     F1_bump_chart        — Bump chart estabilidad de ranks
#     F2_stability_bar     — Barplot de robustez mejorado
#   FIGURAS MULTIPANEL (para manuscrito)
#     FIG1_QC_DE          — Volcano + MA
#     FIG2_pathways       — Hallmarks
#     FIG3_scoring        — Lollipop + Dot matrix
#     FIG4_network        — Scatter logFC vs grado
#     FIG5_sensitivity    — Bump chart + Stability bar
#
# Ambiente: omics-R
# Ejecución:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/17_pub_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(patchwork)
  library(ggrepel)
  library(stringr)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
})

# ── Detectar directorio del proyecto ─────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  proj_dir    <- dirname(dirname(script_path))
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml")))
    setwd(proj_dir)
}

cat("=== 17_pub_figures.R ===\n")
cat("Working directory:", getwd(), "\n")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ── Directorio de salida ──────────────────────────────────────────────────────
pub_dir <- "results/figures/pub"
dir.create(pub_dir, recursive = TRUE, showWarnings = FALSE)

# ── Log ───────────────────────────────────────────────────────────────────────
log_dir  <- "logs"
dir.create(log_dir, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("17_pub_figures_", timestamp, ".log"))
sink(log_file, split = TRUE)
cat("Inicio:", format(Sys.time()), "\n\n")

# =============================================================================
# SISTEMA DE ESTILOS CENTRALIZADO
# =============================================================================

# ── Presets (dimensiones en mm) ───────────────────────────────────────────────
PRESETS <- list(
  double_col = list(w = 180, h = 120, base = 8.5, title = 10,
                    axis = 8.5, leg = 7.5, tick = 7.5, lwd = 0.7, pt = 1.8),
  single_col = list(w = 85,  h = 65,  base = 8.0, title = 9.0,
                    axis = 8.0, leg = 7.0, tick = 7.0, lwd = 0.6, pt = 1.6)
)

# ── Paleta Okabe-Ito (colorblind-safe) ────────────────────────────────────────
OKB       <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")
DE_COLS   <- c(up = "#D55E00", down = "#0072B2", ns = "#CCCCCC")
COND_COLS <- c(Tumor = "#D55E00", "Adjacent Normal" = "#0072B2")
PHASE_COLS <- c(
  "Phase IV (Approved)" = "#009E73", "Phase III" = "#56B4E9",
  "Phase II" = "#E69F00", "Phase I" = "#CC79A7", "Unknown" = "#CCCCCC"
)

# Etiquetas legibles de clase de fármaco (script 10 exporta drug_class = A/B/C/D)
DRUG_CLASS_LABELS <- c(
  A = "HNSCC-approved",
  B = "Other cancer",
  C = "Non-oncology",
  D = "Experimental"
)

# ── Tema centralizado ─────────────────────────────────────────────────────────
theme_pub <- function(preset = "double_col") {
  p <- PRESETS[[preset]]
  theme_classic(base_size = p$base, base_family = "sans") +
    theme(
      plot.title    = element_text(size = p$title, face = "bold", hjust = 0,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(size = p$base, color = "grey40"),
      axis.title    = element_text(size = p$axis),
      axis.text     = element_text(size = p$tick, color = "black"),
      legend.title  = element_text(size = p$leg, face = "bold"),
      legend.text   = element_text(size = p$leg),
      legend.key.size    = unit(3, "mm"),
      legend.background  = element_blank(),
      strip.text         = element_text(size = p$base, face = "bold"),
      strip.background   = element_blank(),
      axis.line          = element_line(linewidth = p$lwd * 0.35, color = "black"),
      axis.ticks         = element_line(linewidth = p$lwd * 0.35, color = "black"),
      axis.ticks.length  = unit(1.5, "mm"),
      plot.margin        = margin(4, 6, 4, 6, "mm")
    )
}

# ── Funciones de guardado ─────────────────────────────────────────────────────
save_pub <- function(p, name, preset = "double_col", w_add = 0, h_add = 0) {
  d <- PRESETS[[preset]]
  w <- d$w + w_add
  h <- d$h + h_add
  pdf_path <- file.path(pub_dir, paste0(name, ".pdf"))
  png_path <- file.path(pub_dir, paste0(name, ".png"))
  ggsave(pdf_path, p, width = w, height = h, units = "mm",
         device = cairo_pdf, limitsize = FALSE)
  ggsave(png_path, p, width = w, height = h, units = "mm",
         dpi = 300, limitsize = FALSE)
  cat(sprintf("  Saved: %s (%g x %g mm)\n", basename(pdf_path), w, h))
  invisible(pdf_path)
}

# Para ComplexHeatmap (no usa ggsave)
save_ch <- function(ht, name, preset = "double_col", w_add = 0, h_add = 0) {
  d <- PRESETS[[preset]]
  w_in <- (d$w + w_add) / 25.4
  h_in <- (d$h + h_add) / 25.4
  pdf_path <- file.path(pub_dir, paste0(name, ".pdf"))
  png_path <- file.path(pub_dir, paste0(name, ".png"))
  pdf(pdf_path, width = w_in, height = h_in)
  draw(ht)
  dev.off()
  png(png_path,
      width  = (d$w + w_add) * 300 / 25.4,
      height = (d$h + h_add) * 300 / 25.4,
      res    = 300)
  draw(ht)
  dev.off()
  cat(sprintf("  Saved: %s (%.1f x %.1f in)\n", basename(pdf_path), w_in, h_in))
  invisible(pdf_path)
}

# ── Helper: asignar categoría a columnas de evidencia por patrón ──────────────
# Robusto ante cambios de idioma/formato en nombres de columna de script 13
assign_ev_category <- function(col_names) {
  sapply(col_names, function(cn) {
    cl <- tolower(iconv(cn, to = "ASCII//TRANSLIT"))
    if      (grepl("proteom|logfc|de_sig|significant|expr", cl)) "Proteomics"
    else if (grepl("multi.*source|n_source|cmap|pathway|fuente|reversor", cl)) "Databases"
    else if (grepl("hub|network|ppi|degree|betweenness", cl))    "Network"
    else if (grepl("phase|trial|clinic|fase|hnscc", cl))         "Clinical"
    else if (grepl("pubmed|paper|literature|literat", cl))       "Literature"
    else if (grepl("driver|cancer|cosmic|oncogen", cl))          "Genomics"
    else                                                          "Other"
  })
}

# =============================================================================
# CARGAR DATOS
# =============================================================================
cat("\n-- Cargando datos...\n")

de        <- read_tsv("results/tables/de_limma/01_TVsS_all_proteins.tsv",
                      show_col_types = FALSE)
hallmarks <- read_tsv("results/tables/pathway_enrichment/03_Hallmarks_GSEA.tsv",
                      show_col_types = FALSE)
drug_sum  <- read_tsv("results/tables/drug_targets/08_drug_summary_per_drug.tsv",
                      show_col_types = FALSE)
net_nodes <- read_tsv("results/tables/network/09_network_node_metrics.tsv",
                      show_col_types = FALSE)
top20     <- read_tsv("results/tables/10_top20_candidates.tsv",
                      show_col_types = FALSE)
evidence  <- read_tsv("results/tables/13_evidence_matrix.tsv",
                      show_col_types = FALSE)
sens      <- read_tsv("results/tables/15_sensitivity_ranks.tsv",
                      show_col_types = FALSE)

cat("  Datos cargados OK\n")

# ── Columna de intensidad media para MA plot ──────────────────────────────────
# limma estándar exporta "AveExpr"; script 01 puede haberla renombrado
avg_col <- if ("avg_intensity_TVsS" %in% colnames(de)) {
  "avg_intensity_TVsS"
} else if ("AveExpr" %in% colnames(de)) {
  cat("  NOTE: usando 'AveExpr' para MA plot (avg_intensity_TVsS no encontrada)\n")
  "AveExpr"
} else {
  stop("Columna de intensidad media no encontrada (avg_intensity_TVsS o AveExpr)")
}

# ── Derivar dirección DE ──────────────────────────────────────────────────────
de <- de %>%
  mutate(
    direction = case_when(
      sig.FDR_TVsS ==  1 ~ "up",
      sig.FDR_TVsS == -1 ~ "down",
      TRUE               ~ "ns"
    ),
    direction = factor(direction, levels = c("up", "down", "ns"))
  )
n_up   <- sum(de$direction == "up")
n_down <- sum(de$direction == "down")
n_ns   <- sum(de$direction == "ns")
cat(sprintf("  DE: %d up, %d down, %d NS\n", n_up, n_down, n_ns))

# ── Unir columnas DE a net_nodes si no están presentes ───────────────────────
# Script 09 exporta métricas de red; logFC y adj.P.Val vienen del TSV de script 01
de_cols_needed <- c("logFC_TVsS", "adj.P.Val_TVsS")
if (!all(de_cols_needed %in% colnames(net_nodes))) {
  cat("  NOTE: uniendo logFC_TVsS/adj.P.Val_TVsS a net_nodes desde tabla DE\n")
  net_nodes <- net_nodes %>%
    left_join(
      de %>% select(gene_symbol, logFC_TVsS, adj.P.Val_TVsS),
      by = "gene_symbol"
    )
}

# ── Matriz de expresión por muestra ──────────────────────────────────────────
sample_cols <- c("M1S","M1T","M2S","M2T","M3S","M3T","M4S","M4T","M5S","M5T",
                 "M6S","M6T","M7S","M7T","M8S","M8T","M9S","M9T","M10S","M10T")
sample_cols <- intersect(sample_cols, colnames(de))  # solo columnas que existen

expr_mat <- de %>%
  select(gene_symbol, all_of(sample_cols)) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

# =============================================================================
# SECCIÓN A — PROTEÓMICA QC / DE
# =============================================================================
cat("\n--- Sección A: QC / DE ---\n")

# ── A1: Volcano plot ──────────────────────────────────────────────────────────
# Top 10 por dirección para garantizar etiquetas en ambos grupos (up Y down)
top_vol <- bind_rows(
  de %>% filter(direction == "up")   %>% arrange(desc(logFC_TVsS)) %>% slice_head(n = 10),
  de %>% filter(direction == "down") %>% arrange(logFC_TVsS)       %>% slice_head(n = 10)
)

p_volcano <- ggplot(de, aes(x = logFC_TVsS, y = -log10(adj.P.Val_TVsS),
                            color = direction)) +
  geom_point(aes(alpha = direction), size = PRESETS$double_col$pt, stroke = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = 0.4, color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             linewidth = 0.35, color = "grey50") +
  geom_label_repel(
    data = top_vol, aes(label = gene_symbol),
    size = 1.8, label.padding = 0.1, label.size = 0.15,
    max.overlaps = 20, segment.size = 0.25,
    min.segment.length = 0.2, show.legend = FALSE
  ) +
  scale_color_manual(
    values = DE_COLS,
    labels = c(up   = sprintf("Up (%d)",   n_up),
               down = sprintf("Down (%d)", n_down),
               ns   = sprintf("NS (%d)",   n_ns))
  ) +
  scale_alpha_manual(values = c(up = 0.85, down = 0.85, ns = 0.25),
                     guide = "none") +
  labs(title = "Differential proteome — Tumor vs. Adjacent Normal",
       x     = expression(log[2]~"Fold Change  (TVsS)"),
       y     = expression(-log[10]~"FDR"),
       color = NULL) +
  theme_pub() +
  theme(legend.position = c(0.82, 0.88))

save_pub(p_volcano, "A1_volcano")
cat("  A1: Volcano — OK\n")

# ── A2: MA plot ───────────────────────────────────────────────────────────────
# Top 8 up + 7 down para garantizar etiquetas en ambos grupos (up Y down)
top_ma <- bind_rows(
  de %>% filter(direction == "up")   %>% arrange(desc(logFC_TVsS)) %>% slice_head(n = 8),
  de %>% filter(direction == "down") %>% arrange(logFC_TVsS)       %>% slice_head(n = 7)
)

p_ma <- ggplot(de, aes(x = .data[[avg_col]], y = logFC_TVsS,
                       color = direction)) +
  geom_point(aes(alpha = direction), size = PRESETS$double_col$pt, stroke = 0) +
  geom_hline(yintercept = 0,       linetype = "solid",
             linewidth = 0.45, color = "grey40") +
  geom_hline(yintercept = c(-1,1), linetype = "dashed",
             linewidth = 0.35, color = "grey60") +
  geom_label_repel(
    data = top_ma, aes(label = gene_symbol),
    size = 1.8, label.padding = 0.1, label.size = 0.15,
    max.overlaps = 15, segment.size = 0.25, show.legend = FALSE
  ) +
  scale_color_manual(
    values = DE_COLS,
    labels = c(up   = sprintf("Up (%d)",   n_up),
               down = sprintf("Down (%d)", n_down),
               ns   = sprintf("NS (%d)",   n_ns))
  ) +
  scale_alpha_manual(values = c(up = 0.85, down = 0.85, ns = 0.2),
                     guide = "none") +
  labs(title = "MA plot — abundance bias assessment",
       x     = expression("Mean"~log[2]~"intensity (A)"),
       y     = expression(log[2]~"Fold Change  (M)"),
       color = NULL) +
  theme_pub() +
  theme(legend.position = c(0.15, 0.15))

save_pub(p_ma, "A2_MA_plot")
cat("  A2: MA plot — OK\n")

# ── A3: PCA biplot ────────────────────────────────────────────────────────────
pca_res <- prcomp(t(expr_mat), scale. = TRUE)
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("sample_id") %>%
  mutate(
    condition = ifelse(str_ends(sample_id, "T"), "Tumor", "Adjacent Normal"),
    patient   = str_remove(sample_id, "[TS]$")
  )

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_line(aes(group = patient), color = "grey70",
            linewidth = 0.4, linetype = "dotted") +
  geom_point(aes(shape = condition), size = 2.4, stroke = 0.5) +
  scale_color_manual(values = COND_COLS) +
  scale_shape_manual(values = c(Tumor = 16, "Adjacent Normal" = 1)) +
  labs(title = "PCA — paired tumor/adjacent normal samples",
       x     = sprintf("PC1  (%.1f%%)", var_exp[1]),
       y     = sprintf("PC2  (%.1f%%)", var_exp[2]),
       color = "Condition", shape = "Condition") +
  guides(color = guide_legend(override.aes = list(size = 2.5)),
         shape = guide_legend(override.aes = list(size = 2.5))) +
  theme_pub() +
  theme(legend.position = c(0.82, 0.15))

save_pub(p_pca, "A3_PCA")
cat("  A3: PCA — OK\n")

# ── A4: Heatmap top 30 proteínas DE ──────────────────────────────────────────
top_up   <- de %>% filter(direction == "up")   %>%
  slice_max(logFC_TVsS, n = 15) %>% pull(gene_symbol)
top_down <- de %>% filter(direction == "down") %>%
  slice_min(logFC_TVsS, n = 15) %>% pull(gene_symbol)
top_genes <- c(top_up, top_down)

heat_mat <- expr_mat[top_genes[top_genes %in% rownames(expr_mat)], , drop = FALSE]
heat_z   <- t(scale(t(heat_mat)))   # Z-score por proteína (gen)

# Ordenar columnas: agrupar Tumor / Adjacent Normal
col_cond  <- ifelse(str_ends(colnames(heat_mat), "T"), "Tumor", "Adjacent Normal")
col_order <- order(col_cond)
heat_z    <- heat_z[, col_order]
col_cond  <- col_cond[col_order]

ht_col_ann <- HeatmapAnnotation(
  Condition = col_cond,
  col = list(Condition = c(Tumor = "#D55E00", "Adjacent Normal" = "#0072B2")),
  annotation_name_gp = gpar(fontsize = 6.5),
  simple_anno_size   = unit(3, "mm")
)

# row_split: up genes primero, down genes después
n_up_heat   <- sum(rownames(heat_z) %in% top_up)
n_down_heat <- sum(rownames(heat_z) %in% top_down)
row_split <- factor(
  c(rep("Up-regulated", n_up_heat), rep("Down-regulated", n_down_heat)),
  levels = c("Up-regulated", "Down-regulated")
)

# Sin anotación derecha de logFC: el row_split y los colores del mapa ya
# diferencian claramente las dos direcciones
ht_de <- Heatmap(
  heat_z,
  name    = "Z-score",
  col     = colorRamp2(c(-2.5, 0, 2.5), c("#0072B2", "white", "#D55E00")),
  top_annotation    = ht_col_ann,
  row_split         = row_split,
  cluster_row_slices = FALSE,
  cluster_rows       = TRUE,
  cluster_columns    = FALSE,
  show_column_names  = TRUE,
  column_names_gp    = gpar(fontsize = 5.5),
  row_names_gp       = gpar(fontsize = 6.5),
  row_title_gp       = gpar(fontsize = 8, fontface = "bold"),
  column_title       = "Top 30 differentially expressed proteins — Tumor vs. Adjacent Normal",
  column_title_gp    = gpar(fontsize = 8, fontface = "bold"),
  rect_gp            = gpar(col = "grey90", lwd = 0.4),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5)
  )
)

save_ch(ht_de, "A4_heatmap_topDE", "double_col", h_add = 50)
cat("  A4: Heatmap top DE — OK\n")

# =============================================================================
# SECCIÓN B — ENRIQUECIMIENTO FUNCIONAL
# =============================================================================
cat("\n--- Sección B: Pathways ---\n")

# ── B1: Hallmarks GSEA — barplot horizontal ───────────────────────────────────
hall_df <- hallmarks %>%
  mutate(
    pathway  = str_remove(Description, "^HALLMARK_") %>%
               str_replace_all("_", " ") %>%
               str_to_title() %>%
               str_trunc(35),
    gsea_dir = factor(
      ifelse(NES < 0, "Downregulated in tumor", "Upregulated in tumor"),
      levels = c("Upregulated in tumor", "Downregulated in tumor")
    ),
    pval_label  = sprintf("p=%.1e", p.adjust),
    label_x     = ifelse(NES > 0, NES + 0.05, NES - 0.05),
    label_hjust = ifelse(NES > 0, 0, 1)
  ) %>%
  arrange(NES)

p_hall <- ggplot(hall_df,
                 aes(x = NES, y = reorder(pathway, NES), fill = gsea_dir)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "black") +
  geom_text(aes(x = label_x, hjust = label_hjust, label = pval_label),
            size = 1.9, color = "grey35") +
  scale_fill_manual(
    values = c("Upregulated in tumor"   = "#D55E00",
               "Downregulated in tumor" = "#0072B2")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.25))) +
  labs(title = "Hallmark gene sets — GSEA  (Tumor vs. Adjacent Normal)",
       x     = "Normalized enrichment score (NES)",
       y     = NULL,
       fill  = NULL) +
  theme_pub() +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

save_pub(p_hall, "B1_hallmarks_barplot", "double_col", h_add = 50)
cat("  B1: Hallmarks GSEA barplot — OK\n")

# =============================================================================
# SECCIÓN C — BASES DE DATOS DE FÁRMACOS
# =============================================================================
cat("\n--- Sección C: Drug databases ---\n")

# ── C1: Candidatos por fuente ─────────────────────────────────────────────────
source_long <- drug_sum %>%
  select(drug_name_norm, sources) %>%
  mutate(
    DGIdb       = str_detect(sources, "DGIdb"),
    ChEMBL      = str_detect(sources, "ChEMBL"),
    OpenTargets = str_detect(sources, "OpenTargets")
  ) %>%
  pivot_longer(c(DGIdb, ChEMBL, OpenTargets),
               names_to = "source", values_to = "present") %>%
  filter(present) %>%
  count(source, name = "n_drugs") %>%
  mutate(source = factor(source, levels = c("DGIdb", "ChEMBL", "OpenTargets")))

p_sources <- ggplot(source_long,
                    aes(x = n_drugs, y = reorder(source, n_drugs), fill = source)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n_drugs), hjust = -0.25, size = 2.6) +
  scale_fill_manual(values = OKB[1:3]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Drug candidates per database source",
       x = "Number of drugs", y = NULL) +
  theme_pub() +
  theme(legend.position = "none",
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_pub(p_sources, "C1_drug_sources_bar")
cat("  C1: Drug sources bar — OK\n")

# ── C2: Distribución de fases clínicas ────────────────────────────────────────
phase_df <- drug_sum %>%
  mutate(
    phase_label = case_when(
      max_phase == 4 ~ "Phase IV (Approved)",
      max_phase == 3 ~ "Phase III",
      max_phase == 2 ~ "Phase II",
      max_phase == 1 ~ "Phase I",
      TRUE           ~ "Unknown"
    ),
    phase_label = factor(phase_label,
      levels = c("Phase IV (Approved)", "Phase III", "Phase II",
                 "Phase I", "Unknown"))
  ) %>%
  count(phase_label) %>%
  filter(!is.na(phase_label))

p_phase <- ggplot(phase_df,
                  aes(x = phase_label, y = n, fill = phase_label)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = PHASE_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Clinical development phase of drug candidates",
       x = NULL, y = "Number of drugs") +
  theme_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

save_pub(p_phase, "C2_drug_phase_dist")
cat("  C2: Phase distribution — OK\n")

# =============================================================================
# SECCIÓN D — RED PPI + SCORING
# =============================================================================
cat("\n--- Sección D: Network + Scoring ---\n")

# ── D1: logFC vs grado PPI ────────────────────────────────────────────────────
net_df <- net_nodes %>%
  mutate(
    de_dir = case_when(
      logFC_TVsS > 0 & adj.P.Val_TVsS < 0.05 ~ "up",
      logFC_TVsS < 0 & adj.P.Val_TVsS < 0.05 ~ "down",
      TRUE ~ "ns"
    ),
    de_dir = factor(de_dir, levels = c("up", "down", "ns"))
  )

hub_labels <- net_df %>%
  filter(de_dir != "ns") %>%
  slice_max(degree, n = 14)

p_deg_fc <- ggplot(net_df, aes(x = degree, y = logFC_TVsS,
                                color = de_dir, size = betweenness_norm)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             linewidth = 0.4, color = "grey55") +
  geom_point(aes(alpha = de_dir), stroke = 0) +
  geom_label_repel(
    data = hub_labels, aes(label = gene_symbol),
    size = 1.8, label.padding = 0.1, label.size = 0.15,
    max.overlaps = 12, segment.size = 0.25, show.legend = FALSE
  ) +
  scale_x_continuous(trans = "log10",
                     labels = label_number(accuracy = 1)) +
  scale_color_manual(
    values = DE_COLS,
    labels = c(up   = "Up-regulated",
               down = "Down-regulated",
               ns   = "Non-significant")
  ) +
  scale_alpha_manual(values = c(up = 0.85, down = 0.85, ns = 0.25),
                     guide = "none") +
  scale_size_continuous(range = c(0.5, 4.5), name = "Betweenness\n(norm.)") +
  labs(
    title = "PPI network centrality vs. proteomic fold change",
    x     = "PPI degree  (log scale)",
    y     = expression(log[2]~"Fold Change  (TVsS)"),
    color = NULL
  ) +
  theme_pub() +
  theme(legend.position = c(0.82, 0.82))

save_pub(p_deg_fc, "D1_degree_vs_logFC")
cat("  D1: Degree vs logFC scatter — OK\n")

# ── D2: Dot matrix de componentes de score ────────────────────────────────────
score_cols   <- c("s_logfc","s_sig","s_clinical","s_cmap","s_pathway","s_network")
score_labels <- c("Proteomics\n(logFC)","Proteomics\n(FDR)","Clinical\nphase",
                  "CMap\nconn.","Pathway\nrelevance","Network\nhub")

scores_long <- top20 %>%
  slice_head(n = 20) %>%
  mutate(drug_label = str_to_title(drug_name_norm) %>% str_trunc(26)) %>%
  select(drug_label, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "component", values_to = "score") %>%
  mutate(
    component  = factor(component, levels = score_cols, labels = score_labels),
    drug_label = factor(drug_label, levels = rev(unique(drug_label)))
  )

p_dot_matrix <- ggplot(scores_long,
                       aes(x = component, y = drug_label,
                           size = score, color = score)) +
  geom_point() +
  scale_size_continuous(range = c(0.3, 5.5), name = "Score\n(0 – 1)",
                        limits = c(0, 1)) +
  scale_color_gradientn(
    colors  = c("#0072B2", "white", "#D55E00"),
    limits  = c(0, 1),
    name    = "Score\n(0 – 1)"
  ) +
  labs(title = "Scoring components — Top 20 candidates",
       x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 7),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    legend.position  = "right"
  )

save_pub(p_dot_matrix, "D2_scoring_components_dot", "double_col", h_add = 60)
cat("  D2: Scoring dot matrix — OK\n")

# ── D3: Top 20 lollipop ───────────────────────────────────────────────────────
# drug_class (A/B/C/D) → etiqueta legible; si no existe la columna se advierte
top20_plot <- top20 %>%
  slice_head(n = 20) %>%
  mutate(
    drug_label  = str_to_title(drug_name_norm) %>% str_trunc(28),
    drug_label  = factor(drug_label, levels = rev(drug_label)),
    class_short = if ("drug_class" %in% colnames(top20)) {
      DRUG_CLASS_LABELS[drug_class]
    } else {
      cat("  WARN: columna 'drug_class' no encontrada en top20 — sin color por clase\n")
      rep("Unknown", n())
    }
  )

p_lollipop <- ggplot(top20_plot,
                     aes(x = final_score, y = drug_label,
                         color = class_short)) +
  geom_segment(aes(xend = 0, yend = drug_label), linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1))) +
  scale_color_manual(values = setNames(OKB[seq_along(DRUG_CLASS_LABELS)],
                                       DRUG_CLASS_LABELS),
                     na.value = "grey50") +
  labs(title = "Top 20 drug candidates — composite score",
       x = "Final composite score",
       y = NULL,
       color = "Drug class") +
  theme_pub() +
  theme(
    legend.position    = "right",
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

save_pub(p_lollipop, "D3_top20_lollipop", "double_col", h_add = 55)
cat("  D3: Top 20 lollipop — OK\n")

# =============================================================================
# SECCIÓN E — EVIDENCIA FINAL
# =============================================================================
cat("\n--- Sección E: Evidencia final ---\n")

# ── E1: Heatmap de evidencia binaria (ComplexHeatmap) ────────────────────────
ev_mat_raw <- evidence %>%
  slice_head(n = 20) %>%
  select(-any_of(c("final_rank", "evidence_level")))

# Detectar columna de nombre de fármaco (puede ser drug_name o drug_name_norm)
ev_drug_col <- intersect(c("drug_name", "drug_name_norm"), colnames(ev_mat_raw))[1]
if (is.na(ev_drug_col))
  stop("Columna de nombre de fármaco no encontrada en 13_evidence_matrix.tsv")

ev_drug_names <- str_to_title(ev_mat_raw[[ev_drug_col]]) %>% str_trunc(28)
ev_mat <- ev_mat_raw %>%
  select(-all_of(ev_drug_col)) %>%
  mutate(across(everything(), as.integer)) %>%
  as.matrix()
rownames(ev_mat) <- ev_drug_names

# Categorías por patrón de nombre de columna (robusto ante cambios de idioma)
ev_cats  <- assign_ev_category(colnames(ev_mat))
cat_cols <- c(Proteomics = "#E69F00", Databases = "#56B4E9",
              Network    = "#009E73", Clinical  = "#D55E00",
              Literature = "#CC79A7", Genomics  = "#0072B2",
              Other      = "#AAAAAA")

ev_level <- if ("evidence_level" %in% colnames(evidence)) {
  evidence %>% slice_head(n = 20) %>% pull(evidence_level)
} else {
  NULL
}

ha_ev_top <- HeatmapAnnotation(
  Category = ev_cats,
  col = list(Category = cat_cols[intersect(names(cat_cols), unique(ev_cats))]),
  annotation_name_gp   = gpar(fontsize = 6.5),
  simple_anno_size     = unit(3, "mm"),
  show_annotation_name = TRUE
)

ha_ev_row <- if (!is.null(ev_level)) {
  rowAnnotation(
    Level = as.character(ev_level),
    col = list(Level = c("1" = "#D55E00", "2" = "#E69F00",
                         "3" = "#56B4E9", "4" = "#009E73")),
    annotation_name_gp = gpar(fontsize = 6),
    simple_anno_size   = unit(3, "mm")
  )
} else {
  NULL
}

ht_ev <- Heatmap(
  ev_mat,
  name     = "Evidence",
  col      = c("0" = "#F5F5F5", "1" = "#0072B2"),
  top_annotation    = ha_ev_top,
  right_annotation  = ha_ev_row,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  row_names_gp      = gpar(fontsize = 7),
  column_names_gp   = gpar(fontsize = 7),
  column_names_rot  = 40,
  rect_gp           = gpar(col = "grey80", lwd = 0.5),
  column_title      = "Evidence matrix — Top 20 candidates",
  column_title_gp   = gpar(fontsize = 8, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5),
    at        = c(0, 1),
    labels    = c("Absent", "Present")
  )
)

save_ch(ht_ev, "E1_evidence_heatmap", "double_col", w_add = 20, h_add = 55)
cat("  E1: Evidence heatmap — OK\n")

# =============================================================================
# SECCIÓN F — ANÁLISIS DE SENSIBILIDAD
# =============================================================================
cat("\n--- Sección F: Sensibilidad ---\n")

config_names <- c(
  baseline        = "Baseline",
  clinical_heavy  = "Clinical\nheavy",
  molecular_heavy = "Molecular\nheavy",
  network_heavy   = "Network\nheavy",
  pathway_heavy   = "Pathway\nheavy",
  equal_weights   = "Equal\nweights"
)

# ── F1: Bump chart ────────────────────────────────────────────────────────────
bump_drugs <- sens %>%
  filter(n_configs_top20 >= 4) %>%
  arrange(baseline) %>%
  pull(drug_name_norm)

n_bump <- length(bump_drugs)
bump_palette <- setNames(colorRampPalette(OKB)(n_bump),
                         str_to_title(bump_drugs) %>% str_trunc(22))

bump_df <- sens %>%
  filter(drug_name_norm %in% bump_drugs) %>%
  select(drug_name_norm, all_of(names(config_names))) %>%
  mutate(drug_label = str_to_title(drug_name_norm) %>% str_trunc(22)) %>%
  pivot_longer(all_of(names(config_names)),
               names_to = "config", values_to = "rank") %>%
  mutate(
    config     = factor(config, levels = names(config_names)),
    config_num = as.numeric(config)
  )

p_bump <- ggplot(bump_df, aes(x = config_num, y = rank,
                               group = drug_label, color = drug_label)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 2.0) +
  geom_text(
    data = bump_df %>% filter(config_num == 1),
    aes(label = drug_label, x = 0.85), hjust = 1, size = 1.85,
    check_overlap = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:6,
    labels = unname(config_names),
    expand = expansion(add = c(2.8, 0.4))
  ) +
  # limits en orden natural (min, max); scale_y_reverse invierte la visualización
  scale_y_reverse(breaks = seq(1, 20, 2), limits = c(1, 20)) +
  scale_color_manual(values = bump_palette) +
  labs(title = "Rank stability across weight configurations",
       x     = "Weight configuration",
       y     = "Rank in Top 20") +
  theme_pub() +
  theme(
    legend.position    = "none",
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    axis.line.y        = element_blank(),
    axis.text.x        = element_text(size = 7)
  )

save_pub(p_bump, "F1_bump_chart", "double_col", w_add = 30, h_add = 10)
cat("  F1: Bump chart — OK\n")

# ── F2: Stability bar ─────────────────────────────────────────────────────────
stab_df <- sens %>%
  mutate(
    drug_label = str_to_title(drug_name_norm) %>% str_trunc(25),
    drug_label = factor(drug_label,
                        levels = drug_label[order(n_configs_top20,
                                                  -baseline, na.last = TRUE)])
  )

p_stab <- ggplot(stab_df, aes(x = n_configs_top20, y = drug_label,
                               fill = factor(n_configs_top20))) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 6, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 5.9, y = 1.5, label = "All configs",
           hjust = 1, size = 2.0, color = "grey40") +
  scale_fill_manual(
    values = c("1" = "#CCCCCC", "2" = "#CC79A7", "3" = "#E69F00",
               "4" = "#56B4E9", "5" = "#009E73", "6" = "#D55E00"),
    name = "# configs\nin Top 20"
  ) +
  scale_x_continuous(breaks = 1:6,
                     expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Candidate robustness — presence across weight configurations",
       x     = "# weight configurations in Top 20  (max = 6)",
       y     = NULL) +
  theme_pub() +
  theme(
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

save_pub(p_stab, "F2_stability_bar", "double_col", h_add = 70)
cat("  F2: Stability bar — OK\n")

# =============================================================================
# FIGURAS MULTIPANEL (para manuscrito)
# =============================================================================
cat("\n--- Figuras multipanel para manuscrito ---\n")

tag_theme <- theme(plot.tag = element_text(size = 10, face = "bold"))

# FIG1: QC/DE — Volcano + MA
p_fig1 <- (p_volcano | p_ma) +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_pub(p_fig1, "FIG1_QC_DE", "double_col")
cat("  FIG1: Volcano + MA — OK\n")

# FIG2: Pathways — Hallmarks
p_fig2 <- p_hall +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_pub(p_fig2, "FIG2_pathways", "double_col", h_add = 50)
cat("  FIG2: Pathways — OK\n")

# FIG3: Scoring — Lollipop + Dot matrix
p_fig3 <- (p_lollipop | p_dot_matrix) +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_pub(p_fig3, "FIG3_scoring", "double_col", h_add = 60)
cat("  FIG3: Scoring — OK\n")

# FIG4: Network — Degree vs logFC
p_fig4 <- p_deg_fc +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_pub(p_fig4, "FIG4_network", "double_col")
cat("  FIG4: Network — OK\n")

# FIG5: Sensibilidad — Bump + Stability
p_fig5 <- (p_bump | p_stab) +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_pub(p_fig5, "FIG5_sensitivity", "double_col", w_add = 20, h_add = 70)
cat("  FIG5: Sensibilidad — OK\n")

# =============================================================================
# RESUMEN FINAL
# =============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("Script 17 COMPLETO\n")
cat("Output:", pub_dir, "\n")

pub_files <- list.files(pub_dir, pattern = "\\.pdf$", full.names = FALSE)
cat(sprintf("PDFs generados: %d\n", length(pub_files)))
cat("  ", paste(pub_files, collapse = "\n  "), "\n")
cat("Fin:", format(Sys.time()), "\n")

sink()
