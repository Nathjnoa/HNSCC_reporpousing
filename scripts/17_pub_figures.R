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
#     A4_heatmap_topDE    — Heatmap top 40 proteínas DE × 20 muestras (20 up + 20 down)
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
  library(igraph)
  library(ggraph)
  library(tidygraph)
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
  "Aprobado"   = "#009E73",
  "Fase III"   = "#56B4E9",
  "Fase II"    = "#E69F00",
  "Fase I"     = "#CC79A7",
  "Sin datos"  = "#CCCCCC"
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
net_nodes    <- read_tsv("results/tables/network/09_network_node_metrics.tsv",
                        show_col_types = FALSE)
net_edges    <- read_tsv("results/tables/network/09_network_edges.tsv",
                        show_col_types = FALSE)
drg_hubs     <- read_tsv("results/tables/network/09_druggable_hubs.tsv",
                        show_col_types = FALSE)
net_modules  <- read_tsv("results/tables/network/09_modules.tsv",
                        show_col_types = FALSE)
top20     <- read_tsv("results/tables/10_top20_candidates.tsv",
                      show_col_types = FALSE)
evidence  <- read_tsv("results/tables/13_evidence_matrix.tsv",
                      show_col_types = FALSE)
sens      <- read_tsv("results/tables/15_sensitivity_ranks.tsv",
                      show_col_types = FALSE)
multi_src <- read_tsv("results/tables/drug_targets/08_multi_source_candidates.tsv",
                      show_col_types = FALSE)

cat("  Datos cargados OK\n")

# ── Cargar metadata (HPV status) ──────────────────────────────────────────────
meta <- read_delim("data/raw/metadata.csv", delim = ";", show_col_types = FALSE) %>%
  mutate(hpv = ifelse(tolower(vph) == "positive", "HPV+", "HPV-"))
hpv_lookup <- setNames(meta$hpv, meta$sample_id)
HPV_COLS   <- c("HPV+" = "#009E73", "HPV-" = "#E69F00")
cat("  Metadata cargada OK —", sum(meta$hpv == "HPV+"), "HPV+,",
    sum(meta$hpv == "HPV-"), "HPV-\n")

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
# Top 8 por dirección + proteínas de interés terapéutico forzadas independientemente del rango
force_label_genes <- c("EGFR", "PSMB5", "PSMB2", "PSMA1", "PSMD11",
                        "DNMT1", "TOP2A", "NDUFA9", "NDUFS3", "NDUFB3")
top_vol <- bind_rows(
  de %>% filter(direction == "up")   %>% arrange(desc(logFC_TVsS)) %>% slice_head(n = 8),
  de %>% filter(direction == "down") %>% arrange(logFC_TVsS)       %>% slice_head(n = 8),
  de %>% filter(gene_symbol %in% force_label_genes, direction != "ns")
) %>% distinct(gene_symbol, .keep_all = TRUE)
force_found <- intersect(force_label_genes, de$gene_symbol[de$direction != "ns"])
cat(sprintf("  A1: %d proteínas forzadas en etiquetas: %s\n",
            length(force_found), paste(force_found, collapse = ", ")))

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
  labs(x     = expression(log[2]~"cambio de expresión (TVsS)"),
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
  labs(       x     = expression("Mean"~log[2]~"intensity (A)"),
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
    patient   = str_remove(sample_id, "[TS]$"),
    hpv       = hpv_lookup[sample_id]
  )

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = hpv)) +
  geom_line(aes(group = patient), color = "grey70",
            linewidth = 0.4, linetype = "dotted") +
  geom_point(size = 2.5, stroke = 0.6) +
  scale_color_manual(values = COND_COLS) +
  scale_shape_manual(values = c("HPV+" = 17, "HPV-" = 16), na.value = 1) +
  labs(       x     = sprintf("PC1  (%.1f%%)", var_exp[1]),
       y     = sprintf("PC2  (%.1f%%)", var_exp[2]),
       color = "Condition", shape = "HPV status") +
  guides(
    color = guide_legend(override.aes = list(size = 2.5, shape = 15)),
    shape = guide_legend(override.aes = list(size = 2.5))
  ) +
  theme_pub() +
  theme(legend.position = c(0.80, 0.12))

save_pub(p_pca, "A3_PCA")
cat("  A3: PCA — OK\n")

# ── A4: Heatmap top 40 proteínas DE con detección completa ───────────────────
# Criterio de selección: top proteínas DE detectadas en las 20/20 muestras
# (complete cases). Con este filtro hay 194 up y 137 down elegibles, más que
# suficientes para elegir 20+20 sin necesidad de imputación.
# Justificación en Methods: "Se seleccionaron las proteínas DE con detección
# en el 100% de las muestras, ordenadas por |log2FC|."
de_detected <- de %>%
  mutate(n_detected = rowSums(!is.na(across(all_of(sample_cols))))) %>%
  filter(n_detected == length(sample_cols))   # 20/20 muestras

top_up   <- de_detected %>% filter(direction == "up")   %>%
  slice_max(logFC_TVsS, n = 20) %>% pull(gene_symbol)
top_down <- de_detected %>% filter(direction == "down") %>%
  slice_min(logFC_TVsS, n = 20) %>% pull(gene_symbol)
top_genes <- c(top_up, top_down)

cat(sprintf("  A4: %d up + %d down elegibles (20/20 detección); seleccionados %d+%d\n",
            sum(de_detected$direction == "up"),
            sum(de_detected$direction == "down"),
            length(top_up), length(top_down)))

# Construir matriz desde expr_mat (ya es complete cases, todas las proteínas aquí están)
heat_mat <- expr_mat[top_genes[top_genes %in% rownames(expr_mat)], , drop = FALSE]
heat_z   <- t(scale(t(heat_mat)))   # Z-score estándar por proteína

# Ordenar columnas: agrupar Tumor / Adjacent Normal
col_cond  <- ifelse(str_ends(colnames(heat_z), "T"), "Tumor", "Adjacent Normal")
col_order <- order(col_cond)
heat_z    <- heat_z[, col_order]
col_cond  <- col_cond[col_order]

hpv_ann_cols <- hpv_lookup[colnames(heat_mat)[col_order]]
ht_col_ann <- HeatmapAnnotation(
  Condition = col_cond,
  HPV       = hpv_ann_cols,
  col = list(
    Condition = c(Tumor = "#D55E00", "Adjacent Normal" = "#0072B2"),
    HPV       = HPV_COLS
  ),
  annotation_name_gp = gpar(fontsize = 6.5),
  simple_anno_size   = unit(3, "mm")
)

# row_split: up genes primero, down genes después
n_up_heat   <- sum(rownames(heat_z) %in% top_up)
n_down_heat <- sum(rownames(heat_z) %in% top_down)
cat(sprintf("  A4: %d up + %d down en heatmap\n", n_up_heat, n_down_heat))
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
  column_title       = "Top 40 proteínas diferencialmente expresadas — Tumor vs. Tejido Normal Adyacente",
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
  labs(       x     = "Normalized enrichment score (NES)",
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
    OpenTargets = str_detect(sources, "OpenTargets"),
    L2S2        = str_detect(sources, "L2S2")
  ) %>%
  pivot_longer(c(DGIdb, ChEMBL, OpenTargets, L2S2),
               names_to = "source", values_to = "present") %>%
  filter(present) %>%
  count(source, name = "n_drugs") %>%
  mutate(source = factor(source, levels = c("DGIdb", "OpenTargets", "L2S2", "ChEMBL")))

p_sources <- ggplot(source_long,
                    aes(x = n_drugs, y = reorder(source, n_drugs), fill = source)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = n_drugs), hjust = -0.25, size = 2.6) +
  scale_fill_manual(values = OKB[1:4]) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Número de fármacos", y = NULL) +
  theme_pub() +
  theme(legend.position = "none",
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_pub(p_sources, "C1_drug_sources_bar")
cat("  C1: Drug sources bar — OK\n")

# ── C2: Distribución de fases clínicas (candidatos multi-fuente) ──────────────
# Usa multi_src en lugar de drug_sum para evitar el dominio de DGIdb-only (91% NA).
# Combina max_phase (ChEMBL) con is_approved (flag DGIdb) para recuperar
# los 222 aprobados cuya fase no está en ChEMBL.
phase_df <- multi_src %>%
  mutate(
    phase_label = case_when(
      max_phase == 4                         ~ "Aprobado",
      is_approved == TRUE & is.na(max_phase) ~ "Aprobado",
      max_phase == 3                         ~ "Fase III",
      max_phase == 2                         ~ "Fase II",
      max_phase == 1                         ~ "Fase I",
      TRUE                                   ~ "Sin datos"
    ),
    phase_label = factor(phase_label, levels = names(PHASE_COLS))
  ) %>%
  count(phase_label) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat(sprintf("  C2: %d candidatos multi-fuente; %d aprobados (%.0f%%)\n",
            nrow(multi_src),
            sum(phase_df$n[phase_df$phase_label == "Aprobado"]),
            sum(phase_df$pct[phase_df$phase_label == "Aprobado"])))

p_phase <- ggplot(phase_df,
                  aes(x = phase_label, y = n, fill = phase_label)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            vjust = -0.3, size = 2.3, lineheight = 0.9) +
  scale_fill_manual(values = PHASE_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    subtitle = sprintf("n = %d candidatos con evidencia en \u22652 bases de datos  \u2022  Aprobado = Fase IV ChEMBL o aprobado en DGIdb",
                       nrow(multi_src)),
    x = NULL, y = "Número de fármacos"
  ) +
  theme_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

save_pub(p_phase, "C2_drug_phase_dist")
cat("  C2: Phase distribution (multi-source) — OK\n")

# ── C3: UpSet plot — overlap entre bases de datos (OE1_FigC) ─────────────────
multi_src <- read_tsv("results/tables/drug_targets/08_multi_source_candidates.tsv",
                      show_col_types = FALSE)

upset_list <- list(
  DGIdb       = multi_src$drug_name_norm[str_detect(multi_src$sources, "DGIdb")],
  ChEMBL      = multi_src$drug_name_norm[str_detect(multi_src$sources, "ChEMBL")],
  OpenTargets = multi_src$drug_name_norm[str_detect(multi_src$sources, "OpenTargets")],
  L2S2        = multi_src$drug_name_norm[str_detect(multi_src$sources, "L2S2")]
)

comb_m <- make_comb_mat(upset_list)

ht_upset <- UpSet(
  comb_m,
  comb_col    = "#0072B2",
  pt_size     = unit(3.5, "mm"),
  lwd         = 1.5,
  set_order   = order(set_size(comb_m), decreasing = TRUE),
  comb_order  = order(comb_size(comb_m), decreasing = TRUE),
  top_annotation = upset_top_annotation(
    comb_m,
    add_numbers      = TRUE,
    numbers_gp       = gpar(fontsize = 7),
    annotation_name_gp = gpar(fontsize = 7)
  ),
  right_annotation = upset_right_annotation(
    comb_m,
    add_numbers      = TRUE,
    numbers_gp       = gpar(fontsize = 7),
    annotation_name_gp = gpar(fontsize = 7)
  ),
  row_names_gp = gpar(fontsize = 8),
  column_title = "Solapamiento de fármacos candidatos entre bases de datos",
  column_title_gp = gpar(fontsize = 8, fontface = "bold")
)

save_ch(ht_upset, "C3_upset_sources", "double_col", w_add = 10, h_add = 20)
cat("  C3: UpSet plot sources — OK\n")

# =============================================================================
# SECCIÓN OE2 — RED PPI + HUBS DRUGGABLES + PANEL FINAL
# =============================================================================
cat("\n--- Sección OE2: Red PPI + panel final ---\n")

# ── OE2_FigA: Red PPI de proteínas DE ─────────────────────────────────────────
# Helper para construir y graficar red PPI (reutilizado en v1 y v2)
build_network_fig <- function(gene_set, edge_tbl, node_meta, hub_genes,
                               drg_hub_genes, label_genes, layout_algo,
                               seed = 42, title_sfx = "") {
  e <- edge_tbl %>% filter(gene_A %in% gene_set, gene_B %in% gene_set)
  n <- node_meta %>%
    filter(gene_symbol %in% gene_set) %>%
    mutate(
      de_dir     = case_when(
        logFC_TVsS > 0 & adj.P.Val_TVsS < 0.05 ~ "Up",
        logFC_TVsS < 0 & adj.P.Val_TVsS < 0.05 ~ "Down",
        TRUE                                     ~ "NS"
      ),
      de_dir     = factor(de_dir, levels = c("Up", "Down", "NS")),
      is_drg_hub = gene_symbol %in% drg_hub_genes,
      nd_size    = log1p(degree),
      nd_alpha   = ifelse(de_dir == "NS", 0.20, 0.82),
      nd_stroke  = ifelse(is_drg_hub, 0.7, 0.0),
      nd_color   = ifelse(is_drg_hub, "black",
                   ifelse(de_dir == "Up",  "#D55E00",
                   ifelse(de_dir == "Down","#0072B2", "#AAAAAA")))
    )

  g <- graph_from_data_frame(e, directed = FALSE, vertices = n$gene_symbol)
  for (col in c("de_dir","is_drg_hub","nd_size","nd_alpha","nd_stroke","nd_color","degree")) {
    vertex_attr(g, col) <- n[[col]][match(V(g)$name, n$gene_symbol)]
  }
  tg <- as_tbl_graph(g)

  set.seed(seed)
  lay <- create_layout(tg, layout = layout_algo)

  lbl_df <- lay %>%
    filter(name %in% label_genes[label_genes %in% gene_set]) %>%
    select(x, y, name, de_dir)

  p <- ggraph(lay) +
    geom_edge_link(color = "grey78", alpha = 0.28, linewidth = 0.14) +
    geom_node_point(
      aes(fill = de_dir, size = nd_size, alpha = nd_alpha,
          color = nd_color, stroke = nd_stroke),
      shape = 21
    ) +
    ggrepel::geom_label_repel(
      data          = lbl_df,
      aes(x = x, y = y, label = name, fill = de_dir),
      color         = "black",
      size          = 2.6,
      fontface      = "plain",
      label.size    = 0.18,
      label.padding = unit(0.11, "lines"),
      box.padding   = unit(0.28, "lines"),
      segment.size  = 0.25,
      segment.color = "grey40",
      max.overlaps  = 30,
      alpha         = 0.88,
      show.legend   = FALSE
    ) +
    scale_fill_manual(
      values = c(Up = "#D55E00", Down = "#0072B2", NS = "#CCCCCC"),
      name   = "DE direction",
      guide  = guide_legend(override.aes = list(size = 3.5, alpha = 0.92,
                                                stroke = 0, shape = 21))
    ) +
    scale_color_identity() +
    scale_size_continuous(range = c(0.7, 5.0), guide = "none") +
    scale_alpha_identity(guide = "none") +
    labs() +
    theme_graph(base_family = "sans", base_size = 8) +
    theme(
      legend.title    = element_text(size = 7.5, face = "bold"),
      legend.text     = element_text(size = 7),
      legend.position = c(0.02, 0.10),
      plot.margin     = margin(4, 4, 4, 4, "mm")
    )
  p
}

# Etiquetas curadas (presentes en ambas versiones si el nodo existe)
hub_label_genes <- c(
  "EGFR", "PSMB3", "PSMA2", "MMP9", "CASP3", "ENO1",
  "NDUFS3", "NDUFS2", "NDUFV1", "ATP5F1C", "SDHA", "UQCRC2", "NDUFA9"
)

# ── v1: red completa filtrada (degree > 8 | hub) — layout stress ──────────────
g_full <- graph_from_data_frame(net_edges, directed = FALSE,
                                vertices = net_nodes$gene_symbol)
comp        <- components(g_full)
giant_genes <- names(comp$membership[comp$membership == which.max(comp$csize)])

core_genes_v1 <- net_nodes %>%
  filter(gene_symbol %in% giant_genes,
         degree > 8 | is_hub == TRUE) %>%
  pull(gene_symbol)

p_net_v1 <- build_network_fig(
  gene_set    = core_genes_v1,
  edge_tbl    = net_edges,
  node_meta   = net_nodes,
  hub_genes   = net_nodes %>% filter(is_hub) %>% pull(gene_symbol),
  drg_hub_genes = drg_hubs$gene_symbol,
  label_genes = hub_label_genes,
  layout_algo = "stress",   # graphlayouts: minimiza diferencias en longitud de arista
  title_sfx   = " — all DE proteins (degree > 8 or hub)"
)
save_pub(p_net_v1, "OE2_FigA_ppi_network", "double_col", w_add = 60, h_add = 90)
cat("  OE2_FigA v1 (stress, degree>8|hub): PPI network — OK\n")

# ── v2: solo hubs (top 10% degree o betweenness) — layout stress ──────────────
hub_only_genes <- net_nodes %>%
  filter(gene_symbol %in% giant_genes, is_hub == TRUE) %>%
  pull(gene_symbol)

p_net_v2 <- build_network_fig(
  gene_set    = hub_only_genes,
  edge_tbl    = net_edges,
  node_meta   = net_nodes,
  hub_genes   = hub_only_genes,
  drg_hub_genes = drg_hubs$gene_symbol,
  label_genes = hub_label_genes,
  layout_algo = "stress",
  title_sfx   = " — hub proteins only (top 10%)"
)
save_pub(p_net_v2, "OE2_FigA_v2_hubs_only", "double_col", w_add = 60, h_add = 90)
cat("  OE2_FigA v2 (stress, hubs only): PPI network — OK\n")

# ── OE2_FigB: Fármacos del panel final × hubs que apuntan ────────────────────
# Diseño:
#   - Filas: fármacos priorizados (top20, ordenados por score)
#   - Columnas: hubs druggables clave (representativos por módulo biológico)
#   - Punto: el fármaco apunta a ese hub (presencia en top_drugs del hub)
#   - Color punto: dirección DE del hub (azul = down, naranja = up)
#   - Anotación columnas: módulo biológico del hub
# Objetivo: conectar el panel de fármacos con la biología de red

MODULE_LABELS <- c(
  "8"  = "OXPHOS",
  "4"  = "Proteasoma",
  "2"  = "EGFR / Señalización",
  "5"  = "Resp. inmune / ECM",
  "14" = "Metab. peq. mol."
)

# Módulos clave para mostrar (los más terapéuticamente relevantes)
KEY_MODULES <- c(8, 2, 4, 5, 14)

# Hubs representativos: top 2 por módulo clave
node_de_dir <- net_nodes %>%
  transmute(gene_symbol,
            de_dir = case_when(
              logFC_TVsS > 0 & adj.P.Val_TVsS < 0.05 ~ "Up",
              logFC_TVsS < 0 & adj.P.Val_TVsS < 0.05 ~ "Down",
              TRUE ~ "NS"
            ))

key_hubs_b <- drg_hubs %>%
  filter(module_id %in% KEY_MODULES) %>%
  left_join(node_de_dir, by = "gene_symbol") %>%
  mutate(module_label = MODULE_LABELS[as.character(module_id)]) %>%
  group_by(module_id, module_label) %>%
  slice_max(degree, n = 2) %>%
  ungroup() %>%
  arrange(de_dir, module_id, desc(degree))

# Cargar master table gen-fármaco
master_tbl <- read_tsv("results/tables/drug_targets/08_drug_target_master_table.tsv",
                       show_col_types = FALSE)

# ── OE2_FigB: Candidatos por módulo biológico (barplot apilado por clase) ─────
# Pregunta: ¿Cuántos fármacos apuntan a cada módulo? ¿Con qué respaldo clínico?
# Datos: todos los candidatos que apuntan a algún hub druggable, por módulo

hub_module_map <- drg_hubs %>%
  filter(module_id %in% KEY_MODULES) %>%
  mutate(module_label = MODULE_LABELS[as.character(module_id)]) %>%
  select(gene_symbol, module_label, module_id)

module_order <- c("EGFR / Señalización", "Proteasoma", "OXPHOS",
                  "Resp. inmune / ECM", "Metab. peq. mol.")

drug_per_module <- master_tbl %>%
  filter(gene_symbol %in% hub_module_map$gene_symbol,
         drug_class %in% c("A", "B", "C")) %>%
  inner_join(hub_module_map, by = "gene_symbol") %>%
  distinct(drug_name_norm, module_label, drug_class) %>%
  count(module_label, drug_class) %>%
  mutate(
    module_label = factor(module_label, levels = module_order),
    drug_class   = factor(drug_class, levels = c("A", "B", "C")),
    class_label  = DRUG_CLASS_LABELS[as.character(drug_class)],
    class_label  = factor(class_label,
                          levels = DRUG_CLASS_LABELS[c("A","B","C")])
  )

p_module_bar <- ggplot(drug_per_module,
                       aes(y = module_label, x = n, fill = class_label)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(
    values = c("HNSCC-approved" = OKB[1],
               "Other cancer"   = OKB[2],
               "Non-oncology"   = OKB[3]),
    name = "Clase de fármaco"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    subtitle = "Candidatos dirigidos a hubs  \u2022  Agrupados por módulo STRING PPI",
    x = "N° fármacos candidatos", y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y     = element_text(size = 7.5),
    legend.position = "right"
  )

save_pub(p_module_bar, "OE2_FigB_module_barplot", "double_col", h_add = 0)
cat("  OE2_FigB: N\u00b0 candidatos por m\u00f3dulo — OK\n")

# ── OE2_FigC: Distribución clase clínica A/B/C/D (candidatos hub-targeting) ───
# Pregunta: ¿Qué proporción de los candidatos son aprobados vs experimentales?

drug_class_counts <- master_tbl %>%
  filter(gene_symbol %in% drg_hubs$gene_symbol) %>%
  distinct(drug_name_norm, drug_class) %>%
  count(drug_class, name = "n_drugs") %>%
  mutate(
    class_label = DRUG_CLASS_LABELS[drug_class],
    class_label = factor(class_label, levels = DRUG_CLASS_LABELS),
    pct         = round(n_drugs / sum(n_drugs) * 100, 1)
  )

CLASS_COLS_OE2 <- setNames(OKB[1:4], DRUG_CLASS_LABELS)

p_class_bar <- ggplot(drug_class_counts,
                      aes(x = class_label, y = n_drugs, fill = class_label)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = paste0(n_drugs, "\n(", pct, "%)")),
            vjust = -0.35, size = 2.8, lineheight = 1.1) +
  scale_fill_manual(values = CLASS_COLS_OE2, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title    = "Clinical classification of hub-targeting candidates",
    subtitle = paste0("N = ", sum(drug_class_counts$n_drugs),
                      " unique candidates targeting druggable hub proteins"),
    x = NULL, y = "N\u00b0 drug candidates"
  ) +
  theme_pub() +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 7.5),
    plot.title   = element_text(size = 8.5),
    plot.subtitle = element_text(size = 6.5)
  )

save_pub(p_class_bar, "OE2_FigC_class_distribution", "double_col", h_add = 10)
cat("  OE2_FigC: Distribuci\u00f3n clase A/B/C/D — OK\n")

# ── OE2_FigD — Opciones para el panel de candidatos priorizados ───────────────
# Opción 1: Dot plot clínico (módulo × candidato, color=clase, tamaño=fuentes)
# Opción 2: Tabla visual como figura (filas=candidatos, columnas=atributos clave)

# --- Datos comunes: 32 candidatos LOD-stable con hub primario y módulo ---
lod_stable <- read_tsv("results/tables/15_lod_stability.tsv",
                       show_col_types = FALSE)

# Etiquetas legibles para todos los módulos de la red
module_labels_full <- net_modules %>%
  distinct(module_id, module_name) %>%
  mutate(
    mod_label = module_name %>%
      str_remove("^M\\d+_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence() %>%
      str_trunc(24)
  )
# Sobrescribir con etiquetas curadas para módulos clave
module_labels_full <- module_labels_full %>%
  mutate(mod_label = case_when(
    module_id == 2  ~ "EGFR / Signaling",
    module_id == 4  ~ "Proteasome",
    module_id == 8  ~ "OXPHOS",
    module_id == 5  ~ "Immune / ECM",
    module_id == 14 ~ "Small mol. metab.",
    TRUE ~ mod_label
  ))

top_lod_ann <- top20 %>%
  inner_join(lod_stable %>% filter(lod_stable) %>% select(drug_name_norm),
             by = "drug_name_norm") %>%
  arrange(desc(final_score)) %>%
  mutate(hub_gene = primary_target) %>%
  left_join(net_modules %>% distinct(gene_symbol, module_id),
            by = c("hub_gene" = "gene_symbol")) %>%
  left_join(module_labels_full %>% select(module_id, mod_label),
            by = "module_id") %>%
  mutate(
    module_label = coalesce(mod_label, "—"),
    drug_label   = str_to_title(drug_name_norm) %>% str_trunc(26),
    drug_label   = factor(drug_label, levels = rev(unique(drug_label))),
    class_label  = factor(DRUG_CLASS_LABELS[drug_class],
                          levels = DRUG_CLASS_LABELS),
    hub_gene     = ifelse(is.na(hub_gene), "—", hub_gene)
  )

# alias para compatibilidad con código siguiente
top20_ann <- top_lod_ann

# --- Opción 1: dot plot clínico ---
p_dot_clinical <- ggplot(top20_ann,
                          aes(x = module_label, y = drug_label)) +
  geom_point(aes(color = class_label, size = n_sources), alpha = 0.90) +
  geom_text(aes(label = hub_gene),
            hjust = -0.30, size = 1.85, color = "grey30",
            fontface = "italic") +
  scale_color_manual(values = CLASS_COLS_OE2, name = "Drug class",
                     drop = FALSE) +
  scale_size_continuous(range = c(2.5, 6.5), name = "N\u00b0 sources",
                        breaks = c(1, 2, 3, 4)) +
  labs(
    title    = "Top prioritized candidates — clinical & biological profile",
    subtitle = "Ordered by composite score (top \u2192 bottom)  \u2022  Color: clinical class  \u2022  Size: data sources  \u2022  Label: hub target",
    x = "Biological module (hub target)", y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x      = element_text(angle = 18, hjust = 1, size = 7),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = "right"
  )

save_pub(p_dot_clinical, "OE2_FigD_opt1_dot_clinical", "double_col", h_add = 60)
cat("  OE2_FigD opt1: dot plot cl\u00ednico — OK\n")

# --- Opción 2: tabla visual como figura ---
# Columnas: Fármaco | Módulo | Hub | Clase | Fuentes | Score (barra)
tbl_df <- top20_ann %>%
  arrange(desc(final_score)) %>%
  mutate(
    rank       = seq_len(n()),
    score_pct  = final_score / max(final_score),
    src_dots   = strrep("\u25cf", n_sources)   # puntos como icono de fuentes
  ) %>%
  select(rank, drug_label, module_label, hub_gene, class_label,
         src_dots, final_score, score_pct, drug_class)

# Posiciones de columnas (en unidades del plot, x de 0 a 1)
COL_X <- c(rank=0.03, drug=0.16, module=0.38, hub=0.57,
            class=0.70, sources=0.83, score=0.94)
N_ROWS  <- nrow(tbl_df)
ROW_Y   <- seq(N_ROWS, 1) / (N_ROWS + 2)  # de arriba a abajo
HEADER_Y <- (N_ROWS + 1.5) / (N_ROWS + 2)

class_bg <- c(A="#FFF3CD", B="#D6EAF8", C="#D5F5E3", D="#F5EEF8")

p_vis_table <- ggplot() +
  # Fondo por fila (color clase)
  geom_rect(
    data = tbl_df,
    aes(xmin = 0, xmax = 1,
        ymin = ROW_Y - 0.5/(N_ROWS+2),
        ymax = ROW_Y + 0.5/(N_ROWS+2),
        fill = drug_class),
    alpha = 0.35
  ) +
  scale_fill_manual(values = class_bg, guide = "none") +
  # Barra de score (fondo gris)
  geom_rect(
    data = tbl_df,
    aes(xmin = COL_X["score"] - 0.02,
        xmax = COL_X["score"] + 0.055,
        ymin = ROW_Y - 0.38/(N_ROWS+2),
        ymax = ROW_Y + 0.38/(N_ROWS+2)),
    fill = "grey88"
  ) +
  # Barra de score (valor)
  geom_rect(
    data = tbl_df,
    aes(xmin = COL_X["score"] - 0.02,
        xmax = COL_X["score"] - 0.02 + score_pct * 0.075,
        ymin = ROW_Y - 0.38/(N_ROWS+2),
        ymax = ROW_Y + 0.38/(N_ROWS+2),
        fill = drug_class),
    alpha = 0.80
  ) +
  # Texto: rank
  geom_text(data=tbl_df, aes(x=COL_X["rank"], y=ROW_Y, label=rank),
            size=2.3, hjust=0.5, color="grey40") +
  # Texto: fármaco (negrita)
  geom_text(data=tbl_df, aes(x=COL_X["drug"], y=ROW_Y, label=drug_label),
            size=2.35, hjust=0, fontface="bold") +
  # Texto: módulo
  geom_text(data=tbl_df, aes(x=COL_X["module"], y=ROW_Y, label=module_label),
            size=2.1, hjust=0, color="grey25") +
  # Texto: hub
  geom_text(data=tbl_df, aes(x=COL_X["hub"], y=ROW_Y, label=hub_gene),
            size=2.2, hjust=0, fontface="italic", color="#333333") +
  # Texto: clase (coloreado)
  geom_text(data=tbl_df, aes(x=COL_X["class"], y=ROW_Y,
                               label=as.character(class_label), color=drug_class),
            size=2.0, hjust=0, fontface="bold") +
  scale_color_manual(values=c(A="#B8860B",B="#1A6EA8",C="#1A8A5A",D="#7B4FA0"),
                     guide="none") +
  # Texto: fuentes (puntos)
  geom_text(data=tbl_df, aes(x=COL_X["sources"], y=ROW_Y, label=src_dots),
            size=2.1, hjust=0.5, color="grey30") +
  # Línea separadora header
  geom_hline(yintercept = HEADER_Y - 0.5/(N_ROWS+2),
             color="grey50", linewidth=0.4) +
  # Headers
  annotate("text", x=COL_X["rank"],    y=HEADER_Y, label="#",
           size=2.5, fontface="bold", hjust=0.5) +
  annotate("text", x=COL_X["drug"],    y=HEADER_Y, label="Fármaco candidato",
           size=2.5, fontface="bold", hjust=0) +
  annotate("text", x=COL_X["module"],  y=HEADER_Y, label="Módulo",
           size=2.5, fontface="bold", hjust=0) +
  annotate("text", x=COL_X["hub"],     y=HEADER_Y, label="Target hub",
           size=2.5, fontface="bold", hjust=0) +
  annotate("text", x=COL_X["class"],   y=HEADER_Y, label="Clase",
           size=2.5, fontface="bold", hjust=0) +
  annotate("text", x=COL_X["sources"], y=HEADER_Y, label="Fuentes",
           size=2.5, fontface="bold", hjust=0.5) +
  annotate("text", x=COL_X["score"],   y=HEADER_Y, label="Puntuación",
           size=2.5, fontface="bold", hjust=0.5) +
  scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
  theme_void(base_family="sans") +
  theme(plot.margin = margin(4,4,4,4,"mm"))

save_pub(p_vis_table, "OE2_FigD_opt2_vis_table", "double_col", h_add = 80)
cat("  OE2_FigD opt2: tabla visual — OK\n")

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
score_cols   <- c("s_pi_stat","s_clinical","s_cmap","s_pathway","s_network")
score_labels <- c("Proteomics\n(pi-stat)","Clinical\nphase",
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
  labs(       x = NULL, y = NULL) +
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
  labs(       x = "Final composite score",
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
  labs(       x     = "Weight configuration",
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
  labs(       x     = "# weight configurations in Top 20  (max = 6)",
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
