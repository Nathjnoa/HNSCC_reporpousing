#!/usr/bin/env Rscript
# =============================================================================
# 17_pub_figures.R
# =============================================================================
# Genera las 7 figuras de publicación incluidas en el reporte revisado
# (HNSCC_DrugRepurposing_Figuras_reviewed.docx).
#
# Figuras generadas (results/figures/pub/main/):
#   Sec0_FigB_volcano          — Volcano plot Tumor vs. Adjacent Normal
#   Sec0_FigC_heatmap_topDE    — Heatmap top 40 proteínas DE × 20 muestras
#   OE1_FigA_drug_sources_bar  — Nº candidatos por fuente de BD
#   OE1_FigB_drug_phase_dist   — Distribución fases clínicas (multi-fuente)
#   OE2_FigA_ppi_network       — Red PPI de proteínas DE (degree > 8 | hub)
#   OE2_FigB_module_barplot    — Candidatos aprobados por módulo PPI
#   OE2_FigC_class_distribution — Clasificación clínico-regulatoria (A/B/C/D)
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
pub_dir <- "results/figures/pub/main"
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

# =============================================================================
# CARGAR DATOS
# =============================================================================
cat("\n-- Cargando datos...\n")

de        <- read_tsv("results/tables/de_limma/01_TVsS_all_proteins.tsv",
                      show_col_types = FALSE)
drug_sum  <- read_tsv("results/tables/drug_targets/08_drug_summary_per_drug.tsv",
                      show_col_types = FALSE)
multi_src <- read_tsv("results/tables/drug_targets/08_multi_source_candidates.tsv",
                      show_col_types = FALSE)
master_tbl <- read_tsv("results/tables/drug_targets/08_drug_target_master_table.tsv",
                       show_col_types = FALSE)
net_nodes    <- read_tsv("results/tables/network/09_network_node_metrics.tsv",
                        show_col_types = FALSE)
net_edges    <- read_tsv("results/tables/network/09_network_edges.tsv",
                        show_col_types = FALSE)
drg_hubs     <- read_tsv("results/tables/network/09_druggable_hubs.tsv",
                        show_col_types = FALSE)
net_modules  <- read_tsv("results/tables/network/09_modules.tsv",
                        show_col_types = FALSE)

cat("  Datos cargados OK\n")

# ── Cargar metadata (HPV status) ──────────────────────────────────────────────
meta <- read_delim("data/raw/metadata.csv", delim = ";", show_col_types = FALSE) %>%
  mutate(hpv = ifelse(tolower(vph) == "positive", "HPV+", "HPV-"))
hpv_lookup <- setNames(meta$hpv, meta$sample_id)
HPV_COLS   <- c("HPV+" = "#009E73", "HPV-" = "#E69F00")
cat("  Metadata cargada OK —", sum(meta$hpv == "HPV+"), "HPV+,",
    sum(meta$hpv == "HPV-"), "HPV-\n")

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
sample_cols <- intersect(sample_cols, colnames(de))

expr_mat <- de %>%
  select(gene_symbol, all_of(sample_cols)) %>%
  filter(complete.cases(.)) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

# =============================================================================
# Sec0_FigB — VOLCANO PLOT (Tumor vs. Adjacent Normal)
# =============================================================================
cat("\n--- Sec0_FigB: Volcano plot ---\n")

force_label_genes <- c("EGFR", "PSMB5", "PSMB2", "PSMA1", "PSMD11",
                        "DNMT1", "TOP2A", "NDUFA9", "NDUFS3", "NDUFB3")
top_vol <- bind_rows(
  de %>% filter(direction == "up")   %>% arrange(desc(logFC_TVsS)) %>% slice_head(n = 8),
  de %>% filter(direction == "down") %>% arrange(logFC_TVsS)       %>% slice_head(n = 8),
  de %>% filter(gene_symbol %in% force_label_genes, direction != "ns")
) %>% distinct(gene_symbol, .keep_all = TRUE)
force_found <- intersect(force_label_genes, de$gene_symbol[de$direction != "ns"])
cat(sprintf("  %d proteínas forzadas en etiquetas: %s\n",
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
  labs(title = "Differential proteome — Tumor vs. Adjacent Normal",
       x     = expression(log[2]~"Fold Change  (TVsS)"),
       y     = expression(-log[10]~"FDR"),
       color = NULL) +
  theme_pub() +
  theme(legend.position = c(0.82, 0.88))

save_pub(p_volcano, "Sec0_FigB_volcano")
cat("  Sec0_FigB: Volcano — OK\n")

# =============================================================================
# Sec0_FigC — HEATMAP TOP 40 PROTEÍNAS DE
# =============================================================================
cat("\n--- Sec0_FigC: Heatmap top DE ---\n")

de_detected <- de %>%
  mutate(n_detected = rowSums(!is.na(across(all_of(sample_cols))))) %>%
  filter(n_detected == length(sample_cols))

top_up   <- de_detected %>% filter(direction == "up")   %>%
  slice_max(logFC_TVsS, n = 20) %>% pull(gene_symbol)
top_down <- de_detected %>% filter(direction == "down") %>%
  slice_min(logFC_TVsS, n = 20) %>% pull(gene_symbol)
top_genes <- c(top_up, top_down)

cat(sprintf("  %d up + %d down elegibles (20/20 detección); seleccionados %d+%d\n",
            sum(de_detected$direction == "up"),
            sum(de_detected$direction == "down"),
            length(top_up), length(top_down)))

heat_mat <- expr_mat[top_genes[top_genes %in% rownames(expr_mat)], , drop = FALSE]
heat_z   <- t(scale(t(heat_mat)))

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

n_up_heat   <- sum(rownames(heat_z) %in% top_up)
n_down_heat <- sum(rownames(heat_z) %in% top_down)
cat(sprintf("  %d up + %d down en heatmap\n", n_up_heat, n_down_heat))
row_split <- factor(
  c(rep("Up-regulated", n_up_heat), rep("Down-regulated", n_down_heat)),
  levels = c("Up-regulated", "Down-regulated")
)

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
  column_title       = "Top 40 differentially expressed proteins — Tumor vs. Adjacent Normal",
  column_title_gp    = gpar(fontsize = 8, fontface = "bold"),
  rect_gp            = gpar(col = "grey90", lwd = 0.4),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5)
  )
)

save_ch(ht_de, "Sec0_FigC_heatmap_topDE", "double_col", h_add = 50)
cat("  Sec0_FigC: Heatmap top DE — OK\n")

# =============================================================================
# OE1_FigA — CANDIDATOS POR FUENTE DE BASE DE DATOS
# =============================================================================
cat("\n--- OE1_FigA: Drug sources bar ---\n")

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
  labs(title = "Drug candidates per database source",
       x = "Number of drugs", y = NULL) +
  theme_pub() +
  theme(legend.position = "none",
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

save_pub(p_sources, "OE1_FigA_drug_sources_bar")
cat("  OE1_FigA: Drug sources bar — OK\n")

# =============================================================================
# OE1_FigB — DISTRIBUCIÓN FASES CLÍNICAS (candidatos multi-fuente)
# =============================================================================
cat("\n--- OE1_FigB: Phase distribution ---\n")

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

cat(sprintf("  %d candidatos multi-fuente; %d aprobados (%.0f%%)\n",
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
    title    = "Clinical phase of multi-source drug candidates",
    subtitle = sprintf("n = %d candidatos con evidencia en \u22652 bases de datos  \u2022  Aprobado = Fase IV ChEMBL o aprobado en DGIdb",
                       nrow(multi_src)),
    x = NULL, y = "N\u00famero de f\u00e1rmacos"
  ) +
  theme_pub() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

save_pub(p_phase, "OE1_FigB_drug_phase_dist")
cat("  OE1_FigB: Phase distribution — OK\n")

# =============================================================================
# OE2_FigA — RED PPI DE PROTEÍNAS DE
# =============================================================================
cat("\n--- OE2_FigA: Red PPI ---\n")

# Helper para construir y graficar red PPI
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

hub_label_genes <- c(
  "EGFR", "PSMB3", "PSMA2", "MMP9", "CASP3", "ENO1",
  "NDUFS3", "NDUFS2", "NDUFV1", "ATP5F1C", "SDHA", "UQCRC2", "NDUFA9"
)

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
  layout_algo = "stress",
  title_sfx   = " — all DE proteins (degree > 8 or hub)"
)
save_pub(p_net_v1, "OE2_FigA_ppi_network", "double_col", w_add = 60, h_add = 90)
cat("  OE2_FigA: PPI network — OK\n")

# =============================================================================
# OE2_FigB — CANDIDATOS APROBADOS POR MÓDULO BIOLÓGICO
# =============================================================================
cat("\n--- OE2_FigB: Module barplot ---\n")

MODULE_LABELS <- c(
  "8"  = "OXPHOS",
  "4"  = "Proteasoma",
  "2"  = "EGFR / Señalización",
  "5"  = "Resp. inmune / ECM",
  "14" = "Metab. peq. mol."
)
KEY_MODULES <- c(8, 2, 4, 5, 14)

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
    name = "Drug class"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Drug candidates by biological module",
    subtitle = "Hub-targeting candidates  \u2022  Grouped by STRING PPI module",
    x = "N\u00b0 drug candidates", y = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y     = element_text(size = 7.5),
    legend.position = "right"
  )

save_pub(p_module_bar, "OE2_FigB_module_barplot", "double_col", h_add = 0)
cat("  OE2_FigB: Module barplot — OK\n")

# =============================================================================
# OE2_FigC — CLASIFICACIÓN CLÍNICO-REGULATORIA (A/B/C/D)
# =============================================================================
cat("\n--- OE2_FigC: Class distribution ---\n")

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
cat("  OE2_FigC: Class distribution — OK\n")

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
