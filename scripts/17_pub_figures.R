#!/usr/bin/env Rscript
# =============================================================================
# 17_pub_figures.R
# =============================================================================
# Genera las 7 figuras de publicación incluidas en el reporte revisado
# (HNSCC_DrugRepurposing_Figuras_reviewed.docx).
#
# Figuras generadas (results/figures/pub/main/):
#   Fig2A_volcano          — Volcano plot Tumor vs. Adjacent Normal
#   Fig2B_heatmap_topDE    — Heatmap top 40 proteínas DE × 20 muestras
#   Fig2C_hallmarks_gsea   — Hallmarks GSEA dotplot (sobrescrito por 17c; versión canónica allí)
#   Fig3A_funnel           — Embudo de filtrado de candidatos (3513 → 458 → 35 → 32)
#   Fig3B_drug_class       — Clase clínico-regulatoria A/B/C/D sobre los 458 multi-fuente
#   Fig4A_ppi_network      — Red PPI de proteínas DE (degree > 8 | hub)
#   Fig4B_module_barplot   — Candidatos aprobados por módulo PPI
#   (supp) FigS_drug_phase_dist — Fase clínica de los 458 multi-fuente (suplementario)
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

# ── Working directory (raíz del proyecto vía scripts/_setup.R) ───────────────
cat("=== 17_pub_figures.R ===\n")
source(here::here("scripts", "_setup.R"))
setup_project()
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
# SISTEMA DE ESTILOS CENTRALIZADO (fuente única: scripts/_fig_style.R)
# =============================================================================
# Provee PRESETS, OKB, DE_COLS, COND_COLS, HPV_COLS, PHASE_COLS,
# DRUG_CLASS_LABELS, theme_pub(), save_pub(), save_ch(), save_tiff(),
# save_panel_obj().  save_pub/save_ch escriben en results/figures/pub/main
# por defecto (out_dir).
source(here::here("scripts", "_fig_style.R"))

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
top_ranked   <- read_tsv("results/tables/10_top20_candidates.tsv",
                        show_col_types = FALSE)
lod_stab     <- read_tsv("results/tables/15_lod_stability.tsv",
                        show_col_types = FALSE)

cat("  Datos cargados OK\n")

# ── Cargar metadata (HPV status) ──────────────────────────────────────────────
meta <- read_delim("data/raw/metadata.csv", delim = ";", show_col_types = FALSE) %>%
  mutate(hpv = ifelse(tolower(vph) == "positive", "HPV+", "HPV-"))
hpv_lookup <- setNames(meta$hpv, meta$sample_id)
# HPV_COLS proviene de _fig_style.R
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
# Fig2A — VOLCANO PLOT (Tumor vs. Adjacent Normal)
# =============================================================================
cat("\n--- Fig2A: Volcano plot ---\n")

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
  geom_text_repel(
    data = top_vol, aes(label = gene_symbol),
    size = 2.3, max.overlaps = 20, segment.size = 0.25,
    min.segment.length = 0.2, show.legend = FALSE,
    bg.color = "white", bg.r = 0.12          # halo en vez de recuadro
  ) +
  scale_color_manual(
    values = DE_COLS,
    labels = c(up   = sprintf("Up (%d)",   n_up),
               down = sprintf("Down (%d)", n_down),
               ns   = sprintf("NS (%d)",   n_ns))
  ) +
  scale_alpha_manual(values = c(up = 0.85, down = 0.85, ns = 0.25),
                     guide = "none") +
  labs(title = NULL,
       x     = expression(log[2]~"fold change"),
       y     = expression(-log[10]~"FDR"),
       color = NULL) +
  theme_pub() +
  theme(legend.position  = c(0.22, 0.15),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(2.2, "mm"),
        legend.background = element_rect(fill = alpha("white", 0.6), color = NA))

save_pub(p_volcano, "Fig2A_volcano")
save_panel_obj(p_volcano, "Fig2A_volcano")   # insumo multipanel
cat("  Fig2A: Volcano — OK\n")

# =============================================================================
# Fig2B — HEATMAP TOP 40 PROTEÍNAS DE
# =============================================================================
cat("\n--- Fig2B: Heatmap top DE ---\n")

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
  annotation_name_gp = gpar(fontsize = 8.5),
  simple_anno_size   = unit(3, "mm"),
  annotation_legend_param = list(
    Condition = list(title_gp = gpar(fontsize = 9, fontface = "bold"),
                     labels_gp = gpar(fontsize = 8.5)),
    HPV       = list(title_gp = gpar(fontsize = 9, fontface = "bold"),
                     labels_gp = gpar(fontsize = 8.5))
  )
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
  column_names_gp    = gpar(fontsize = 7.2),
  row_names_gp       = gpar(fontsize = 8.5),
  row_title_gp       = gpar(fontsize = 10.5, fontface = "bold"),
  column_title       = NULL,
  rect_gp            = gpar(col = "grey90", lwd = 0.4),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8.5)
  )
)

save_ch(ht_de, "Fig2B_heatmap_topDE", "double_col", h_add = 50)
save_panel_obj(ht_de, "Fig2B_heatmap_topDE")   # insumo multipanel
cat("  Fig2B: Heatmap top DE — OK\n")

# =============================================================================
# Fig3A — FUNNEL DE FILTRADO DE CANDIDATOS
# =============================================================================
cat("\n--- Fig3A: Selection funnel ---\n")

n_total <- nrow(drug_sum)
n_multi <- nrow(multi_src)
n_top   <- nrow(top_ranked)
n_lod   <- sum(lod_stab$lod_stable == TRUE, na.rm = TRUE)
cat(sprintf("  Funnel: %d -> %d -> %d -> %d\n", n_total, n_multi, n_top, n_lod))

funnel_df <- tibble(
  stage = c("Candidate drugs", "Multi-source candidates",
            "Top-ranked", "LOD-stable panel"),
  n     = c(n_total, n_multi, n_top, n_lod),
  note  = c("across 4 databases",
            "≥2 DB or approved", "composite score", "final panel")
)
nstg <- nrow(funnel_df)

# Anchos decorativos (decrecientes pero suaves: el texto va DENTRO de cada banda;
# bandas anchas para que las etiquetas no se desborden al escalar en el multipanel)
w  <- c(1.00, 0.88, 0.78, 0.72)
wn <- c(w, w[nstg])               # ancho inferior del último = su propio ancho

# Polígonos (trapecios encadenados → silueta de embudo)
funnel_poly <- do.call(rbind, lapply(seq_len(nstg), function(k) {
  ytop <- nstg - (k - 1); ybot <- nstg - k
  data.frame(
    id    = k,
    stage = funnel_df$stage[k],
    x     = c(-w[k] / 2, w[k] / 2, wn[k + 1] / 2, -wn[k + 1] / 2),
    y     = c(ytop, ytop, ybot, ybot)
  )
}))
funnel_poly$stage <- factor(funnel_poly$stage, levels = funnel_df$stage)

funnel_lab <- funnel_df %>%
  mutate(ymid = nstg - row_number() + 0.5,
         # texto uniforme negro: la paleta se aclaró para que el contraste sea correcto
         txtcol = "grey15")

FUNNEL_COLS <- c(
  "Candidate drugs"         = "#D6E6F2",
  "Multi-source candidates" = "#A9CCE3",
  "Top-ranked"              = "#7FB3D5",
  "LOD-stable panel"        = "#F5B36B"   # acento naranja claro: resalta el panel final
)                                          # paleta aclarada → texto negro legible en las 4

p_funnel <- ggplot() +
  geom_polygon(data = funnel_poly,
               aes(x = x, y = y, group = id, fill = stage),
               color = "white", linewidth = 0.5) +
  # Etiquetas DENTRO de cada banda (stage · conteo · nota), centradas.
  # family="sans" explícito para igualar la tipografía de theme_pub() (otros paneles).
  # color = txtcol → contraste adaptativo (scale_color_identity).
  geom_text(data = funnel_lab, aes(x = 0, y = ymid + 0.19, label = stage, color = txtcol),
            fontface = "bold", family = "sans", size = 3.1) +
  geom_text(data = funnel_lab, aes(x = 0, y = ymid - 0.07, label = scales::comma(n), color = txtcol),
            fontface = "bold", family = "sans", size = 4.4) +
  geom_text(data = funnel_lab, aes(x = 0, y = ymid - 0.32, label = note, color = txtcol),
            family = "sans", size = 2.4) +
  scale_fill_manual(values = FUNNEL_COLS) +
  scale_color_identity() +
  coord_cartesian(xlim = c(-0.55, 0.55), clip = "off") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_void(base_family = "sans") +
  theme(legend.position = "none",
        plot.margin = margin(6, 8, 6, 8, "mm"))

save_pub(p_funnel, "Fig3A_funnel")
save_panel_obj(p_funnel, "Fig3A_funnel")
cat("  Fig3A: Selection funnel — OK\n")

# =============================================================================
# Fig3B — DISTRIBUCIÓN FASES CLÍNICAS (candidatos multi-fuente)
# =============================================================================
cat("\n--- Fig3B: Phase distribution ---\n")

phase_df <- multi_src %>%
  mutate(
    phase_label = case_when(
      max_phase == 4                         ~ "Approved",
      is_approved == TRUE & is.na(max_phase) ~ "Approved",
      max_phase == 3                         ~ "Phase III",
      max_phase == 2                         ~ "Phase II",
      max_phase == 1                         ~ "Phase I",
      TRUE                                   ~ "No data"
    ),
    phase_label = factor(phase_label, levels = names(PHASE_COLS))
  ) %>%
  count(phase_label) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat(sprintf("  %d multi-source candidates; %d approved (%.0f%%)\n",
            nrow(multi_src),
            sum(phase_df$n[phase_df$phase_label == "Approved"]),
            sum(phase_df$pct[phase_df$phase_label == "Approved"])))

p_phase <- ggplot(phase_df,
                  aes(x = phase_label, y = n, fill = phase_label)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            vjust = -0.3, size = 2.6, lineheight = 0.9) +
  scale_fill_manual(values = PHASE_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(
    title    = NULL,
    x = NULL, y = "Number of drugs"
  ) +
  theme_pub() +
  theme(legend.position = "none")   # labels de eje horizontales (estándar uniforme)

save_pub(p_phase, "FigS_drug_phase_dist", out_dir = "results/figures/pub/supp")
cat("  FigS (supp): Phase distribution — OK\n")

# =============================================================================
# Fig3B — CLASE CLINICO-REGULATORIA (A/B/C/D) sobre los 458 multi-fuente
# =============================================================================
cat("\n--- Fig3C: Drug class distribution (multi-source) ---\n")

class_df <- multi_src %>%
  count(drug_class, name = "n_drugs") %>%
  mutate(
    class_label = DRUG_CLASS_LABELS[drug_class],
    class_label = factor(class_label, levels = DRUG_CLASS_LABELS),
    pct         = round(n_drugs / sum(n_drugs) * 100, 1)
  )

cat(sprintf("  Clases sobre %d multi-fuente: %s\n", nrow(multi_src),
            paste(sprintf("%s=%d(%.0f%%)", class_df$drug_class,
                          class_df$n_drugs, class_df$pct), collapse = " ")))

CLASS_COLS <- setNames(OKB[1:4], DRUG_CLASS_LABELS)

# Barras horizontales: categorías largas caben en el eje Y (sin solape, labels
# horizontales como el panel C → orientación de eje uniforme).
class_df <- class_df %>%
  mutate(class_label = factor(class_label, levels = rev(DRUG_CLASS_LABELS)))

p_class <- ggplot(class_df,
                  aes(x = n_drugs, y = class_label, fill = class_label)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", n_drugs, pct)),
            hjust = -0.15, size = 2.6) +
  scale_fill_manual(values = CLASS_COLS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(title = NULL, x = "Number of drugs", y = NULL) +
  theme_pub()

save_pub(p_class, "Fig3B_drug_class")
save_panel_obj(p_class, "Fig3B_drug_class")
cat("  Fig3B: Drug class distribution — OK\n")

# =============================================================================
# Fig3C — UpSet: solapamiento entre las 4 bases de datos (justifica filtro >=2)
# =============================================================================
cat("\n--- Fig3C: UpSet source overlap ---\n")

# Construido sobre los 458 candidatos multi-fuente (mismo set que el funnel),
# para consistencia entre paneles. Los 3 singletons son aprobados (clase A/B)
# conservados por la regla "n_sources>=2 OR aprobado" del script 08.
src_sets <- list(
  DGIdb       = multi_src$drug_name_norm[str_detect(multi_src$sources, "DGIdb")],
  ChEMBL      = multi_src$drug_name_norm[str_detect(multi_src$sources, "ChEMBL")],
  OpenTargets = multi_src$drug_name_norm[str_detect(multi_src$sources, "OpenTargets")],
  L2S2        = multi_src$drug_name_norm[str_detect(multi_src$sources, "L2S2")]
)
m <- make_comb_mat(src_sets)
cat(sprintf("  Combinaciones: %d | drogas (suma = multi_src): %d\n",
            length(comb_size(m)), sum(comb_size(m))))

ht_upset <- UpSet(
  m,
  comb_order = order(-comb_size(m)),
  set_order  = order(-set_size(m)),
  pt_size    = unit(2.6, "mm"), lwd = 1.3,
  comb_col   = "#0072B2",
  # Top annotation manual: nombre "Intersection size" SIN \n (upset_top_annotation
  # lo hardcodea en dos líneas), rotado 90° como título de eje Y de una sola línea.
  top_annotation = HeatmapAnnotation(
    "Intersection size" = anno_barplot(
      comb_size(m), gp = gpar(fill = "#0072B2", col = NA),
      border = FALSE,                       # sin recuadro alrededor del barplot
      height = unit(42, "mm"),
      axis_param = list(gp = gpar(fontsize = 7))
    ),
    annotation_name_gp   = gpar(fontsize = 11),   # = título de eje theme_pub
    annotation_name_side = "left",
    annotation_name_rot  = 90
  ),
  right_annotation = upset_right_annotation(
    m, gp = gpar(fill = "grey55", col = NA),
    width = unit(16, "mm"),
    annotation_name_gp = gpar(fontsize = 11)
  ),
  row_names_gp = gpar(fontsize = 8.5)
)

save_upset(ht_upset, m, "Fig3C_upset_overlap", "double_col", h_add = -30)
save_panel_obj(list(ht = ht_upset, m = m), "Fig3C_upset_overlap")  # ht + m para números en multipanel
cat("  Fig3C: UpSet source overlap — OK\n")

# =============================================================================
# Fig4A — RED PPI DE PROTEÍNAS DE
# =============================================================================
cat("\n--- Fig4A: Red PPI ---\n")

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
save_pub(p_net_v1, "Fig4A_ppi_network", "double_col", w_add = 60, h_add = 90)
cat("  Fig4A: PPI network — OK\n")

# =============================================================================
# Fig4B — CANDIDATOS APROBADOS POR MÓDULO BIOLÓGICO
# =============================================================================
cat("\n--- Fig4B: Module barplot ---\n")

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
  theme(legend.position = "right")   # fuente de ejes uniforme (theme_pub)

save_pub(p_module_bar, "Fig4B_module_barplot", "double_col", h_add = 0)
cat("  Fig4B: Module barplot — OK\n")

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
