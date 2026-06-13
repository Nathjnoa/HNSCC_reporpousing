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
#   Fig4A_hub_subnetwork   — Subred de hubs druggables coloreada por módulo
#   Fig4B_module_barplot   — Candidatos por módulo × clase (labels data-driven)
#   Fig4C_module_enrichment — Enriquecimiento GO BP por módulo (valida nombres)
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
# NOTA: clusterProfiler/org.Hs.eg.db NO se cargan con library() — entran en conflicto
# con el dibujo de ComplexHeatmap (C stack). Se usan con `::` solo en Fig4C (al final).

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
# FigS — FUNNEL DE SELECCIÓN (suplementario; overview del pipeline completo)
# =============================================================================
cat("\n--- FigS: Selection funnel (supp) ---\n")

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

# El funnel completo (3513→458→35→32) ya NO es panel de Fig3: con un solo filtro
# hasta multi-source no amerita embudo, y las etapas top-ranked/LOD pertenecen a la
# priorización (Fig5). Se conserva como suplementario (overview del pipeline).
save_pub(p_funnel, "FigS_selection_funnel", out_dir = "results/figures/pub/supp")
cat("  FigS (supp): Selection funnel — OK\n")

# =============================================================================
# Fig3A — DISTRIBUCIÓN DE FASE CLÍNICA (candidatos multi-fuente, n=458)
# =============================================================================
# Reemplaza al funnel como panel A: describe el MISMO set de 458 (igual que B y C),
# y prepara la accionabilidad clínica sin adelantar la priorización (Fig5).
cat("\n--- Fig3A: Phase distribution ---\n")

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

save_pub(p_phase, "Fig3A_phase")
save_panel_obj(p_phase, "Fig3A_phase")
cat("  Fig3A: Phase distribution — OK\n")

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
  scale_x_continuous(expand = expansion(mult = c(0, 0.42))) +   # aire para labels "n (xx%)"
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
# CLASIFICACIÓN DE MÓDULOS POR DROGABILIDAD (data-driven; reemplaza selección manual)
# =============================================================================
# Cada módulo Louvain del giant component se clasifica por evidencia drogable:
#   approved  : ≥1 fármaco APROBADO dirigido a un nodo del módulo  → color propio
#   hubs_only : 0 aprobados pero ≥2 hubs druggables                → color flag (slate)
#   none      : sin potencial drogable                             → gris (contexto)
# Nombres de display: 5 curados + término GO BP top para el resto.
# M2 NO es "EGFR/Signaling" (su GO es adhesión/membrana); EGFR es un hub interno.
.gA  <- igraph::graph_from_data_frame(net_edges, directed = FALSE,
                                      vertices = net_nodes$gene_symbol)
.cmp <- igraph::components(.gA)
giant_set <- names(.cmp$membership[.cmp$membership == which.max(.cmp$csize)])

.mod_appr <- master_tbl %>%
  filter(gene_symbol %in% giant_set) %>%
  left_join(net_nodes %>% select(gene_symbol, module_id), by = "gene_symbol") %>%
  group_by(module_id) %>%
  summarise(n_approved = n_distinct(drug_name_norm[is_approved %in% c(TRUE, "True", "TRUE")]),
            .groups = "drop")
.mod_hubs <- drg_hubs %>% count(module_id, name = "n_drug_hubs")

MODULE_NAMES <- c(
  "8" = "OXPHOS", "4" = "Proteasome", "5" = "Immune response",
  "14" = "Amino acid metabolism", "2" = "Cell adhesion / membrane",
  "10" = "mRNA processing", "11" = "Muscle contraction",
  "6" = "Carboxylic acid metabolism", "7" = "Carbohydrate catabolism",
  "1" = "ECM / development", "3" = "Translation",
  "17" = "Chromatin remodeling", "13" = "Microtubule process"
)

module_class <- net_nodes %>% filter(gene_symbol %in% giant_set) %>%
  count(module_id, name = "n_nodes") %>%
  left_join(.mod_appr, by = "module_id") %>%
  left_join(.mod_hubs, by = "module_id") %>%
  mutate(
    n_approved  = tidyr::replace_na(n_approved, 0L),
    n_drug_hubs = tidyr::replace_na(n_drug_hubs, 0L),
    tier = dplyr::case_when(n_approved >= 1   ~ "approved",
                            n_drug_hubs >= 2  ~ "hubs_only",
                            TRUE              ~ "none"),
    name = MODULE_NAMES[as.character(module_id)]
  ) %>%
  arrange(dplyr::desc(tier == "approved"), dplyr::desc(n_nodes))

# Paleta colorblind-safe (Carto "Safe", 11 colores) para módulos con aprobados;
# hubs_only = slate (flag visual); none/contexto = gris.
SAFE11 <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499",
            "#44AA99","#999933","#882255","#661100","#6699CC")
.appr <- module_class %>% filter(tier == "approved")  %>% pull(name)
.huo  <- module_class %>% filter(tier == "hubs_only") %>% pull(name)
MODULE_COLS <- setNames(SAFE11[seq_along(.appr)], .appr)
for (nm in .huo) MODULE_COLS[nm] <- "#7D8CA3"   # slate: hubs druggables, sin aprobado
MODULE_COLS["Other"] <- "#CFCFCF"               # none + contexto

# Etiqueta de módulo para plotear: name si drogable, "Other" si tier == none
module_lab_of <- function(mid) {
  t <- module_class$tier[match(mid, module_class$module_id)]
  ifelse(is.na(t) | t == "none", "Other",
         module_class$name[match(mid, module_class$module_id)])
}
# Orden B/C y niveles de factor: aprobados (por tamaño) → hubs_only ; none excluido
MODULE_ORDER    <- module_class %>% filter(tier != "none") %>% pull(name)
COLORED_MODULES <- module_class %>% filter(tier != "none") %>% pull(module_id)
KEY_MODULES     <- COLORED_MODULES   # compat con paneles B/C
cat(sprintf("  Módulos drogables coloreados: %d | hubs_only: %s | gris(none): %s\n",
            length(COLORED_MODULES), paste(.huo, collapse=","),
            paste(module_class$name[module_class$tier=="none"], collapse=",")))

# =============================================================================
# Fig4A — RED COMPLETA (giant component) COLOREADA POR MÓDULO
# =============================================================================
# Layout sobre la red donde se DEFINIERON los módulos (Louvain, script 09) → las
# comunidades aparecen como clusters espaciales reales. No-hubs en gris/tenue;
# hubs druggables resaltados (tamaño + borde) y etiquetados (top-3 por degree).
cat("\n--- Fig4A: Full network colored by module ---\n")

hub_set <- drg_hubs$gene_symbol
g_all <- graph_from_data_frame(net_edges, directed = FALSE,
                               vertices = net_nodes$gene_symbol)
comp   <- components(g_all)
giant  <- names(comp$membership[comp$membership == which.max(comp$csize)])

n_all <- net_nodes %>%
  filter(gene_symbol %in% giant) %>%
  mutate(
    module_lab  = module_lab_of(module_id),
    module_lab  = factor(module_lab, levels = c(MODULE_ORDER, "Other")),
    is_drug_hub = gene_symbol %in% hub_set,
    nd_size     = ifelse(is_drug_hub, log1p(degree) + 1, log1p(degree))
  )
e_all <- net_edges %>% filter(gene_A %in% giant, gene_B %in% giant)
cat(sprintf("  Giant component: %d nodos, %d aristas | %d hubs druggables\n",
            nrow(n_all), nrow(e_all), sum(n_all$is_drug_hub)))

# Labels: top-2 druggable hubs por módulo coloreado, por DEGREE (con 12 módulos,
# 3 labels c/u satura; 2 mantiene legibilidad)
hub_label_genes <- drg_hubs %>%
  filter(module_id %in% COLORED_MODULES) %>%
  group_by(module_id) %>%
  slice_max(degree, n = 2, with_ties = FALSE) %>%
  ungroup() %>%
  pull(gene_symbol)

g_net <- graph_from_data_frame(e_all, directed = FALSE, vertices = n_all$gene_symbol)
for (col in c("module_lab", "is_drug_hub", "nd_size", "degree")) {
  vertex_attr(g_net, col) <- n_all[[col]][match(V(g_net)$name, n_all$gene_symbol)]
}
tg_net <- as_tbl_graph(g_net)

# Layout stress con PESOS intra-módulo reforzados: las aristas dentro de un mismo
# módulo Louvain pesan más → cada comunidad se compacta y se separa de las otras.
# Sin esto, los 2 módulos más grandes y centrales (Immune, Cell adhesion) se
# superponen en el núcleo denso. No altera la topología; enfatiza la estructura
# que Louvain ya definió. (graphlayouts::layout_as_backbone requeriría 'oaqc'.)
el_idx   <- igraph::as_edgelist(g_net, names = FALSE)
mod_vec  <- n_all$module_id[match(V(g_net)$name, n_all$gene_symbol)]
same_mod <- mod_vec[el_idx[, 1]] == mod_vec[el_idx[, 2]]
net_w    <- ifelse(same_mod, 8, 1)   # intra-módulo 8x → separa comunidades

set.seed(42)
lay_net <- create_layout(tg_net, layout = "stress", weights = net_w)

# Estiramiento x suave para llenar el panel ancho (full-width arriba en multipanel).
NET_XSTRETCH <- 1.2
lay_net$x <- lay_net$x * NET_XSTRETCH

lbl_net <- lay_net %>%
  filter(name %in% hub_label_genes[hub_label_genes %in% giant]) %>%
  select(x, y, name)

p_net <- ggraph(lay_net) +
  geom_edge_link(color = "grey60", alpha = 0.30, linewidth = 0.14) +
  # Tamaño por CATEGORÍA (no por degree): brecha CHICA. El contexto gris
  # periférico sube de tamaño; los coloreados (módulos clave) quedan ~original.
  # contexto (módulos no-clave): gris VISIBLE, ahora MÁS GRANDE (menos brecha)
  # OJO: shape 21 requiere `color` no-NA; con color=NA ggplot ELIMINA las filas
  # (borra el nodo). Usar "transparent" = borde invisible pero fila conservada.
  geom_node_point(aes(fill = module_lab,
                      filter = !is_drug_hub & module_lab == "Other"),
                  shape = 21, color = "transparent", alpha = 0.72, size = 3.5) +
  # miembros no-hub de módulos clave: color saturado
  geom_node_point(aes(fill = module_lab,
                      filter = !is_drug_hub & module_lab != "Other"),
                  shape = 21, color = "transparent", alpha = 0.90, size = 4.5) +
  # hubs druggables: protagonistas (borde + opacos)
  geom_node_point(aes(fill = module_lab, filter = is_drug_hub),
                  shape = 21, color = "grey20", stroke = 0.4, alpha = 0.95, size = 5.5) +
  ggrepel::geom_text_repel(
    data = lbl_net, aes(x = x, y = y, label = name),
    size = 3.3, fontface = "bold", family = "sans",
    bg.color = "white", bg.r = 0.12, max.overlaps = 40,
    segment.size = 0.25, segment.color = "grey40", min.segment.length = 0
  ) +
  # Leyenda de módulos VISIBLE: con 12 comunidades drogables, nombrarlas en la red.
  scale_fill_manual(values = MODULE_COLS, name = "Module (Louvain)",
                    breaks = c(MODULE_ORDER, "Other"),
                    labels = c(MODULE_ORDER, "Other (below threshold)"),
                    guide = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  theme_graph(base_family = "sans", base_size = 8) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 11),
        legend.key.height = unit(5.5, "mm"),
        plot.margin = margin(4, 4, 4, 4, "mm"))

save_pub(p_net, "Fig4A_network_modules", "double_col", w_add = 60, h_add = 80)
save_panel_obj(p_net, "Fig4A_network_modules")
cat("  Fig4A: Full network by module — OK\n")

# =============================================================================
# Fig4B — CANDIDATOS POR MÓDULO BIOLÓGICO (× clase)
# =============================================================================
cat("\n--- Fig4B: Module barplot ---\n")

# Métrica ESTRUCTURAL (Fig4 es pre-priorización): nº de hubs druggables por módulo,
# coloreado POR MÓDULO (coherente con A). La dimensión clínica (clase de fármaco) y
# la priorización se reservan para Fig5 — aquí solo se cuantifica la estructura.
hubs_per_module <- drg_hubs %>%
  filter(module_id %in% COLORED_MODULES) %>%
  mutate(module_label = MODULE_NAMES[as.character(module_id)]) %>%
  distinct(gene_symbol, module_label) %>%
  count(module_label, name = "n_hubs") %>%
  mutate(module_label = factor(module_label, levels = rev(MODULE_ORDER)))

p_module_bar <- ggplot(hubs_per_module,
                       aes(y = module_label, x = n_hubs, fill = module_label)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = MODULE_COLS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = NULL, x = "Druggable network hubs", y = NULL) +
  theme_pub() +
  theme(legend.position = "none")

save_pub(p_module_bar, "Fig4B_module_barplot", "double_col", h_add = 0)
save_panel_obj(p_module_bar, "Fig4B_module_barplot")
cat("  Fig4B: Module barplot — OK\n")

# =============================================================================
# Fig4C — ENRIQUECIMIENTO FUNCIONAL POR MÓDULO (valida los nombres de módulo)
# =============================================================================
cat("\n--- Fig4C: Per-module functional enrichment ---\n")

orgdb <- org.Hs.eg.db::org.Hs.eg.db   # `::` evita attach (no rompe ComplexHeatmap)
mod_tbl <- read_tsv("results/tables/network/09_modules.tsv", show_col_types = FALSE)
universe_entrez <- clusterProfiler::bitr(unique(mod_tbl$gene_symbol), "SYMBOL",
                                         "ENTREZID", OrgDb = orgdb)$ENTREZID

enr_list <- lapply(KEY_MODULES, function(mid) {
  genes <- mod_tbl$gene_symbol[mod_tbl$module_id == mid]
  eg <- clusterProfiler::bitr(genes, "SYMBOL", "ENTREZID", OrgDb = orgdb)$ENTREZID
  ego <- clusterProfiler::enrichGO(eg, OrgDb = orgdb, ont = "BP",
                  universe = universe_entrez, pAdjustMethod = "BH",
                  qvalueCutoff = 0.2, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  df_raw <- as.data.frame(ego)   # pre-simplify (para fallback por keyword)
  # Reducción de redundancia semántica (GO BP es jerárquico): conserva el
  # representante más significativo de cada grupo de términos similares.
  if (nrow(df_raw) > 1) {
    ego <- clusterProfiler::simplify(ego, cutoff = 0.6,
                                     by = "p.adjust", select_fun = min)
  }
  df <- as.data.frame(ego)[order(as.data.frame(ego)$p.adjust), ]
  # Validación de nombre: para el módulo curado de aminoácidos, mostrar el término
  # ESPECÍFICO de aminoácidos (no el padre genérico que sobrevive a simplify).
  if (MODULE_NAMES[as.character(mid)] == "Amino acid metabolism") {
    aa <- df_raw[order(df_raw$p.adjust), ]
    aa <- aa[grepl("amino acid", aa$Description, ignore.case = TRUE), ]
    if (nrow(aa) > 0) df <- aa
  }
  # Top-1: el término más significativo = el que da NOMBRE al módulo (valida el label)
  df <- df[seq_len(min(1, nrow(df))), ]
  df$module_label <- MODULE_NAMES[as.character(mid)]
  df$Count <- as.integer(df$Count)
  df$GeneRatio <- df$Count / length(genes)   # fracción del módulo en el término
  df[, c("module_label", "Description", "p.adjust", "Count", "GeneRatio")]
})
enr_df <- bind_rows(enr_list) %>%
  mutate(
    module_label = factor(module_label, levels = rev(MODULE_ORDER)),
    Description  = stringr::str_wrap(stringr::str_to_sentence(Description), 30),
    neglog10FDR  = -log10(p.adjust)
  ) %>%
  arrange(module_label, GeneRatio) %>%
  mutate(term = factor(Description, levels = unique(Description)))

write_tsv(enr_df, "results/tables/network/17_module_enrichment_top3.tsv")
cat(sprintf("  %d terminos (top-1 x %d modulos) exportados\n",
            nrow(enr_df), length(KEY_MODULES)))

p_enr <- ggplot(enr_df, aes(x = GeneRatio, y = term,
                            color = module_label, size = Count)) +
  geom_point() +
  scale_color_manual(values = MODULE_COLS, name = "Module",
                     guide = guide_legend(override.aes = list(size = 3.5))) +
  scale_size_continuous(name = "Genes", range = c(2, 6)) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = NULL, x = "Gene ratio (term ∩ module / module)", y = NULL) +
  theme_pub() +
  theme(legend.position = "right")

save_pub(p_enr, "Fig4C_module_enrichment", "double_col", h_add = 0)
save_panel_obj(p_enr, "Fig4C_module_enrichment")
cat("  Fig4C: Per-module enrichment — OK\n")

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
