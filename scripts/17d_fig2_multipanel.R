#!/usr/bin/env Rscript
# =============================================================================
# 17d_fig2_multipanel.R — Ensambla la Figura 2 multipanel (publicación)
# =============================================================================
# Combina los tres paneles ya generados y cacheados como objetos por:
#   - 17_pub_figures.R   → Fig2A_volcano (ggplot), Fig2B_heatmap_topDE (Heatmap)
#   - 17c_pub_figD_gsea.R → Fig2C_hallmarks_gsea (ggplot, 11 hallmarks)
#
# Lee los .rds de results/figures/pub/.objects/ (NO re-grafica: evita
# divergencia de contenido). El estilo proviene de scripts/_fig_style.R.
#
# Layout (orden de lectura A→B→C; el heatmap a ancho completo para que su
# cuerpo no se colapse, y con más altura para que las etiquetas respiren):
#   A (volcano)  |  B (GSEA)
#   C (heatmap, full width)
#
# Output (results/figures/pub/main/):
#   Fig2_multipanel.tif  — TIFF 600 dpi LZW (submission)
#   Fig2_multipanel.pdf  — PDF vectorial (referencia)
#
# Ambiente: omics-R
# Ejecución (correr DESPUÉS de 17 y 17c):
#   conda activate omics-R
#   Rscript scripts/17_pub_figures.R && Rscript scripts/17c_pub_figD_gsea.R
#   Rscript scripts/17d_fig2_multipanel.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ComplexHeatmap)
})

cat("=== 17d_fig2_multipanel.R ===\n")
source(here::here("scripts", "_setup.R"))
setup_project()
source(here::here("scripts", "_fig_style.R"))

# ── Cargar objetos de panel cacheados ─────────────────────────────────────────
obj_dir <- PANEL_OBJ_DIR
need <- c("Fig2A_volcano", "Fig2B_heatmap_topDE", "Fig2C_hallmarks_gsea")
paths <- file.path(obj_dir, paste0(need, ".rds"))
missing <- need[!file.exists(paths)]
if (length(missing)) {
  stop("Faltan objetos de panel: ", paste(missing, collapse = ", "),
       "\n  Corré antes:  Rscript scripts/17_pub_figures.R && Rscript scripts/17c_pub_figD_gsea.R")
}

p_volcano <- readRDS(file.path(obj_dir, "Fig2A_volcano.rds"))
ht_de     <- readRDS(file.path(obj_dir, "Fig2B_heatmap_topDE.rds"))
p_gsea    <- readRDS(file.path(obj_dir, "Fig2C_hallmarks_gsea.rds"))
cat("  Objetos de panel cargados OK\n")

# ── Convertir el ComplexHeatmap a grob (única vía de componer con ggplot) ──────
# A ancho completo el cuerpo del heatmap se renderiza correctamente; leyendas a
# la derecha con merge_legend compacto.
heat_grob <- grid.grabExpr(
  draw(ht_de, merge_legend = TRUE,
       heatmap_legend_side = "bottom", annotation_legend_side = "bottom"),
  width = unit(205, "mm"), height = unit(152, "mm")
)
panel_C <- wrap_elements(full = heat_grob)

# ── Ensamblar (patchwork) ─────────────────────────────────────────────────────
# Orden de adición fija el tag: p_volcano→A, p_gsea→B, panel_C→C (heatmap).
# Fila 1: A (volcano, más ancho) | B (GSEA)   ·   Fila 2: C (heatmap, ancho completo, más alto)
design <- "AB
CC"

fig2 <- p_volcano + p_gsea + panel_C +
  plot_layout(design = design, widths = c(1.3, 1), heights = c(1, 1.75)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# ── Export: TIFF 600 dpi LZW + PDF ────────────────────────────────────────────
save_tiff(fig2, "Fig2_multipanel", width_mm = 215, height_mm = 255)

# PNG de revisión rápida (no es submission, pero útil en Claude Code)
ggsave("results/figures/pub/main/Fig2_multipanel.png", fig2,
       width = 215, height = 255, units = "mm", dpi = 300, limitsize = FALSE)
cat("  PNG de revisión: Fig2_multipanel.png\n")

cat("\nFig2_multipanel — OK\n")
cat("Fin:", format(Sys.time()), "\n")
