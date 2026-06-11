#!/usr/bin/env Rscript
# =============================================================================
# 17e_fig3_multipanel.R — Ensambla la Figura 3 multipanel (drug landscape)
# =============================================================================
# Combina los tres paneles ya cacheados como objetos por 17_pub_figures.R:
#   Fig3A_funnel         (ggplot)              — embudo de filtrado
#   Fig3B_drug_class     (ggplot)              — clase A/B/C/D sobre los 458
#   Fig3C_upset_overlap  (ComplexHeatmap/UpSet) — solapamiento de BD
#
# Lee los .rds de results/figures/pub/.objects/ (NO re-grafica). Estilo de
# scripts/_fig_style.R. El UpSet (grid) se captura como grob vía grid.grabExpr.
#
# Layout:
#   A (funnel)  |  B (class)
#   C (UpSet, full width)
#
# Output (results/figures/pub/main/):
#   Fig3_multipanel.tif  — TIFF 600 dpi LZW (submission)
#   Fig3_multipanel.pdf / .png
#
# Ambiente: omics-R
# Ejecución (DESPUÉS de 17):
#   Rscript scripts/17_pub_figures.R && Rscript scripts/17e_fig3_multipanel.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ComplexHeatmap)
})

cat("=== 17e_fig3_multipanel.R ===\n")
source(here::here("scripts", "_setup.R"))
setup_project()
source(here::here("scripts", "_fig_style.R"))

# ── Cargar objetos de panel cacheados ─────────────────────────────────────────
obj_dir <- PANEL_OBJ_DIR
need <- c("Fig3A_funnel", "Fig3B_drug_class", "Fig3C_upset_overlap")
paths <- file.path(obj_dir, paste0(need, ".rds"))
missing <- need[!file.exists(paths)]
if (length(missing)) {
  stop("Faltan objetos de panel: ", paste(missing, collapse = ", "),
       "\n  Corré antes:  Rscript scripts/17_pub_figures.R")
}

p_funnel <- readRDS(file.path(obj_dir, "Fig3A_funnel.rds"))
p_class  <- readRDS(file.path(obj_dir, "Fig3B_drug_class.rds"))
upset    <- readRDS(file.path(obj_dir, "Fig3C_upset_overlap.rds"))   # list(ht, m)
cat("  Objetos de panel cargados OK\n")

# ── UpSet (grid) → grob para componer con ggplot ──────────────────────────────
# draw_upset_panel (de _fig_style.R) dibuja el UpSet + números sobre las barras.
upset_grob <- grid.grabExpr(draw_upset_panel(upset$ht, upset$m),
                            width = unit(178, "mm"), height = unit(95, "mm"))
panel_C <- wrap_elements(full = upset_grob)

# ── Ensamblar (patchwork) ─────────────────────────────────────────────────────
# Orden de adición fija el tag: p_funnel→A, p_class→B, panel_C→C (UpSet).
design <- "AB
CC"

fig3 <- p_funnel + p_class + panel_C +
  plot_layout(design = design, widths = c(1.25, 1), heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# ── Export: TIFF 600 dpi LZW + PDF + PNG de revisión ──────────────────────────
save_tiff(fig3, "Fig3_multipanel", width_mm = 205, height_mm = 190)
ggsave("results/figures/pub/main/Fig3_multipanel.png", fig3,
       width = 205, height = 190, units = "mm", dpi = 300, limitsize = FALSE)
cat("  PNG de revisión: Fig3_multipanel.png\n")

cat("\nFig3_multipanel — OK\n")
cat("Fin:", format(Sys.time()), "\n")
