#!/usr/bin/env Rscript
# =============================================================================
# 17f_fig4_multipanel.R — Ensambla la Figura 4 multipanel (network biology)
# =============================================================================
# Combina los tres paneles ggplot cacheados por 17_pub_figures.R:
#   Fig4A_hub_subnetwork    — subred de hubs druggables coloreada por módulo
#   Fig4B_module_barplot    — candidatos por módulo × clase
#   Fig4C_module_enrichment — enriquecimiento GO BP por módulo (top-3, simplify)
#
# Los 3 son ggplot → composición directa con patchwork (sin grid.grabExpr).
# Estilo de scripts/_fig_style.R.
#
# Layout narrativo de descubrimiento:
#   A (red, ancho completo): módulos coloreados SIN nombre (aún no se sabe qué son)
#   B (enrichment, abajo-izq): revela el nombre de cada color (leyenda Module + términos)
#   C (drogas, abajo-der): payoff de accionabilidad por módulo ya nombrado
#
# Output (results/figures/pub/main/):
#   Fig4_multipanel.tif  — TIFF 600 dpi LZW (submission) + .pdf + .png
#
# Ambiente: omics-R
# Ejecución (DESPUÉS de 17):
#   Rscript scripts/17_pub_figures.R && Rscript scripts/17f_fig4_multipanel.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

cat("=== 17f_fig4_multipanel.R ===\n")
source(here::here("scripts", "_setup.R"))
setup_project()
source(here::here("scripts", "_fig_style.R"))

# ── Cargar objetos de panel cacheados ─────────────────────────────────────────
obj_dir <- PANEL_OBJ_DIR
need <- c("Fig4A_network_modules", "Fig4B_module_barplot", "Fig4C_module_enrichment")
paths <- file.path(obj_dir, paste0(need, ".rds"))
missing <- need[!file.exists(paths)]
if (length(missing)) {
  stop("Faltan objetos de panel: ", paste(missing, collapse = ", "),
       "\n  Corré antes:  Rscript scripts/17_pub_figures.R")
}

p_net   <- readRDS(file.path(obj_dir, "Fig4A_network_modules.rds"))
p_bar   <- readRDS(file.path(obj_dir, "Fig4B_module_barplot.rds"))
p_enr   <- readRDS(file.path(obj_dir, "Fig4C_module_enrichment.rds"))
cat("  Objetos de panel cargados OK\n")

# Leyenda Module aparece en A (red) y C (enrichment) → la dejo SOLO en la red
# (mismos colores). El enrichment conserva solo la leyenda de tamaño (Genes).
# La leyenda de la red va abajo (horizontal) para no robarle ancho al panel.
# Narrativa de descubrimiento: la red (A) ya viene SIN leyenda de módulos (script 17);
# el enrichment (B) es el que nombra cada color.
p_enr <- p_enr + theme(legend.position = "right")    # aquí se nombran los módulos

# ── Ensamblar (patchwork) ─────────────────────────────────────────────────────
# Orden de adición fija el tag: p_net→A, p_enr→B, p_bar→C.
# A red ANCHO COMPLETO arriba (no se corta); B enrichment abajo-izq, C drogas abajo-der.
design <- "AA
BC"

fig4 <- p_net + p_enr + p_bar +
  plot_layout(design = design, widths = c(1, 1), heights = c(1.6, 0.75)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# ── Export: TIFF 600 dpi LZW + PDF + PNG de revisión ──────────────────────────
# Figura ancha (B/C respiran) y alta (la red no se aplasta vertical).
save_tiff(fig4, "Fig4_multipanel", width_mm = 340, height_mm = 235)
ggsave("results/figures/pub/main/Fig4_multipanel.png", fig4,
       width = 340, height = 235, units = "mm", dpi = 300, limitsize = FALSE)
cat("  PNG de revisión: Fig4_multipanel.png\n")

cat("\nFig4_multipanel — OK\n")
cat("Fin:", format(Sys.time()), "\n")
