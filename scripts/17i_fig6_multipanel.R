#!/usr/bin/env Rscript
# =============================================================================
# 17i_fig6_multipanel.R — Ensambla la Figura 6 multipanel (validación externa)
# =============================================================================
# Combina los paneles ggplot cacheados por 16_external_validation.R:
#   Fig6A_concordance   — concordancia global proteoma DIA vs TCGA RNA-seq
#   Fig6B_targets_tcga  — dianas priorizadas (shortlist) validadas en TCGA
#
# Layout: A izquierda (concordancia global) | B derecha (dianas priorizadas).
# Ancho relativo 1:1.2 para dar más espacio al panel B (lollipop + barras).
#
# Output (results/figures/pub/main/):
#   Fig6_multipanel.tif — TIFF 600 dpi LZW (submission) + .pdf + .png
#
# Narrativa de cierre: Fig2 (biología desregulada) → Fig5 (dianas priorizadas)
#   → Fig6 (esas dianas son reales y reproducibles en TCGA).
#
# Ambiente: omics-R
# Ejecución (DESPUÉS de 16_external_validation.R):
#   Rscript scripts/16_external_validation.R && \
#   Rscript scripts/17i_fig6_multipanel.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

cat("=== 17i_fig6_multipanel.R ===\n")
source(here::here("scripts", "_setup.R")); setup_project()
source(here::here("scripts", "_fig_style.R"))

obj_dir <- PANEL_OBJ_DIR
need    <- c("Fig6A_concordance", "Fig6B_targets_tcga")
paths   <- file.path(obj_dir, paste0(need, ".rds"))
missing <- need[!file.exists(paths)]
if (length(missing))
  stop("Faltan objetos de panel: ", paste(missing, collapse = ", "),
       "\n  Corré antes:  Rscript scripts/16_external_validation.R")

p_concord <- readRDS(file.path(obj_dir, "Fig6A_concordance.rds"))
p_targets <- readRDS(file.path(obj_dir, "Fig6B_targets_tcga.rds"))
cat("  Objetos de panel cargados OK\n")

# A (concordancia global) | B (dianas priorizadas)
# wrap_elements() evita que el patchwork interno p_targets propague sub-tags (C, D...)
fig6 <- p_concord + patchwork::wrap_elements(p_targets) +
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Dimensiones: más alta que Fig5 porque B tiene ~14 filas (lollipop)
save_tiff(fig6, "Fig6_multipanel", width_mm = 340, height_mm = 180)
ggsave("results/figures/pub/main/Fig6_multipanel.png", fig6,
       width = 340, height = 180, units = "mm", dpi = 300, limitsize = FALSE)
cat("  PNG de revisión: Fig6_multipanel.png\n")

cat("\nFig6_multipanel — OK\n")
cat("Fin:", format(Sys.time()), "\n")
