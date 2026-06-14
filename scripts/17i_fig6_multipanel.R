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
need    <- c("Fig6A_cptac_concordance", "Fig6B_tcga_concordance", "Fig6C_targets_unified")
paths   <- file.path(obj_dir, paste0(need, ".rds"))
missing <- need[!file.exists(paths)]
if (length(missing))
  stop("Faltan objetos de panel: ", paste(missing, collapse = ", "),
       "\n  Corré antes:  Rscript scripts/16_external_validation.R")

p_cptac   <- readRDS(file.path(obj_dir, "Fig6A_cptac_concordance.rds"))
p_tcga    <- readRDS(file.path(obj_dir, "Fig6B_tcga_concordance.rds"))
p_targets <- readRDS(file.path(obj_dir, "Fig6C_targets_unified.rds"))
cat("  Objetos de panel cargados OK\n")

# Layout: A+B en fila superior | C en fila inferior (ancho completo)
# patchwork design: A ocupa col1 fila1, B col2 fila1, C col1+col2 fila2.
# wrap_elements() en C para tag "C" visible y sin propagación de sub-tags.
design <- "AB
           CC"

fig6 <- p_cptac + p_tcga + patchwork::wrap_elements(p_targets) +
  plot_layout(design = design, heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# 340mm ancho × 280mm alto — A/B cuadradas arriba, C panorámica abajo
save_tiff(fig6, "Fig6_multipanel", width_mm = 340, height_mm = 280)
ggsave("results/figures/pub/main/Fig6_multipanel.png", fig6,
       width = 340, height = 280, units = "mm", dpi = 300, limitsize = FALSE)
cat("  PNG de revisión: Fig6_multipanel.png\n")

cat("\nFig6_multipanel — OK\n")
cat("Fin:", format(Sys.time()), "\n")
