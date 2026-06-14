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

# A (CPTAC scatter) | B (TCGA scatter) | C (targets unificado)
# Sin wrap_elements(): C se estira al alto completo de la figura.
# tag_level='new' en el layout interno de p_targets (script 16) evita que las
# sub-etiquetas se propaguen al patchwork externo.
fig6 <- p_cptac + p_tcga + p_targets +
  plot_layout(widths = c(1, 1, 2.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# 540mm ancho; 175mm alto — el eje discreto de C rellena el espacio vertical
save_tiff(fig6, "Fig6_multipanel", width_mm = 560, height_mm = 175)
ggsave("results/figures/pub/main/Fig6_multipanel.png", fig6,
       width = 560, height = 175, units = "mm", dpi = 300, limitsize = FALSE)
cat("  PNG de revisión: Fig6_multipanel.png\n")

cat("\nFig6_multipanel — OK\n")
cat("Fin:", format(Sys.time()), "\n")
