#!/usr/bin/env Rscript
# =============================================================================
# 17g_fig5_multipanel.R — Ensambla la Figura 5 multipanel (priorización)
# =============================================================================
# Combina los paneles ggplot cacheados por 17b_fig5_panels.R:
#   Fig5A_shortlist   — shortlist rankeado + descomposición TP/DV (faceta por tier)
#   Fig5B_tpdv_space  — espacio bifactor TargetPriority × DrugViability
#
# Layout: A ancho a la izquierda (la historia: candidatos priorizados por módulo),
#         B a la derecha (el espacio bifactor que explica la estructura del score).
#
# Output (results/figures/pub/main/):
#   Fig5_multipanel.tif — TIFF 600 dpi LZW (submission) + .pdf + .png
#
# Ambiente: omics-R
# Ejecución (DESPUÉS de 17b):
#   Rscript scripts/17b_fig5_panels.R && Rscript scripts/17g_fig5_multipanel.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork)
})

cat("=== 17g_fig5_multipanel.R ===\n")
source(here::here("scripts", "_setup.R")); setup_project()
source(here::here("scripts", "_fig_style.R"))

obj_dir <- PANEL_OBJ_DIR
need  <- c("Fig5A_shortlist", "Fig5B_tpdv_space")
paths <- file.path(obj_dir, paste0(need, ".rds"))
missing <- need[!file.exists(paths)]
if (length(missing))
  stop("Faltan objetos de panel: ", paste(missing, collapse = ", "),
       "\n  Corré antes:  Rscript scripts/17b_fig5_panels.R")

p_short <- readRDS(file.path(obj_dir, "Fig5A_shortlist.rds"))
p_space <- readRDS(file.path(obj_dir, "Fig5B_tpdv_space.rds"))
cat("  Objetos de panel cargados OK\n")

# A (shortlist) más ancho que B (espacio bifactor)
fig5 <- p_short + p_space +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

save_tiff(fig5, "Fig5_multipanel", width_mm = 340, height_mm = 165)
ggsave("results/figures/pub/main/Fig5_multipanel.png", fig5,
       width = 340, height = 165, units = "mm", dpi = 300, limitsize = FALSE)
cat("  PNG de revisión: Fig5_multipanel.png\n")

cat("\nFig5_multipanel — OK\n")
cat("Fin:", format(Sys.time()), "\n")
