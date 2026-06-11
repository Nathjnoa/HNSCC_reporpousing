#!/usr/bin/env Rscript
# =============================================================================
# 17c_pub_figD_gsea.R — Fig2C: Hallmarks GSEA dotplot (publication figure)
# =============================================================================
# Standalone script: only requires 03_Hallmarks_GSEA.tsv from script 03.
# Does NOT require intermediate pipeline files (network, drug_targets).
#
# Input:  results/tables/pathway_enrichment/03_Hallmarks_GSEA.tsv
# Output: results/figures/pub/main/Fig2C_hallmarks_gsea.pdf
#         results/figures/pub/main/Fig2C_hallmarks_gsea.png
#
# Ambiente: omics-R
# Ejecución:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/17c_pub_figD_gsea.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tools)
})

# ── Working directory (raíz del proyecto vía scripts/_setup.R) ───────────────
cat("=== 17c_pub_figD_gsea.R ===\n")
source(here::here("scripts", "_setup.R"))
setup_project()

# ── Directorios ───────────────────────────────────────────────────────────────
pub_dir <- "results/figures/pub/main"
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

# ── Estilos (fuente única: scripts/_fig_style.R) ──────────────────────────────
source(here::here("scripts", "_fig_style.R"))

# ── Cargar datos ──────────────────────────────────────────────────────────────
hallmarks_file <- "results/tables/pathway_enrichment/03_Hallmarks_GSEA.tsv"
if (!file.exists(hallmarks_file)) {
  stop("Missing: ", hallmarks_file, " — run scripts/03_pathway_enrichment.R first")
}

gsea_h <- read.delim(hallmarks_file, stringsAsFactors = FALSE)
n_sets <- nrow(gsea_h)
cat(sprintf("Hallmarks GSEA: %d gene sets significativos (FDR < 0.05)\n", n_sets))
cat("NES range:", round(min(gsea_h$NES), 2), "to", round(max(gsea_h$NES), 2), "\n\n")

# ── Preparar datos para figura ────────────────────────────────────────────────
# Limpiar nombres: quitar prefijo HALLMARK_ y convertir a Title Case
gsea_h <- gsea_h %>%
  mutate(
    Description = gsub("^HALLMARK_", "", Description),
    Description = gsub("_", " ", Description),
    Description = toTitleCase(tolower(Description))
  )

# Con <= 20 gene sets significativos se grafican TODOS (evita ocultar hallmarks
# como EMT); si hubiera muchos, se toman top 10 activados + 10 suprimidos por p.adjust.
if (n_sets <= 20) {
  plot_df <- gsea_h %>% arrange(NES)
} else {
  top_act  <- gsea_h %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = 10)
  top_supp <- gsea_h %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = 10)
  plot_df  <- bind_rows(top_act, top_supp) %>% arrange(NES)
}
plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)

cat(sprintf("Graficados: %d activados + %d suprimidos (total %d)\n",
            sum(plot_df$NES > 0), sum(plot_df$NES < 0), nrow(plot_df)))

# ── Figura ────────────────────────────────────────────────────────────────────
p_gsea <- ggplot(plot_df,
                 aes(x = NES, y = Description,
                     size = setSize, color = p.adjust)) +
  geom_point() +
  geom_vline(xintercept = 0, linewidth = 0.4,
             color = "grey50", linetype = "dashed") +
  scale_color_gradient(
    low  = "#D55E00", high = "#56B4E9",
    name = "FDR", trans = "log10",
    guide = guide_colorbar(reverse = TRUE,
                           barwidth  = unit(3, "mm"),
                           barheight = unit(18, "mm"))
  ) +
  scale_size_continuous(
    name = "Set size", range = c(1.8, 5.5),
    guide = guide_legend(override.aes = list(color = "grey50"))
  ) +
  labs(
    title    = NULL,
    subtitle = NULL,
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y     = element_text(size = 9.7),
    legend.position = "right"
  )

n_rows <- nrow(plot_df)
h_extra <- max(0, (n_rows - 10) * 6)
save_pub(p_gsea, "Fig2C_hallmarks_gsea", h_add = h_extra)
save_panel_obj(p_gsea, "Fig2C_hallmarks_gsea")   # insumo multipanel (canónico, 11 sets)

cat("\nFig2C_hallmarks_gsea — OK\n")
cat("Fin:", format(Sys.time()), "\n")
