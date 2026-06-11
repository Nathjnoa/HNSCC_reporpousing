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

# ── Estilos (coincide con 17_pub_figures.R) ───────────────────────────────────
PRESETS <- list(
  double_col = list(w = 180, h = 120, base = 8.5, title = 10,
                    axis = 8.5, leg = 7.5, tick = 7.5, lwd = 0.7, pt = 1.8)
)

theme_pub <- function() {
  p <- PRESETS$double_col
  theme_classic(base_size = p$base, base_family = "sans") +
    theme(
      plot.title    = element_text(size = p$title, face = "bold", hjust = 0,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(size = p$base, color = "grey40"),
      axis.title    = element_text(size = p$axis),
      axis.text     = element_text(size = p$tick, color = "black"),
      legend.title  = element_text(size = p$leg, face = "bold"),
      legend.text   = element_text(size = p$leg),
      legend.key.size    = unit(3, "mm"),
      legend.background  = element_blank(),
      axis.line          = element_line(linewidth = p$lwd * 0.35, color = "black"),
      axis.ticks         = element_line(linewidth = p$lwd * 0.35, color = "black"),
      axis.ticks.length  = unit(1.5, "mm"),
      plot.margin        = margin(4, 6, 4, 6, "mm")
    )
}

save_pub <- function(p, name, w_add = 0, h_add = 0) {
  d <- PRESETS$double_col
  w <- d$w + w_add
  h <- d$h + h_add
  pdf_path <- file.path(pub_dir, paste0(name, ".pdf"))
  png_path <- file.path(pub_dir, paste0(name, ".png"))
  ggsave(pdf_path, p, width = w, height = h, units = "mm",
         device = cairo_pdf, limitsize = FALSE)
  ggsave(png_path, p, width = w, height = h, units = "mm",
         dpi = 300, limitsize = FALSE)
  cat(sprintf("  Saved: %s (%g x %g mm)\n", basename(pdf_path), w, h))
  invisible(pdf_path)
}

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

# Selección balanceada: hasta 10 activados (NES > 0) + 10 suprimidos (NES < 0)
n_each <- min(10, floor(n_sets / 2))
top_act  <- gsea_h %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = n_each)
top_supp <- gsea_h %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = n_each)
plot_df  <- bind_rows(top_act, top_supp) %>% arrange(NES)
plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)

cat(sprintf("Seleccionados: %d activados + %d suprimidos\n",
            nrow(top_act), nrow(top_supp)))

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
    title    = "MSigDB Hallmarks — GSEA (Tumor vs. Adjacent Normal)",
    subtitle = sprintf("Ranked by \u03c0-statistic (sign(logFC) \u00d7 |logFC| \u00d7 -log\u2081\u2080(FDR)) | %d gene sets FDR < 0.05",
                       n_sets),
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_pub() +
  theme(
    axis.text.y     = element_text(size = 7.5),
    legend.position = "right"
  )

n_rows <- nrow(plot_df)
h_extra <- max(0, (n_rows - 10) * 6)
save_pub(p_gsea, "Fig2C_hallmarks_gsea", h_add = h_extra)

cat("\nFig2C_hallmarks_gsea — OK\n")
cat("Fin:", format(Sys.time()), "\n")
