#!/usr/bin/env Rscript
# =============================================================================
# 17b_volcano_oxphos.R
# Volcán con proteínas OXPHOS/Complejo I destacadas para póster
# Output: results/figures/pub/main/Sec0_FigB_volcano_oxphos.png/.pdf
# Ambiente: omics-R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggrepel)
})

# Directorio del proyecto
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  proj_dir <- dirname(dirname(normalizePath(sub("^--file=", "", script_flag))))
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml")))
    setwd(proj_dir)
}
cat("Working directory:", getwd(), "\n")

# ── Colores y tema ────────────────────────────────────────────────────────────
# Paleta base (misma que script 17)
DE_COLS <- c(up = "#D55E00", down = "#0072B2", ns = "#CCCCCC", oxphos = "#0072B2")

theme_pub <- function() {
  theme_classic(base_size = 8.5, base_family = "sans") +
    theme(
      plot.title   = element_text(size = 10, face = "bold", hjust = 0,
                                  margin = margin(b = 2)),
      axis.title   = element_text(size = 8.5),
      axis.text    = element_text(size = 7.5, color = "black"),
      legend.title = element_text(size = 7.5, face = "bold"),
      legend.text  = element_text(size = 7.5),
      legend.key.size   = unit(3, "mm"),
      legend.background = element_blank(),
      axis.line  = element_line(linewidth = 0.25, color = "black"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      plot.margin = margin(4, 6, 4, 6, "mm")
    )
}

# ── Cargar datos ──────────────────────────────────────────────────────────────
de <- read_tsv("results/tables/de_limma/01_TVsS_all_proteins.tsv",
               show_col_types = FALSE) %>%
  mutate(
    direction = case_when(
      sig.FDR_TVsS ==  1 ~ "up",
      sig.FDR_TVsS == -1 ~ "down",
      TRUE               ~ "ns"
    ),
    direction = factor(direction, levels = c("up", "down", "ns", "oxphos"))
  )

n_up   <- sum(de$direction == "up")
n_down <- sum(de$direction == "down")
n_ns   <- sum(de$direction == "ns")

# ── Proteínas OXPHOS a destacar ───────────────────────────────────────────────
# Representantes de cada complejo, seleccionadas por logFC más negativo
oxphos_highlight <- c(
  # Complejo I
  "NDUFS3", "NDUFS2", "NDUFB4", "NDUFB1", "NDUFA2",
  # Complejo II
  "SDHB",
  # Complejo III
  "CYC1", "UQCRC1",
  # Complejo IV
  "COX4I1",
  # Complejo V
  "ATP5F1B", "ATP5F1C"
)

# Reclasificar las OXPHOS como grupo propio (solo las DE significativas)
de <- de %>%
  mutate(
    group = case_when(
      gene_symbol %in% oxphos_highlight & direction %in% c("down", "up") ~ "oxphos",
      TRUE ~ as.character(direction)
    ),
    group = factor(group, levels = c("up", "down", "oxphos", "ns"))
  )

# Subsets para capas
de_bg    <- de %>% filter(group != "oxphos")
de_ox    <- de %>% filter(group == "oxphos")


# ── Figura ────────────────────────────────────────────────────────────────────
p <- ggplot(de_bg, aes(x = logFC_TVsS, y = -log10(adj.P.Val_TVsS))) +

  # 1. Fondo: todos los puntos no-OXPHOS
  geom_point(aes(color = group, alpha = group),
             size = 1.8, stroke = 0) +

  # 2. Líneas de corte
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = 0.4, color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             linewidth = 0.35, color = "grey50") +

  # 3. Puntos OXPHOS encima (borde negro para destacar sobre los otros azules)
  geom_point(data = de_ox,
             aes(x = logFC_TVsS, y = -log10(adj.P.Val_TVsS)),
             color = "#0072B2", size = 2.8, stroke = 0) +

  # 4. Etiquetas OXPHOS
  geom_label_repel(
    data = de_ox,
    aes(label = gene_symbol),
    color = "#0072B2", fill = "white",
    size = 1.8, label.padding = 0.1, label.size = 0.25,
    max.overlaps = 20, segment.size = 0.3,
    min.segment.length = 0.1, show.legend = FALSE,
    fontface = "bold"
  ) +

  # 6. Escalas
  scale_color_manual(
    values = DE_COLS,
    labels = c(
      up   = sprintf("Up (%d)",   n_up),
      down = sprintf("Down (%d)", n_down),
      ns   = sprintf("NS (%d)",   n_ns)
    ),
    breaks = c("up", "down", "ns")
  ) +
  scale_alpha_manual(
    values = c(up = 0.85, down = 0.85, oxphos = 1, ns = 0.25),
    guide  = "none"
  ) +

  labs(
    title = "Differential proteome — Tumor vs. Adjacent Normal",
    x     = expression(log[2]~"Fold Change  (Tumor vs. Normal)"),
    y     = expression(-log[10]~"FDR"),
    color = NULL
  ) +
  theme_pub() +
  theme(legend.position = c(0.82, 0.88)) +
  guides(color = guide_legend(override.aes = list(size = 2.5, stroke = 0)))

# ── Guardar ───────────────────────────────────────────────────────────────────
out_dir <- "results/figures/pub/main"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "Sec0_FigB_volcano_oxphos.pdf"),
       p, width = 180, height = 120, units = "mm",
       device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(out_dir, "Sec0_FigB_volcano_oxphos.png"),
       p, width = 180, height = 120, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("Guardado: Sec0_FigB_volcano_oxphos.png/.pdf\n")
