#!/usr/bin/env Rscript
# fig1_workflow_draft.R — Bosquejo Fig1: workflow drug repurposing HNSCC
# Ambiente: omics-R

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ── Paleta por etapa ────────────────────────────────────────────────────────
COL <- list(
  input   = list(fill = "#DBEAFE", border = "#3B82F6", text = "#1E3A5F"),
  de      = list(fill = "#DCFCE7", border = "#16A34A", text = "#14532D"),
  db      = list(fill = "#FEF9C3", border = "#CA8A04", text = "#713F12"),
  score   = list(fill = "#EDE9FE", border = "#7C3AED", text = "#3B0764"),
  valid   = list(fill = "#FFE4E6", border = "#E11D48", text = "#881337")
)

# ── Helper: caja con texto ───────────────────────────────────────────────────
box <- function(xmin, xmax, ymin, ymax, label, sublabel = NULL, col,
                label_size = 3.2, sub_size = 2.5) {
  out <- list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = col$fill, color = col$border, linewidth = 0.7),
    annotate("text", x = (xmin + xmax) / 2, y = if (is.null(sublabel)) (ymin + ymax) / 2
             else (ymin + ymax) / 2 + 0.18,
             label = label, size = label_size, fontface = "bold",
             color = col$text, hjust = 0.5, vjust = 0.5)
  )
  if (!is.null(sublabel)) {
    out <- c(out, list(
      annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2 - 0.22,
               label = sublabel, size = sub_size, color = col$text,
               hjust = 0.5, vjust = 0.5, fontface = "plain")
    ))
  }
  out
}

# ── Helper: flecha vertical ──────────────────────────────────────────────────
arrow_v <- function(x, y_from, y_to, label = NULL) {
  out <- list(
    annotate("segment", x = x, xend = x, y = y_from, yend = y_to + 0.05,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
             color = "#6B7280", linewidth = 0.6)
  )
  if (!is.null(label)) {
    out <- c(out, list(
      annotate("text", x = x + 0.08, y = (y_from + y_to) / 2,
               label = label, size = 2.3, color = "#6B7280",
               hjust = 0, fontface = "italic")
    ))
  }
  out
}

arrow_h <- function(x_from, x_to, y, label = NULL) {
  out <- list(
    annotate("segment", x = x_from, xend = x_to - 0.05, y = y, yend = y,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             color = "#6B7280", linewidth = 0.5)
  )
  if (!is.null(label)) {
    out <- c(out, list(
      annotate("text", x = (x_from + x_to) / 2, y = y + 0.1,
               label = label, size = 2.2, color = "#6B7280", hjust = 0.5)
    ))
  }
  out
}

# ── Canvas ───────────────────────────────────────────────────────────────────
p <- ggplot() +
  xlim(0, 10) + ylim(0, 11) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 10, 10, 10))

# ── ETAPA 1: Input (y = 9.5-10.5) ────────────────────────────────────────────
p <- p +
  box(2, 8, 9.4, 10.4,
      label    = "Proteómica DIA-MS  ·  10 pares tumor/normal (HNSCC)",
      sublabel = "MaxQuant  ->  398 proteínas cuantificadas",
      col = COL$input)

# ── ETAPA 2: Análisis diferencial (y = 7.6-8.6) ──────────────────────────────
p <- p +
  arrow_v(5, 9.4, 8.6) +
  box(1.2, 4.6, 7.6, 8.6,
      label    = "Expresión diferencial",
      sublabel = "limma · adj.P<0.05 · |log2FC|>0.5",
      col = COL$de) +
  box(5.4, 8.8, 7.6, 8.6,
      label    = "Enriquecimiento de vías",
      sublabel = "GSEA · Hallmarks · KEGG · GO",
      col = COL$de)

# ── ETAPA 3: Consulta DBs (y = 5.8-6.8) ──────────────────────────────────────
# flechas desde etapa 2 hacia etapa 3
p <- p +
  arrow_v(2.9, 7.6, 6.8) +
  arrow_v(7.1, 7.6, 6.8) +
  # 4 cajas de DBs
  box(0.3, 2.2, 5.8, 6.8, label = "DGIdb",        col = COL$db) +
  box(2.5, 4.4, 5.8, 6.8, label = "ChEMBL",       col = COL$db) +
  box(4.7, 6.6, 5.8, 6.8, label = "OpenTargets",  col = COL$db) +
  box(6.9, 8.8, 5.8, 6.8, label = "L2S2",         col = COL$db) +
  # etiqueta de sección
  annotate("text", x = 5, y = 7.25, label = "Consulta de bases de datos de interacciones droga-gen",
           size = 2.6, color = "#92400E", hjust = 0.5, fontface = "italic")

# ── línea de convergencia de las 4 DBs ───────────────────────────────────────
p <- p +
  annotate("segment", x = 1.25, xend = 7.75, y = 5.8, yend = 5.8,
           color = "#CA8A04", linewidth = 0.5, linetype = "dashed") +
  arrow_v(5, 5.8, 5.1)

# ── ETAPA 4: Scoring + Red (y = 4.1-5.1) ────────────────────────────────────
p <- p +
  box(0.8, 4.5, 4.1, 5.1,
      label    = "Score compuesto",
      sublabel = "5 dimensiones: clínica · molecular\nred · vías · CMap  ·  N=398",
      col = COL$score) +
  box(5.5, 9.2, 4.1, 5.1,
      label    = "Red PPI STRING",
      sublabel = "Detección de módulos\n(Louvain · n=4 módulos)",
      col = COL$score) +
  # flecha horizontal entre los dos
  arrow_h(4.5, 5.5, 4.6)

# ── análisis LOD (sub-etapa, y = 2.9-3.8) ────────────────────────────────────
p <- p +
  arrow_v(2.65, 4.1, 3.8) +
  box(0.8, 4.5, 2.9, 3.8,
      label    = "Análisis de sensibilidad LOD",
      sublabel = "6 esquemas de pesos  ·  32 candidatos estables",
      col = COL$score)

# ── ETAPA 5: Validación externa (y = 1.0-2.1) ────────────────────────────────
p <- p +
  arrow_v(2.65, 2.9, 2.1) +
  arrow_v(7.35, 4.1, 2.1) +
  box(0.5, 4.8, 1.0, 2.1,
      label    = "Validación TCGA-HNSC",
      sublabel = "RNA-seq · 528 tumores\nConcordancia proteoma-transcriptoma",
      col = COL$valid) +
  box(5.2, 9.5, 1.0, 2.1,
      label    = "Validación CPTAC",
      sublabel = "Proteómica independiente\nSupervivencia global (KM)",
      col = COL$valid)

# ── Leyenda de etapas ────────────────────────────────────────────────────────
legend_y <- 0.55
legend_items <- data.frame(
  x     = c(0.5, 2.5, 4.5, 6.5, 8.5),
  label = c("Input", "Análisis", "Bases de datos", "Integración", "Validación"),
  fill  = c(COL$input$fill, COL$de$fill, COL$db$fill, COL$score$fill, COL$valid$fill),
  border = c(COL$input$border, COL$de$border, COL$db$border, COL$score$border, COL$valid$border)
)

for (i in seq_len(nrow(legend_items))) {
  p <- p +
    annotate("rect",
             xmin = legend_items$x[i] - 0.35, xmax = legend_items$x[i] + 0.35,
             ymin = legend_y - 0.22,           ymax = legend_y + 0.22,
             fill = legend_items$fill[i], color = legend_items$border[i], linewidth = 0.5) +
    annotate("text", x = legend_items$x[i], y = legend_y,
             label = legend_items$label[i], size = 2.2,
             color = "#374151", hjust = 0.5, fontface = "bold")
}

# ── Export ───────────────────────────────────────────────────────────────────
out_base <- "results/figures/pub/main/Fig1_workflow"

ggsave(paste0(out_base, ".pdf"), p, width = 10, height = 11, device = "pdf")
ggsave(paste0(out_base, ".png"), p, width = 10, height = 11, dpi = 300)

cat("Saved:", out_base, ".pdf/.png\n")
