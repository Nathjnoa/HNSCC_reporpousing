#!/usr/bin/env Rscript
# fig1_workflow.R — Figura 1: diseño del estudio y flujo analítico (HNSCC drug repurposing)
# Ambiente: omics-R
# Layout: horizontal con carriles, 5 etapas en pastel (Okabe-Ito, colorblind-safe).
# Números anclados al manuscrito post-auditoría (figures.md, results/tables/pub/*).

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(here)
})
source(here::here("scripts", "_fig_style.R"))

FAM <- "sans"   # misma familia que theme_pub() → consistencia con Fig2–6

# ── Paleta por etapa: pastel (fill) + acento (border/header) derivados de OKB ──
# OKB: blue #56B4E9, green #009E73, amber #E69F00, purple #CC79A7, verm #D55E00
STAGE <- list(
  input    = list(fill = "#DCEEFA", acc = "#1F78A8", lab = "INPUT"),
  signature= list(fill = "#D2EDE4", acc = "#00795A", lab = "SIGNATURE"),
  repurpose= list(fill = "#FBEBCB", acc = "#B97D00", lab = "REPURPOSING"),
  prioritize=list(fill = "#F2DEE9", acc = "#A85786", lab = "PRIORITIZATION"),
  validate = list(fill = "#FBDCC9", acc = "#A84A00", lab = "VALIDATION")
)
TXT  <- "#1F2937"   # texto cuerpo
GREY <- "#6B7280"   # flechas

# ── Helpers ───────────────────────────────────────────────────────────────────
stage_box <- function(xmin, xmax, ymin, ymax, st) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = st$fill, color = st$acc, linewidth = 0.7),
    annotate("text", x = (xmin + xmax) / 2, y = ymax + 2.4,
             label = st$lab, size = 3.0, fontface = "bold", family = FAM,
             color = st$acc, hjust = 0.5)
  )
}

ctext <- function(x, y, label, size = 2.8, face = "plain", col = TXT) {
  annotate("text", x = x, y = y, label = label, size = size, family = FAM,
           fontface = face, color = col, hjust = 0.5, lineheight = 0.95)
}

mini_box <- function(xmin, xmax, ymin, ymax, label, st, size = 2.6) {
  list(
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = "white", color = st$acc, linewidth = 0.45),
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, size = size, fontface = "bold", family = FAM,
             color = st$acc, hjust = 0.5, lineheight = 0.9)
  )
}

arrow_h <- function(x_from, x_to, y) {
  annotate("segment", x = x_from, xend = x_to, y = y, yend = y,
           arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
           color = GREY, linewidth = 0.7)
}

# ── Geometría del canvas ──────────────────────────────────────────────────────
# x: 0–100 ; y: 0–50.  Cajas y 9–40 (centro 24.5).  Etiqueta de etapa en y ~42.4.
# Prioritization (S4) y Validation (S5) ensanchadas para mejor distribución.
yb <- c(9, 40); yc <- mean(yb)
S1 <- c(2.0, 15.5); S2 <- c(18.0, 34.0); S3 <- c(36.5, 60.0)
S4 <- c(62.5, 81.0); S5 <- c(83.5, 99.8)

p <- ggplot() + xlim(0, 100) + ylim(0, 47) + theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(6, 6, 6, 6))

# ── Flechas entre etapas ──────────────────────────────────────────────────────
p <- p +
  arrow_h(S1[2], S2[1], yc) + arrow_h(S2[2], S3[1], yc) +
  arrow_h(S3[2], S4[1], yc) + arrow_h(S4[2], S5[1], yc)

# ── Etapa 1: INPUT ────────────────────────────────────────────────────────────
p <- p + stage_box(S1[1], S1[2], yb[1], yb[2], STAGE$input) +
  ctext(mean(S1), 34, "DIA-MS\nproteomics", 3.0, "bold") +
  ctext(mean(S1), 24, "10 paired\ntumor / normal", 2.8) +
  ctext(mean(S1), 14, "3,352 proteins\nquantified", 2.8, col = STAGE$input$acc)

# ── Etapa 2: SIGNATURE ────────────────────────────────────────────────────────
p <- p + stage_box(S2[1], S2[2], yb[1], yb[2], STAGE$signature) +
  ctext(mean(S2), 34, "Differential\nabundance\n(limma)", 3.0, "bold") +
  annotate("text", x = mean(S2), y = 25.5, parse = TRUE, family = FAM, size = 2.8,
           color = TXT, hjust = 0.5, label = "group('|', log[2]*FC, '|') > 1") +
  annotate("text", x = mean(S2), y = 21.8, family = FAM, size = 2.8,
           color = TXT, hjust = 0.5, label = "FDR < 0.05") +
  ctext(mean(S2), 14, "666 DE proteins\n(329 up / 337 dn)", 2.8, col = STAGE$signature$acc)

# ── Etapa 3: REPURPOSING (4 bases de datos, 2x2) ──────────────────────────────
p <- p + stage_box(S3[1], S3[2], yb[1], yb[2], STAGE$repurpose) +
  ctext(mean(S3), 36, "Drug mapping\n(4 resources)", 2.9, "bold")
cx1 <- 43.6; cx2 <- 52.9; mbw <- 8.8      # centros de columna + ancho mini-caja
mb <- function(cx, cy, lab) mini_box(cx - mbw/2, cx + mbw/2, cy - 2.6, cy + 2.6, lab, STAGE$repurpose, 2.5)
p <- p +
  mb(cx1, 27.5, "DGIdb") + mb(cx2, 27.5, "ChEMBL") +
  mb(cx1, 21.0, "Open\nTargets") + mb(cx2, 21.0, "L2S2") +
  ctext(mean(S3), 13, "458 candidate drugs", 2.8, "bold", STAGE$repurpose$acc)

# ── Etapa 4: PRIORITIZATION ───────────────────────────────────────────────────
p <- p + stage_box(S4[1], S4[2], yb[1], yb[2], STAGE$prioritize) +
  ctext(mean(S4), 34, "STRING PPI\nLouvain modules", 3.0, "bold") +
  ctext(mean(S4), 24, "Composite score", 2.8) +
  ctext(mean(S4), 14.5, "0.60 · TP\n+ 0.40 · DV", 2.9, "bold", STAGE$prioritize$acc)

# ── Etapa 5: VALIDATION (robustez + 2 cohortes) ───────────────────────────────
p <- p + stage_box(S5[1], S5[2], yb[1], yb[2], STAGE$validate) +
  ctext(mean(S5), 36, "Robustness\n(weights · LOD)", 2.9, "bold")
p <- p +
  mini_box(S5[1] + 1.6, S5[2] - 1.6, 24.0, 30.0, "CPTAC\nproteome", STAGE$validate, 2.5) +
  mini_box(S5[1] + 1.6, S5[2] - 1.6, 15.0, 21.0, "TCGA-HNSC\nRNA-seq", STAGE$validate, 2.5)

# ── Export ────────────────────────────────────────────────────────────────────
out_dir <- here::here("results", "figures", "pub", "main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "Fig1_workflow.png"), p,
       width = 190, height = 92, units = "mm", dpi = 300)
save_tiff(p, "Fig1_workflow", width_mm = 190, height_mm = 92, out_dir = out_dir)

cat("Fig1 workflow saved to", out_dir, "\n")
