# =============================================================================
# _fig_style.R — Sistema de estilo centralizado para figuras de publicación
# =============================================================================
# Fuente única de verdad para paleta, presets, tema y funciones de guardado.
# TODOS los scripts de figuras deben:  source(here::here("scripts","_fig_style.R"))
#
# Provee:
#   PRESETS            — dimensiones (mm) y tamaños de fuente por formato
#   OKB                — paleta Okabe-Ito (colorblind-safe)
#   DE_COLS / COND_COLS / HPV_COLS / PHASE_COLS / DRUG_CLASS_LABELS
#   theme_pub()        — tema ggplot estandarizado
#   save_pub()         — guarda ggplot en PNG (300 dpi) + PDF vectorial
#   save_ch()          — guarda ComplexHeatmap en PNG + PDF
#   save_tiff()        — guarda figura ENSAMBLADA en TIFF 600 dpi LZW + PDF
#
# Política de exportación del proyecto:
#   - Paneles individuales  → PNG (revisión) + PDF (vector)   [save_pub/save_ch]
#   - Figura multipanel final → TIFF 600 dpi LZW + PDF        [save_tiff]
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# ── Presets (dimensiones en mm + tamaños de fuente en pt) ─────────────────────
PRESETS <- list(
  double_col = list(w = 180, h = 120, base = 11, title = 13,
                    axis = 11, leg = 9.7, tick = 9.7, lwd = 0.7, pt = 1.8),
  single_col = list(w = 85,  h = 65,  base = 8.0, title = 9.0,
                    axis = 8.0, leg = 7.0, tick = 7.0, lwd = 0.6, pt = 1.6)
)

# ── Paleta Okabe-Ito (colorblind-safe) ────────────────────────────────────────
OKB <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
         "#0072B2", "#D55E00", "#CC79A7", "#000000")

# ── Constantes de color/etiqueta compartidas ──────────────────────────────────
DE_COLS    <- c(up = "#D55E00", down = "#0072B2", ns = "#CCCCCC")
COND_COLS  <- c(Tumor = "#D55E00", "Adjacent Normal" = "#0072B2")
HPV_COLS   <- c("HPV+" = "#009E73", "HPV-" = "#E69F00")
PHASE_COLS <- c(
  "Approved"  = "#009E73",
  "Phase III" = "#56B4E9",
  "Phase II"  = "#E69F00",
  "Phase I"   = "#CC79A7",
  "No data"   = "#CCCCCC"
)
# Clase clínico-regulatoria (script 10 exporta drug_class = A/B/C/D)
DRUG_CLASS_LABELS <- c(
  A = "HNSCC-approved",
  B = "Other cancer",
  C = "Non-oncology",
  D = "Experimental"
)

# ── Tema centralizado ─────────────────────────────────────────────────────────
theme_pub <- function(preset = "double_col") {
  p <- PRESETS[[preset]]
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
      strip.text         = element_text(size = p$base, face = "bold"),
      strip.background   = element_blank(),
      axis.line          = element_line(linewidth = p$lwd * 0.35, color = "black"),
      axis.ticks         = element_line(linewidth = p$lwd * 0.35, color = "black"),
      axis.ticks.length  = unit(1.5, "mm"),
      plot.margin        = margin(4, 6, 4, 6, "mm")
    )
}

# ── Guardado de paneles individuales: PNG (300 dpi) + PDF vectorial ────────────
save_pub <- function(p, name, preset = "double_col", w_add = 0, h_add = 0,
                     out_dir = "results/figures/pub/main") {
  d <- PRESETS[[preset]]
  w <- d$w + w_add
  h <- d$h + h_add
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))
  ggsave(pdf_path, p, width = w, height = h, units = "mm",
         device = cairo_pdf, limitsize = FALSE)
  ggsave(png_path, p, width = w, height = h, units = "mm",
         dpi = 300, limitsize = FALSE)
  cat(sprintf("  Saved: %s (%g x %g mm)\n", basename(pdf_path), w, h))
  invisible(pdf_path)
}

# ── Guardado de ComplexHeatmap (grid, no ggsave): PNG + PDF ────────────────────
save_ch <- function(ht, name, preset = "double_col", w_add = 0, h_add = 0,
                    out_dir = "results/figures/pub/main") {
  d <- PRESETS[[preset]]
  w_in <- (d$w + w_add) / 25.4
  h_in <- (d$h + h_add) / 25.4
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))
  pdf(pdf_path, width = w_in, height = h_in)
  ComplexHeatmap::draw(ht)
  dev.off()
  png(png_path,
      width  = (d$w + w_add) * 300 / 25.4,
      height = (d$h + h_add) * 300 / 25.4,
      res    = 300)
  ComplexHeatmap::draw(ht)
  dev.off()
  cat(sprintf("  Saved: %s (%.1f x %.1f in)\n", basename(pdf_path), w_in, h_in))
  invisible(pdf_path)
}

# ── Guardado de figura MULTIPANEL final: TIFF 600 dpi LZW + PDF ────────────────
# Acepta un objeto graficable (ggplot/patchwork) o un grob/gTree.
save_tiff <- function(plot_obj, name, width_mm, height_mm,
                      out_dir = "results/figures/pub/main", dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tif_path <- file.path(out_dir, paste0(name, ".tif"))
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  is_grid <- inherits(plot_obj, c("grob", "gTree", "gList"))

  # TIFF 600 dpi con compresión LZW (submission a revista)
  tiff(tif_path, width = width_mm, height = height_mm, units = "mm",
       res = dpi, compression = "lzw")
  if (is_grid) grid.draw(plot_obj) else print(plot_obj)
  dev.off()

  # PDF vectorial (referencia/edición)
  cairo_pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4)
  if (is_grid) grid.draw(plot_obj) else print(plot_obj)
  dev.off()

  cat(sprintf("  Saved: %s + .pdf  (%g x %g mm, %d dpi TIFF/LZW)\n",
              basename(tif_path), width_mm, height_mm, dpi))
  invisible(tif_path)
}

# ── UpSet: dibujo con números sobre barras (decoración) ───────────────────────
# El top annotation debe llamarse "Intersection size" (anno_barplot manual).
# Se usa tanto en el guardado standalone como en el ensamblado multipanel, para
# que los números aparezcan en ambos.
draw_upset_panel <- function(ht, m) {
  ComplexHeatmap::draw(ht)
  cs  <- ComplexHeatmap::comb_size(m)
  ord <- order(-cs)                       # mismo orden que comb_order del UpSet
  cs_ord <- cs[ord]
  ComplexHeatmap::decorate_annotation("Intersection size", {
    n <- length(cs_ord)
    grid::grid.text(cs_ord, x = seq_len(n),
                    y = grid::unit(cs_ord, "native") + grid::unit(1, "mm"),
                    just = "bottom", default.units = "native",
                    gp = grid::gpar(fontsize = 6.5))
  })
}

# Guardado de panel UpSet (PNG + PDF) con números, vía draw_upset_panel.
save_upset <- function(ht, m, name, preset = "double_col", w_add = 0, h_add = 0,
                       out_dir = "results/figures/pub/main") {
  d <- PRESETS[[preset]]
  w_in <- (d$w + w_add) / 25.4
  h_in <- (d$h + h_add) / 25.4
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))
  cairo_pdf(pdf_path, width = w_in, height = h_in); draw_upset_panel(ht, m); dev.off()
  png(png_path, width = (d$w + w_add) * 300 / 25.4,
      height = (d$h + h_add) * 300 / 25.4, res = 300)
  draw_upset_panel(ht, m); dev.off()
  cat(sprintf("  Saved: %s (%.1f x %.1f in)\n", basename(pdf_path), w_in, h_in))
  invisible(pdf_path)
}

# ── Directorio para serializar objetos de panel (insumo de multipaneles) ──────
PANEL_OBJ_DIR <- "results/figures/pub/.objects"
save_panel_obj <- function(obj, name, dir = PANEL_OBJ_DIR) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".rds"))
  saveRDS(obj, path)
  cat(sprintf("  Panel object cached: %s\n", basename(path)))
  invisible(path)
}
