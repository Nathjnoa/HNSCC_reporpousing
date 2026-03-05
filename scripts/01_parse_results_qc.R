# ============================================================
# Script: 01_parse_results_qc.R
# Proyecto: HNSCC Drug Repurposing
# Propósito: Parsear resultados de proteómica (proteoDA/limma),
#            generar QC plots y exportar proteínas DE del
#            contraste TVsS (Tumor vs Sano) para análisis
#            downstream.
#
# Input:
#   data/raw/results_proteomica.tsv    — output de write_limma_tables()
#   data/raw/metadata.csv             — metadatos de muestras (condition, vph)
#
# Output:
#   results/tables/de_limma/
#     01_TVsS_all_proteins.tsv         — todas las proteínas con estadísticos TVsS
#     01_TVsS_significant.tsv          — sig.FDR != 0 (up + down)
#     01_TVsS_upregulated.tsv          — sig.FDR == 1
#     01_TVsS_downregulated.tsv        — sig.FDR == -1
#   results/figures/qc/
#     01_boxplot_intensidades.pdf      — distribución log2 por muestra
#     01_PCA_muestras.pdf              — PCA coloreado por condición y VPH
#     01_volcano_TVsS.pdf              — volcano plot mejorado
#     01_resumen_DE.pdf                — conteos up/down
#
# Ambiente: omics-R
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(viridis)
})

# ── Rutas ────────────────────────────────────────────────────
proj_dir  <- getwd()   # ejecutar desde la raíz del proyecto
raw_file  <- file.path(proj_dir, "data/raw/results_proteomica.tsv")
meta_file <- file.path(proj_dir, "data/raw/metadata.csv")

out_tables  <- file.path(proj_dir, "results/tables/de_limma")
out_figures <- file.path(proj_dir, "results/figures/qc")

dir.create(out_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_figures, showWarnings = FALSE, recursive = TRUE)

# ── Log ──────────────────────────────────────────────────────
log_dir  <- file.path(proj_dir, "logs")
dir.create(log_dir, showWarnings = FALSE)
log_file <- file.path(log_dir,
  paste0("01_parse_qc_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("=== 01_parse_results_qc.R ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# ============================================================
# 1. CARGAR DATOS
# ============================================================
# El TSV tiene 3 filas de encabezados fusionados (Excel-style)
# Row 4 contiene los nombres reales de columnas
cat("-- Cargando results_proteomica.tsv...\n")

raw <- read.table(
  file             = raw_file,
  sep              = "\t",
  dec              = ",",          # separador decimal europeo
  skip             = 3,            # saltar filas 1-3 (headers fusionados)
  header           = TRUE,
  stringsAsFactors = FALSE,
  check.names      = TRUE,
  fill             = TRUE,
  quote            = ""
)

cat("  Dimensiones brutas:", nrow(raw), "proteínas ×", ncol(raw), "columnas\n")

# ── Renombrar columnas por comparación ──────────────────────
# Posiciones (1-indexed en R):
#  1-3   : gene, UniProt.ID, entryname
#  4-23  : intensidades log2 (M1S, M1T, ..., M10S, M10T)
#  24-33 : TVsS   estadísticos
#  34-43 : TpVsTn estadísticos
#  44-53 : TpVsSp estadísticos

stat_cols <- c("logFC", "CI.L", "CI.R", "avg_intensity",
               "t_stat", "B_stat", "P.Value", "adj.P.Val",
               "sig.PVal", "sig.FDR")

colnames(raw)[24:33] <- paste0(stat_cols, "_TVsS")
colnames(raw)[34:43] <- paste0(stat_cols, "_TpVsTn")
colnames(raw)[44:53] <- paste0(stat_cols, "_TpVsSp")

# Renombrar anotación para mayor claridad
colnames(raw)[1] <- "gene_symbol"
colnames(raw)[2] <- "uniprot_id"
colnames(raw)[3] <- "entry_name"

# ── Columnas de intensidades ─────────────────────────────────
intensity_cols <- colnames(raw)[4:23]   # M1S, M1T, ..., M10S, M10T

# Convertir intensidades a numérico (pueden quedar como character)
raw[, intensity_cols] <- lapply(raw[, intensity_cols], as.numeric)

# Convertir estadísticos TVsS a numérico
raw[, paste0(stat_cols, "_TVsS")] <-
  lapply(raw[, paste0(stat_cols, "_TVsS")], function(x) as.numeric(as.character(x)))

cat("  Columnas de intensidades:", paste(intensity_cols, collapse = ", "), "\n\n")

# ============================================================
# 2. CARGAR METADATA
# ============================================================
meta <- read.csv2(meta_file, stringsAsFactors = FALSE)
# meta tiene: sample_id, condition, vph, group
cat("-- Metadata:\n")
print(table(meta$condition, meta$vph))
cat("\n")

# ============================================================
# 3. RESUMEN DE PROTEÍNAS SIGNIFICATIVAS (TVsS)
# ============================================================
cat("-- Resumen TVsS (Tumor vs Sano):\n")
cat("  Total proteínas cuantificadas:", nrow(raw), "\n")

# Nota: proteoDA puede producir valores no enteros (ej. 0.73, 1.71) en sig.FDR
# para proteínas borderline. Solo incluimos valores exactos 1 y -1.
n_up   <- sum(raw$sig.FDR_TVsS ==  1, na.rm = TRUE)
n_down <- sum(raw$sig.FDR_TVsS == -1, na.rm = TRUE)
n_ns   <- nrow(raw) - n_up - n_down   # resto = NS (incluye NA y borderline)

cat("  Upregulated   (sig.FDR = 1) :", n_up, "\n")
cat("  Downregulated (sig.FDR = -1):", n_down, "\n")
cat("  No significativas           :", n_ns, "\n")
cat("  Total significativas        :", n_up + n_down, "\n\n")

# Rango de logFC en significativas
sig_fc <- raw$logFC_TVsS[raw$sig.FDR_TVsS != 0]
cat(sprintf("  logFC rango: %.2f a %.2f\n\n",
            min(sig_fc, na.rm = TRUE), max(sig_fc, na.rm = TRUE)))

# ============================================================
# 4. EXPORTAR TABLAS DE PROTEÍNAS DE
# ============================================================
cat("-- Exportando tablas...\n")

# Columnas a incluir en los outputs
tvs_stat_cols <- paste0(stat_cols, "_TVsS")
out_cols <- c("gene_symbol", "uniprot_id", "entry_name",
              tvs_stat_cols, intensity_cols)

# Todas las proteínas con estadísticos TVsS
df_all <- raw[, out_cols]
write.table(df_all,
  file.path(out_tables, "01_TVsS_all_proteins.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Significativas (up + down) — filtro estricto: solo 1 o -1
df_sig <- raw[!is.na(raw$sig.FDR_TVsS) & raw$sig.FDR_TVsS %in% c(1, -1), out_cols]
write.table(df_sig,
  file.path(out_tables, "01_TVsS_significant.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Solo upreguladas
df_up <- raw[!is.na(raw$sig.FDR_TVsS) & raw$sig.FDR_TVsS == 1, out_cols]
write.table(df_up,
  file.path(out_tables, "01_TVsS_upregulated.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# Solo downreguladas
df_down <- raw[!is.na(raw$sig.FDR_TVsS) & raw$sig.FDR_TVsS == -1, out_cols]
write.table(df_down,
  file.path(out_tables, "01_TVsS_downregulated.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("  01_TVsS_all_proteins.tsv   : %d proteínas\n", nrow(df_all)))
cat(sprintf("  01_TVsS_significant.tsv    : %d proteínas\n", nrow(df_sig)))
cat(sprintf("  01_TVsS_upregulated.tsv    : %d proteínas\n", nrow(df_up)))
cat(sprintf("  01_TVsS_downregulated.tsv  : %d proteínas\n\n", nrow(df_down)))

# ============================================================
# 5. FIGURA QC-1: BOXPLOT DE INTENSIDADES POR MUESTRA
# ============================================================
cat("-- Generando figura QC-1: boxplot de intensidades...\n")

# Formato largo para ggplot
intens_long <- raw[, c("gene_symbol", intensity_cols)] %>%
  pivot_longer(cols = all_of(intensity_cols),
               names_to  = "sample_id",
               values_to = "log2_intensity") %>%
  filter(!is.na(log2_intensity)) %>%
  left_join(meta, by = "sample_id")

p_box <- ggplot(intens_long,
                aes(x = sample_id, y = log2_intensity,
                    fill = condition, color = vph)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = c("NORMAL" = "#4393C3", "TUMORAL" = "#D6604D"),
                    name = "Condición") +
  scale_color_manual(values = c("NEGATIVE" = "grey40", "POSITIVE" = "#F4A582"),
                     name = "VPH") +
  labs(title = "Distribución de intensidades por muestra",
       subtitle = "Log2 Cyclic Loess Normalized | azul = Normal, rojo = Tumor | borde = VPH",
       x = NULL, y = "Log2 Intensidad") +
  theme_classic(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        plot.title   = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave(file.path(out_figures, "01_boxplot_intensidades.pdf"),
       p_box, width = 10, height = 5)
cat("  Guardado: 01_boxplot_intensidades.pdf\n")

# ============================================================
# 6. FIGURA QC-2: PCA DE MUESTRAS
# ============================================================
cat("-- Generando figura QC-2: PCA...\n")

# Matriz de intensidades transpuesta: muestras × proteínas
mat_int <- as.matrix(raw[, intensity_cols])
rownames(mat_int) <- raw$gene_symbol

# Filtrar proteínas con NA en cualquier muestra para PCA limpio
mat_complete <- mat_int[complete.cases(mat_int), ]
cat(sprintf("  Proteínas sin NA para PCA: %d\n", nrow(mat_complete)))

pca_res <- prcomp(t(mat_complete), scale. = TRUE, center = TRUE)

# Varianza explicada
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample_id <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = "sample_id")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                             color = condition, shape = vph,
                             label = sample_id)) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(size = 2.8, max.overlaps = 15,
                  box.padding = 0.3, point.padding = 0.2) +
  scale_color_manual(values = c("NORMAL" = "#4393C3", "TUMORAL" = "#D6604D"),
                     name = "Condición") +
  scale_shape_manual(values = c("NEGATIVE" = 16, "POSITIVE" = 17),
                     name = "VPH") +
  labs(title = "PCA — Muestras proteómica HNSCC",
       subtitle = "Cada punto = 1 muestra | círculo = VPH−, triángulo = VPH+",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave(file.path(out_figures, "01_PCA_muestras.pdf"),
       p_pca, width = 6.5, height = 5)
cat("  Guardado: 01_PCA_muestras.pdf\n")

# ============================================================
# 7. FIGURA QC-3: VOLCANO PLOT TVsS
# ============================================================
cat("-- Generando figura QC-3: volcano plot TVsS...\n")

# Preparar dataframe de volcáno con categoría DE
df_volc <- raw %>%
  select(gene_symbol, logFC_TVsS, adj.P.Val_TVsS, sig.FDR_TVsS) %>%
  filter(!is.na(logFC_TVsS), !is.na(adj.P.Val_TVsS)) %>%
  mutate(
    neg_log10_p = -log10(adj.P.Val_TVsS + 1e-10),  # +epsilon para evitar -Inf
    DE_status = case_when(
      sig.FDR_TVsS ==  1 & !is.na(sig.FDR_TVsS) ~ "Up",
      sig.FDR_TVsS == -1 & !is.na(sig.FDR_TVsS) ~ "Down",
      TRUE                                        ~ "NS"
    ),
    # Etiquetar top genes (|logFC| alto Y p-value bajo)
    label = ifelse(
      (logFC_TVsS > 3.5 | logFC_TVsS < -3.5) &
        adj.P.Val_TVsS < 0.01,
      gene_symbol, ""
    )
  )

# Umbrales para líneas guía
fc_thresh  <- 1.0
pv_thresh  <- -log10(0.05)

p_volc <- ggplot(df_volc,
                 aes(x = logFC_TVsS, y = neg_log10_p,
                     color = DE_status, label = label)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = pv_thresh,
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh),
             linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_text_repel(size = 2.5, max.overlaps = 20,
                  box.padding = 0.4, show.legend = FALSE) +
  scale_color_manual(
    values = c("Up" = "#D6604D", "Down" = "#4393C3", "NS" = "grey70"),
    labels = c(
      paste0("Up (n=", n_up, ")"),
      paste0("Down (n=", n_down, ")"),
      paste0("NS (n=", n_ns, ")")
    ),
    name = NULL
  ) +
  labs(
    title    = "Volcáno plot — Tumor vs Sano (TVsS)",
    subtitle = paste0("Criterio: adj.P.Val < 0.05 & |logFC| > 1 (BH) | n=",
                      n_up + n_down, " significativas"),
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"(adj.P.Val)")
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "top")

ggsave(file.path(out_figures, "01_volcano_TVsS.pdf"),
       p_volc, width = 6.5, height = 6)
cat("  Guardado: 01_volcano_TVsS.pdf\n")

# ============================================================
# 8. FIGURA QC-4: RESUMEN CONTEOS DE/BARRAS
# ============================================================
cat("-- Generando figura QC-4: resumen conteos DE...\n")

df_counts <- data.frame(
  Direction = factor(c("Upregulated", "Downregulated"),
                     levels = c("Upregulated", "Downregulated")),
  Count     = c(n_up, n_down),
  Color     = c("#D6604D", "#4393C3")
)

p_counts <- ggplot(df_counts, aes(x = Direction, y = Count, fill = Direction)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = Count), vjust = -0.4, fontface = "bold", size = 4.5) +
  scale_fill_manual(values = c("Upregulated" = "#D6604D",
                                "Downregulated" = "#4393C3")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title    = "Proteínas DE — TVsS",
       subtitle = "adj.P.Val < 0.05 & |logFC| > 1",
       x = NULL, y = "Número de proteínas") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

ggsave(file.path(out_figures, "01_resumen_DE.pdf"),
       p_counts, width = 4, height = 4.5)
cat("  Guardado: 01_resumen_DE.pdf\n\n")

# ============================================================
# 9. RESUMEN FINAL
# ============================================================
cat("=== RESUMEN FINAL ===\n")
cat("Proteínas totales cuantificadas :", nrow(raw), "\n")
cat("Proteínas DE significativas     :", n_up + n_down, "\n")
cat("  Upregulated                  :", n_up, "\n")
cat("  Downregulated                :", n_down, "\n")
cat("logFC mínimo (sig)             :", round(min(sig_fc, na.rm = TRUE), 2), "\n")
cat("logFC máximo (sig)             :", round(max(sig_fc, na.rm = TRUE), 2), "\n")
cat("\nTablas exportadas en:", out_tables, "\n")
cat("Figuras exportadas en:", out_figures, "\n")
cat("\nFin:", format(Sys.time()), "\n")
sink()

message("Script 01 completado. Ver log en: ", log_file)
