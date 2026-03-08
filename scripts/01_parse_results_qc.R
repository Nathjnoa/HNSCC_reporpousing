# ============================================================
# Script: 01_parse_results_qc.R
# Proyecto: HNSCC Drug Repurposing
# Propósito: Reanalizar proteómica con modelo limma HPV-ajustado
#            (TVsS con HPV como covariable).
#
# Modelo: ~ condición (Tumor/Normal) + estado_VPH
#         bloqueo intra-paciente vía duplicateCorrelation()
#
# Input:
#   data/raw/results_proteomica.tsv    — matrix procesada (log2, MinProb, median-norm)
#   data/raw/metadata.csv             — metadatos (sample_id, condition, vph)
#
# Output:
#   results/tables/de_limma/
#     01_TVsS_all_proteins.tsv         — todas las proteínas con estadísticos
#     01_TVsS_significant.tsv          — adj.P.Val < 0.05 & |logFC| > 1
#     01_TVsS_upregulated.tsv
#     01_TVsS_downregulated.tsv
#     01_missingness_analysis.tsv      — % missing por proteína/grupo (bias check)
#   results/figures/qc/
#     01_boxplot_intensidades.pdf
#     01_PCA_muestras.pdf
#     01_volcano_TVsS.pdf
#     01_resumen_DE.pdf
#     01_missingness_vs_direction.pdf  — NUEVO: sesgo MinProb imputation
#
# Ambiente: omics-R
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(ggrepel)
  library(patchwork)
  library(viridis)
  library(yaml)
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
on.exit(sink(), add = TRUE)
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

stopifnot("Número de columnas inesperado en results_proteomica.tsv (esperado >= 53)" = ncol(raw) >= 53)
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

cat("  Columnas de intensidades:", paste(intensity_cols, collapse = ", "), "\n\n")

# ============================================================
# 2. CARGAR METADATA — incluir patient_id para diseño pareado
# ============================================================
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
# Extraer patient_id desde nombre de muestra: M1S/M1T -> patient "M1"
meta <- meta %>%
  mutate(patient_id = sub("[ST]$", "", sample_id))

cat("-- Metadata:\n")
print(table(meta$condition, meta$vph))
cat("\nPacientes únicos:", n_distinct(meta$patient_id), "\n")
cat("(esperado: 10 pacientes, 2 muestras c/u = 20 total)\n\n")

# Verificar alineación muestra-columna
stopifnot("Muestras en metadata no coinciden con columnas de intensidades" =
          all(intensity_cols %in% meta$sample_id))
meta <- meta[match(intensity_cols, meta$sample_id), ]  # alinear orden

# ============================================================
# 3. ANÁLISIS DE MISSINGNESS — verificar sesgo de imputación MinProb
# ============================================================
cat("-- Análisis de missingness por grupo...\n")

expr_raw_preimputation <- as.matrix(raw[, intensity_cols])
total_na <- sum(is.na(expr_raw_preimputation))
cat(sprintf("  Total NAs en matriz pre-imputación: %d (de %d valores)\n",
            total_na, length(expr_raw_preimputation)))
if (total_na == 0) {
  cat("  ADVERTENCIA: No se encontraron NAs — los datos parecen ya estar imputados.\n")
  cat("  El análisis de missingness puede ser no informativo.\n")
}
rownames(expr_raw_preimputation) <- raw$gene_symbol

normal_cols  <- meta$sample_id[meta$condition == "NORMAL"]
tumor_cols   <- meta$sample_id[meta$condition == "TUMORAL"]

missing_df <- data.frame(
  gene_symbol   = raw$gene_symbol,
  pct_missing_normal = rowMeans(is.na(expr_raw_preimputation[, normal_cols])) * 100,
  pct_missing_tumor  = rowMeans(is.na(expr_raw_preimputation[, tumor_cols]))  * 100
) %>%
  mutate(
    missing_diff = pct_missing_tumor - pct_missing_normal,
    bias_flag    = abs(missing_diff) > 30  # diferencia >30% entre grupos = sospechoso
  )

n_biased <- sum(missing_df$bias_flag, na.rm = TRUE)
cat(sprintf("  Proteínas con diferencial missingness >30%% entre grupos: %d / %d\n",
            n_biased, nrow(missing_df)))

write.table(missing_df,
  file.path(out_tables, "01_missingness_analysis.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# 4. MODELO LIMMA HPV-AJUSTADO
#    ~ condition + vph_status
#    bloqueo intra-paciente (duplicateCorrelation)
# ============================================================
cat("-- Construyendo modelo limma con covariable VPH...\n")

params <- yaml::read_yaml(file.path(proj_dir, "config/analysis_params.yaml"))
fc_thresh  <- params$de$log2fc_threshold      # 1.0
pval_thr   <- params$de$adj_pval_threshold    # 0.05

# Matriz de expresión: proteínas × muestras (ya log2, imputado, normalizado)
expr_mat <- as.matrix(raw[, intensity_cols])
rownames(expr_mat) <- raw$gene_symbol

# Factores del diseño
condition <- factor(meta$condition, levels = c("NORMAL", "TUMORAL"))
vph       <- factor(meta$vph,       levels = c("NEGATIVE", "POSITIVE"))
patient   <- factor(meta$patient_id)

# Diseño sin intercepto para facilitar makeContrasts
design <- model.matrix(~ 0 + condition + vph, data = data.frame(condition, vph))
cat("  Columnas del diseño (antes de renombrar):", paste(colnames(design), collapse = ", "), "\n")
colnames(design) <- c("NORMAL", "TUMORAL", "VPH_POSITIVE")
stopifnot("Diseño tiene número inesperado de columnas (esperado: 3)" = ncol(design) == 3)
cat(sprintf("  Diseño: %d muestras × %d coeficientes\n", nrow(design), ncol(design)))
cat("  Coeficientes:", paste(colnames(design), collapse = ", "), "\n")

# Correlación inter-paciente (diseño pareado)
cat("  Estimando correlación inter-paciente (duplicateCorrelation)...\n")
corfit <- duplicateCorrelation(expr_mat, design, block = patient)
cat(sprintf("  Correlación intra-paciente (consensus): %.3f\n",
            corfit$consensus.correlation))

# Ajuste del modelo con bloqueo por paciente
fit <- lmFit(expr_mat, design,
             block       = patient,
             correlation = corfit$consensus.correlation)

# Contraste: Tumoral vs Normal (ajustado por VPH)
contrasts <- makeContrasts(TumoralVsNormal = TUMORAL - NORMAL, levels = design)
fit2      <- contrasts.fit(fit, contrasts)
fit2      <- eBayes(fit2, trend = TRUE)  # trend=TRUE para proteómica (heterocedasticidad)

# Extraer resultados completos
results_all <- topTable(fit2, coef = "TumoralVsNormal",
                        n = Inf, adjust.method = "BH",
                        sort.by = "none")
results_all$gene_symbol <- rownames(results_all)

# Unir con anotaciones originales (uniprot_id, entry_name)
annot_cols <- raw %>% select(gene_symbol, uniprot_id, entry_name)
results_all <- results_all %>%
  left_join(annot_cols, by = "gene_symbol") %>%
  # Renombrar columnas para compatibilidad downstream
  # Nota: topTable() no devuelve CI.L / CI.R — se omiten
  rename(
    logFC_TVsS         = logFC,
    avg_intensity_TVsS = AveExpr,
    t_stat_TVsS        = t,
    B_stat_TVsS        = B,
    P.Value_TVsS       = P.Value,
    adj.P.Val_TVsS     = adj.P.Val
  ) %>%
  mutate(
    sig.FDR_TVsS = case_when(
      adj.P.Val_TVsS < pval_thr & logFC_TVsS >  fc_thresh ~  1L,
      adj.P.Val_TVsS < pval_thr & logFC_TVsS < -fc_thresh ~ -1L,
      TRUE                                                ~  0L
    )
  ) %>%
  select(gene_symbol, uniprot_id, entry_name,
         logFC_TVsS, avg_intensity_TVsS,
         t_stat_TVsS, B_stat_TVsS, P.Value_TVsS, adj.P.Val_TVsS,
         sig.FDR_TVsS)

# Añadir columnas de intensidades para downstream (script 02 las espera)
results_all <- results_all %>%
  left_join(raw %>% select(gene_symbol, all_of(intensity_cols)), by = "gene_symbol")

# Resumen
n_up   <- sum(results_all$sig.FDR_TVsS ==  1, na.rm = TRUE)
n_down <- sum(results_all$sig.FDR_TVsS == -1, na.rm = TRUE)
n_ns   <- sum(results_all$sig.FDR_TVsS == 0, na.rm = TRUE)

cat("\n=== Resultados limma HPV-ajustado ===\n")
cat(sprintf("  Total proteínas:    %d\n", nrow(results_all)))
cat(sprintf("  Upregulated:        %d\n", n_up))
cat(sprintf("  Downregulated:      %d\n", n_down))
cat(sprintf("  No significativas:  %d\n", n_ns))
cat(sprintf("  Criterio: adj.P.Val < %.2f & |logFC| > %.1f (BH)\n", pval_thr, fc_thresh))

sig_fc <- results_all$logFC_TVsS[!is.na(results_all$sig.FDR_TVsS) & results_all$sig.FDR_TVsS != 0]
if (length(sig_fc) > 0) {
  cat(sprintf("  logFC rango sig: %.2f a %.2f\n\n",
              min(sig_fc, na.rm = TRUE), max(sig_fc, na.rm = TRUE)))
} else {
  cat("  logFC rango sig: ninguna proteína significativa\n\n")
}

# ============================================================
# 5. EXPORTAR TABLAS
# ============================================================
cat("-- Exportando tablas...\n")

out_cols <- c("gene_symbol", "uniprot_id", "entry_name",
              "logFC_TVsS", "avg_intensity_TVsS",
              "t_stat_TVsS", "B_stat_TVsS", "P.Value_TVsS", "adj.P.Val_TVsS",
              "sig.FDR_TVsS", intensity_cols)

df_all  <- results_all
df_sig  <- results_all %>% filter(sig.FDR_TVsS %in% c(1L, -1L))
df_up   <- results_all %>% filter(sig.FDR_TVsS ==  1L)
df_down <- results_all %>% filter(sig.FDR_TVsS == -1L)

for (pair in list(
  list(df_all,  "01_TVsS_all_proteins.tsv"),
  list(df_sig,  "01_TVsS_significant.tsv"),
  list(df_up,   "01_TVsS_upregulated.tsv"),
  list(df_down, "01_TVsS_downregulated.tsv")
)) {
  write.table(pair[[1]], file.path(out_tables, pair[[2]]),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  %s: %d proteínas\n", pair[[2]], nrow(pair[[1]])))
}

# ============================================================
# 6. FIGURA QC-1: BOXPLOT DE INTENSIDADES POR MUESTRA
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
# 7. FIGURA QC-2: PCA DE MUESTRAS
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
# 8. FIGURA QC-3: VOLCANO PLOT TVsS
# ============================================================
cat("-- Generando figura QC-3: volcano plot TVsS...\n")

# Preparar dataframe de volcáno con categoría DE
df_volc <- results_all %>%
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
# 9. FIGURA QC-4: RESUMEN CONTEOS DE/BARRAS
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
# FIGURA QC-5: Missingness vs dirección DE (sesgo MinProb)
# ============================================================
cat("-- Generando figura QC-5: missingness vs dirección DE...\n")

miss_de <- missing_df %>%
  left_join(results_all %>% select(gene_symbol, logFC_TVsS, sig.FDR_TVsS),
            by = "gene_symbol") %>%
  mutate(
    DE_status = case_when(
      sig.FDR_TVsS ==  1 ~ "Up",
      sig.FDR_TVsS == -1 ~ "Down",
      TRUE               ~ "NS"
    )
  )

p_miss <- ggplot(miss_de,
                 aes(x = pct_missing_tumor, y = pct_missing_normal,
                     color = DE_status, alpha = DE_status)) +
  geom_point(size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Up" = "#D6604D", "Down" = "#4393C3", "NS" = "grey60")) +
  scale_alpha_manual(values = c("Up" = 0.9, "Down" = 0.9, "NS" = 0.2)) +
  labs(
    title    = "Missingness por grupo — verificación sesgo MinProb",
    subtitle = "Puntos sobre la diagonal = más missing en tumor (posible artefacto down-regulation)",
    x = "% Missing en TUMORAL", y = "% Missing en NORMAL",
    color = "DE status"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  guides(alpha = "none")

ggsave(file.path(out_figures, "01_missingness_vs_direction.pdf"),
       p_miss, width = 6.5, height = 6)
cat("  Guardado: 01_missingness_vs_direction.pdf\n")

# ============================================================
# 10. RESUMEN FINAL
# ============================================================
cat("=== RESUMEN FINAL ===\n")
cat("Proteínas totales cuantificadas :", nrow(results_all), "\n")
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
