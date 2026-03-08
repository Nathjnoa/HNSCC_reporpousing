#!/usr/bin/env Rscript
# =============================================================================
# 13_evidence_summary.R
# =============================================================================
# Compila toda la evidencia y genera el reporte final de candidatos.
# Integra: scoring (script 10) + ensayos clínicos (script 11) +
#          cancer drivers (script 12) + clase de reposicionamiento (script 08)
#
# Outputs:
#   results/tables/13_FINAL_drug_candidates.xlsx  — tabla maestra con 5 hojas
#   results/tables/13_evidence_matrix.tsv         — matriz de evidencia binaria
#   results/figures/final/13_evidence_heatmap.pdf — heatmap de evidencia
#   results/figures/final/13_final_ranking.pdf    — lollipop ranking final
#   results/figures/final/13_evidence_radar.pdf   — radar de top 5
#   results/figures/final/13_multipanel_summary.pdf — figura multipanel
#
# Ejecución:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/13_evidence_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(openxlsx)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(scales)
  library(stringr)
})

# ── Detectar directorio del proyecto ─────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  proj_dir    <- dirname(dirname(script_path))
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml")))
    setwd(proj_dir)
}

cat("Working directory:", getwd(), "\n")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file  <- file.path("logs", paste0("13_evidence_summary_", timestamp, ".log"))
dir.create("logs", showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "[%H:%M:%S]"), " ", ..., "\n")
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
}

log_msg("=" , strrep("=", 58))
log_msg("Script 13: Evidence summary y reporte final")
log_msg("=" , strrep("=", 58))

# ── Crear directorios de output ───────────────────────────────────────────────
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/final", recursive = TRUE, showWarnings = FALSE)

# ── Cargar todos los inputs ───────────────────────────────────────────────────
log_msg("Cargando inputs...")

# Script 10: scoring y top 20
top20 <- tryCatch(
  read_tsv("results/tables/10_top20_candidates.tsv", show_col_types = FALSE),
  error = function(e) { log_msg("ERROR: no encontrado 10_top20_candidates.tsv"); stop(e) }
)
all_scored <- tryCatch(
  read_tsv("results/tables/10_all_candidates_scored.tsv", show_col_types = FALSE),
  error = function(e) { log_msg("WARN: 10_all_candidates_scored.tsv no encontrado"); top20 }
)

# Script 11: evidencia clínica
clinical <- tryCatch(
  read_tsv("results/tables/evidence/11_clinical_evidence.tsv", show_col_types = FALSE),
  error = function(e) { log_msg("WARN: 11_clinical_evidence.tsv no encontrado"); NULL }
)

# Script 12: cancer drivers
drivers <- tryCatch(
  read_tsv("results/tables/evidence/12_top20_with_drivers.tsv", show_col_types = FALSE),
  error = function(e) { log_msg("WARN: 12_top20_with_drivers.tsv no encontrado"); NULL }
)

log_msg("Top 20 candidatos: ", nrow(top20), " filas")
log_msg("Columnas top20: ", paste(names(top20), collapse = ", "))

# Identificar columna de nombre de fármaco
name_col <- intersect(c("drug_name_norm", "drug_name", "name"), names(top20))[1]
log_msg("Columna nombre fármaco: ", name_col)

# ── Integrar evidencia en tabla maestra ───────────────────────────────────────
log_msg("Integrando evidencia...")

df_master <- top20 %>%
  rename(drug_name = all_of(name_col))

# Añadir evidencia clínica
if (!is.null(clinical)) {
  clin_cols <- c("drug_name_norm", "ct_hnscc_trials", "ct_active_trials",
                 "pubmed_hnscc", "evidence_score")
  clin_cols <- intersect(clin_cols, names(clinical))
  clin_sub  <- clinical %>%
    select(all_of(clin_cols)) %>%
    rename(drug_name = drug_name_norm)

  df_master <- df_master %>%
    left_join(clin_sub, by = "drug_name")
  log_msg("Evidencia clínica integrada: ", nrow(clin_sub), " filas")
}

# Añadir info de cancer drivers
if (!is.null(drivers)) {
  drv_cols <- c("drug_name_norm", "n_target_drivers", "has_driver_target")
  drv_cols <- intersect(drv_cols, names(drivers))
  if (length(drv_cols) >= 2) {
    drv_sub  <- drivers %>%
      select(all_of(drv_cols)) %>%
      rename(drug_name = if ("drug_name_norm" %in% drv_cols) "drug_name_norm" else drv_cols[1])
    df_master <- df_master %>%
      left_join(drv_sub, by = "drug_name")
    log_msg("Cancer driver info integrada")
  }
}

# Rellenar NAs con 0
df_master <- df_master %>%
  mutate(
    ct_hnscc_trials   = replace_na(as.integer(ct_hnscc_trials), 0L),
    ct_active_trials  = replace_na(as.integer(ct_active_trials), 0L),
    pubmed_hnscc      = replace_na(as.integer(pubmed_hnscc), 0L),
    evidence_score    = replace_na(as.numeric(evidence_score), 0),
    n_target_drivers  = replace_na(as.integer(n_target_drivers), 0L),
    has_driver_target = replace_na(as.logical(has_driver_target), FALSE)
  )

# ── Calcular nivel de evidencia ───────────────────────────────────────────────
log_msg("Calculando niveles de evidencia...")

# Clase de reposicionamiento (A/B/C/D)
class_col <- intersect(c("drug_class", "repurposing_class", "class", "clase"),
                       names(df_master))[1]
log_msg("Columna clase: ", class_col)

df_master <- df_master %>%
  mutate(
    repurposing_class = if (!is.na(class_col)) .data[[class_col]] else "D",
    has_hnscc_trial   = ct_hnscc_trials > 0,
    has_active_trial  = ct_active_trials > 0,
    has_pubmed        = pubmed_hnscc > 0,
    # Score compuesto de final_score
    scoring_rank      = rank(-final_score, ties.method = "min"),
    # Nivel de evidencia 1-4
    evidence_level = case_when(
      repurposing_class == "A"                               ~ 1L,
      repurposing_class == "B" & has_hnscc_trial             ~ 2L,
      repurposing_class == "B"                               ~ 2L,
      repurposing_class == "C" & has_hnscc_trial & has_pubmed ~ 3L,
      repurposing_class == "C"                               ~ 3L,
      TRUE                                                   ~ 4L
    ),
    evidence_label = case_when(
      evidence_level == 1L ~ "Nivel 1: Aprobado HNSCC",
      evidence_level == 2L ~ "Nivel 2: Aprobado otro cáncer",
      evidence_level == 3L ~ "Nivel 3: Aprobado no-oncológico",
      evidence_level == 4L ~ "Nivel 4: Experimental",
    )
  )

log_msg("Distribución de niveles de evidencia:")
log_msg(capture.output(print(table(df_master$evidence_level))) |> paste(collapse = "\n"))

# Ranking final combinando scoring con evidencia clínica
df_master <- df_master %>%
  mutate(
    combined_rank_score = final_score * 0.6 + evidence_score * 0.4
  ) %>%
  arrange(desc(combined_rank_score)) %>%
  mutate(final_rank = row_number())

# ── Matriz de evidencia binaria ───────────────────────────────────────────────
log_msg("Construyendo matriz de evidencia...")

evidence_matrix <- df_master %>%
  select(drug_name, final_rank, evidence_level) %>%
  mutate(
    `Proteómica DE`       = TRUE,
    `Multi-fuente DB`     = if ("n_sources" %in% names(df_master))
                              df_master$n_sources >= 2 else TRUE,
    `CMap reversor`       = if ("s_cmap" %in% names(df_master))
                              df_master$s_cmap > 0 else FALSE,
    `Pathway relevante`   = if ("s_pathway" %in% names(df_master))
                              df_master$s_pathway > 0 else FALSE,
    `Hub PPI`             = if ("s_network" %in% names(df_master))
                              df_master$s_network > 0.3 else FALSE,
    `Fase clínica ≥ 3`    = if ("max_phase" %in% names(df_master))
                              (df_master$max_phase >= 3) else FALSE,
    `Trial HNSCC`         = df_master$has_hnscc_trial,
    `Trial HNSCC activo`  = df_master$has_active_trial,
    `PubMed HNSCC`        = df_master$has_pubmed,
    `Target cancer driver`= df_master$has_driver_target,
  )

write_tsv(evidence_matrix, "results/tables/13_evidence_matrix.tsv")
log_msg("Matriz de evidencia guardada")

# ── Heatmap de evidencia ──────────────────────────────────────────────────────
log_msg("Generando heatmap de evidencia...")

bin_cols <- c("Proteómica DE", "Multi-fuente DB", "CMap reversor",
              "Pathway relevante", "Hub PPI", "Fase clínica ≥ 3",
              "Trial HNSCC", "Trial HNSCC activo", "PubMed HNSCC",
              "Target cancer driver")

mat_bin <- evidence_matrix %>%
  select(all_of(bin_cols)) %>%
  mutate(across(everything(), as.integer)) %>%
  as.matrix()
rownames(mat_bin) <- str_trunc(df_master$drug_name, 22)

# Colores por nivel de evidencia
level_colors <- c(
  "1" = "#d62728", "2" = "#ff7f0e", "3" = "#2ca02c", "4" = "#aec7e8"
)
level_labels <- c("1" = "Nivel 1", "2" = "Nivel 2",
                  "3" = "Nivel 3", "4" = "Nivel 4")

row_ha <- rowAnnotation(
  Nivel = as.character(df_master$evidence_level),
  col = list(Nivel = level_colors),
  annotation_legend_param = list(
    Nivel = list(labels = level_labels)
  )
)

col_fun <- colorRamp2(c(0, 1), c("white", "#1a6db3"))

pdf("results/figures/final/13_evidence_heatmap.pdf", width = 11, height = 8)
ht <- Heatmap(
  mat_bin,
  name = "Evidencia",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 40,
  right_annotation = row_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (mat_bin[i, j] == 1)
      grid.text("✓", x, y, gp = gpar(fontsize = 9, col = "white", fontface = "bold"))
  },
  column_title = "Dimensiones de evidencia",
  row_title = "Top 20 candidatos (ordenados por ranking final)",
  heatmap_legend_param = list(at = c(0, 1), labels = c("No", "Sí")),
)
draw(ht)
dev.off()
log_msg("Heatmap guardado: results/figures/final/13_evidence_heatmap.pdf")

# ── Lollipop ranking final ────────────────────────────────────────────────────
log_msg("Generando lollipop de ranking final...")

level_pal <- c("1" = "#d62728", "2" = "#ff7f0e", "3" = "#2ca02c", "4" = "#aec7e8")

p_lollipop <- df_master %>%
  mutate(
    drug_label = str_trunc(drug_name, 25),
    drug_label = factor(drug_label, levels = rev(drug_label)),
    nivel_f    = factor(as.character(evidence_level), levels = c("1","2","3","4"))
  ) %>%
  ggplot(aes(x = combined_rank_score, y = drug_label)) +
  geom_segment(aes(x = 0, xend = combined_rank_score, yend = drug_label),
               color = "grey70", linewidth = 0.7) +
  geom_point(aes(color = nivel_f, size = final_score), alpha = 0.9) +
  geom_text(aes(label = sprintf("%.3f", combined_rank_score)),
            hjust = -0.2, size = 2.8, color = "grey30") +
  scale_color_manual(
    values = level_pal,
    labels = c("1" = "Nivel 1: HNSCC", "2" = "Nivel 2: Otro cáncer",
               "3" = "Nivel 3: No-oncológico", "4" = "Nivel 4: Experimental"),
    name = "Nivel evidencia"
  ) +
  scale_size_continuous(range = c(3, 8), name = "Scoring\nmulti-criterio") +
  scale_x_continuous(limits = c(0, max(df_master$combined_rank_score) * 1.15)) +
  labs(
    title    = "Ranking final: Top 20 candidatos HNSCC",
    subtitle = "Score combinado = 0.60 × scoring multi-criterio + 0.40 × evidencia clínica",
    x        = "Score combinado final",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.box         = "horizontal",
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey40", size = 9)
  )

ggsave("results/figures/final/13_final_ranking.pdf",
       p_lollipop, width = 10, height = 8, device = "pdf")
log_msg("Lollipop guardado: results/figures/final/13_final_ranking.pdf")

# ── Radar chart de top 5 ─────────────────────────────────────────────────────
log_msg("Generando radar chart de top 5...")

# Columnas de score parcial
score_cols <- c("s_logfc", "s_sig", "s_clinical", "s_cmap", "s_pathway", "s_network")
score_cols <- intersect(score_cols, names(df_master))
score_labels <- c(
  s_logfc    = "|logFC|",
  s_sig      = "Significancia",
  s_clinical = "Fase clínica",
  s_cmap     = "L2S2",
  s_pathway  = "Pathway",
  s_network  = "Red PPI"
)

if (length(score_cols) >= 3) {
  top5 <- df_master %>% slice_head(n = 5)

  radar_data <- top5 %>%
    select(drug_name, all_of(score_cols)) %>%
    pivot_longer(-drug_name, names_to = "dimension", values_to = "score") %>%
    mutate(
      drug_label = str_trunc(drug_name, 22),
      dim_label  = dplyr::coalesce(score_labels[dimension], dimension)
    )

  p_radar <- ggplot(radar_data, aes(x = dim_label, y = score,
                                    fill = drug_label)) +
    geom_col(position = "dodge", color = "white", linewidth = 0.3) +
    scale_fill_brewer(palette = "Dark2", name = "Fármaco") +
    scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
    labs(
      title    = "Perfil de scoring — Top 5 candidatos",
      subtitle = "Valores normalizados 0-1 por dimensión",
      x        = NULL, y = "Score"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(face = "bold", size = 12),
      axis.text.x     = element_text(size = 9)
    )

  ggsave("results/figures/final/13_evidence_radar.pdf",
         p_radar, width = 9, height = 5, device = "pdf")
  log_msg("Scoring bars guardado: results/figures/final/13_evidence_radar.pdf")
} else {
  log_msg("WARN: columnas de score parcial no encontradas; omitiendo figura de scores")
}

# ── Figura resumen de estadísticas ────────────────────────────────────────────
log_msg("Generando figura estadísticas de candidatos...")

# Panel A: distribución de niveles de evidencia
p_levels <- df_master %>%
  count(evidence_level, evidence_label) %>%
  mutate(
    nivel_f = factor(as.character(evidence_level), levels = c("1","2","3","4"))
  ) %>%
  ggplot(aes(x = nivel_f, y = n, fill = nivel_f)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = level_pal, guide = "none") +
  labs(title = "Niveles de evidencia (Top 20)",
       x = "Nivel", y = "N candidatos") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

# Panel B: trials HNSCC vs PubMed HNSCC
p_evidence <- df_master %>%
  mutate(drug_label = str_trunc(drug_name, 18)) %>%
  ggplot(aes(x = pubmed_hnscc, y = ct_hnscc_trials,
             color = factor(as.character(evidence_level)),
             size = final_score)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = level_pal, name = "Nivel") +
  scale_size_continuous(range = c(3, 10), name = "Scoring") +
  scale_x_log10(labels = comma, breaks = c(1, 10, 100, 1000, 10000)) +
  labs(
    title = "Ensayos HNSCC vs publicaciones",
    x     = "PubMed HNSCC (log10)",
    y     = "Ensayos ClinicalTrials (HNSCC)"
  ) +
  theme_minimal(base_size = 11) +
  annotation_logticks(sides = "b")

# Panel C: clase de reposicionamiento
p_class <- df_master %>%
  count(repurposing_class) %>%
  ggplot(aes(x = repurposing_class, y = n,
             fill = repurposing_class)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("A"="#d62728","B"="#ff7f0e","C"="#2ca02c","D"="#aec7e8"),
                    guide = "none") +
  labs(title = "Clase de reposicionamiento",
       x = "Clase", y = "N candidatos") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank())

p_multi <- (p_levels | p_class) / p_evidence +
  plot_annotation(
    title   = "Resumen de evidencia — Top 20 candidatos HNSCC",
    subtitle = paste0("Pipeline: proteómica DIA → 4 DBs → red PPI → scoring multi-criterio → evidencia clínica"),
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("results/figures/final/13_multipanel_summary.pdf",
       p_multi, width = 12, height = 9, device = "pdf")
log_msg("Multipanel guardado: results/figures/final/13_multipanel_summary.pdf")

# ── Exportar Excel con 5 hojas ────────────────────────────────────────────────
log_msg("Exportando Excel...")

wb <- createWorkbook()
addWorksheet(wb, "Top20_Final")
addWorksheet(wb, "Todos_Candidatos")
addWorksheet(wb, "Evidencia_Matriz")
addWorksheet(wb, "Ensayos_Clinicos")
addWorksheet(wb, "Metodologia")

# Estilos
header_style <- createStyle(
  fontColour = "white", bgFill = "#1a4a7a",
  textDecoration = "Bold", halign = "center",
  border = "Bottom", borderColour = "white"
)
level1_style <- createStyle(bgFill = "#ffd7d7")
level2_style <- createStyle(bgFill = "#ffe4cc")
level3_style <- createStyle(bgFill = "#d5f0d5")
level4_style <- createStyle(bgFill = "#e8f0fe")

# Hoja 1: Top 20 final
writeData(wb, "Top20_Final", df_master %>% arrange(final_rank))
addStyle(wb, "Top20_Final", header_style, rows = 1,
         cols = 1:ncol(df_master), gridExpand = TRUE)

# Colorear por nivel
for (i in 1:nrow(df_master)) {
  lvl <- df_master$evidence_level[i]
  sty <- switch(as.character(lvl),
    "1" = level1_style, "2" = level2_style,
    "3" = level3_style, "4" = level4_style, level4_style
  )
  addStyle(wb, "Top20_Final", sty, rows = i + 1,
           cols = 1:ncol(df_master), gridExpand = TRUE, stack = TRUE)
}
setColWidths(wb, "Top20_Final", cols = 1:ncol(df_master), widths = "auto")

# Hoja 2: Todos los candidatos scored
all_out <- if (nrow(all_scored) > 0) all_scored else df_master
writeData(wb, "Todos_Candidatos", all_out)
addStyle(wb, "Todos_Candidatos", header_style, rows = 1,
         cols = 1:ncol(all_out), gridExpand = TRUE)
setColWidths(wb, "Todos_Candidatos", cols = 1:ncol(all_out), widths = "auto")

# Hoja 3: Matriz de evidencia
writeData(wb, "Evidencia_Matriz", evidence_matrix)
addStyle(wb, "Evidencia_Matriz", header_style, rows = 1,
         cols = 1:ncol(evidence_matrix), gridExpand = TRUE)
setColWidths(wb, "Evidencia_Matriz", cols = 1:ncol(evidence_matrix), widths = "auto")

# Hoja 4: Ensayos clínicos
if (!is.null(clinical)) {
  writeData(wb, "Ensayos_Clinicos", clinical)
  addStyle(wb, "Ensayos_Clinicos", header_style, rows = 1,
           cols = 1:ncol(clinical), gridExpand = TRUE)
  setColWidths(wb, "Ensayos_Clinicos", cols = 1:ncol(clinical), widths = "auto")
} else {
  writeData(wb, "Ensayos_Clinicos",
            data.frame(nota = "Datos no disponibles — ejecutar script 11"))
}

# Hoja 5: Metodología
metodo <- data.frame(
  Fase        = c("1", "2", "3", "3", "3", "3", "3", "4", "4", "5", "5", "5"),
  Script      = c("01-02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13"),
  Descripcion = c(
    "Parseo proteómica, QC, mapeo ID",
    "Enriquecimiento GO/KEGG/Reactome/Hallmarks",
    "DGIdb: interacciones gen-fármaco",
    "ChEMBL: fármacos fase ≥ 3",
    "Open Targets: evidencia HNSCC",
    "L2S2: reversores firma tumoral",
    "Integración 4 fuentes, clasificación A/B/C/D",
    "Red PPI STRING (score≥700), hubs druggables",
    "Scoring multi-criterio, Top 20 candidatos",
    "ClinicalTrials API v2 + PubMed",
    "COSMIC CGC / NCG7 / literatura HNSCC",
    "Reporte final (este archivo)"
  ),
  N_Resultados = c(
    "3,352 proteínas; 520 DE",
    "GO:77, KEGG:16, Reactome:14, Hallmarks:15",
    "2,846 interacciones; 2,252 fármacos",
    "90 pares gen-fármaco fase≥3",
    "354/520 genes con evidencia HNSCC",
    "1,044 drugs FDA-aprobados",
    "2,421 fármacos; 187 multi-fuente",
    "403 nodos; 2,001 aristas; 41 hubs",
    "187 candidatos; Top 20 de 9 targets",
    "8/20 con trials HNSCC; 6 activos",
    "2 genes driver DE; EGFR validado",
    "Top 20 con evidencia integrada"
  )
)
writeData(wb, "Metodologia", metodo)
addStyle(wb, "Metodologia", header_style, rows = 1, cols = 1:4, gridExpand = TRUE)
setColWidths(wb, "Metodologia", cols = 1:4, widths = c(8, 10, 50, 35))

saveWorkbook(wb, "results/tables/13_FINAL_drug_candidates.xlsx", overwrite = TRUE)
log_msg("Excel guardado: results/tables/13_FINAL_drug_candidates.xlsx")

# ── Resumen final ─────────────────────────────────────────────────────────────
log_msg("")
log_msg("=" , strrep("=", 58))
log_msg("RESUMEN FINAL")
log_msg("=" , strrep("=", 58))
log_msg("Candidatos Top 20 integrados: ", nrow(df_master))
log_msg("Distribución por nivel de evidencia:")
tbl <- table(df_master$evidence_level)
for (n in names(tbl)) log_msg("  Nivel ", n, ": ", tbl[[n]])
log_msg("")
log_msg("Top 10 candidatos (score combinado):")
print_cols <- intersect(
  c("drug_name", "evidence_level", "combined_rank_score",
    "final_score", "ct_hnscc_trials", "pubmed_hnscc"),
  names(df_master)
)
top10_print <- df_master %>% slice_head(n = 10) %>% select(all_of(print_cols))
log_msg(capture.output(print(as.data.frame(top10_print), row.names = FALSE)) |>
        paste(collapse = "\n"))
log_msg("")
log_msg("Outputs:")
log_msg("  results/tables/13_FINAL_drug_candidates.xlsx")
log_msg("  results/tables/13_evidence_matrix.tsv")
log_msg("  results/figures/final/13_evidence_heatmap.pdf")
log_msg("  results/figures/final/13_final_ranking.pdf")
log_msg("  results/figures/final/13_multipanel_summary.pdf")
log_msg("Script 13 completado.")
