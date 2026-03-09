#!/usr/bin/env Rscript
# =============================================================================
# Script 15: Análisis de sensibilidad de pesos del scoring
# HNSCC Drug Repurposing — Validación de robustez
# =============================================================================
# Objetivo: Verificar que el ranking del Top 20 es robusto ante cambios en los
#           pesos del scoring compuesto. Reutiliza los 6 componentes individuales
#           de score ya calculados en script 10, sin re-ejecutar el pipeline.
#
# Input:  results/tables/10_all_candidates_scored.tsv  (tiene s_logfc..s_network)
#         config/analysis_params.yaml                  (para re-aplicar exclusiones)
#
# Output: results/tables/15_sensitivity_ranks.tsv
#         results/figures/15_rank_heatmap.pdf
#         results/figures/15_stability_bar.pdf
#         results/figures/15_score_distribution.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(yaml)
})

# --- Working directory -------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  proj_dir    <- dirname(dirname(script_path))
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml")))
    setwd(proj_dir)
}
cat("Working dir:", getwd(), "\n")

# --- Log setup ---------------------------------------------------------------
dir.create("logs",   showWarnings = FALSE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

log_file <- paste0("logs/15_sensitivity_analysis_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
con_log <- file(log_file, open = "wt")
sink(con_log, type = "output")
sink(con_log, type = "message", append = TRUE)

cat("=== 15_sensitivity_analysis.R ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# =============================================================================
# CONFIGURACIONES DE PESOS
# =============================================================================
weight_configs <- list(
  baseline       = c(logfc=0.20, sig=0.15, clinical=0.20, cmap=0.15, pathway=0.15, network=0.15),
  clinical_heavy = c(logfc=0.10, sig=0.10, clinical=0.45, cmap=0.10, pathway=0.15, network=0.10),
  molecular_heavy= c(logfc=0.35, sig=0.25, clinical=0.15, cmap=0.10, pathway=0.10, network=0.05),
  network_heavy  = c(logfc=0.15, sig=0.10, clinical=0.20, cmap=0.10, pathway=0.10, network=0.35),
  pathway_heavy  = c(logfc=0.15, sig=0.10, clinical=0.20, cmap=0.10, pathway=0.35, network=0.10),
  equal_weights  = c(logfc=1/6,  sig=1/6,  clinical=1/6,  cmap=1/6,  pathway=1/6,  network=1/6)
)

cat("Configuraciones de pesos:\n")
for (nm in names(weight_configs)) {
  w <- weight_configs[[nm]]
  cat(sprintf("  %-15s: logFC=%.2f sig=%.2f clin=%.2f cmap=%.2f path=%.2f net=%.2f (sum=%.2f)\n",
              nm, w["logfc"], w["sig"], w["clinical"], w["cmap"], w["pathway"], w["network"],
              sum(w)))
}

# =============================================================================
# CARGAR DATOS
# =============================================================================
cat("\n--- Cargando datos ---\n")

params  <- yaml::read_yaml("config/analysis_params.yaml")
top_n   <- params$candidates$top_n

all_scored <- read.delim("results/tables/10_all_candidates_scored.tsv",
                         stringsAsFactors = FALSE)
cat(sprintf("Candidatos cargados: %d\n", nrow(all_scored)))

# Re-aplicar exclusiones (igual que script 10)
excluded_drugs  <- toupper(unlist(params$exclusions$drugs))
excl_targets    <- toupper(unlist(params$exclusions$exclude_single_target_genes))

# Verificar que columna 'bonus' existe (calculada y exportada por script 10)
if (!"bonus" %in% colnames(all_scored)) {
  cat("WARN: columna 'bonus' no encontrada en 10_all_candidates_scored.tsv — usando 0\n")
  all_scored$bonus <- 0
}

all_scored <- all_scored %>%
  # drug_name_norm está en minúsculas; excluded_drugs en mayúsculas → toupper antes de comparar
  filter(!toupper(drug_name_norm) %in% excluded_drugs) %>%
  # de_genes puede contener múltiples genes separados por "|"; comparar gen a gen
  filter(!sapply(de_genes, function(g) {
    if (is.na(g) || g == "") return(FALSE)
    any(toupper(str_split(g, "\\|")[[1]]) %in% excl_targets)
  }))

cat(sprintf("Candidatos tras exclusiones: %d\n", nrow(all_scored)))

# Extraer gen diana primario (primer gen en de_genes)
all_scored <- all_scored %>%
  mutate(primary_target = sapply(de_genes, function(g) {
    if (is.na(g) || g == "") return("unknown")
    str_split(g, "\\|")[[1]][1]
  }))

# =============================================================================
# HELPER: Top N con diversidad (max 3 por target primario)
# =============================================================================
get_top_diverse <- function(df, score_col, n = top_n) {
  df <- df %>% arrange(desc(.data[[score_col]]))
  df %>%
    group_by(primary_target) %>%
    slice_head(n = 3) %>%
    ungroup() %>%
    arrange(desc(.data[[score_col]])) %>%
    slice_head(n = n)
}

# =============================================================================
# CALCULAR RANKING PARA CADA CONFIGURACION
# =============================================================================
cat("\n--- Calculando rankings para 6 configuraciones ---\n")

score_cols <- c("s_logfc", "s_sig", "s_clinical", "s_cmap", "s_pathway", "s_network")

rank_list <- list()

for (config_name in names(weight_configs)) {
  w <- weight_configs[[config_name]]

  scored_cfg <- all_scored %>%
    mutate(
      composite_score_cfg = w["logfc"]   * s_logfc   +
                            w["sig"]     * s_sig     +
                            w["clinical"]* s_clinical +
                            w["cmap"]    * s_cmap    +
                            w["pathway"] * s_pathway +
                            w["network"] * s_network,
      final_score_cfg = pmin(composite_score_cfg + bonus, 1.0)
    )

  top20_cfg <- get_top_diverse(scored_cfg, "final_score_cfg", top_n)

  # Guardar rank para todos los candidatos que quedaron en top 20
  ranks_cfg <- top20_cfg %>%
    arrange(desc(final_score_cfg)) %>%
    mutate(rank = row_number(),
           config = config_name,
           score  = final_score_cfg) %>%
    select(drug_name_norm, config, rank, score)

  rank_list[[config_name]] <- ranks_cfg

  cat(sprintf("  %s: top1=%s (%.4f), top5 incluye: %s\n",
              config_name,
              ranks_cfg$drug_name_norm[1],
              ranks_cfg$score[1],
              paste(ranks_cfg$drug_name_norm[1:min(5, nrow(ranks_cfg))], collapse=", ")))
}

# Tabla consolidada larga
rank_long <- bind_rows(rank_list)

# Tabla ancha: candidato × config → rank
all_candidates_ranked <- all_scored$drug_name_norm
rank_wide <- rank_long %>%
  pivot_wider(id_cols = drug_name_norm,
              names_from = config,
              values_from = rank)

# Candidatos en alguna config top 20
candidates_in_any <- unique(rank_long$drug_name_norm)

rank_wide_filtered <- rank_wide %>%
  filter(drug_name_norm %in% candidates_in_any) %>%
  mutate(n_configs_top20 = rowSums(!is.na(select(., -drug_name_norm)))) %>%
  arrange(desc(n_configs_top20), baseline)

cat(sprintf("\nCandidatos en al menos 1 configuracion top 20: %d\n", nrow(rank_wide_filtered)))

# =============================================================================
# ANÁLISIS DE ROBUSTEZ
# =============================================================================
cat("\n=== ROBUSTEZ DEL RANKING ===\n")
cat(sprintf("  %-35s %4s  %s\n", "Candidato", "N/6", "Configs en top20"))
cat(paste(rep("-", 80), collapse=""), "\n")

stable_candidates <- rank_wide_filtered %>%
  filter(n_configs_top20 == 6)
robust_candidates <- rank_wide_filtered %>%
  filter(n_configs_top20 >= 5)

for (i in seq_len(min(25, nrow(rank_wide_filtered)))) {
  r <- rank_wide_filtered[i, ]
  configs_present <- names(r)[!is.na(r) & names(r) != "drug_name_norm" &
                                names(r) != "n_configs_top20"]
  cat(sprintf("  %-35s %4d  %s\n",
              substr(r$drug_name_norm, 1, 35),
              r$n_configs_top20,
              paste(configs_present, collapse=", ")))
}

cat(sprintf("\nCandidatos en TODAS las 6 configs (altamente estables): %d\n",
            nrow(stable_candidates)))
if (nrow(stable_candidates) > 0) {
  cat(paste("  ", stable_candidates$drug_name_norm, collapse="\n"), "\n")
}
cat(sprintf("Candidatos en >= 5 configs (robustos): %d\n", nrow(robust_candidates)))

# =============================================================================
# FIGURAS
# =============================================================================
cat("\n--- Generando figuras ---\n")

safe_pdf <- function(file, expr, w = 12, h = 8) {
  tryCatch({
    pdf(file, width = w, height = h)
    force(expr)
    dev.off()
    cat(sprintf("  %s\n", file))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("  WARN %s: %s\n", file, e$message))
  })
}

config_labels <- c(
  baseline        = "Baseline\n(actual)",
  clinical_heavy  = "Clinical\nheavy",
  molecular_heavy = "Molecular\nheavy",
  network_heavy   = "Network\nheavy",
  pathway_heavy   = "Pathway\nheavy",
  equal_weights   = "Equal\nweights"
)

# 1. Heatmap de ranks: candidato (ordenado por estabilidad) × config
safe_pdf("results/figures/15_rank_heatmap.pdf", w = 13, h = 10, {
  # Usar solo candidatos en >= 1 config top 20
  plot_data <- rank_long %>%
    # Factor de candidato ordenado por n_configs_top20 + rank baseline
    left_join(rank_wide_filtered %>% select(drug_name_norm, n_configs_top20, baseline_rank = baseline),
              by = "drug_name_norm") %>%
    mutate(
      drug_label = str_to_title(str_to_lower(drug_name_norm)),
      drug_label = factor(drug_label,
                          levels = rev(str_to_title(str_to_lower(
                            rank_wide_filtered$drug_name_norm)))),
      config_label = factor(config_labels[config],
                            levels = config_labels),
      # Categoria de rank para color
      rank_cat = case_when(
        rank <= 5  ~ "Top 5",
        rank <= 10 ~ "Top 6-10",
        rank <= 20 ~ "Top 11-20",
        TRUE       ~ "Fuera top 20"
      ),
      rank_cat = factor(rank_cat, levels = c("Top 5", "Top 6-10", "Top 11-20", "Fuera top 20"))
    )

  p <- ggplot(plot_data, aes(x = config_label, y = drug_label, fill = rank_cat)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = rank), size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c(
      "Top 5"        = "#2ca02c",
      "Top 6-10"     = "#98df8a",
      "Top 11-20"    = "#d4e6c3",
      "Fuera top 20" = "white"
    ), na.value = "white", name = "Posicion") +
    labs(
      title   = "Estabilidad del ranking — Analisis de sensibilidad de pesos",
      subtitle = sprintf("Cada celda = rank en esa configuracion (verde = top 5) | %d candidatos en >= 1 config top 20",
                         length(candidates_in_any)),
      x = "Configuracion de pesos",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 9),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  print(p)
})

# 2. Barplot de estabilidad: cuántas configs aparece cada candidato
safe_pdf("results/figures/15_stability_bar.pdf", w = 12, h = 8, {
  plot_data <- rank_wide_filtered %>%
    mutate(drug_label = str_to_title(str_to_lower(drug_name_norm)),
           drug_label = factor(drug_label,
                               levels = drug_label[order(n_configs_top20,
                                                          -coalesce(baseline, 21L))]),
           stability = case_when(
             n_configs_top20 == 6 ~ "6/6 (muy estable)",
             n_configs_top20 >= 5 ~ "5/6 (robusto)",
             n_configs_top20 >= 4 ~ "4/6 (moderado)",
             TRUE                 ~ "<= 3/6 (variable)"
           ),
           stability = factor(stability, levels = c("6/6 (muy estable)", "5/6 (robusto)",
                                                     "4/6 (moderado)", "<= 3/6 (variable)")))

  p <- ggplot(plot_data, aes(x = drug_label, y = n_configs_top20, fill = stability)) +
    geom_col(width = 0.75, alpha = 0.9) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "gray40", linewidth = 0.7) +
    annotate("text", x = nrow(plot_data) * 0.9, y = 5.2,
             label = "umbral robusto (5/6)", color = "gray40", size = 3.5) +
    scale_fill_manual(values = c(
      "6/6 (muy estable)" = "#006d2c",
      "5/6 (robusto)"     = "#31a354",
      "4/6 (moderado)"    = "#fd8d3c",
      "<= 3/6 (variable)" = "#bdbdbd"
    ), name = "Estabilidad") +
    scale_y_continuous(breaks = 0:6, limits = c(0, 6.5)) +
    coord_flip() +
    labs(
      title    = "Estabilidad de candidatos a traves de 6 configuraciones de pesos",
      subtitle = "N configuraciones en las que el candidato aparece en el Top 20",
      x = NULL,
      y = "N configuraciones en top 20 (max = 6)"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "bottom")
  print(p)
})

# 3. Boxplot de scores: variacion de final_score por candidato a traves de configs
safe_pdf("results/figures/15_score_distribution.pdf", w = 12, h = 7, {
  # Solo candidatos robustos (>= 4 configs)
  robust_names <- rank_wide_filtered %>%
    filter(n_configs_top20 >= 4) %>%
    pull(drug_name_norm)

  # Calcular score en cada config para cada candidato robusto
  score_long <- lapply(names(weight_configs), function(config_name) {
    w <- weight_configs[[config_name]]
    all_scored %>%
      filter(drug_name_norm %in% robust_names) %>%
      mutate(
        final_score_cfg = pmin(
          w["logfc"]    * s_logfc   +
          w["sig"]      * s_sig     +
          w["clinical"] * s_clinical +
          w["cmap"]     * s_cmap    +
          w["pathway"]  * s_pathway +
          w["network"]  * s_network + bonus, 1.0),
        config = config_name
      ) %>%
      select(drug_name_norm, config, final_score_cfg)
  }) %>% bind_rows()

  # Ordenar por score mediano
  med_order <- score_long %>%
    group_by(drug_name_norm) %>%
    summarise(med = median(final_score_cfg)) %>%
    arrange(desc(med))

  score_long <- score_long %>%
    mutate(drug_label = str_to_title(str_to_lower(drug_name_norm)),
           drug_label = factor(drug_label,
                               levels = rev(str_to_title(str_to_lower(med_order$drug_name_norm)))))

  p <- ggplot(score_long, aes(x = drug_label, y = final_score_cfg)) +
    geom_boxplot(fill = "#a6cee3", alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = config), width = 0.15, size = 2.5, alpha = 0.8) +
    scale_color_brewer(palette = "Dark2", name = "Config") +
    coord_flip() +
    labs(
      title    = "Distribucion de scores a traves de configuraciones (candidatos robustos)",
      subtitle = "Cada punto = score en una configuracion | Boxplot = rango inter-configuraciones",
      x = NULL,
      y = "Final score (0-1)"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "bottom")
  print(p)
})

# =============================================================================
# EXPORTAR TABLA
# =============================================================================
cat("\n--- Exportando ---\n")

# =============================================================================
# ANÁLISIS 2: Drop-one-database (leave-one-out por fuente)
# =============================================================================
cat("\n=== ANÁLISIS 2: Drop-one-database ===\n")

min_db         <- params$candidates$min_databases   # 2 por defecto según config
sources_to_drop <- c("DGIdb", "ChEMBL", "OpenTargets", "L2S2")

lod_results <- lapply(sources_to_drop, function(drop_src) {
  cat(sprintf("  Excluyendo fuente: %s\n", drop_src))

  # Eliminar la fuente del string de fuentes y recontarlas
  cands_filtered <- all_scored %>%
    mutate(
      sources_adj   = str_remove_all(sources, drop_src),
      # Limpiar pipes dobles o iniciales/finales que puedan quedar tras la eliminación
      sources_adj   = str_replace_all(sources_adj, "\\|{2,}", "|"),
      sources_adj   = str_remove_all(sources_adj, "^\\||\\|$"),
      # Contar fuentes restantes: si vacío = 0, si no vacío contar separadores + 1
      n_sources_adj = ifelse(sources_adj == "", 0L,
                             str_count(sources_adj, "\\|") + 1L)
    ) %>%
    filter(n_sources_adj >= min_db | drug_class %in% c("A", "B"))

  if (nrow(cands_filtered) == 0) return(NULL)

  scored_lod <- all_scored %>%
    filter(drug_name_norm %in% cands_filtered$drug_name_norm) %>%
    arrange(desc(final_score)) %>%
    slice_head(n = top_n)

  data.frame(
    drug_name_norm = scored_lod$drug_name_norm,
    rank_lod       = seq_len(nrow(scored_lod)),
    drop_source    = drop_src,
    stringsAsFactors = FALSE
  )
})

df_lod <- bind_rows(lod_results)

lod_stability <- df_lod %>%
  group_by(drug_name_norm) %>%
  summarise(n_lod_topN = n(), .groups = "drop") %>%
  mutate(lod_stable = n_lod_topN == length(sources_to_drop))

cat(sprintf("Drugs estables en top %d con cualquier fuente excluida: %d\n",
            top_n, sum(lod_stability$lod_stable)))

# =============================================================================
# ANÁLISIS 3: Permutation test (distribución nula del composite score)
# =============================================================================
cat("\n=== ANÁLISIS 3: Permutation test (n=1000) ===\n")

# Cargar tabla DE (proteómica) para obtener pi_stat por proteína
de_file <- "results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv"
if (!file.exists(de_file)) {
  cat("WARN: archivo DE no encontrado:", de_file, "— omitiendo permutation test\n")
} else {
  gene_de <- read.delim(de_file, stringsAsFactors = FALSE)

  # Calcular pi_stat si no existe
  if (!"pi_stat" %in% colnames(gene_de)) {
    gene_de <- gene_de %>%
      mutate(pi_stat = sign(logFC_TVsS) * abs(logFC_TVsS) *
               (-log10(adj.P.Val_TVsS + 1e-300)))
  }

  max_pi_real <- max(abs(gene_de$pi_stat), na.rm = TRUE)

  # Pesos del scoring — usar config baseline del script 15
  w_base <- weight_configs[["baseline"]]
  # Baseline: logfc=0.20, sig=0.15, clinical=0.20, cmap=0.15, pathway=0.15, network=0.15
  # Para el permutation test colapsamos logfc+sig en la componente molecular (pi_stat)
  w_molec <- w_base["logfc"] + w_base["sig"]   # 0.35 combinado
  w_clin  <- w_base["clinical"]                # 0.20
  w_cmap  <- w_base["cmap"]                   # 0.15
  w_path  <- w_base["pathway"]                # 0.15
  w_net   <- w_base["network"]               # 0.15
  # Renormalizar para que sumen 1 (ya suman 1 por construcción)

  set.seed(2026)
  n_perm <- 1000

  perm_score_fn <- function() {
    gene_de_perm <- gene_de %>%
      mutate(pi_stat_perm = sample(pi_stat))
    max_pi_perm <- max(abs(gene_de_perm$pi_stat_perm), na.rm = TRUE)
    if (max_pi_perm == 0) max_pi_perm <- 1  # evitar division por cero

    # Calcular s_pi_perm por candidato a partir de sus genes diana
    s_pi_vec <- sapply(all_scored$de_genes, function(g) {
      genes   <- if (is.na(g) || g == "") character(0) else
                 strsplit(g, "\\|")[[1]]
      pi_vals <- abs(gene_de_perm$pi_stat_perm[gene_de_perm$symbol_org %in% genes])
      if (length(pi_vals) == 0) 0 else mean(pi_vals, na.rm = TRUE) / max_pi_perm
    })

    # Composite score permutado (estructura análoga a baseline)
    s_comp_perm <- w_molec * s_pi_vec +
                   w_clin  * all_scored$s_clinical +
                   w_cmap  * all_scored$s_cmap +
                   w_path  * all_scored$s_pathway +
                   w_net   * all_scored$s_network

    max(s_comp_perm, na.rm = TRUE)
  }

  perm_null       <- replicate(n_perm, perm_score_fn())
  true_top_score  <- max(all_scored$composite_score, na.rm = TRUE)
  perm_pval       <- mean(perm_null >= true_top_score)

  cat(sprintf("Score real top-1 (composite_score): %.4f\n", true_top_score))
  cat(sprintf("Permutation p-value (top score):     %.4f\n", perm_pval))
  cat(sprintf("Distribucion nula: media=%.4f, SD=%.4f, p95=%.4f\n",
              mean(perm_null), sd(perm_null), quantile(perm_null, 0.95)))

  # Exportar resultados adicionales
  write.table(lod_stability,
              "results/tables/15_lod_stability.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)

  write.table(data.frame(
    analysis       = "permutation_null",
    n_permutations = n_perm,
    true_top_score = true_top_score,
    perm_pval      = perm_pval,
    perm_mean      = mean(perm_null),
    perm_sd        = sd(perm_null),
    perm_95pct     = as.numeric(quantile(perm_null, 0.95))
  ), "results/tables/15_permutation_test.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Exportados: 15_lod_stability.tsv y 15_permutation_test.tsv\n")
}

write.table(rank_wide_filtered,
            "results/tables/15_sensitivity_ranks.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Exportado: 15_sensitivity_ranks.tsv (%d candidatos)\n",
            nrow(rank_wide_filtered)))

# =============================================================================
# RESUMEN FINAL
# =============================================================================
cat("\n=== RESUMEN: ROBUSTEZ DEL SCORING ===\n")
cat(sprintf("  Configuraciones evaluadas: %d\n", length(weight_configs)))
cat(sprintf("  Candidatos en top 20 en alguna config: %d\n", length(candidates_in_any)))
cat(sprintf("\n  Candidatos altamente estables (6/6 configs):\n"))
for (nm in stable_candidates$drug_name_norm) {
  r <- rank_wide_filtered %>% filter(drug_name_norm == nm)
  cat(sprintf("    %s (rank baseline=%s)\n", nm,
              ifelse(is.na(r$baseline), "fuera", r$baseline)))
}
cat(sprintf("\n  Candidatos robustos (5/6 configs):\n"))
for (nm in setdiff(robust_candidates$drug_name_norm, stable_candidates$drug_name_norm)) {
  cat(sprintf("    %s\n", nm))
}

cat(sprintf("\nCONCLUSION: El Top 5 del ranking final "))
top5_baseline <- rank_wide_filtered %>%
  filter(!is.na(baseline)) %>%
  arrange(baseline) %>%
  slice_head(n = 5)
stability_top5 <- mean(top5_baseline$n_configs_top20, na.rm = TRUE)
if (stability_top5 >= 5.5) {
  cat("es ALTAMENTE ROBUSTO — aparece en todas o casi todas las configuraciones.\n")
} else if (stability_top5 >= 4.5) {
  cat("es ROBUSTO — se mantiene en la mayoria de configuraciones.\n")
} else {
  cat("muestra variabilidad moderada — considerar analisis adicional.\n")
}

cat("\nFin:", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con_log)
cat("Script completado. Log en:", log_file, "\n")
