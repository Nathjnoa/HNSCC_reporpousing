#!/usr/bin/env Rscript
# =============================================================================
# Script 15: Robustez del scoring v3 (TargetPriority × DrugViability)
# HNSCC Drug Repurposing — Validación de robustez (→ figura SUPLEMENTARIA)
# =============================================================================
# Reescrito para el composite de dos factores del script 10 v3:
#   composite = w_target·TP + w_drug·DV
#     TP = w_cent·centrality + w_pi·pi_directional         (nivel-diana)
#     DV = w_rev·reversal + w_reg·regulatory + w_brd·breadth (nivel-droga)
#
# Dos análisis de robustez (ambos reutilizan las componentes ya exportadas por
# script 10, sin re-ejecutar el pipeline):
#   1. Sensibilidad a pesos — 6 configuraciones que varían el balance
#      target/drug y los sub-pesos; ¿se mantiene el ranking?
#   2. Leave-one-database — quita cada fuente (DGIdb/ChEMBL/OpenTargets/L2S2),
#      recalcula breadth (y reversal si cae L2S2) y re-rankea.
#
# NOTA: el test de permutación se eliminó (2026-06). Para un score de PRIORIZACIÓN
# por agregación de evidencia (no un estadístico de descubrimiento), la validez se
# sostiene en recuperación de controles positivos (eje EGFR aprobado) + robustez
# (sensibilidad a pesos + LOD) + validación externa (Fig6), no en un p de
# permutación. Permutar solo la abundancia diferencial testeaba un sub-componente
# menor del composite v3 y resultaba engañoso.
#
# Input:  results/tables/10_all_candidates_scored.tsv
#         config/analysis_params.yaml (scoring_v3)
# Output: results/tables/15_sensitivity_ranks.tsv
#         results/tables/15_lod_stability.tsv
#   (figura suplementaria: scripts/17h_figS_robustness.R)
# Ambiente: omics-R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(yaml)
})

source(here::here("scripts", "_setup.R")); setup_project()
dir.create("logs", showWarnings = FALSE)
log_file <- paste0("logs/15_sensitivity_analysis_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
con_log <- file(log_file, open = "wt")
sink(con_log, type = "output"); sink(con_log, type = "message", append = TRUE)

cat("=== 15_sensitivity_analysis.R (v3: TP × DV) ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

params <- yaml::read_yaml("config/analysis_params.yaml")
top_n  <- params$candidates$top_n
sc     <- params$scoring_v3

scored <- read.delim("results/tables/10_all_candidates_scored.tsv",
                     stringsAsFactors = FALSE)
cat(sprintf("Candidatos cargados: %d\n", nrow(scored)))

# Componente de composite a partir de las columnas crudas del script 10
composite_from <- function(df, wt, wd, wc, wp, wr, wg, wb,
                           cent = df$primary_target_centrality, pi = df$primary_pi,
                           rev = df$s_reversal, reg = df$s_regulatory, brd = df$s_breadth) {
  TP <- wc * cent + wp * pi
  DV <- wr * rev + wg * reg + wb * brd
  wt * TP + wd * DV
}

# =============================================================================
# 1. SENSIBILIDAD A PESOS — 6 configuraciones
# =============================================================================
# Cada config: target, drug (composite); cent, pi (TP); rev, reg, brd (DV).
# target+drug=1 ; cent+pi=1 ; rev+reg+brd=1.
weight_configs <- list(
  baseline         = list(t=sc$weight_target, d=sc$weight_drug,
                          c=sc$weight_centrality, p=sc$weight_pi,
                          r=sc$weight_reversal, g=sc$weight_regulatory, b=sc$weight_breadth),
  target_heavy     = list(t=0.75, d=0.25, c=0.55, p=0.45, r=0.40, g=0.40, b=0.20),
  drug_heavy       = list(t=0.40, d=0.60, c=0.55, p=0.45, r=0.40, g=0.40, b=0.20),
  centrality_heavy = list(t=0.60, d=0.40, c=0.75, p=0.25, r=0.40, g=0.40, b=0.20),
  abundance_heavy  = list(t=0.60, d=0.40, c=0.30, p=0.70, r=0.40, g=0.40, b=0.20),
  equal_weights    = list(t=0.50, d=0.50, c=0.50, p=0.50, r=1/3,  g=1/3,  b=1/3)
)

cat("\n--- 1. Sensibilidad a pesos (6 configs) ---\n")
rank_list <- list()
for (nm in names(weight_configs)) {
  w <- weight_configs[[nm]]
  cfg <- scored %>%
    mutate(score_cfg = composite_from(., w$t, w$d, w$c, w$p, w$r, w$g, w$b)) %>%
    arrange(desc(score_cfg)) %>%
    slice_head(n = top_n) %>%
    mutate(rank = row_number(), config = nm) %>%
    select(drug_name_norm, config, rank)
  rank_list[[nm]] <- cfg
  cat(sprintf("  %-16s top1=%s\n", nm, cfg$drug_name_norm[1]))
}
rank_long <- bind_rows(rank_list)
rank_wide <- rank_long %>%
  pivot_wider(id_cols = drug_name_norm, names_from = config, values_from = rank) %>%
  mutate(n_configs_topN = rowSums(!is.na(select(., -drug_name_norm)))) %>%
  arrange(desc(n_configs_topN), baseline)
cat(sprintf("Candidatos en >=1 config top %d: %d | en las 6: %d\n",
            top_n, nrow(rank_wide), sum(rank_wide$n_configs_topN == 6)))

# =============================================================================
# 2. LEAVE-ONE-DATABASE
# =============================================================================
cat("\n--- 2. Leave-one-database ---\n")
min_db <- params$candidates$min_databases
srcs   <- c("DGIdb", "ChEMBL", "OpenTargets", "L2S2")
wb <- weight_configs$baseline

lod_results <- lapply(srcs, function(drop_src) {
  df <- scored %>%
    mutate(
      sources_adj   = str_remove_all(sources, drop_src),
      sources_adj   = str_replace_all(sources_adj, "\\|{2,}", "|"),
      sources_adj   = str_remove_all(sources_adj, "^\\||\\|$"),
      n_sources_adj = ifelse(sources_adj == "", 0L, str_count(sources_adj, "\\|") + 1L),
      s_breadth_adj = pmin(n_sources_adj / 4, 1.0),
      # Si se elimina L2S2, la reversión transcriptómica deja de aportar
      s_reversal_adj = if (drop_src == "L2S2") 0 else s_reversal,
      score_adj = composite_from(., wb$t, wb$d, wb$c, wb$p, wb$r, wb$g, wb$b,
                                  rev = s_reversal_adj, brd = s_breadth_adj)
    ) %>%
    filter(n_sources_adj >= min_db | drug_class %in% c("A", "B")) %>%
    arrange(desc(score_adj)) %>%
    slice_head(n = top_n)
  data.frame(drug_name_norm = df$drug_name_norm, drop_source = drop_src)
})
df_lod <- bind_rows(lod_results)
lod_stability <- df_lod %>%
  group_by(drug_name_norm) %>%
  summarise(n_lod_topN = n(), .groups = "drop") %>%
  mutate(lod_stable = n_lod_topN == length(srcs))
cat(sprintf("Estables en top %d con cualquier fuente excluida: %d\n",
            top_n, sum(lod_stability$lod_stable)))

# =============================================================================
# EXPORTAR
# =============================================================================
write.table(rank_wide, "results/tables/15_sensitivity_ranks.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lod_stability, "results/tables/15_lod_stability.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nExportados: 15_sensitivity_ranks / 15_lod_stability\n")

cat("\nFin:", format(Sys.time()), "\n")
sink(type = "message"); sink(); close(con_log)
cat("Script completado. Log en:", log_file, "\n")
