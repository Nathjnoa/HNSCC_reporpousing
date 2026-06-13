#!/usr/bin/env Rscript
# =============================================================================
# 17h_figS_robustness.R — Figura SUPLEMENTARIA: robustez del scoring v3
# =============================================================================
# Heatmap de estabilidad del ranking a través de 6 configuraciones de pesos,
# anotado con estabilidad leave-one-database (LOD). Reemplaza los paneles base-R
# legacy del antiguo script 15. Estilo de scripts/_fig_style.R.
#
# Input:  results/tables/15_sensitivity_ranks.tsv
#         results/tables/15_lod_stability.tsv
# Output: results/figures/pub/supp/FigS_robustness.{tif,pdf,png}
# Ambiente: omics-R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(readr)
})

cat("=== 17h_figS_robustness.R ===\n")
source(here::here("scripts", "_setup.R")); setup_project()
source(here::here("scripts", "_fig_style.R"))

sens <- read_tsv("results/tables/15_sensitivity_ranks.tsv", show_col_types = FALSE)
lod  <- read_tsv("results/tables/15_lod_stability.tsv",     show_col_types = FALSE)

config_order <- c("baseline","target_heavy","drug_heavy",
                  "centrality_heavy","abundance_heavy","equal_weights")
config_lab <- c(baseline="Baseline", target_heavy="Target\nheavy",
                drug_heavy="Drug\nheavy", centrality_heavy="Centrality\nheavy",
                abundance_heavy="Abundance\nheavy", equal_weights="Equal\nweights")

# Orden de fármacos por estabilidad (n_configs) y rank baseline
ord <- sens %>% arrange(desc(n_configs_topN), baseline)
lod_lab <- lod %>%
  mutate(lod_txt = ifelse(lod_stable, "LOD-stable", sprintf("%d/4", n_lod_topN))) %>%
  select(drug_name_norm, lod_txt, lod_stable)

long <- sens %>%
  pivot_longer(all_of(config_order), names_to = "config", values_to = "rank") %>%
  left_join(lod_lab, by = "drug_name_norm") %>%
  mutate(
    drug_label = str_to_title(str_to_lower(drug_name_norm)),
    drug_label = factor(drug_label,
                        levels = rev(str_to_title(str_to_lower(ord$drug_name_norm)))),
    config = factor(config_lab[config], levels = config_lab),
    rank_cat = factor(case_when(
      is.na(rank)  ~ "Outside top",
      rank <= 5    ~ "Top 5",
      rank <= 15   ~ "Top 6-15",
      TRUE         ~ "Top 16-35"
    ), levels = c("Top 5","Top 6-15","Top 16-35","Outside top"))
  )

RANK_COLS <- c("Top 5"="#08519c","Top 6-15"="#6baed6",
               "Top 16-35"="#c6dbef","Outside top"="grey92")

p <- ggplot(long, aes(config, drug_label, fill = rank_cat)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = rank), size = 2.2, color = "grey15", na.rm = TRUE) +
  scale_fill_manual(values = RANK_COLS, name = "Rank en configuración") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

n_drug <- nlevels(long$drug_label)
save_pub(p, "FigS_robustness", out_dir = "results/figures/pub/supp",
         h_add = max(0, n_drug * 3 - 120))
save_tiff(p, "FigS_robustness", width_mm = 180, height_mm = max(120, n_drug * 5 + 40),
          out_dir = "results/figures/pub/supp")

cat(sprintf("  FigS_robustness — %d candidatos × 6 configs\n", n_drug))
cat("Fin:", format(Sys.time()), "\n")
