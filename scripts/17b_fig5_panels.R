#!/usr/bin/env Rscript
# =============================================================================
# 17b_fig5_panels.R — Paneles de la Figura 5 (priorización hub-céntrica)
# =============================================================================
# Figura-clímax: aterriza el resultado del método. Organizada por la estructura
# de la red (módulo → hub drogable → fármaco), con el composite de dos factores
# (TargetPriority × DrugViability) del script 10 v3.
#
# Narrativa de dos niveles:
#   - hub_central     : fármaco anclado a un HUB de red (EGFR=validación;
#                       OXPHOS/metformina, proteasoma/carfilzomib = novedad)
#   - peripheral_diff : diana en red diferencialmente abundante pero no-hub
#                       (DNMT1/epigenético, antimetabolitos, repurposing clásico)
#
# Paneles (cacheados como .rds para 17g_fig5_multipanel.R):
#   Fig5A_shortlist   — shortlist rankeado; barras apiladas = descomposición del
#                       composite en aporte de diana (0.6·TP) vs droga (0.4·DV);
#                       faceteado por tier. Muestra QUÉ rankea y POR QUÉ.
#   Fig5B_tpdv_space  — espacio bifactor TargetPriority × DrugViability de todos
#                       los candidatos; shortlist etiquetado. Muestra candidatos
#                       anclados-en-diana vs anclados-en-droga.
#
# Input:  results/tables/10_all_candidates_scored.tsv
#         config/analysis_params.yaml  (pesos w_target / w_drug)
# Ambiente: omics-R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(ggrepel); library(readr); library(yaml)
})

cat("=== 17b_fig5_panels.R ===\n")
source(here::here("scripts", "_setup.R")); setup_project()
source(here::here("scripts", "_fig_style.R"))

N_SHORTLIST <- 14   # nº de hubs-ancla a mostrar (top por composite, tier != off_network)

# ── Pesos del composite (para descomponer en aporte diana vs droga) ───────────
params   <- yaml::read_yaml("config/analysis_params.yaml")
w_target <- params$scoring_v3$weight_target
w_drug   <- params$scoring_v3$weight_drug

# ── Datos ─────────────────────────────────────────────────────────────────────
scored <- read_tsv("results/tables/10_all_candidates_scored.tsv", show_col_types = FALSE)

# Etiqueta de módulo legible (quita prefijo M#_, normaliza)
mod_label <- function(x) {
  x <- str_replace(x, "^M[0-9]+_", "")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace(x, "_?$", "")
  str_to_sentence(str_trim(x))
}
TIER_LAB <- c(hub_central     = "Network-central hub",
              peripheral_diff = "Peripheral, differentially abundant")
TIER_COLS <- c("Network-central hub" = "#0072B2",
               "Peripheral, differentially abundant" = "#E69F00")
COMP_COLS <- c("Target priority"  = "#009E73",
               "Drug viability"   = "#CC79A7")

cand <- scored %>%
  filter(tier != "off_network") %>%
  mutate(tier_lab = factor(TIER_LAB[tier], levels = TIER_LAB),
         module_lab = mod_label(module_name),
         drug_title = str_to_title(str_to_lower(drug_name_norm)))

# Shortlist: top fármaco por hub-ancla, luego top N por composite ──────────────
short <- cand %>%
  group_by(primary_target) %>%
  slice_max(composite_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  slice_max(composite_score, n = N_SHORTLIST) %>%
  mutate(label = sprintf("%s  (%s)", drug_title, primary_target))

cat(sprintf("Shortlist: %d hubs-ancla | hub_central=%d peripheral=%d\n",
            nrow(short), sum(short$tier == "hub_central"),
            sum(short$tier == "peripheral_diff")))

# =============================================================================
# PANEL A — Shortlist rankeado + descomposición del composite (TP vs DV)
# =============================================================================
short_long <- short %>%
  mutate(`Target priority` = w_target * TargetPriority,
         `Drug viability`   = w_drug   * DrugViability) %>%
  select(label, tier_lab, module_lab, composite_score,
         `Target priority`, `Drug viability`) %>%
  pivot_longer(c(`Target priority`, `Drug viability`),
               names_to = "component", values_to = "contrib") %>%
  mutate(component = factor(component, levels = names(COMP_COLS)),
         label = reorder(label, composite_score))

p_short <- ggplot(short_long, aes(contrib, label, fill = component)) +
  geom_col(width = 0.72) +
  geom_text(data = short %>% mutate(label = reorder(label, composite_score)),
            aes(x = composite_score, y = label, label = sprintf("%.2f", composite_score)),
            inherit.aes = FALSE, hjust = -0.18, size = 2.6, color = "grey20") +
  facet_grid(tier_lab ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = COMP_COLS, name = "Composite component") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = NULL, x = "Composite score", y = NULL) +
  theme_pub() +
  theme(legend.position = "top",
        panel.grid.major.x = element_line(linewidth = 0.25, color = "grey88"),
        strip.text.y = element_text(angle = 90, size = 9))

# =============================================================================
# PANEL B — Espacio bifactor TargetPriority × DrugViability
# =============================================================================
p_space <- ggplot(cand, aes(TargetPriority, DrugViability)) +
  # iso-líneas de composite (referencia)
  geom_abline(slope = -w_target / w_drug,
              intercept = c(0.45, 0.55, 0.65) / w_drug,
              linetype = "dotted", color = "grey80", linewidth = 0.3) +
  geom_point(aes(color = tier_lab), alpha = 0.55, size = 1.6) +
  geom_point(data = short, aes(color = tier_lab), size = 2.4) +
  geom_text_repel(data = short, aes(label = drug_title, color = tier_lab),
                  size = 2.3, max.overlaps = 20, min.segment.length = 0,
                  box.padding = 0.3, show.legend = FALSE) +
  scale_color_manual(values = TIER_COLS, name = NULL) +
  labs(title = NULL, x = "Target priority", y = "Drug viability") +
  theme_pub() +
  theme(legend.position = "top")

# ── Guardar paneles individuales + cache .rds ─────────────────────────────────
save_pub(p_short, "Fig5A_shortlist", h_add = 10)
save_panel_obj(p_short, "Fig5A_shortlist")
save_pub(p_space, "Fig5B_tpdv_space")
save_panel_obj(p_space, "Fig5B_tpdv_space")

cat("\n17b_fig5_panels — OK\n")
cat("Fin:", format(Sys.time()), "\n")
