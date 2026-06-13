#!/usr/bin/env Rscript
# =============================================================================
# 18_pub_tables.R
# =============================================================================
# Genera las tablas formateadas para publicación.
# Sección 0 (base metodológica), OE1 (fuentes farmacológicas) y OE2 (priorización).
#
# Tablas generadas (results/tables/pub/):
#   main/
#     Tab1_resumen_DE.tsv          — Resumen estadístico del análisis DE
#     Tab2_top_proteinas.tsv       — Top proteínas DE (20 up + 20 down)
#     Tab3_top_candidatos.tsv       — Top 30 candidatos por n_fuentes
#     Tab4_EGFR_LOD_stable.tsv      — Candidatos EGFR LOD-stable
#     Tab5_noEGFR_LOD_stable.tsv    — Candidatos no-EGFR LOD-stable
#   supp/
#     TabS1_candidatos_extendidos_noEGFR.tsv — Candidatos robustos a pesos, no LOD-stable
#
# Ambiente: omics-R
# Ejecución:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/18_pub_tables.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

# ── Working directory (raíz del proyecto vía scripts/_setup.R) ───────────────
cat("=== 18_pub_tables.R ===\n")
source(here::here("scripts", "_setup.R"))
setup_project()
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ── Lista curada: aprobación regulatoria real en HNSCC ───────────────────────
# Fuente: FDA/EMA. Solo incluir si hay aprobación regulatoria específica para
# HNSCC (no off-label ni solo ensayos clínicos).
HNSCC_APPROVED_CURATED <- toupper(c(
  "cetuximab",       # FDA 2006, HNSCC platino-refractario y primera línea
  "afatinib",        # Aprobado HNSCC (ErbB family blocker)
  "pembrolizumab",   # FDA 2019, HNSCC primera línea y R/M
  "nivolumab",       # FDA 2016, HNSCC platino-refractario
  "capecitabine"     # Aprobado HNSCC en algunas indicaciones
))

# ── Directorios de salida ─────────────────────────────────────────────────────
out_main <- "results/tables/pub/main"
out_supp <- "results/tables/pub/supp"
dir.create(out_main, recursive = TRUE, showWarnings = FALSE)
dir.create(out_supp, recursive = TRUE, showWarnings = FALSE)

# ── Log ───────────────────────────────────────────────────────────────────────
log_dir  <- "logs"
log_file <- file.path(log_dir, paste0("18_pub_tables_", timestamp, ".log"))
sink(log_file, split = TRUE)
cat("Inicio:", format(Sys.time()), "\n\n")

# ── Función de guardado ───────────────────────────────────────────────────────
save_tsv <- function(df, name, dir) {
  path <- file.path(dir, paste0(name, ".tsv"))
  write_tsv(df, path)
  cat(sprintf("  Saved: %s (%d rows x %d cols)\n", basename(path), nrow(df), ncol(df)))
  invisible(path)
}

# =============================================================================
# CARGAR DATOS
# =============================================================================
cat("\n-- Cargando datos...\n")

de_all   <- read_tsv("results/tables/de_limma/01_TVsS_all_proteins.tsv",
                     show_col_types = FALSE)
de_sig   <- read_tsv("results/tables/de_limma/01_TVsS_significant.tsv",
                     show_col_types = FALSE)
# up/down se derivan de de_all (los splits separados son intermedios redundantes,
# eliminados en la limpieza); filtro idéntico al que aplicaba el script 01.
de_up    <- de_all %>% filter(sig.FDR_TVsS ==  1)
de_down  <- de_all %>% filter(sig.FDR_TVsS == -1)
multi    <- read_tsv("results/tables/drug_targets/08_multi_source_candidates.tsv",
                     show_col_types = FALSE)

cat("  Datos cargados OK\n")

# =============================================================================
# SECCIÓN 0 — TABLAS
# =============================================================================
cat("\n--- Sección 0: Tablas proteómica diferencial ---\n")

# ── Tab1: Resumen del análisis DE ───────────────────────────────────────
sec0_tab1 <- tibble(
  Parámetro   = c(
    "Proteínas cuantificadas (total)",
    "Proteínas mapeadas a Entrez ID",
    "Umbral log2 Fold Change",
    "Umbral FDR (adj. p-valor)",
    "Proteínas DE significativas",
    "  — Sobreexpresadas en tumor (up)",
    "  — Subexpresadas en tumor (down)",
    "HPV+ pacientes / pares",
    "HPV− pacientes / pares"
  ),
  Valor = c(
    nrow(de_all),
    sum(!is.na(de_all$gene_symbol)),
    "> 1.0  (cambio > 2× en escala lineal)",
    "< 0.05",
    nrow(de_sig),
    nrow(de_up),
    nrow(de_down),
    "6 pacientes / 12 muestras",
    "4 pacientes / 8 muestras"
  ),
  Nota = c(
    "Proteómica DIA, MaxQuant, 10 pares tumor/normal",
    "Mapeo con org.Hs.eg.db (Bioconductor)",
    "Aplicado sobre contrastes limma TVsS",
    "Benjamini-Hochberg (limma/proteoDA)",
    paste0(round(nrow(de_sig)/nrow(de_all)*100, 1), "% del proteoma cuantificado"),
    paste0(round(nrow(de_up)/nrow(de_sig)*100, 1), "% del total DE"),
    paste0(round(nrow(de_down)/nrow(de_sig)*100, 1), "% del total DE"),
    "M3, M5, M6, M7, M9, M10",
    "M1, M2, M4, M8"
  )
)
save_tsv(sec0_tab1, "Tab1_resumen_DE", out_main)

# ── Tab2: Top proteínas DE ──────────────────────────────────────────────
top_up_tab <- de_up %>%
  arrange(desc(logFC_TVsS)) %>%
  slice_head(n = 20) %>%
  mutate(Dirección = "Sobreexpresada")

top_down_tab <- de_down %>%
  arrange(logFC_TVsS) %>%
  slice_head(n = 20) %>%
  mutate(Dirección = "Subexpresada")

sec0_tab2 <- bind_rows(top_up_tab, top_down_tab) %>%
  select(
    `Gen`           = gene_symbol,
    `UniProt ID`    = uniprot_id,
    `log2FC`        = logFC_TVsS,
    `FDR (adj.p)`   = adj.P.Val_TVsS,
    `Dirección`     = Dirección
  ) %>%
  mutate(
    `log2FC`      = round(`log2FC`, 3),
    `FDR (adj.p)` = signif(`FDR (adj.p)`, 3)
  )
save_tsv(sec0_tab2, "Tab2_top_proteinas", out_main)

# =============================================================================
# OE1 — TABLAS
# =============================================================================
cat("\n--- OE1: Tablas bases de datos de fármacos ---\n")

# ── Tab3: Top candidatos por número de fuentes ───────────────────────────
multi_cols_querer <- c("drug_name_norm", "n_sources", "sources", "max_phase",
                        "is_approved", "hnscc_indication", "has_cancer_indication",
                        "chembl_id")
multi_cols_exist  <- intersect(multi_cols_querer, colnames(multi))

oe1_tab2 <- multi %>%
  select(all_of(multi_cols_exist)) %>%
  arrange(desc(n_sources), desc(is_approved)) %>%
  slice_head(n = 30) %>%
  mutate(
    `Fármaco`             = str_to_title(drug_name_norm),
    `N fuentes`           = n_sources,
    `Fuentes`             = sources,
    `Fase clínica máx.`  = case_when(
      max_phase == 4 ~ "Aprobado (Fase IV)",
      max_phase == 3 ~ "Fase III",
      max_phase == 2 ~ "Fase II",
      max_phase == 1 ~ "Fase I",
      TRUE           ~ "No reportada"
    ),
    `Indicación HNSCC`     = ifelse(drug_name_norm %in% HNSCC_APPROVED_CURATED, "Sí", "No"),
    `Indicación oncológica` = ifelse(has_cancer_indication, "Sí", "No")
  ) %>%
  select(`Fármaco`, `N fuentes`, `Fuentes`, `Fase clínica máx.`,
         `Indicación HNSCC`, `Indicación oncológica`,
         `ChEMBL ID` = any_of("chembl_id"))
save_tsv(oe1_tab2, "Tab3_top_candidatos", out_main)

# =============================================================================
# OE2 — TABLAS: Priorización hub-céntrica (scoring v3, dos niveles)
# =============================================================================
# Tab4 = eje EGFR (validación del método); Tab5 = candidatos novedosos
# priorizados por módulo (top fármaco por hub-ancla); TabS1 = lista extendida.
# Definidas por TIER/MÓDULO (no por regex EGFR ni por lod_stable, que ahora está
# dominado por el clúster de ~30 fármacos anti-EGFR). lod_stable y robustez a
# pesos se reportan como COLUMNAS de anotación.
cat("\n--- OE2: Tablas priorización hub-céntrica ---\n")

N_TAB5 <- 12  # nº de hubs-ancla no-EGFR en la tabla principal (top por composite)

scored     <- read_tsv("results/tables/10_all_candidates_scored.tsv", show_col_types = FALSE)
lod        <- read_tsv("results/tables/15_lod_stability.tsv",         show_col_types = FALSE)
sens_ranks <- read_tsv("results/tables/15_sensitivity_ranks.tsv",     show_col_types = FALSE)

# Anotaciones de robustez por fármaco
rob <- lod %>%
  select(drug_name_norm, lod_stable) %>%
  full_join(sens_ranks %>% select(drug_name_norm, n_configs_topN), by = "drug_name_norm")

mod_label <- function(x) {
  x <- str_replace(x, "^M[0-9]+_", ""); x <- str_replace_all(x, "_", " ")
  str_to_sentence(str_trim(x))
}
phase_lab <- function(p) case_when(
  p == 4 ~ "Aprobado (Fase IV)", p == 3 ~ "Fase III", p == 2 ~ "Fase II",
  p == 1 ~ "Fase I", TRUE ~ "No reportada")

base <- scored %>%
  filter(tier != "off_network") %>%
  left_join(rob, by = "drug_name_norm") %>%
  mutate(
    Fármaco           = str_to_title(str_to_lower(drug_name_norm)),
    `Módulo funcional` = mod_label(module_name),
    `Fase clínica máx.` = phase_lab(max_phase),
    `LOD-stable`      = ifelse(coalesce(lod_stable, FALSE), "Sí", "No"),
    `Robustez pesos`  = sprintf("%d/6", coalesce(n_configs_topN, 0L)),
    `Composite score` = round(composite_score, 3),
    TP                = round(TargetPriority, 3),
    DV                = round(DrugViability, 3)
  )

# Subclases EGFR (para Tab4)
tki_gen1_2 <- toupper(c("gefitinib","erlotinib","afatinib","lapatinib","icotinib",
                        "neratinib","dacomitinib","canertinib dihydrochloride"))
tki_gen3   <- toupper(c("osimertinib","lazertinib","olmutinib","abivertinib",
                        "aumolertinib","firmonertinib","rociletinib","mobocertinib"))
mab_egfr   <- toupper(c("cetuximab","cetuximab sarotalocan","necitumumab","nimotuzumab",
                        "panitumumab","amivantamab","depatuxizumab mafodotin"))

# ── Tab4: eje EGFR — validación del método ────────────────────────────────────
oe2_tab1 <- base %>%
  filter(primary_target == "EGFR") %>%
  mutate(
    hnscc_approved_curated = drug_name_norm %in% HNSCC_APPROVED_CURATED,
    `Subclase` = case_when(
      drug_name_norm %in% tki_gen1_2 ~ "TKI EGFR 1ª-2ª generación",
      drug_name_norm %in% tki_gen3   ~ "TKI EGFR 3ª generación",
      drug_name_norm %in% mab_egfr   ~ "Anticuerpo/ADC anti-EGFR",
      TRUE                           ~ "Inhibidor multi-quinasa (EGFR+)"
    ),
    `Aprobado HNSCC` = ifelse(hnscc_approved_curated, "Sí", "No")
  ) %>%
  arrange(`Subclase`, desc(composite_score)) %>%
  select(Fármaco, `Subclase`, `Target primario` = primary_target,
         `Fase clínica máx.`, `Aprobado HNSCC`, `Composite score`, TP, DV,
         `N fuentes` = n_sources, `LOD-stable`, `Robustez pesos`)
save_tsv(oe2_tab1, "Tab4_EGFR_validation", out_main)

# ── Tab5: candidatos novedosos priorizados por módulo (top fármaco por hub) ────
# Un representante (mejor composite) por hub-ancla no-EGFR; top N por composite.
tier_lab <- c(hub_central = "Hub central de red",
              peripheral_diff = "Periférico diferencial")
hub_top <- base %>%
  filter(primary_target != "EGFR") %>%
  group_by(primary_target) %>%
  slice_max(composite_score, n = 1, with_ties = FALSE) %>%
  ungroup()

oe2_tab2 <- hub_top %>%
  slice_max(composite_score, n = N_TAB5) %>%
  mutate(Tier = tier_lab[tier]) %>%
  arrange(factor(tier, levels = c("hub_central","peripheral_diff")), desc(composite_score)) %>%
  select(Fármaco, Tier, `Módulo funcional`, `Hub / diana` = primary_target,
         `Fase clínica máx.`, `Composite score`, TP, DV,
         `N fuentes` = n_sources, `LOD-stable`, `Robustez pesos`)
save_tsv(oe2_tab2, "Tab5_novel_candidates_by_module", out_main)

# ── TabS1: lista extendida (todos los hubs-ancla no-EGFR) ──────────────────────
oe2_tabs1 <- hub_top %>%
  mutate(Tier = tier_lab[tier]) %>%
  arrange(factor(tier, levels = c("hub_central","peripheral_diff")), desc(composite_score)) %>%
  select(Fármaco, Tier, `Módulo funcional`, `Hub / diana` = primary_target,
         `Fase clínica máx.`, `Composite score`, TP, DV,
         `N fuentes` = n_sources, `LOD-stable`, `Robustez pesos`)
save_tsv(oe2_tabs1, "TabS1_extended_candidates_by_module", out_supp)

cat(sprintf("  Tab4: %d fármacos eje EGFR (validación)\n", nrow(oe2_tab1)))
cat(sprintf("  Tab5: %d hubs-ancla novedosos (top por módulo)\n", nrow(oe2_tab2)))
cat(sprintf("  TabS1: %d hubs-ancla no-EGFR (lista extendida)\n", nrow(oe2_tabs1)))

# =============================================================================
# RESUMEN
# =============================================================================
cat("\n============================================================\n")
cat("Script 18 COMPLETO\n")
cat("Output main:  results/tables/pub/main/\n")
cat("Output supp:  results/tables/pub/supp/\n")
cat(sprintf("Fin: %s\n", format(Sys.time())))
sink()
