#!/usr/bin/env Rscript
# =============================================================================
# 18_pub_tables.R
# =============================================================================
# Genera las tablas formateadas para publicación.
# Sección 0 (base metodológica), OE1 (fuentes farmacológicas) y OE2 (priorización).
#
# Tablas generadas (results/tables/pub/):
#   main/
#     Sec0_Tab1_resumen_DE.tsv          — Resumen estadístico del análisis DE
#     Sec0_Tab2_top_proteinas.tsv       — Top proteínas DE (20 up + 20 down)
#     OE1_Tab2_top_candidatos.tsv       — Top 30 candidatos por n_fuentes
#     OE2_Tab1_EGFR_LOD_stable.tsv      — Candidatos EGFR LOD-stable
#     OE2_Tab2_noEGFR_LOD_stable.tsv    — Candidatos no-EGFR LOD-stable
#   supp/
#     OE2_TabS1_candidatos_extendidos_noEGFR.tsv — Candidatos robustos a pesos, no LOD-stable
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

# ── Directorio del proyecto ───────────────────────────────────────────────────
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  proj_dir    <- dirname(dirname(script_path))
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml")))
    setwd(proj_dir)
}

cat("=== 18_pub_tables.R ===\n")
cat("Working directory:", getwd(), "\n")
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
de_up    <- read_tsv("results/tables/de_limma/01_TVsS_upregulated.tsv",
                     show_col_types = FALSE)
de_down  <- read_tsv("results/tables/de_limma/01_TVsS_downregulated.tsv",
                     show_col_types = FALSE)
multi    <- read_tsv("results/tables/drug_targets/08_multi_source_candidates.tsv",
                     show_col_types = FALSE)

cat("  Datos cargados OK\n")

# =============================================================================
# SECCIÓN 0 — TABLAS
# =============================================================================
cat("\n--- Sección 0: Tablas proteómica diferencial ---\n")

# ── Sec0_Tab1: Resumen del análisis DE ───────────────────────────────────────
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
save_tsv(sec0_tab1, "Sec0_Tab1_resumen_DE", out_main)

# ── Sec0_Tab2: Top proteínas DE ──────────────────────────────────────────────
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
save_tsv(sec0_tab2, "Sec0_Tab2_top_proteinas", out_main)

# =============================================================================
# OE1 — TABLAS
# =============================================================================
cat("\n--- OE1: Tablas bases de datos de fármacos ---\n")

# ── OE1_Tab2: Top candidatos por número de fuentes ───────────────────────────
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
save_tsv(oe1_tab2, "OE1_Tab2_top_candidatos", out_main)

# =============================================================================
# OE2 — TABLAS: Candidatos LOD-stable
# =============================================================================
cat("\n--- OE2: Tablas candidatos LOD-stable ---\n")

scored     <- read_tsv("results/tables/10_all_candidates_scored.tsv", show_col_types = FALSE)
lod        <- read_tsv("results/tables/15_lod_stability.tsv",         show_col_types = FALSE)
sens_ranks <- read_tsv("results/tables/15_sensitivity_ranks.tsv",     show_col_types = FALSE)

lod_stable_set <- lod %>% filter(lod_stable == TRUE) %>% pull(drug_name_norm)

# Subclases EGFR para clasificación
tki_gen1_2 <- toupper(c(
  "gefitinib", "erlotinib", "afatinib", "lapatinib", "icotinib",
  "neratinib", "dacomitinib", "canertinib dihydrochloride"
))
tki_gen3 <- toupper(c(
  "osimertinib", "lazertinib", "olmutinib", "abivertinib",
  "aumolertinib", "firmonertinib", "rociletinib", "mobocertinib"
))
mab_egfr <- toupper(c(
  "cetuximab", "cetuximab sarotalocan", "necitumumab", "nimotuzumab",
  "panitumumab", "amivantamab", "depatuxizumab mafodotin"
))

oe2_base <- scored %>%
  filter(drug_name_norm %in% lod_stable_set) %>%
  mutate(
    hnscc_approved_curated = drug_name_norm %in% HNSCC_APPROVED_CURATED,
    is_egfr_primary = str_detect(primary_target, "EGFR|ERBB|HER"),
    egfr_subclass = case_when(
      drug_name_norm %in% tki_gen1_2 ~ "TKI EGFR 1ª-2ª generación",
      drug_name_norm %in% tki_gen3   ~ "TKI EGFR 3ª generación",
      drug_name_norm %in% mab_egfr   ~ "Anticuerpo/ADC anti-EGFR",
      is_egfr_primary                ~ "Inhibidor multi-quinasa (EGFR+)",
      TRUE                           ~ NA_character_
    ),
    `Fase clínica máx.` = case_when(
      max_phase == 4 ~ "Aprobado (Fase IV)",
      max_phase == 3 ~ "Fase III",
      max_phase == 2 ~ "Fase II",
      max_phase == 1 ~ "Fase I",
      TRUE           ~ "No reportada"
    )
  )

# ── OE2_Tab1: LOD-stable EGFR — validación del método ───────────────────────
oe2_tab1 <- oe2_base %>%
  filter(!is.na(egfr_subclass)) %>%
  arrange(egfr_subclass, desc(composite_score)) %>%
  mutate(
    Fármaco           = str_to_title(drug_name_norm),
    `Aprobado HNSCC`  = ifelse(hnscc_approved_curated, "Sí", "No"),
    `Composite score` = round(composite_score, 3)
  ) %>%
  select(
    Fármaco,
    `Subclase`          = egfr_subclass,
    `Target primario`   = primary_target,
    `Fase clínica máx.`,
    `Aprobado HNSCC`,
    `Composite score`,
    `N fuentes`         = n_sources
  )
save_tsv(oe2_tab1, "OE2_Tab1_EGFR_LOD_stable", out_main)

# ── OE2_Tab2: LOD-stable no-EGFR — candidatos de repurposing ─────────────────
mecanismo_map <- c(
  "DECITABINE"      = "Inhibidor DNMT",
  "AZACITIDINE"     = "Inhibidor DNMT",
  "CEDAZURIDINE"    = "Inhibidor desaminasa de citidina (potenciador DNMT)",
  "CARFILZOMIB"     = "Inhibidor de proteasoma",
  "MITAPIVAT"       = "Activador piruvato quinasa (PK-R)",
  "TRANYLCYPROMINE" = "Inhibidor LSD1/MAO (epigenético)"
)
indicacion_map <- c(
  "DECITABINE"      = "Síndrome mielodisplásico / LMA",
  "AZACITIDINE"     = "Síndrome mielodisplásico / LMA",
  "CEDAZURIDINE"    = "Síndrome mielodisplásico (combinado con decitabine)",
  "CARFILZOMIB"     = "Mieloma múltiple",
  "MITAPIVAT"       = "Anemia hemolítica (deficiencia piruvato quinasa)",
  "TRANYLCYPROMINE" = "Depresión / investigacional en oncología"
)

oe2_tab2 <- oe2_base %>%
  filter(is.na(egfr_subclass)) %>%
  arrange(desc(composite_score)) %>%
  mutate(
    Fármaco           = str_to_title(drug_name_norm),
    Mecanismo         = mecanismo_map[drug_name_norm],
    `Indicación actual` = indicacion_map[drug_name_norm],
    `Composite score` = round(composite_score, 3)
  ) %>%
  select(
    Fármaco,
    Mecanismo,
    `Target primario`   = primary_target,
    `Indicación actual`,
    `Fase clínica máx.`,
    `Composite score`,
    `N fuentes`         = n_sources
  )
save_tsv(oe2_tab2, "OE2_Tab2_noEGFR_LOD_stable", out_main)

# ── OE2_TabS1: Candidatos extendidos no-EGFR (robustos a pesos, no LOD) ──────
# Drogas que aparecen en ≥5/6 configuraciones de pesos (top-35)
# pero NO son LOD-stable (sensibles al tamaño del pool de corte)
oe2_tabs1 <- sens_ranks %>%
  filter(
    !drug_name_norm %in% lod_stable_set,
    n_configs_top20 >= 5
  ) %>%
  left_join(
    scored %>% select(drug_name_norm, primary_target, composite_score,
                      max_phase, n_sources),
    by = "drug_name_norm"
  ) %>%
  filter(!str_detect(primary_target, "EGFR|ERBB|HER")) %>%
  arrange(desc(n_configs_top20), desc(composite_score)) %>%
  mutate(
    Fármaco             = str_to_title(drug_name_norm),
    `N configs (de 6)`  = n_configs_top20,
    `Fase clínica máx.` = case_when(
      max_phase == 4 ~ "Aprobado (Fase IV)",
      max_phase == 3 ~ "Fase III",
      max_phase == 2 ~ "Fase II",
      max_phase == 1 ~ "Fase I",
      TRUE           ~ "No reportada"
    ),
    `Composite score`   = round(composite_score, 3)
  ) %>%
  select(
    Fármaco,
    `Target primario`   = primary_target,
    `N configs (de 6)`,
    `Fase clínica máx.`,
    `Composite score`,
    `N fuentes`         = n_sources
  )
save_tsv(oe2_tabs1, "OE2_TabS1_candidatos_extendidos_noEGFR", out_supp)

cat(sprintf("  OE2_Tab1: %d drugs EGFR LOD-stable\n",  nrow(oe2_tab1)))
cat(sprintf("  OE2_Tab2: %d drugs no-EGFR LOD-stable\n", nrow(oe2_tab2)))
cat(sprintf("  OE2_TabS1: %d candidatos extendidos no-EGFR\n", nrow(oe2_tabs1)))

# =============================================================================
# RESUMEN
# =============================================================================
cat("\n============================================================\n")
cat("Script 18 COMPLETO\n")
cat("Output main:  results/tables/pub/main/\n")
cat("Output supp:  results/tables/pub/supp/\n")
cat(sprintf("Fin: %s\n", format(Sys.time())))
sink()
