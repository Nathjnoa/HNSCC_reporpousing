#!/usr/bin/env Rscript
# =============================================================================
# 18_pub_tables.R
# =============================================================================
# Genera las tablas formateadas para tesis/artículo.
# Sección 0 (base metodológica) y OE1 (relacionar DE con fármacos).
#
# Tablas generadas (results/tables/pub/):
#   main/
#     Sec0_Tab1_resumen_DE.tsv          — Resumen estadístico del análisis DE
#     Sec0_Tab2_top_proteinas.tsv       — Top proteínas DE (20 up + 20 down)
#     OE1_Tab1_resumen_bases_datos.tsv  — Resumen de consulta por base de datos
#     OE1_Tab2_top_candidatos.tsv       — Top candidatos por n_fuentes
#   supp/
#     Sec0_TabS1_proteinas_DE_completo.tsv   — 666 proteínas DE completas
#     Sec0_TabS2_hallmarks_gsea.tsv          — GSEA Hallmarks completa
#     OE1_TabS1_dgidb_pares.tsv             — Pares gen-fármaco DGIdb
#     OE1_TabS2_chembl_drugs.tsv            — Fármacos ChEMBL fase>=3
#     OE1_TabS3_opentargets_drugs.tsv       — Fármacos OpenTargets
#
# También genera: pub_tables.xlsx con todas las tablas en hojas separadas.
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
  library(openxlsx)
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
hallmarks <- read_tsv("results/tables/pathway_enrichment/03_Hallmarks_GSEA.tsv",
                      show_col_types = FALSE)
dgidb    <- read_tsv("results/tables/drug_targets/04_dgidb_raw.tsv",
                     show_col_types = FALSE)
chembl   <- read_tsv("results/tables/drug_targets/05_chembl_drugs.tsv",
                     show_col_types = FALSE)
ot_drugs <- read_tsv("results/tables/drug_targets/06_opentargets_gene_drugs.tsv",
                     show_col_types = FALSE)
l2s2     <- read_tsv("results/tables/drug_targets/07_l2s2_results.tsv",
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

# ── Sec0_TabS1: Tabla completa 666 proteínas DE ──────────────────────────────
sec0_tabs1 <- de_sig %>%
  mutate(
    Dirección = ifelse(logFC_TVsS > 0, "Up", "Down"),
    `log2FC`  = round(logFC_TVsS, 3),
    `FDR`     = signif(adj.P.Val_TVsS, 3)
  ) %>%
  select(
    Gen          = gene_symbol,
    `UniProt ID` = uniprot_id,
    `log2FC`,
    `FDR`,
    Dirección
  ) %>%
  arrange(desc(`log2FC`))
save_tsv(sec0_tabs1, "Sec0_TabS1_proteinas_DE_completo", out_supp)

# ── Sec0_TabS2: GSEA Hallmarks ───────────────────────────────────────────────
# Seleccionar columnas relevantes si existen
hall_cols_querer <- c("Description", "NES", "pvalue", "p.adjust",
                       "setSize", "enrichmentScore", "core_enrichment")
hall_cols_exist  <- intersect(hall_cols_querer, colnames(hallmarks))

sec0_tabs2 <- hallmarks %>%
  select(all_of(hall_cols_exist)) %>%
  mutate(
    Dirección   = ifelse(NES > 0, "Activado en tumor", "Reprimido en tumor"),
    Pathway     = str_remove(Description, "^HALLMARK_") %>%
                  str_replace_all("_", " ") %>% str_to_title(),
    NES         = round(NES, 3),
    pvalue      = signif(pvalue, 3),
    p.adjust    = signif(p.adjust, 3)
  ) %>%
  arrange(NES) %>%
  select(Pathway, NES, `p-valor` = pvalue, `FDR` = p.adjust,
         `Genes en set` = setSize, Dirección, everything(),
         -Description, -any_of("enrichmentScore"))
save_tsv(sec0_tabs2, "Sec0_TabS2_hallmarks_gsea", out_supp)

# =============================================================================
# OE1 — TABLAS
# =============================================================================
cat("\n--- OE1: Tablas bases de datos de fármacos ---\n")

# ── OE1_Tab1: Resumen por base de datos ──────────────────────────────────────
# Calcular métricas desde archivos reales
n_de_genes <- nrow(de_sig)

# DGIdb
dgi_gene_col <- grep("query_gene|gene_symbol|symbol", names(dgidb),
                     ignore.case = TRUE, value = TRUE)[1]
dgi_drug_col <- grep("drug_name|drug", names(dgidb),
                     ignore.case = TRUE, value = TRUE)[1]
dgi_n_genes  <- n_distinct(dgidb[[dgi_gene_col]])
dgi_n_drugs  <- n_distinct(dgidb[[dgi_drug_col]])
dgi_n_pairs  <- nrow(dgidb)

# ChEMBL
chmb_gene_col   <- grep("gene_symbol|symbol", names(chembl), ignore.case=TRUE, value=TRUE)[1]
chmb_n_genes    <- n_distinct(chembl[[chmb_gene_col]])
chmb_n_approved <- sum(chembl$max_phase >= 4, na.rm = TRUE)
chmb_n_phase3   <- sum(chembl$max_phase == 3, na.rm = TRUE)

# OpenTargets
ot_gene_col  <- grep("gene_symbol|symbol", names(ot_drugs), ignore.case=TRUE, value=TRUE)[1]
ot_drug_col  <- grep("drug.*name|pref_name|drug.*label", names(ot_drugs),
                     ignore.case=TRUE, value=TRUE)[1]
ot_n_genes   <- n_distinct(ot_drugs[[ot_gene_col]])
ot_n_drugs   <- n_distinct(ot_drugs[[ot_drug_col]])

# L2S2
l2s2_pert_col <- grep("drug_name|pert$|pert_name", names(l2s2),
                      ignore.case=TRUE, value=TRUE)[1]
l2s2_n_drugs  <- n_distinct(l2s2[[l2s2_pert_col]])

oe1_tab1 <- tibble(
  `Base de datos`             = c("DGIdb v5", "ChEMBL 34", "Open Targets",  "L2S2"),
  `Descripción`               = c(
    "Drug-Gene Interaction Database — interacciones gen-fármaco consolidadas de ~30 fuentes",
    "European Bioinformatics Institute — compuestos con mecanismo de acción curado",
    "EBI/Sanger — evidencia multi-tipo para asociaciones gen-enfermedad",
    "Library of Signatures and Similarities — conectividad transcriptómica inversa"
  ),
  `Genes DE consultados`      = n_de_genes,
  `Genes con hits`            = c(dgi_n_genes, chmb_n_genes, ot_n_genes, NA_integer_),
  `Fármacos/pert. únicos`     = c(dgi_n_drugs, NA_integer_, ot_n_drugs, l2s2_n_drugs),
  `Pares gen-fármaco`         = c(dgi_n_pairs, nrow(chembl), nrow(ot_drugs), nrow(l2s2)),
  `Criterio de filtro`        = c(
    "Cualquier interacción reportada",
    "Fase clínica >= 3 (Phase III o aprobado)",
    "Cualquier fármaco conocido para el gen",
    "Conectividad negativa (reversión de firma)"
  ),
  `Aprobados (Fase IV)`       = c(NA_integer_, chmb_n_approved, NA_integer_, NA_integer_),
  `Fase III`                  = c(NA_integer_, chmb_n_phase3,   NA_integer_, NA_integer_)
)
save_tsv(oe1_tab1, "OE1_Tab1_resumen_bases_datos", out_main)

# ── OE1_Tab2: Top candidatos por número de fuentes ───────────────────────────
# Columnas útiles de la tabla multi_source
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

# ── OE1_TabS1: Pares gen-fármaco DGIdb ───────────────────────────────────────
oe1_tabs1 <- dgidb %>%
  select(
    Gen               = any_of(c("query_gene", "gene_symbol")),
    `Fármaco`         = any_of(c("drug_name")),
    `Tipo interacción` = any_of(c("interaction_types")),
    `Score`           = any_of(c("interaction_score")),
    `log2FC`          = any_of(c("logFC", "logFC_TVsS")),
    `FDR`             = any_of(c("adj_pval", "adj.P.Val_TVsS")),
    `Dirección DE`    = any_of(c("direction"))
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  arrange(desc(abs(`log2FC`)))
save_tsv(oe1_tabs1, "OE1_TabS1_dgidb_pares", out_supp)

# ── OE1_TabS2: Fármacos ChEMBL ───────────────────────────────────────────────
chembl_cols_querer <- c("gene_symbol", "pref_name", "max_phase", "phase_label",
                         "first_approval", "molecule_type", "mechanisms",
                         "n_up", "n_down")
chembl_cols_exist  <- intersect(chembl_cols_querer, colnames(chembl))
oe1_tabs2 <- chembl %>%
  select(all_of(chembl_cols_exist)) %>%
  rename_with(~c("Gen", "Fármaco", "Fase", "Fase (etiqueta)",
                  "Año aprobación", "Tipo molécula", "Mecanismo",
                  "Genes up", "Genes down")[seq_along(.)]) %>%
  arrange(desc(Fase))
save_tsv(oe1_tabs2, "OE1_TabS2_chembl_drugs", out_supp)

# ── OE1_TabS3: Fármacos OpenTargets ──────────────────────────────────────────
ot_cols_querer <- c("gene_symbol", "drug_name", "pref_name", "max_phase",
                     "drug_type", "action_type", "indication", "logFC_TVsS",
                     "direction")
ot_cols_exist  <- intersect(ot_cols_querer, colnames(ot_drugs))
oe1_tabs3 <- ot_drugs %>%
  select(all_of(ot_cols_exist)) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  arrange(desc(abs(logFC_TVsS)))
save_tsv(oe1_tabs3, "OE1_TabS3_opentargets_drugs", out_supp)

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
    `CMap score`      = round(cmap_score, 3),
    `Composite score` = round(composite_score, 3)
  ) %>%
  select(
    Fármaco,
    Mecanismo,
    `Target primario`   = primary_target,
    `Indicación actual`,
    `Fase clínica máx.`,
    `CMap score`,
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
# EXCEL UNIFICADO
# =============================================================================
cat("\n-- Generando Excel unificado...\n")

wb <- createWorkbook()

# Función auxiliar para añadir hoja con formato básico
add_sheet <- function(wb, df, sheet_name, title = NULL) {
  addWorksheet(wb, sheet_name)
  if (!is.null(title)) {
    writeData(wb, sheet_name, x = title, startRow = 1, startCol = 1)
    addStyle(wb, sheet_name,
             style = createStyle(fontSize = 11, textDecoration = "bold"),
             rows = 1, cols = 1)
    writeDataTable(wb, sheet_name, x = df, startRow = 3, tableStyle = "TableStyleLight2")
    setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
  } else {
    writeDataTable(wb, sheet_name, x = df, startRow = 1, tableStyle = "TableStyleLight2")
    setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
  }
}

add_sheet(wb, sec0_tab1,  "Sec0_Tab1",  "Resumen análisis DE — HNSCC Tumor vs. Tejido Normal")
add_sheet(wb, sec0_tab2,  "Sec0_Tab2",  "Top 40 proteínas DE (20 sobreexpresadas + 20 subexpresadas)")
add_sheet(wb, sec0_tabs1, "Sec0_TabS1", "Tabla completa — 666 proteínas DE significativas")
add_sheet(wb, sec0_tabs2, "Sec0_TabS2", "GSEA Hallmarks MSigDB — todos los gene sets significativos")
add_sheet(wb, oe1_tab1,   "OE1_Tab1",   "Resumen consulta a bases de datos farmacológicas")
add_sheet(wb, oe1_tab2,   "OE1_Tab2",   "Top 30 candidatos por número de fuentes de respaldo")
add_sheet(wb, oe1_tabs1,  "OE1_TabS1",  "Pares gen-fármaco DGIdb completo")
add_sheet(wb, oe1_tabs2,  "OE1_TabS2",  "Fármacos ChEMBL fase >= 3")
add_sheet(wb, oe1_tabs3,  "OE1_TabS3",  "Fármacos Open Targets")
add_sheet(wb, oe2_tab1,   "OE2_Tab1",   "Candidatos EGFR LOD-stable — validación del método (26 fármacos)")
add_sheet(wb, oe2_tab2,   "OE2_Tab2",   "Candidatos no-EGFR LOD-stable — repurposing real (6 fármacos)")
add_sheet(wb, oe2_tabs1,  "OE2_TabS1",  "Candidatos extendidos no-EGFR — robustos a pesos, no LOD-stable")

excel_path <- "results/tables/pub/HNSCC_DrugRepurposing_Tables_Sec0_OE1_OE2.xlsx"
saveWorkbook(wb, excel_path, overwrite = TRUE)
cat(sprintf("  Excel guardado: %s\n", excel_path))

# =============================================================================
# RESUMEN
# =============================================================================
cat("\n============================================================\n")
cat("Script 18 COMPLETO\n")
cat("Output main:  results/tables/pub/main/\n")
cat("Output supp:  results/tables/pub/supp/\n")
cat("Excel:        results/tables/pub/HNSCC_DrugRepurposing_Tables_Sec0_OE1_OE2.xlsx\n")
cat(sprintf("Fin: %s\n", format(Sys.time())))
sink()
