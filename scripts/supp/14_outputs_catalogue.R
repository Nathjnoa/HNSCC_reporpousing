#!/usr/bin/env Rscript
# =============================================================================
# 14_outputs_catalogue.R — Genera docs/OUTPUTS.md (catálogo de archivos)
# =============================================================================
# Reemplaza al antiguo 14_methods_summary.R como generador de OUTPUTS.md.
#
# NOTA IMPORTANTE: docs/METHODS.md YA NO se autogenera — se mantiene a mano
# (contiene Phase 6 dos-cohortes, Figure Generation, etc., más ricos que cualquier
# plantilla). Este script SOLO escribe el catálogo mecánico de outputs.
#
# Output: docs/OUTPUTS.md — catálogo completo de archivos generados con tamaños
#
# Ambiente: omics-R
# Ejecución (desde raíz del proyecto, tras correr el pipeline):
#   Rscript scripts/supp/14_outputs_catalogue.R
# =============================================================================

suppressPackageStartupMessages(library(fs))
setwd(here::here())

cat_files <- function(dir_path, indent = "") {
  if (!dir_exists(dir_path)) return(character(0))
  files <- dir_ls(dir_path, recurse = TRUE, type = "file")
  rel_paths <- path_rel(files, start = getwd())
  sizes <- file_size(files)
  paste0(indent, "- `", rel_paths, "` (", format(sizes, units = "auto"), ")")
}

tables_prioritization <- {
  all_table_files <- cat_files("results/tables")
  grep("/(10_|15_)", all_table_files, value = TRUE)
}

outputs_text <- paste0(
"# Outputs Catalogue — HNSCC Drug Repurposing

*Generated: ", Sys.Date(), "*

---

## Tables

### DE Analysis (scripts 01-02)
", paste(cat_files("results/tables/de_limma"), collapse = "\n"), "

### Pathway Enrichment (script 03)
", paste(cat_files("results/tables/pathway_enrichment"), collapse = "\n"), "

### Drug Targets (scripts 04-08)
", paste(cat_files("results/tables/drug_targets"), collapse = "\n"), "

### Network (script 09)
", paste(cat_files("results/tables/network"), collapse = "\n"), "

### Prioritization & Sensitivity (scripts 10, 15)
", paste(tables_prioritization, collapse = "\n"), "

### Publication Tables (scripts 18, 16)
", paste(c(cat_files("results/tables/pub/main"), cat_files("results/tables/pub/supp")), collapse = "\n"), "

---

## Figures

### Pathway Figures (script 03)
", paste(cat_files("results/figures/pathway_enrichment"), collapse = "\n"), "

### Publication Figures — main (scripts 17, 17b-17i, 16)
", paste(cat_files("results/figures/pub/main"), collapse = "\n"), "

### Publication Figures — supplementary (scripts 17h, 16)
", paste(cat_files("results/figures/pub/supp"), collapse = "\n"), "

---

## Key Output Files

- `results/tables/10_top20_candidates.tsv` — Top 20 prioritized candidates with composite scores
- `results/tables/10_module_hub_candidates.tsv` — Module -> hub -> drug prioritization summary
- `results/tables/15_lod_stability.tsv` — LOD-stability panel (leave-one-database robustness)
- `results/tables/pub/main/` — Publication tables Tab1-Tab6 (DE summary, top proteins/candidates, EGFR axis, candidates by module, CPTAC+TCGA concordance)
- `results/figures/pub/main/` — Publication figures: panels (PDF + PNG) + `Fig{2,3,4,5,6}_multipanel.tif` (600-dpi TIFF)
- `results/figures/pub/supp/` — Supplementary figures (robustness heatmap, selection funnel, OS survival)

---

*To regenerate all outputs, follow `docs/RUNBOOK.md` in order.*
")

writeLines(outputs_text, "docs/OUTPUTS.md")
cat("[14] OUTPUTS.md regenerado.\n")
