#!/usr/bin/env Rscript
# =============================================================================
# 14_methods_summary.R
# =============================================================================
# Genera la sección METHODS del manuscrito con todos los parámetros exactos
# usados en el pipeline, y el catálogo de outputs (OUTPUTS.md).
#
# Outputs:
#   docs/METHODS.md   — sección Methods lista para insertar en manuscrito
#   docs/OUTPUTS.md   — catálogo completo de archivos generados
#
# Ejecución:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/14_methods_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(fs)
})

# ── Working directory (raíz del proyecto vía scripts/_setup.R) ───────────────
source(here::here("scripts", "_setup.R"))
setup_project()
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_dir   <- "logs"
dir.create(log_dir, showWarnings = FALSE)

cat("[14] Generando METHODS.md y OUTPUTS.md...\n")

# Cargar parámetros de config
cfg <- tryCatch(
  yaml::read_yaml("config/analysis_params.yaml"),
  error = function(e) list()
)

# ── Versiones de paquetes clave ───────────────────────────────────────────────
r_packages <- c("limma", "clusterProfiler", "ReactomePA",
                 "igraph", "ggraph", "openxlsx", "ComplexHeatmap")
pkg_versions <- sapply(r_packages, function(p) {
  tryCatch(as.character(packageVersion(p)), error = function(e) "N/A")
})

py_script_pkgs <- c("requests", "pandas", "numpy", "matplotlib")

# Operador %||% — definido antes de cualquier uso en paste0()
`%||%` <- function(a, b) if (is.null(a)) b else a

# ── Construir METHODS.md ──────────────────────────────────────────────────────
methods_text <- paste0(
'# Methods — HNSCC Drug Repurposing Pipeline

*Generated automatically by 14_methods_summary.R on ', Sys.Date(), '*

---

## Data

Quantitative proteomics data (Data-Independent Acquisition, DIA) from 10 paired
tumor/normal tissue samples from head and neck squamous cell carcinoma (HNSCC)
patients (n=20 total; 6 HPV-positive, 4 HPV-negative pairs) were processed
with MaxQuant and differential expression analysis performed using proteoDA
(a limma wrapper). The primary comparison was Tumor vs. Adjacent Normal (TVsS), averaging
over HPV status (TVsS = Tumor vs. Sano/Adjacent Normal, column suffix in all
output files). Proteins were considered significantly differentially expressed
at |log2FC| > ', cfg$de$log2fc_threshold %||% 1.0, ' and adjusted P-value
(Benjamini-Hochberg) < ', cfg$de$adj_pval_threshold %||% 0.05, '.

---

## Phase 1: Data Preparation

### Identifier mapping (Script 02)

UniProt IDs were mapped to Entrez gene IDs and HGNC gene symbols using
`org.Hs.eg.db` via `bitr()` in the clusterProfiler package.

---

## Phase 2: Functional Enrichment Analysis (Script 03)

Gene Set Enrichment Analysis (GSEA) was performed using `clusterProfiler`
(v', pkg_versions["clusterProfiler"], ') with proteins ranked by a pi-statistic
(sign(log2FC) x |log2FC| x -log10(adjusted P)). Gene set collections: MSigDB
Hallmarks, Gene Ontology Biological Process, KEGG pathways, and Reactome
pathways (ReactomePA). Parameters: minGSSize=', cfg$pathway$min_gs_size %||% 10,
', maxGSSize=', cfg$pathway$max_gs_size %||% 500, ', P-adjusted < ',
cfg$pathway$pval_cutoff %||% 0.05, ' (Benjamini-Hochberg). A fixed random seed
(set.seed(42), seed=TRUE) was used for reproducible permutation results.

The Hallmarks GSEA result is presented as the narrative enrichment panel of
Figure 2 (GSEA panel). The pooled leading-edge genes of the significant
GO BP + KEGG + Reactome GSEA sets feed the pathway-relevance dimension of the
prioritization score (Script 10). Over-representation analysis (ORA) was not used.

---

## Phase 3: Drug-Target Database Queries (Scripts 04-07)

- **DGIdb**: GraphQL API v5, all 520 DE genes queried.
- **ChEMBL**: REST API v33, drugs in clinical phase >= 3.
- **Open Targets**: GraphQL API, HNSCC association (EFO_0000181) + known drugs.
- **L2S2 (LINCS L1000 Signature Search)**: API GraphQL publica (l2s2.maayanlab.cloud);
  enriquecimiento bidireccional: top 150 proteinas UP vs firmas DOWN del farmaco + top 150 DOWN vs firmas UP;
  filtro filterFda=TRUE; pvalue < 0.001; 248 lineas celulares; 1,044 drugs FDA-aprobados evaluados;
  reversal_score = -log10(best_p) x mean_log2(OR) x log2(n_sigs+1), normalizado a [-1, 0].

---

## Phase 3: Drug Integration (Script 08)

Four sources unified; drugs classified A (HNSCC-approved), B (other cancer),
C (non-oncology), D (experimental). Multi-source candidates (>= 2 databases)
prioritized for scoring.

---

## Phase 4: PPI Network (Script 09)

STRING v12 REST API queried for DE proteins; combined_score >= ',
cfg$network$string_score_threshold %||% 700, ' (high confidence).
Metrics: degree, betweenness, eigenvector centrality (igraph v',
pkg_versions["igraph"], '). Hubs: top ',
100 - (cfg$network$hub_percentile %||% 90), '% by degree.

---

## Phase 4: Multi-Criteria Scoring (Script 10)

Six-dimension composite score:
- |log2FC| (w=', cfg$scoring$weight_log2fc %||% 0.20, ')
- Significance (w=', cfg$scoring$weight_significance %||% 0.15, ')
- Clinical phase (w=', cfg$scoring$weight_clinical_phase %||% 0.20, ')
- L2S2 reversal (w=', cfg$scoring$weight_cmap_connectivity %||% 0.15, ')
- Pathway relevance (w=', cfg$scoring$weight_pathway_relevance %||% 0.15, '): fraction of drug DE targets in the pooled GSEA leading-edge (GO BP + KEGG + Reactome)
- Network centrality (w=', cfg$scoring$weight_network_centrality %||% 0.15, ')

Top 20 selected with diversity filter (max 3 per primary target gene).

---

## Phase 5: Sensitivity & LOD Stability (Script 15)

Weight-perturbation sensitivity analysis and limit-of-detection (LOD) stability
of the candidate ranking; candidates stable across LOD scenarios define the final
panel (`15_lod_stability.tsv`).

---

## Phase 6: External Validation (Script 16)

Differential expression concordance against TCGA-HNSC (TCGAbiolinks, DESeq2):
tumor-vs-normal log2FC correlation between the proteomic cohort and TCGA-HNSC.

---

## Figure Generation (Scripts 17, 17c, 17d)

Publication figures were produced in R (ggplot2 v3.x, ComplexHeatmap v',
pkg_versions["ComplexHeatmap"], ') under a single centralized style module
(`scripts/_fig_style.R`) sourced by every figure script, ensuring consistent
typography, palette and export settings. A colorblind-safe Okabe-Ito palette was
used throughout (up-regulated #D55E00, down-regulated #0072B2). Embedded plot
titles were omitted; panel descriptions are carried in the figure legends.
Individual panels (Script 17: volcano, top-40 differentially expressed protein
heatmap, drug-source and clinical-phase bars, STRING PPI network, biological-module
and drug-class bars; Script 17c: MSigDB Hallmarks GSEA dotplot) were exported as
300-dpi PNG (review) and vector PDF. Multipanel composite figures (Script 17d)
were assembled with patchwork — ComplexHeatmap panels captured as grobs via
grid.grabExpr() — and exported as 600-dpi LZW-compressed TIFF (journal submission)
plus vector PDF. Figure 2 combines the differential-proteome volcano (panel A),
the MSigDB Hallmarks GSEA dotplot (panel B) and the top-40 DE heatmap annotated by
condition and HPV status (panel C); GSEA gene sets are ordered by the pi-statistic
(sign(log2FC) x |log2FC| x -log10(FDR)).

---

## Software

**R** (', R.version$version.string, ')

Key R packages:
', paste(sprintf("- %s v%s", names(pkg_versions), pkg_versions), collapse = "\n"), '

**Python 3** (scripts 04, 05, 06, 11, 12)
Key packages: requests, pandas, numpy, matplotlib

**Databases:** DGIdb v5, ChEMBL v33, Open Targets (Mar 2026),
L2S2 LINCS L1000 (l2s2.maayanlab.cloud, Evangelista et al. 2025), STRING v12, ClinicalTrials.gov API v2,
NCBI PubMed, MSigDB Hallmark, NCG7

**Analysis date:** ', Sys.Date(), '

**Reproducibility:** Parameters in `config/analysis_params.yaml`; figure style
centralized in `scripts/_fig_style.R`; execution order in `docs/RUNBOOK.md`.
')

# Escribir METHODS.md
writeLines(methods_text, "docs/METHODS.md")
cat("[14] METHODS.md generado\n")

# ── Construir OUTPUTS.md ──────────────────────────────────────────────────────
# Catalogar todos los archivos en results/
cat_files <- function(dir_path, indent = "") {
  if (!dir_exists(dir_path)) return(character(0))
  files <- dir_ls(dir_path, recurse = TRUE, type = "file")
  rel_paths <- path_rel(files, start = getwd())
  sizes <- file_size(files)
  paste0(indent, "- `", rel_paths, "` (", format(sizes, units = "auto"), ")")
}

# Archivos de priorización (script 10) y sensibilidad (script 15) en results/tables raíz
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

### Publication Tables (script 18)
", paste(c(cat_files("results/tables/pub/main"), cat_files("results/tables/pub/supp")), collapse = "\n"), "

---

## Figures

### QC Figures (script 01)
", paste(cat_files("results/figures/qc"), collapse = "\n"), "

### Pathway Figures (script 03)
", paste(cat_files("results/figures/pathway_enrichment"), collapse = "\n"), "

### Drug Target Figures (script 08)
", paste(cat_files("results/figures/drug_targets"), collapse = "\n"), "

### Network Figures (script 09)
", paste(cat_files("results/figures/network"), collapse = "\n"), "

### Publication Figures (scripts 16-19)
", paste(cat_files("results/figures/pub/main"), collapse = "\n"), "

---

## Key Output Files

- `results/tables/10_top20_candidates.tsv` — Top 20 prioritized candidates with composite scores
- `results/tables/15_lod_stability.tsv` — LOD-stability panel (sensitivity analysis)
- `results/tables/pub/main/` — Publication tables (DE summary, top candidates, LOD-stable panels, TCGA concordance)
- `results/figures/pub/main/` — Publication figures (volcano, GSEA, candidates, TCGA validation)

---

*To regenerate all outputs, follow `docs/RUNBOOK.md` in order.*
")

writeLines(outputs_text, "docs/OUTPUTS.md")
cat("[14] OUTPUTS.md generado\n")

cat("\n")
cat("[14] Outputs:\n")
cat("  docs/METHODS.md\n")
cat("  docs/OUTPUTS.md\n")
cat("[14] Script 14 completado.\n")
