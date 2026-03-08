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

Over-representation analysis (ORA) was performed using `clusterProfiler`
(v', pkg_versions["clusterProfiler"], ') for Gene Ontology (BP, MF, CC),
KEGG pathways, and Reactome pathways (ReactomePA). Background: all quantified
proteins with Entrez IDs (n=3,262). Significance thresholds: P-adjusted < ',
cfg$pathway$pval_cutoff %||% 0.05, ', Q-value < ', cfg$pathway$qval_cutoff %||% 0.20, '.
GO BP terms were simplified using `simplify()` (semantic similarity cutoff=',
cfg$pathway$simplify_cutoff %||% 0.70, ').

GSEA was performed against GO BP terms and MSigDB Hallmark gene sets.
Parameters: minGSSize=', cfg$pathway$min_gs_size %||% 10, ', maxGSSize=',
cfg$pathway$max_gs_size %||% 500, ', permutations=1000.

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
- Pathway relevance (w=', cfg$scoring$weight_pathway_relevance %||% 0.15, ')
- Network centrality (w=', cfg$scoring$weight_network_centrality %||% 0.15, ')

Top 20 selected with diversity filter (max 3 per primary target gene).

---

## Phase 5: Evidence Validation (Scripts 11-12)

- **ClinicalTrials.gov API v2**: drug + HNSCC keyword search per candidate.
- **PubMed E-utilities**: drug name (salt-simplified) + HNSCC MeSH query.
- **Cancer driver overlap**: COSMIC CGC (if available), NCG7, HNSCC literature.

---

## Phase 5: Final Integration (Script 13)

Combined rank score = 0.60 x multi-criteria score + 0.40 x clinical evidence score.
Evidence levels 1-4 based on drug approval status and HNSCC indication.

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

**Reproducibility:** Parameters in `config/analysis_params.yaml`;
scripts 01-14 in `scripts/`; execution order in `docs/RUNBOOK.md`.
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

# Archivos de priorización y evidencia (scripts 10-13) en results/tables raíz
tables_prioritization <- {
  all_table_files <- cat_files("results/tables")
  grep("/(10_|11_|12_|13_)", all_table_files, value = TRUE)
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

### Prioritization & Evidence (scripts 10-13)
", paste(c(cat_files("results/tables/evidence"), tables_prioritization), collapse = "\n"), "

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

### Evidence Figures (scripts 11-12)
", paste(cat_files("results/figures/evidence"), collapse = "\n"), "

### Final Figures (script 13)
", paste(cat_files("results/figures/final"), collapse = "\n"), "

---

## Key Output File

`results/tables/13_FINAL_drug_candidates.xlsx` — Master results file with 5 sheets:
1. **Top20_Final**: Top 20 candidates with all evidence integrated
2. **Todos_Candidatos**: All 187 multi-source candidates with scores
3. **Evidencia_Matriz**: Binary evidence matrix (10 dimensions x 20 candidates)
4. **Ensayos_Clinicos**: ClinicalTrials.gov results per candidate
5. **Metodologia**: Pipeline summary with result counts per step

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
