# Methods — HNSCC Drug Repurposing Pipeline

*Generated automatically by 14_methods_summary.R on 2026-06-11*

---

## Data

Quantitative proteomics data (Data-Independent Acquisition, DIA) from 10 paired
tumor/normal tissue samples from head and neck squamous cell carcinoma (HNSCC)
patients (n=20 total; 6 HPV-positive, 4 HPV-negative pairs) were processed
with MaxQuant and differential expression analysis performed using proteoDA
(a limma wrapper). The primary comparison was Tumor vs. Adjacent Normal (TVsS), averaging
over HPV status (TVsS = Tumor vs. Sano/Adjacent Normal, column suffix in all
output files). Proteins were considered significantly differentially expressed
at |log2FC| > 1 and adjusted P-value
(Benjamini-Hochberg) < 0.05.

---

## Phase 1: Data Preparation

### Identifier mapping (Script 02)

UniProt IDs were mapped to Entrez gene IDs and HGNC gene symbols using
`org.Hs.eg.db` via `bitr()` in the clusterProfiler package.

---

## Phase 2: Functional Enrichment Analysis (Script 03)

Gene Set Enrichment Analysis (GSEA) was performed using `clusterProfiler`
(v4.14.6) with proteins ranked by a pi-statistic
(sign(log2FC) x |log2FC| x -log10(adjusted P)). Gene set collections: MSigDB
Hallmarks, Gene Ontology Biological Process, KEGG pathways, and Reactome
pathways (ReactomePA). Parameters: minGSSize=10, maxGSSize=500, P-adjusted < 0.05 (Benjamini-Hochberg). A fixed random seed
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

STRING v12 REST API queried for DE proteins; combined_score >= 700 (high confidence).
Metrics: degree, betweenness, eigenvector centrality (igraph v2.2.2). Hubs: top 10% by degree.

---

## Phase 4: Multi-Criteria Scoring (Script 10)

Six-dimension composite score:
- |log2FC| (w=0.2)
- Significance (w=0.15)
- Clinical phase (w=0.2)
- L2S2 reversal (w=0.1)
- Pathway relevance (w=0.15): fraction of drug DE targets in the pooled GSEA leading-edge (GO BP + KEGG + Reactome)
- Network centrality (w=0.15)

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

Publication figures were produced in R (ggplot2 v3.x, ComplexHeatmap v2.22.0) under a single centralized style module
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

**R** (R version 4.4.3 (2025-02-28))

Key R packages:
- limma v3.62.2
- clusterProfiler v4.14.6
- ReactomePA v1.50.0
- igraph v2.2.2
- ggraph v2.2.2
- openxlsx v4.2.8.1
- ComplexHeatmap v2.22.0

**Python 3** (scripts 04, 05, 06, 11, 12)
Key packages: requests, pandas, numpy, matplotlib

**Databases:** DGIdb v5, ChEMBL v33, Open Targets (Mar 2026),
L2S2 LINCS L1000 (l2s2.maayanlab.cloud, Evangelista et al. 2025), STRING v12, ClinicalTrials.gov API v2,
NCBI PubMed, MSigDB Hallmark, NCG7

**Analysis date:** 2026-06-11

**Reproducibility:** Parameters in `config/analysis_params.yaml`; figure style
centralized in `scripts/_fig_style.R`; execution order in `docs/RUNBOOK.md`.

