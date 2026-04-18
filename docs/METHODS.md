# Methods — HNSCC Drug Repurposing Pipeline

---

## Data

Quantitative proteomics data (Data-Independent Acquisition, DIA) from 10 paired
tumor/adjacent normal tissue samples from head and neck squamous cell carcinoma (HNSCC)
patients (n=20 total; 6 HPV-positive, 4 HPV-negative pairs) were processed
with MaxQuant and differential expression analysis performed using proteoDA
(a limma wrapper). The primary comparison was Tumor vs. Adjacent Normal (TVsS).
Proteins were considered significantly differentially expressed at |log2FC| > 1
and adjusted P-value (Benjamini-Hochberg) < 0.05.

---

## Phase 1: Data Preparation

### Differential expression model (Script 01)

The proteoDA output (limma TVsS contrast, HPV-adjusted with `duplicateCorrelation()`
for the paired tumor/normal design) was parsed to extract proteins in three tiers:
all quantified proteins, significantly DE proteins, and directional subsets
(upregulated / downregulated in tumor). Total: 666 significant proteins
(329 up, 337 down) from 3,352 quantified.

### Identifier mapping (Script 02)

UniProt IDs were mapped to Entrez gene IDs and HGNC gene symbols using
`org.Hs.eg.db` via `bitr()` in the clusterProfiler package. 3,262 of 3,352
proteins were mapped (97.3%).

---

## Phase 2: Drug-Target Database Queries (Scripts 04–07)

Four independent pharmacological databases were queried against all 666
significantly DE proteins:

- **DGIdb v5**: GraphQL API; all reported drug-gene interactions.
- **ChEMBL 34**: REST API; drugs in clinical phase ≥ 3 (approved or Phase III).
- **Open Targets Platform**: GraphQL API (api.platform.opentargets.org); drug
  candidates via `drugAndClinicalCandidates` field (API v4, April 2026);
  approval detected by `maximumClinicalStage = "APPROVAL"`.
- **L2S2 (LINCS L1000 Signature Search)**: public GraphQL API
  (l2s2.maayanlab.cloud); bidirectional enrichment (top 150 proteins UP vs.
  drug DOWN signatures, and vice versa); `filterFda=TRUE`, P < 0.001;
  1,044 FDA-approved drugs across 248 cell lines. Reversal score =
  −log10(best_p) × mean_log2(OR) × log2(n_sigs+1), normalized to [−1, 0].

---

## Phase 3: Drug Integration (Script 08)

Outputs from all four databases were unified into a master table and each drug
classified as: A (approved for HNSCC), B (approved for another cancer),
C (approved non-oncology), or D (experimental). Candidates present in ≥ 2
databases were prioritized as multi-source candidates.

---

## Phase 4: PPI Network Analysis (Script 09)

A protein-protein interaction (PPI) network was constructed using the STRING v12
REST API (combined_score ≥ 700, high confidence). Community detection was
performed with the Louvain algorithm. Hub proteins were defined as the top 10%
by betweenness centrality within each module. Network metrics (degree,
betweenness, eigenvector centrality) were computed with igraph v2.2.2.

---

## Phase 5: Multi-Criteria Scoring (Script 10)

A composite score was computed for each multi-source drug candidate using
five dimensions with the following weights (all configurable in
`config/analysis_params.yaml`):

| Dimension | Weight | Description |
| --- | --- | --- |
| π-statistic | 0.325 | sign(logFC) × \|logFC\| × −log10(adj.P), directional |
| Clinical phase | 0.200 | Approved=1, Phase III=0.75, Phase II=0.5, Phase I=0.25 |
| Network centrality | 0.195 | Betweenness centrality of primary target (normalized) |
| Pathway relevance | 0.150 | Target present in a significantly enriched pathway |
| L2S2 reversal | 0.130 | Connectivity score (more negative = better reversal) |

The top 35 candidates by composite score were retained as the candidate pool
for sensitivity analysis. The final panel was defined by LOD stability
(see Phase 6).

---

## Phase 6: Sensitivity Analysis (Script 15)

Robustness of the drug ranking was assessed through three approaches:

1. **Weight sensitivity**: six alternative weight configurations spanning
   plausible ranges for each scoring dimension.
2. **Leave-One-Database (LOD) analysis**: the scoring procedure was repeated
   five times, each time excluding one of the five scoring sources. A drug was
   classified as LOD-stable if it remained in the top 35 in all five LOD runs.
   LOD stability was the primary criterion for inclusion in the final candidate
   panel (32 LOD-stable candidates).
3. **Permutation test**: composite scores were recomputed 1,000 times with
   randomly permuted drug-gene assignments to derive an empirical null
   distribution. Observed scores above the 95th percentile of the null were
   considered non-random.

---

## Software

**R** (v4.4.3, 2025-02-28)

Key packages: limma v3.62.2, clusterProfiler v4.14.6, igraph v2.2.2,
ggraph v2.2.2, ComplexHeatmap v2.22.0, httr2, jsonlite, openxlsx v4.2.8.1,
dplyr, ggplot2.

**Python 3** (scripts 04–07)

Key packages: requests, pandas, numpy.

**Databases:** DGIdb v5, ChEMBL 34, Open Targets Platform (April 2026),
L2S2 LINCS L1000 (l2s2.maayanlab.cloud), STRING v12.

**Analysis date:** 2026-04-11

**Reproducibility:** All parameters in `config/analysis_params.yaml`.
Scripts 01–02, 04–10, 15, 17–18 in `scripts/`; execution order in
`docs/RUNBOOK.md`.
