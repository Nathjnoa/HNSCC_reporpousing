# Methods — HNSCC Drug Repurposing Pipeline

*Generated automatically by 14_methods_summary.R on 2026-03-04*

---

## Data

Quantitative proteomics data (Data-Independent Acquisition, DIA) from 10 paired
tumor/normal tissue samples from head and neck squamous cell carcinoma (HNSCC)
patients (n=20 total; 6 HPV-positive, 4 HPV-negative pairs) were processed
with MaxQuant and differential expression analysis performed using proteoDA
(a limma wrapper). The primary comparison was Tumor vs. Normal (TVsS), averaging
over HPV status. Proteins were considered significantly differentially expressed
at |log2FC| > 1 and adjusted P-value
(Benjamini-Hochberg) < 0.05.

---

## Phase 1: Data Preparation

### Identifier mapping (Script 02)

UniProt IDs were mapped to Entrez gene IDs and HGNC gene symbols using
`org.Hs.eg.db` via `bitr()` in the clusterProfiler package.

---

## Phase 2: Functional Enrichment Analysis (Script 03)

Over-representation analysis (ORA) was performed using `clusterProfiler`
(v4.14.6) for Gene Ontology (BP, MF, CC),
KEGG pathways, and Reactome pathways (ReactomePA). Background: all quantified
proteins with Entrez IDs (n=3,262). Significance thresholds: P-adjusted < 0.05, Q-value < 0.2.
GO BP terms were simplified using `simplify()` (semantic similarity cutoff=0.7).

GSEA was performed against GO BP terms and MSigDB Hallmark gene sets.
Parameters: minGSSize=10, maxGSSize=500, permutations=1000.

---

## Phase 3: Drug-Target Database Queries (Scripts 04-07)

- **DGIdb**: GraphQL API v5, all 520 DE genes queried.
- **ChEMBL**: REST API v33, drugs in clinical phase >= 3.
- **Open Targets**: GraphQL API, HNSCC association (EFO_0000181) + known drugs.
- **CMap2**: `gess_cmap()` from signatureSearch; query = top 150 up + 150 down proteins
  by log2FC; reversors = bottom 5% scaled_score (threshold = -0.718).

---

## Phase 3: Drug Integration (Script 08)

Four sources unified; drugs classified A (HNSCC-approved), B (other cancer),
C (non-oncology), D (experimental). Multi-source candidates (>= 2 databases)
prioritized for scoring.

---

## Phase 4: PPI Network (Script 09)

STRING v12 REST API queried for DE proteins; combined_score >= 700 (high confidence).
Metrics: degree, betweenness, eigenvector centrality (igraph v2.2.1). Hubs: top 10% by degree.

---

## Phase 4: Multi-Criteria Scoring (Script 10)

Six-dimension composite score:
- |log2FC| (w=0.2)
- Significance (w=0.15)
- Clinical phase (w=0.2)
- CMap reversal (w=0.15)
- Pathway relevance (w=0.15)
- Network centrality (w=0.15)

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

**R** (R version 4.4.3 (2025-02-28))

Key R packages:
- limma v3.62.2
- clusterProfiler v4.14.6
- ReactomePA v1.50.0
- signatureSearch v1.20.0
- igraph v2.2.1
- ggraph v2.2.2
- openxlsx v4.2.8.1
- ComplexHeatmap v2.22.0

**Python 3** (scripts 04, 05, 06, 11, 12)
Key packages: requests, pandas, numpy, matplotlib

**Databases:** DGIdb v5, ChEMBL v33, Open Targets (Mar 2026),
CMap2 (EH3223), STRING v12, ClinicalTrials.gov API v2,
NCBI PubMed, MSigDB Hallmark, NCG7

**Analysis date:** 2026-03-04

**Reproducibility:** Parameters in `config/analysis_params.yaml`;
scripts 01-14 in `scripts/`; execution order in `docs/RUNBOOK.md`.

