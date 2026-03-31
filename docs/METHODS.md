# Methods — HNSCC Drug Repurposing Pipeline

*Última actualización: 2026-03-30 (v2: scoring recalibrado, deduplicación, reclasificación, criterio LOD-stability)*

---

## Data

Quantitative proteomics data (Data-Independent Acquisition, DIA) from 10 paired
tumor/normal tissue samples from head and neck squamous cell carcinoma (HNSCC)
patients (n=20 total; 6 HPV-positive, 4 HPV-negative pairs) were processed
with MaxQuant. A total of 3,352 proteins were quantified across all samples.

---

## Phase 1: Differential Expression Analysis (Script 01)

Differential expression was performed using **limma** with a HPV-adjusted model
to account for the paired patient design and HPV status as a confounding variable.

Model: `~ 0 + condition + vph`

Within-patient correlation for paired samples (tumor/normal from the same patient)
was estimated using `duplicateCorrelation()` with `patient_id` as blocking factor,
then incorporated into `lmFit()`. Contrast: Tumor − Adjacent Normal (TVsS).
Moderated t-statistics were computed with `eBayes(trend=TRUE)`.

Proteins were considered significantly differentially expressed at |log2FC| > 1
and adjusted P-value (Benjamini-Hochberg) < 0.05. Result: **3,352 proteins
quantified; 666 significant (329 up, 337 down)**.

---

## Phase 1: Identifier Mapping (Script 02)

UniProt IDs were mapped to Entrez gene IDs and HGNC gene symbols using
`org.Hs.eg.db` via `bitr()` in the clusterProfiler package.

---

## Phase 2: Functional Enrichment Analysis (Script 03)

Over-representation analysis (ORA) was performed using `clusterProfiler`
(v4.14.6) for Gene Ontology (BP, MF, CC), KEGG pathways, and Reactome pathways
(ReactomePA). Background: all quantified proteins with Entrez IDs (n=3,262).
Significance thresholds: P-adjusted < 0.05, Q-value < 0.2. GO BP terms were
simplified using `simplify()` (semantic similarity cutoff=0.7).

GSEA was performed against GO BP terms and MSigDB Hallmark gene sets using a
**pi-statistic ranked list**: `π = sign(logFC) × |logFC| × −log10(adj.P.Val + 1e-300)`.
This combined effect size and significance into a single ranking metric, avoiding
loss of information from using logFC alone. Parameters: minGSSize=10,
maxGSSize=500, permutations=1000.

---

## Phase 3: Drug-Target Database Queries (Scripts 04–07)

All queries include a `query_date` column recording the date of data retrieval
to enable reproducibility assessment over time.

- **DGIdb**: GraphQL API v5, all 666 DE proteins queried. Result: 3,542 interactions,
  2,697 unique drugs.
- **ChEMBL**: REST API v33, drugs in clinical phase ≥ 3. Result: 113 unique drugs,
  74 approved + 39 Phase III.
- **Open Targets**: GraphQL API, HNSCC association (EFO_0000181) with
  `min_hnscc_score ≥ 0.2` filter applied. Known drugs queried per DE gene.
  Result: 1,165 gene-drug pairs, 173 unique drugs.
- **L2S2 (LINCS L1000 Signature Search)**: Public GraphQL API
  (l2s2.maayanlab.cloud); bidirectional enrichment: top 150 UP proteins vs.
  drug DOWN signatures + top 150 DOWN vs. drug UP signatures; filterFda=TRUE;
  p-value < 0.001; 248 cell lines; 1,044 FDA-approved drugs evaluated.
  `reversal_score = −log10(best_p) × mean_log2(OR) × log2(n_sigs+1)`,
  normalized to [−1, 0]. **Limitation:** L2S2 uses transcriptomic (mRNA)
  signatures; this analysis uses proteomics data. Connectivity scores reflect
  mRNA-level perturbation similarity, which may not fully capture
  protein-level effects.

---

## Phase 3: Drug Integration (Script 08)

Four sources unified into a master drug-target table. Prior to aggregation, drug name
normalization was applied to collapse pharmaceutical salt forms and equivalent
formulations into canonical names (e.g., divalproex sodium / valproate sodium /
sodium valproate → valproic acid; acetyldigitoxin / deslanoside → digitoxin;
afatinib dimaleate → afatinib; lapatinib ditosylate → lapatinib).

Drugs classified by indication:

- **Class A**: HNSCC-approved — used as pipeline validation (positive controls);
  not included in repurposing candidates. Cetuximab, pembrolizumab, and nivolumab
  included via manual HNSCC-approved override to correct EFO annotation gaps.
- **Class B**: Approved for other cancers (intra-oncology repurposing). Two drugs
  not captured by Open Targets — neratinib (FDA 2017, HER2+ breast cancer) and
  lazertinib (FDA 2024, NSCLC EGFR+) — were manually corrected from Class C to B.
- **Class C**: Approved for non-oncology indications (classical repurposing)
- **Class D**: Experimental/investigational (phase ≤ 3)

Non-antitumoral agents excluded from candidates: diagnostic imaging agents
(technetium Tc 99m succimer), topical enzymes (collagenase clostridium histolyticum),
ophthalmic proteases (ocriplasmin), anticoagulants, anti-amyloid antibodies,
gabapentinoids, and alcohol dehydrogenase inhibitors (see `config/analysis_params.yaml`).

Multi-source candidates (≥ 2 databases) prioritized for scoring.

---

## Phase 4: PPI Network Analysis (Script 09)

STRING v12 REST API queried for DE proteins; combined_score ≥ 700 (high confidence).
Network constructed with igraph v2.2.1.

**Community detection**: Louvain algorithm (`cluster_louvain()`) applied to the
giant component to identify functional protein modules. 22 modules detected;
17 annotated with GO BP terms via ORA (clusterProfiler).

**Hub definition**: Within each Louvain module, proteins in the top 10%
by betweenness centrality (normalized) were classified as module hubs. This
module-level definition avoids bias toward globally central proteins (e.g.,
ribosomal proteins) that may not be druggable targets.
Result: 74 module hubs across 22 modules; 498 nodes, 2,698 edges in giant component.

Standard centrality metrics (degree, betweenness, eigenvector) computed for
all nodes.

---

## Phase 4: Multi-Criteria Scoring (Script 10)

Five-dimension composite score (weights sum to 1.00):

| Dimension | Weight | Description |
|-----------|--------|-------------|
| π-statistic (proteomics) | 0.325 | `sign(logFC) × abs(logFC) × -log10(adj.P.Val + ε)` |
| Clinical phase | 0.200 | Normalized max development phase across databases |
| Network centrality | 0.195 | Module diversity score in STRING PPI network |
| Pathway relevance | 0.150 | Overlap of drug targets with enriched GO/Hallmark pathways |
| L2S2 connectivity | 0.130 | Transcriptomic reversal score (L2S2/LINCS) |

The clinical literature evidence score (ClinicalTrials.gov + PubMed hit count) was
excluded from the composite to avoid publication bias (well-studied drugs would score
higher regardless of proteomics signal) and circularity with the manual literature
review performed by the research group. It is retained as a descriptive column.

Class A drugs (HNSCC-approved) received no scoring bonus and are reported separately
as positive controls validating pipeline performance.

Top 35 candidates selected by composite score without per-target diversity limits.
Scores normalized [0, 1] per dimension before weighting.

---

## Phase 5: Evidence Validation (Scripts 11–12)

- **ClinicalTrials.gov API v2**: drug + HNSCC keyword search per candidate.
- **PubMed E-utilities**: drug name (salt-simplified) + HNSCC MeSH query.
- **Cancer driver overlap**: COSMIC CGC (if available), NCG7, HNSCC literature.

---

## Phase 5: Final Integration (Script 13)

Combined rank score = 0.60 × multi-criteria score + 0.40 × clinical evidence score.
Evidence levels 1–4 based on drug approval status and HNSCC indication.

---

## Phase 6: Sensitivity Analysis (Script 15)

Robustness evaluated through three complementary analyses:

1. **Weight sensitivity**: 6 alternative weight configurations (clinical-heavy,
   molecular-heavy, network-heavy, pathway-heavy, equal weights) applied to
   assess rank stability across the top 35.
2. **Leave-one-database (LOD)**: Each of the 4 drug databases excluded in turn;
   candidates appearing in the top 35 across all 4 exclusion scenarios classified
   as `lod_stable = TRUE`. This criterion defines the final candidate panel,
   as it does not depend on arbitrary score thresholds or weight choices.
3. **Permutation test** (n=1,000): π-statistics randomly shuffled across proteins;
   empirical p-value computed for top-1 composite score (p = 0.058).

Result: **26 LOD-stable candidates** identified (25 repurposing candidates +
Cetuximab as Class A positive control). Of these, 13 belong to the EGFR inhibitor
cluster (Class B); remaining 12 span OXPHOS, cardiac glycosides, DNMT inhibitors,
MAO inhibitors, HDAC modulators, and other mechanisms.

---

## Software

**R** (R version 4.4.3 (2025-02-28))

Key R packages:
- limma v3.62.2
- clusterProfiler v4.14.6
- ReactomePA v1.50.0
- igraph v2.2.1
- ggraph v2.2.2
- openxlsx v4.2.8.1
- ComplexHeatmap v2.22.0

**Python 3** (scripts 04, 05, 06, 07, 11, 12)
Key packages: requests, pandas, numpy, matplotlib

**Databases:** DGIdb v5, ChEMBL v33, Open Targets (Mar 2026),
L2S2 LINCS L1000 (l2s2.maayanlab.cloud, Evangelista et al. 2025), STRING v12,
ClinicalTrials.gov API v2, NCBI PubMed, MSigDB Hallmark, NCG7

**Analysis date:** 2026-03-08

**Reproducibility:** Parameters in `config/analysis_params.yaml`;
scripts 01–17 in `scripts/`; execution order in `docs/RUNBOOK.md`.
