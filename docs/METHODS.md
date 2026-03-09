# Methods — HNSCC Drug Repurposing Pipeline

*Generated automatically by 14_methods_summary.R on 2026-03-08*

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
Metrics: degree, betweenness, eigenvector centrality (igraph v2.2.1). Hubs: top 10% by degree.

---

## Phase 4: Multi-Criteria Scoring (Script 10)

Six-dimension composite score:
- |log2FC| (w=0.2)
- Significance (w=0.15)
- Clinical phase (w=0.2)
- L2S2 reversal (w=0.15)
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
- igraph v2.2.1
- ggraph v2.2.2
- openxlsx v4.2.8.1
- ComplexHeatmap v2.22.0

**Python 3** (scripts 04, 05, 06, 11, 12)
Key packages: requests, pandas, numpy, matplotlib

**Databases:** DGIdb v5, ChEMBL v33, Open Targets (Mar 2026),
L2S2 LINCS L1000 (l2s2.maayanlab.cloud, Evangelista et al. 2025), STRING v12, ClinicalTrials.gov API v2,
NCBI PubMed, MSigDB Hallmark, NCG7

**Analysis date:** 2026-03-08

**Reproducibility:** Parameters in `config/analysis_params.yaml`;
scripts 01-14 in `scripts/`; execution order in `docs/RUNBOOK.md`.

---

## Limitaciones Metodológicas

### 1. Asunción transcriptómica en L2S2
El análisis de conectividad de firmas (L2S2/LINCS L1000) se basa en respuestas
transcriptómicas (mRNA) de líneas celulares tratadas. Los inputs al análisis
(top 150 proteínas UP/DOWN por logFC proteómico) asumen que los cambios
proteómicos tienen correspondencia directa con cambios en mRNA. La correlación
mRNA-proteína en tumores sólidos es típicamente r ≈ 0.40–0.60 [PMID: 34425047].
El peso de esta dimensión se redujo a 0.10 (vs. 0.15 original) para reflejar
esta limitación. Los resultados de L2S2 deben interpretarse como evidencia
de soporte, no primaria.

### 2. Sesgo de imputación MinProb
Los datos proteómicos fueron imputados con el método MinProb (valores
ausentes → distribución en percentil 2.5). Este método asume que los valores
ausentes son "below detection threshold" (MAR). Si las proteínas mitocondriales
tienen más valores faltantes en tumor que en normal por razones biológicas reales
(downregulation genuina), MinProb amplificará esta señal. El archivo
`results/figures/qc/01_missingness_vs_direction.pdf` y
`results/tables/de_limma/01_missingness_analysis.tsv` documentan el diferencial
de missingness para cada proteína.

### 3. Confounding HPV en el diseño
Los 10 pares tumor/normal incluyen 6 pacientes HPV+ y 4 HPV-. HPV+ y HPV-
tienen perfiles moleculares distintos (vías E6/E7, inmuno-infiltración,
pronóstico diferencial). El modelo limma incluye HPV como covariable aditiva
(~ condición + vph_status), absorbiendo la varianza HPV-dependiente. Sin
embargo, la potencia estadística para detectar efectos de interacción
(proteínas DE específicamente en HPV+ o HPV-) es limitada con n=6/4 por grupo.
Los candidatos del pipeline representan biología compartida entre ambos subtipos.

### 4. Tamaño muestral
n=10 pares tumor/normal es adecuado para efectos de gran magnitud (|logFC| > 1)
pero subóptimo para efectos moderados (0.5–1.0). Ver análisis de potencia en
`docs/FUTURE_WORK.md`.

### 5. Hubs de red: complejos multi-subunidad
La definición de hubs por betweenness centrality (percentil 90 por módulo)
puede identificar múltiples subunidades del mismo complejo físico (e.g., NADH
dehydrogenase Complex I: NDUFA*, NDUFB*, NDUFS*) como hubs independientes.
El script 09 usa detección de módulos Louvain y redefine hubs a nivel de
módulo (intra-módulo betweenness top 10%) para mitigar este efecto. El score
de red en script 10 incluye `module_diversity` (número de módulos distintos
targetados) para penalizar fármacos que solo actúan sobre un complejo.
