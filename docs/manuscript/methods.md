# Methods (draft, Step 3 · WRITING_PLAN.md)

> Condensed/translated English draft from `docs/METHODS.md` (production-ready), integrating
> the fine network numbers removed from Results R3, the directional π note (M1), the regulatory
> class overrides (m4), and "corroboration" framing (M4). Numbers verified against
> `docs/METHODS.md`, `results/tables/network/*`, and `results/audit/*`. Word count ≈ 1,020.

---

## Methods

### Proteomic data and differential abundance

Data-independent acquisition (DIA) proteomics was performed on ten paired tumour and
adjacent-normal specimens from patients with head and neck squamous cell carcinoma (HNSCC; n = 20
samples; six HPV-positive and four HPV-negative pairs). Spectra were processed with MaxQuant, and
differential abundance was assessed with proteoDA, a limma wrapper, for the primary tumour-versus-
adjacent-normal contrast (averaging over HPV status). Proteins were considered differentially
abundant at |log₂FC| > 1 and Benjamini–Hochberg-adjusted P < 0.05. UniProt accessions were mapped
to Entrez gene identifiers and HGNC symbols with `org.Hs.eg.db` (clusterProfiler `bitr`).

### Functional enrichment

Gene-set enrichment analysis (GSEA) was performed with clusterProfiler (v4.14.6), ranking proteins
by a π-statistic [sign(log₂FC) × |log₂FC| × −log₁₀(adjusted P)] against the MSigDB Hallmark, Gene
Ontology Biological Process, KEGG, and Reactome collections (minGSSize = 10, maxGSSize = 500,
adjusted P < 0.05; fixed seed for reproducible permutations). The Hallmark result is shown in
Figure 2C.

### Drug–target mapping and integration

Differentially abundant genes were queried against four pharmacological resources: DGIdb (GraphQL
API v5), ChEMBL (REST API v33, clinical phase ≥ 3), Open Targets (GraphQL API; HNSCC association
EFO_0000181 and known drugs), and L2S2 (LINCS L1000 Signature Search; bidirectional enrichment of
the top-150 up-regulated proteins against drug down-signatures and vice versa, FDA-approved filter,
P < 0.001; reversal score normalized to [−1, 0]). Sources were unified and each drug assigned an
ordinal regulatory class encoding regulatory maturity (not pharmacological equivalence): **A** =
approved with HNSCC evidence, **B** = approved in another cancer, **C** = approved
non-oncological, **D** = experimental (ordinal scores 1.00/0.75/0.50/0.25). Automated indication
matching was manually corrected for known mismatches: three agents were reassigned to class A
(cetuximab, pembrolizumab, nivolumab: approved in HNSCC but missed by EFO matching), two from C to
B (neratinib, lazertinib), and two from B to C (digoxin, digitoxin: cardiac glycosides with
oncology trials but no oncological approval); all overrides are documented in the analysis code.
The candidate set entering prioritization (n = 458) comprised drugs supported by ≥ 2 databases or
already approved (class A/B retained even if single-source), after removal of a curated exclusion
list.

### Protein–protein interaction network

A high-confidence interaction network was built from the differentially abundant proteins using
STRING v12 (combined score ≥ 700). Its giant component (521 proteins, 2,714 interactions) was
characterized by degree, betweenness, and eigenvector centrality (igraph v2.2.2); hubs were defined
as the top decile by degree. Louvain community detection partitioned the giant component into 13
modules, each named by its top-enriched GO Biological Process term. Druggable hubs (n = 99) were
hubs intersecting the drug-target set; module drugability tier was assigned data-driven (containing
an approved-drug target; druggable hubs only; or below threshold). Grey/below-threshold modules
still contained druggable hubs (36 of 99) and prioritized targets (e.g., DNMT1).

### Network-anchored two-factor prioritization

Each candidate was scored with a composite that anchors prioritization in network structure
(network → module → druggable hub → drug):

`composite = 0.60 · TargetPriority + 0.40 · DrugViability`

**TargetPriority** (target-level) = 0.55 · network centrality (mean of min–max-normalized degree,
betweenness, eigenvector) + 0.45 · directional differential abundance (π = sign(log₂FC)·|log₂FC|·
−log₁₀FDR, direction-aware min–max). This π term is signed and rewards inhibition of over-abundant
targets; consequently, under-abundant targets (including potential loss-of-function/tumour-
suppressor nodes) are structurally down-weighted by this term, and metabolic targets such as
Complex I subunits are prioritized through direction-independent centrality and external evidence
rather than through π. **DrugViability** (drug-level; breaks ties among drugs sharing a target) =
0.40 · L2S2 transcriptomic reversal + 0.40 · regulatory class + 0.20 · evidence breadth (databases
/ 4). Each drug was attributed to its most-central credible target, where a drug→target edge is
credible only if curated by ChEMBL or Open Targets (DGIdb-only edges were excluded as attribution
anchors because DGIdb interaction scores are promiscuous and do not reliably reflect mechanism of
action). A target that is a network hub defines the hub-central tier; a credible non-hub target the
peripheral-differential tier. Anchor genes lacking antitumoral plausibility or reflecting sample
contamination were excluded via a curated configuration list. The composite is strongly correlated
with network centrality alone (Spearman ρ ≈ 0.78); network topology is the dominant driver by
design, with DrugViability acting primarily as a clinical/pharmacological tie-breaker among agents
sharing a target.

### Robustness

Ranking stability was assessed by weight-perturbation sensitivity (six configurations varying the
target/drug balance and the sub-weights) and by leave-one-database (LOD) stability. No permutation
null was used: for an evidence-aggregation prioritization score (rather than a de novo discovery
statistic), validity rests on positive-control recovery (approved EGFR-axis agents), robustness to
weighting/database choices, and external corroboration, not on a permutation P-value.

### External corroboration in two independent cohorts

Target directionality was corroborated, not tested for efficacy or rank order, in two independent
cohorts. The CPTAC-HNSCC TMT proteome (Huang et al., 2021) was retrieved with the `cptac` Python
package and analysed by the same method as discovery (paired limma with duplicateCorrelation
blocking on patient; 116 tumour / 66 normal, 66 pairs); CPTAC log₂FC correlated with the DIA log₂FC
across 636 shared genes (Pearson r = 0.789, P = 4.9 × 10⁻¹³⁶; 86.3% directional concordance). The
TCGA-HNSC transcriptome (STAR-Counts via TCGAbiolinks; 520 tumour / 44 normal) was analysed with
DESeq2 (v1.44); TCGA log₂FC correlated with the DIA log₂FC across 663 shared genes (r = 0.601,
P = 2.7 × 10⁻⁶⁶; 76.2% concordance). For the 14 shortlist anchor targets, per-gene log₂FC and
adjusted P were compared against the DIA proteome in both cohorts (concordance and FDR < 0.05
annotated per target). As a supplementary analysis, overall survival was stratified by median
expression of four module-representative genes (EGFR, PSMB10, DNMT1, NDUFS3) in TCGA-HNSC (476
patients with survival data, log-rank); all comparisons were non-significant (P = 0.098–0.422),
consistent with these targets being therapeutic vulnerabilities rather than prognostic biomarkers.

### Software and reproducibility

Analyses used R 4.4.3 (limma v3.62.2, clusterProfiler v4.14.6, ReactomePA v1.50.0, igraph v2.2.2,
ggraph v2.2.2, ComplexHeatmap v2.22.0, DESeq2 v1.44) and Python 3 (requests, pandas, numpy,
matplotlib; `cptac`). External resources: DGIdb v5, ChEMBL v33, Open Targets, L2S2 LINCS L1000
(Evangelista et al., 2025), STRING v12, and MSigDB Hallmark. Analysis parameters are version-
controlled in `config/analysis_params.yaml`, figure styling is centralized, and the execution order
is documented in the project runbook.
