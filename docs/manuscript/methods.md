# Methods (draft, Step 3 · WRITING_PLAN.md)

> Condensed/translated English draft from `docs/METHODS.md` (production-ready). This revision
> strips outcome numbers that belong in Results (correlation coefficients, concordance
> percentages, candidate/hub/module counts, P-values) and keeps only procedure, parameters and
> thresholds. Study/cohort sizes are retained as data description. The directional π note (M1),
> regulatory-class overrides (m4) and "corroboration" framing (M4) are preserved as method
> rationale. Parameters verified against `docs/METHODS.md`, `config/analysis_params.yaml`.

---

## Methods

### Proteomic data and differential abundance

Data-independent acquisition (DIA) proteomics was performed on ten paired tumor and
adjacent-normal specimens from patients with head and neck squamous cell carcinoma (HNSCC; n = 20
samples; six HPV-positive and four HPV-negative pairs). Spectra were processed with MaxQuant, and
differential abundance was assessed with proteoDA, a limma wrapper, for the primary tumor-versus-
adjacent-normal contrast (averaging over HPV status). Proteins were considered differentially
abundant at |log₂FC| > 1 and Benjamini–Hochberg-adjusted P < 0.05. UniProt accessions were mapped
to Entrez gene identifiers and HGNC symbols with `org.Hs.eg.db` (clusterProfiler `bitr`). To check
that differential abundance was not driven by the muscle- and mucosa-rich composition of the
adjacent-normal tongue and larynx specimens, the enrichment and network analyses were repeated
after excluding a curated set of tissue-composition markers (skeletal/cardiac muscle and
airway/secretory proteins); this left the prioritized targets and the Hallmark results essentially
unchanged (Table S2).

### Functional enrichment

Gene-set enrichment analysis (GSEA) was performed with clusterProfiler (v4.14.6), ranking proteins
by a π-statistic [sign(log₂FC) × |log₂FC| × −log₁₀(adjusted P)] against the MSigDB Hallmark, Gene
Ontology Biological Process, KEGG, and Reactome collections (minGSSize = 10, maxGSSize = 500,
adjusted P < 0.05; fixed seed for reproducible permutations).

### Drug–target mapping and integration

Differentially abundant genes were queried against four pharmacological resources: DGIdb (GraphQL
API v5), ChEMBL (REST API v33, clinical phase ≥ 3), Open Targets (GraphQL API; HNSCC association
EFO_0000181 and known drugs), and L2S2 (LINCS L1000 Signature Search; bidirectional enrichment of
the top-150 up-regulated proteins against drug down-signatures and vice versa, FDA-approved filter,
P < 0.001; reversal score normalized to [−1, 0]). Sources were unified and each drug assigned an
ordinal regulatory class encoding regulatory maturity (not pharmacological equivalence): **A** =
approved with HNSCC evidence, **B** = approved in another cancer, **C** = approved
non-oncological, **D** = experimental (ordinal scores 1.00/0.75/0.50/0.25). Automated indication
matching was manually corrected for seven known mismatches (for example, cetuximab, pembrolizumab
and nivolumab were reassigned to class A as agents approved in HNSCC but missed by EFO matching);
the complete list of overrides is given in Table S6. Drugs entering prioritization were those
supported by ≥ 2 databases or already approved (class A/B retained even if single-source), after
removal of a curated exclusion list.

### Protein–protein interaction network

A high-confidence interaction network was built from the differentially abundant proteins using
STRING v12 (combined score ≥ 700). The largest connected component of the network was characterized
by degree, betweenness, and eigenvector centrality (igraph v2.2.2); hubs were defined as the top
decile by degree. Louvain community detection partitioned this component into modules, each named
by its top-enriched GO Biological Process term. Druggable hubs were the hubs intersecting the drug-target set. Each module
was assigned a data-driven drugability tier (containing an approved-drug target; containing
druggable hubs only; or below threshold).

### Network-anchored two-factor prioritization

Each candidate was scored with a composite that anchors prioritization in network structure
(network → module → druggable hub → drug):

`composite = 0.60 · TargetPriority + 0.40 · DrugViability`

**TargetPriority** (target-level) = 0.55 · network centrality (mean of min–max-normalized degree,
betweenness, eigenvector) + 0.45 · directional differential abundance (π = sign(log₂FC)·|log₂FC|·
−log₁₀FDR, direction-aware min–max). This π term is signed and rewards inhibition of over-abundant
targets; as a result, under-abundant targets (including potential loss-of-function/tumor-
suppressor nodes) are down-weighted by this term, and metabolic targets such as Complex I subunits
are prioritized through direction-independent centrality and external evidence rather than
through π. **DrugViability** (drug-level; breaks ties among drugs sharing a target) =
0.40 · L2S2 transcriptomic reversal + 0.40 · regulatory class + 0.20 · evidence breadth (databases
/ 4). Each drug was assigned to its most central credible target, where a drug→target edge is
credible only if curated by ChEMBL or Open Targets (DGIdb-only edges were excluded as attribution
anchors because DGIdb interaction scores are promiscuous and do not reliably reflect mechanism of
action). A target that is a network hub defines the hub-central tier; a credible non-hub target the
peripheral-differential tier. Anchor genes lacking antitumoral plausibility or reflecting sample
contamination were excluded via a curated configuration list. By design, network topology is the
dominant driver of the composite, with DrugViability acting primarily as a clinical/pharmacological
tie-breaker among agents sharing a target.

### Robustness

Ranking stability was assessed by weight-perturbation sensitivity (six configurations varying the
target/drug balance and the sub-weights) and by leave-one-database (LOD) stability. No permutation
null was used: for an evidence-aggregation prioritization score (rather than a de novo discovery
statistic), its validity comes from recovering positive controls (approved EGFR-axis agents), from
robustness to the weighting and database choices, and from external corroboration, not from a
permutation P-value.

### External corroboration in two independent cohorts

Target directionality was corroborated, not tested for efficacy or rank order, in two independent
cohorts. The CPTAC-HNSCC TMT proteome (Huang et al., 2021) was retrieved with the `cptac` Python
package and analyzed by the same method as discovery (paired limma with duplicateCorrelation
blocking on patient; 116 tumor / 66 normal, 66 pairs). The TCGA-HNSC transcriptome (STAR-Counts
via TCGAbiolinks; 520 tumor / 44 normal) was analyzed with DESeq2 (v1.44). For each cohort, the
per-gene log₂FC was correlated with the DIA log₂FC over the set of shared genes (Pearson
correlation and directional-concordance fraction), and, for the shortlist anchor targets, per-gene
log₂FC and adjusted P were compared against the DIA proteome (directional concordance and FDR < 0.05
annotated per target). As a supplementary analysis, overall survival was stratified by median
expression of four module-representative genes (EGFR, PSMB10, DNMT1, NDUFS3) in TCGA-HNSC (476
patients with survival data) and compared by the log-rank test.

### Software and reproducibility

Analyses used R 4.4.3 (limma v3.62.2, clusterProfiler v4.14.6, ReactomePA v1.50.0, igraph v2.2.2,
ggraph v2.2.2, ComplexHeatmap v2.22.0, DESeq2 v1.44) and Python 3 (requests, pandas, numpy,
matplotlib; `cptac`). External resources: DGIdb v5, ChEMBL v33, Open Targets, L2S2 LINCS L1000
(Evangelista et al., 2025), STRING v12, and MSigDB Hallmark. A Data and Code Availability
statement will be provided at submission.
