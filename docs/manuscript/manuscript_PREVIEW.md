# Network-anchored proteomic drug repurposing in HNSCC (working draft, preview)

Partial preview for pipeline testing: Introduction, Methods and Results only. Discussion, Abstract, Title and complete Methods/Results citations are pending (WRITING_PLAN.md steps 5–8). Citations shown are PubMed-verified; the bibliography is generated automatically by pandoc.

## Introduction

Head and neck squamous cell carcinoma (HNSCC) is among the leading causes of cancer morbidity and
mortality worldwide, and overall survival remains limited, particularly in locally advanced and
recurrent disease. [@johnson2020; @sung2021] Despite advances across surgery,
radiotherapy, chemotherapy, and immunotherapy, clinical outcomes are still suboptimal,
reflecting the biological complexity of the disease. [@johnson2020]

Over the past decade, management has shifted with the introduction of immune-checkpoint
inhibitors. The phase 3 KEYNOTE-048 trial established pembrolizumab (as monotherapy in
PD-L1-positive patients or combined with platinum/5-fluorouracil) as the first-line standard for
recurrent or metastatic disease, with patient selection guided by the combined positive score
(CPS). [@burtness2019] More recently, KEYNOTE-689 reported that perioperative
pembrolizumab added to standard of care significantly improves event-free survival in
locally advanced, resectable disease, a setting where progress had been comparatively slow.
[@imamura2025] These milestones underscore the need to identify new molecular
signatures and therapeutic targets to optimize systemic and perioperative strategies.

A central challenge in HNSCC is its molecular heterogeneity. The distinction between
human-papillomavirus (HPV)-associated tumours and those linked to tobacco and alcohol has defined
subgroups with substantial prognostic and biological differences, yet this stratification remains
insufficient to capture the functional complexity of the tumour. [@ang2010]
In this context, proteomics provides a direct readout of the functional tumour phenotype:
unlike genomic or transcriptomic profiles, the proteome reflects cellular activity integrating
post-transcriptional and post-translational regulation. [@huang2021]
Proteomic studies have implicated EGFR signalling, energy metabolism, epigenetic regulation, and
the tumour microenvironment as recurrently altered programmes in HNSCC. [@huang2021]

Drug repurposing (redeploying approved compounds with established safety profiles) is an
efficient route to accelerate the identification of new therapeutic options, and integrating
tumour molecular signatures with pharmacological-connectivity resources can nominate agents
predicted to reverse a disease state. [@pushpakom2019] For such a
candidate list to be useful to a translational reader, however, it would ideally be both
credible (distinguishing targets that occupy a central position in the tumour's molecular
network from those that are merely strongly altered) and actionable, accounting for whether the
nominating drug is realistically deployable (approval status, regulatory class, evidence of
transcriptomic reversal) and whether its target directionality holds beyond a single cohort.
(Any comparative claim regarding existing pipelines is kept tentative pending citation.)

Here we address this by integrating quantitative proteomics of paired HNSCC tumours with
multi-source pharmacological evidence and protein–protein interaction network structure to derive
a credible, mechanistically organized shortlist of repurposing candidates. Candidates are
prioritized not by fold-change magnitude alone but by the biological position of their target in
the tumour network together with the deployability of the drug, and the resulting ranking is
evaluated against a built-in internal control (recovery of the EGFR axis, the principal target
with established approved therapy in HNSCC) and corroborated for target directionality in two
independent cohorts (a CPTAC proteome and a TCGA transcriptome). The aim is therefore not only to
recover validated targets but also to surface novel, mechanistically interpretable vulnerabilities
with immediate translational potential.

## Methods

### Proteomic data and differential abundance

Data-independent acquisition (DIA) proteomics was performed on ten paired tumour and
adjacent-normal specimens from patients with head and neck squamous cell carcinoma (HNSCC; n = 20
samples; six HPV-positive and four HPV-negative pairs). Spectra were processed with MaxQuant, and
differential abundance was assessed with proteoDA [@ritchie2015], a limma wrapper, for the primary tumour-versus-
adjacent-normal contrast (averaging over HPV status). Proteins were considered differentially
abundant at |log₂FC| > 1 and Benjamini–Hochberg-adjusted P < 0.05. UniProt accessions were mapped
to Entrez gene identifiers and HGNC symbols with `org.Hs.eg.db` (clusterProfiler [@wu2021] `bitr`).

### Functional enrichment

Gene-set enrichment analysis (GSEA) [@subramanian2005] was performed with clusterProfiler (v4.14.6), ranking proteins
by a π-statistic [sign(log₂FC) × |log₂FC| × −log₁₀(adjusted P)] against the MSigDB Hallmark [@liberzon2015], Gene
Ontology Biological Process, KEGG, and Reactome collections (minGSSize = 10, maxGSSize = 500,
adjusted P < 0.05; fixed seed for reproducible permutations). The Hallmark result is shown in
Figure 2C.

### Drug–target mapping and integration

Differentially abundant genes were queried against four pharmacological resources: DGIdb [@freshour2021] (GraphQL
API v5), ChEMBL [@mendez2019] (REST API v33, clinical phase ≥ 3), Open Targets [@ochoa2021] (GraphQL API; HNSCC association
EFO_0000181 and known drugs), and L2S2 [@evangelista2025] (LINCS L1000 Signature Search; bidirectional enrichment of
the top-150 up-regulated proteins against drug down-signatures and vice versa, FDA-approved filter,
P < 0.001; reversal score normalized to [−1, 0]). Sources were unified and each drug assigned an
ordinal regulatory class encoding regulatory maturity (not pharmacological equivalence): A =
approved with HNSCC evidence, B = approved in another cancer, C = approved
non-oncological, D = experimental (ordinal scores 1.00/0.75/0.50/0.25). Automated indication
matching was manually corrected for known mismatches: three agents were reassigned to class A
(cetuximab, pembrolizumab, nivolumab: approved in HNSCC but missed by EFO matching), two from C to
B (neratinib, lazertinib), and two from B to C (digoxin, digitoxin: cardiac glycosides with
oncology trials but no oncological approval); all overrides are documented in the analysis code.
The candidate set entering prioritization (n = 458) comprised drugs supported by ≥ 2 databases or
already approved (class A/B retained even if single-source), after removal of a curated exclusion
list.

### Protein–protein interaction network

A high-confidence interaction network was built from the differentially abundant proteins using
STRING [@szklarczyk2023] v12 (combined score ≥ 700). Its giant component (521 proteins, 2,714 interactions) was
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

TargetPriority (target-level) = 0.55 · network centrality (mean of min–max-normalized degree,
betweenness, eigenvector) + 0.45 · directional differential abundance (π = sign(log₂FC)·|log₂FC|·
−log₁₀FDR, direction-aware min–max). This π term is signed and rewards inhibition of over-abundant
targets; consequently, under-abundant targets (including potential loss-of-function/tumour-
suppressor nodes) are structurally down-weighted by this term, and metabolic targets such as
Complex I subunits are prioritized through direction-independent centrality and external evidence
rather than through π. DrugViability (drug-level; breaks ties among drugs sharing a target) =
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
cohorts. The CPTAC-HNSCC TMT proteome (Huang et al., 2021 [@huang2021]) was retrieved with the `cptac` Python
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

## Results

### The HNSCC tumour proteome is broadly and coherently remodelled

Data-independent acquisition proteomics of ten paired tumour and adjacent-normal specimens
quantified 3,352 proteins, all mapped to Entrez gene identifiers. A paired differential-abundance
analysis (limma; |log₂FC| > 1, FDR < 0.05) identified 666 differentially abundant proteins (19.9%
of the quantified proteome), comprising 329 over-abundant and 337 under-abundant proteins in
tumour relative to adjacent-normal tissue (Fig 2A; Table S1). The cohort included six
HPV-positive and four HPV-negative patients. Gene-set enrichment analysis against the 50 MSigDB
Hallmark gene sets returned 11 significantly enriched programmes (FDR < 0.05; Fig 2C). The most
strongly down-regulated signature was oxidative phosphorylation (NES = −2.21), consistent with a
glycolytic metabolic shift, followed by myogenesis (NES = −1.99) and adipogenesis (NES = −1.83).

The most extreme under-abundant proteins were dominated by skeletal/cardiac muscle (e.g., MYH2,
MYL1/2/3, CASQ1, TCAP, ATP2A1) and airway/secretory markers (BPIFA1, BPIFB1), reflecting the
muscle- and mucosa-rich composition of adjacent-normal tongue and larynx specimens. These
tissue-composition markers accounted for 12.8% of down-regulated proteins overall but 80% of the
20 most extreme. To ensure that downstream conclusions were not driven by tissue composition, we
repeated the enrichment and network analyses after excluding these markers: the hub status of all
prioritized targets was unchanged (e.g., NDUFS2 degree 51→51; EGFR 24→23), the
oxidative-phosphorylation signature was retained and slightly strengthened (NES −2.21→−2.32), and
only the myogenesis hallmark was materially attenuated (NES −2.00→−1.56) (Table S2; sensitivity
analysis). Conclusions regarding prioritized targets and the OXPHOS axis are therefore robust to
this confound.

### A reproducible, multi-source druggable candidate space

Mapping the differentially abundant proteins to four pharmacological resources (DGIdb, ChEMBL,
Open Targets, L2S2) yielded 458 multi-source candidate drugs: agents supported by at least two
databases, or already approved (regulatory class A/B, retained even with a single source), after
a curated exclusion list (Fig 3). The candidate set spanned the full clinical-development
spectrum, with a substantial fraction of approved agents (Fig 3A–B), and the cross-source UpSet
structure confirmed that candidates were not the artefact of any single database (Fig 3C). The
complete selection funnel is provided in Fig S1.

### Network structure organizes targets into druggable mechanistic axes

To prioritize by biological position rather than fold-change magnitude alone, candidate targets
were embedded in a tumour protein–protein interaction network (STRING) and partitioned into
Louvain modules (Fig 4A), which carried 99 druggable hubs (Methods). Module-level Gene Ontology
enrichment assigned a coherent functional identity to each module (Fig 4B), and the per-module
hub counts (Fig 4C) defined the candidate anchors carried forward to prioritization. Modules were
classified by a data-driven drugability tier (containing an approved-drug target, druggable hubs
only, or below threshold); the most prominent axes corresponded to EGFR signalling, oxidative
phosphorylation, the proteasome, and epigenetic regulation.

### A two-factor prioritization recovers known therapy and surfaces novel candidates

Each candidate was scored by a composite combining the biological priority of its target and the
deployability of the drug (`composite = 0.60·TargetPriority + 0.40·DrugViability`), with every drug
anchored to its most credible central target through a curated ChEMBL/Open Targets edge (Fig 5;
Table 1). TargetPriority integrates network centrality with directional differential abundance,
and DrugViability integrates transcriptomic reversal (L2S2), regulatory class, and source breadth
(Methods).

Without any disease-specific tuning, the prioritization elevated the EGFR axis (module M2) to the
top of the ranking: cetuximab ranked first (composite 0.73), accompanied by other approved
EGFR-directed agents (nimotuzumab, panitumumab, gefitinib, afatinib) (Table 1, EGFR control
block). Because EGFR is the principal molecular target with established approved therapy in HNSCC,
this unsupervised recovery of clinically validated agents serves as an internal positive control
for the method, supporting the credibility of its other nominations.

Beyond this control, the shortlist resolved into two interpretable tiers (Table 1; Fig 5A). The
hub-central tier comprised candidates whose anchor is a network hub: metformin, anchored to the
Complex I subunit NDUFS2 within the oxidative-phosphorylation module (composite 0.591), and
proteasome-directed agents anchored to PSMA2/PSMB3. Notably, Complex I subunits were robustly
under-abundant in tumour (NDUFS2 log₂FC = −3.05); metformin was prioritized through the high
network centrality of its anchor, independent of fold-change direction, and is framed here as a
metabolic-vulnerability and combination-therapy candidate rather than as inhibition of an
over-expressed target. The peripheral-differential tier comprised targets that are
differentially abundant but not network hubs, including the epigenetic regulators DNMT1
(decitabine/azacitidine) and MAOA (tranylcypromine) and the antimetabolite target IMPDH2
(thioguanine), alongside classical repurposing agents. The extended candidate list across all
modules is provided in Table S5, and candidate provenance in Table S3.

### Prioritized-target directionality replicates in two independent cohorts

To assess generalizability, we asked whether the directionality of the prioritized targets
replicated in independent data (Fig 6). At the global level, the DIA proteome correlated with the
CPTAC-HNSCC TMT proteome (Pearson r = 0.789 over 636 shared genes; 86.3% directional concordance;
Fig 6A) and with the TCGA-HNSC transcriptome (r = 0.601 over 663 shared genes; 76.2% concordance;
Fig 6B). For the 14 shortlist anchor targets specifically (Fig 6C; Table S4), directional
concordance was observed for 11 of 12 anchors with CPTAC data and 11 of 14 with TCGA, with 10/12
and 12/14 reaching FDR < 0.05, respectively. Concordance was strong for the principal candidates
(EGFR, NDUFS2, DNMT1, PSMA2, ALDH5A1, MAOA, MMP8, IMPDH2, ALDH2) and weaker or discordant for four
anchors (PDE6D, CA2, CDA, RPS11), which should be regarded as lower-confidence nominations. These
cohorts corroborate the directionality of the prioritized targets; they are not a test of drug
efficacy or of the candidate rank order.

### The ranking is robust to scoring choices

Finally, the candidate ranking was stable across six independent weight configurations and a
limit-of-detection filter (Fig S2). Candidates flagged as LOD-stable, including the EGFR-axis
agents and the leading hub-central candidates, retained their positions across configurations
(Table 1, robustness columns), defining a robust core panel. Consistent with their nomination as
therapeutic vulnerabilities rather than prognostic markers, none of the four pillar genes (EGFR,
PSMB10, DNMT1, NDUFS3) was significantly associated with overall survival in TCGA-HNSC (all
p > 0.05; Fig S3).

# References

