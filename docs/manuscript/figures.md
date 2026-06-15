# Figure legends: HNSCC Drug Repurposing

> Step 1 of `WRITING_PLAN.md`. Self-contained, publication-style legends in English,
> anchored to post-audit numbers (commit `938d133`). Figures are **final** and not
> regenerated. Panel letters follow the assembled multipanels in `results/figures/pub/`.
> Numbers must match `results/tables/pub/*` and `results/audit/*`.

---

## Figure 1. Study design and analytical workflow.

Overview of the drug-repurposing pipeline for head and neck squamous cell carcinoma (HNSCC).
Ten paired tumor / adjacent-normal specimens were profiled by data-independent acquisition
(DIA) mass spectrometry (MaxQuant), yielding 3,352 quantified proteins. Differential abundance
(limma, paired contrast; |log₂FC| > 1, FDR < 0.05) defined the tumor proteomic signature, which
was mapped to candidate drugs across four resources (DGIdb, ChEMBL, Open Targets, L2S2). Targets
were placed in a tumor protein–protein interaction (PPI) network (STRING), partitioned into
Louvain modules, and each candidate was scored by a two-factor composite
(`0.60·TargetPriority + 0.40·DrugViability`) anchoring every drug to its most credible central
target. Robustness was assessed across weight configurations and a limit-of-detection (LOD)
filter, and the directionality of prioritized targets was corroborated in two independent
cohorts (CPTAC-HNSCC proteome; TCGA-HNSC transcriptome).

## Figure 2. The HNSCC tumor proteome is broadly and coherently remodeled.

(**A**) Volcano plot of differential protein abundance in tumor versus adjacent-normal tissue
(limma paired contrast; n = 10 pairs). Of 3,352 quantified proteins, 666 were differentially
abundant (FDR < 0.05, |log₂FC| > 1): 329 over- and 337 under-abundant in tumor. Labels mark
proteins with |log₂FC| > 3.5 and FDR < 0.01 plus a fixed set of mechanistically relevant targets
(e.g., EGFR, proteasome and Complex I subunits, DNMT1); labeling is illustrative and does not
affect the analysis. (**B**) Heatmap of the 20 most over- and 20 most under-abundant proteins
detected in all 20 samples (100% completeness, no imputation); rows are proteins, columns are
samples, color shows z-scored abundance. (**C**) Hallmark gene-set enrichment analysis (GSEA,
MSigDB Hallmark, all 50 sets tested); the gene sets reaching FDR < 0.05 are shown, ordered by the
π-statistic [sign(log₂FC) × |log₂FC| × −log₁₀(FDR)]. Bars encode normalized enrichment score
(NES); the oxidative-phosphorylation signature is among the strongest down-regulated programs
(NES ≈ −2.2).

## Figure 3. A reproducible, multi-source druggable candidate space.

Characterization of the 458 multi-source candidate drugs: agents with support in ≥ 2 databases
**or** already approved (regulatory class A/B, retained even with a single source), minus a
curated exclusion list. All three panels describe this same set. (**A**) Distribution of maximum
clinical phase. (**B**) Distribution of ordinal regulatory class (A = approved with HNSCC
evidence; B = approved in another cancer; C = approved non-oncological; D = experimental). (**C**)
UpSet plot of candidate overlap across the four sources (DGIdb, ChEMBL, Open Targets, L2S2),
showing the intersection structure that defines multi-source support.

## Figure 4. Network modules organize targets into druggable mechanistic axes.

(**A**) Tumor PPI network (STRING) arranged with a stress layout that emphasizes intra-module
edges, colored by Louvain module and **data-driven drugability tier**: *approved* modules (≥ 1
approved drug targeting a module node; 11 modules, individual colors), *hubs-only* modules (0
approved but ≥ 2 druggable hubs; slate), and *below-threshold* modules (gray). Gray modules still
contain druggable hubs (36/99 of the total) and prioritized targets (e.g., DNMT1 → decitabine in
the chromatin-remodeling module); color marks the mechanistic axes, not exclusivity of targets. (**B**) Representative enriched
Gene Ontology term per module (one term/module). (**C**) Count of druggable hubs per module. Node
centrality in this network is an input to the composite score (Fig. 5).

## Figure 5. Two-factor prioritization ranks candidates and decomposes the score.

(**A**) Prioritized shortlist faceted by tier, with each candidate's composite score
decomposed into its TargetPriority (TP) and DrugViability (DV) contributions. *Composite* =
0.60·TP + 0.40·DV. *TargetPriority* = 0.55·network centrality (degree/betweenness/eigenvector,
min–max) + 0.45·directional differential abundance (π = sign(log₂FC)·|log₂FC|·−log₁₀FDR,
direction-aware min–max, rewarding inhibition of up-regulated targets). *DrugViability* =
0.40·transcriptomic reversal (L2S2) + 0.40·regulatory class (A–D) + 0.20·number of sources / 4.
Each drug is anchored to its most credible central target via a curated ChEMBL/Open Targets edge
(DGIdb excluded for edge anchoring because of noise); *tier* = hub-central (target is a network hub)
versus peripheral-differential. (**B**) Bifactor space of TargetPriority × DrugViability; dashed
lines are iso-contours of constant composite score. No permutation test is reported: for an
evidence-aggregation score, validity derives from recovery of approved-drug controls, robustness,
and external corroboration rather than from a permutation p-value.

## Figure 6. Target directionality replicates across two independent cohorts.

External corroboration in two independent HNSCC cohorts. (**A**) Proteome-versus-proteome
concordance: log₂FC from the DIA proteome (tumor vs adjacent-normal, n = 10 pairs) against
CPTAC-HNSCC TMT proteome (Huang et al., 2021; paired limma, 116 tumor / 66 normal, 66 pairs);
Pearson r = 0.789 over n = 636 shared genes, 86.3% directional concordance; four pillar genes
labeled. (**B**) Proteome-versus-transcriptome concordance: DIA versus TCGA-HNSC RNA-seq
(STAR-Counts, DESeq2; 520 tumor / 44 normal); r = 0.601 over n = 663 shared genes, 76.2%
directional concordance; four pillar genes labeled. Grayscale in A/B marks background
directional concordance/discordance. (**C**) The 14 shortlist anchor targets with their log₂FC in
CPTAC and TCGA (filled symbol = concordant with the DIA proteome; asterisks = FDR < 0.05 / 0.01 /
0.001) next to a composite-score bar split into TargetPriority and DrugViability, ordered
by descending composite. CPTAC: 12/14 with data (11 concordant, 10 FDR < 0.05); TCGA: 14/14 (11
concordant, 11 FDR < 0.05). Color marks the therapeutic pillar (EGFR / Proteasome / Epigenetic /
OXPHOS). Cohorts corroborate target **directionality**; they do not test drug efficacy or the rank
order of candidates.

---

## Supplementary figures

## Figure S1. Ranking robustness.

Heatmap of candidate ranking stability across six weight configurations plus the limit-of-detection
(LOD) filter. Candidates flagged `lod_stable` retain their position across configurations,
defining the robust panel.

