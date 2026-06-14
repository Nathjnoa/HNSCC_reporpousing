# Table legends — HNSCC Drug Repurposing

> Step 1 of `WRITING_PLAN.md`. Self-contained, publication-style legends in English,
> anchored to the actual column structure of `results/tables/pub/*.tsv` (post-audit,
> commit `938d133`). Tables are **final** and not regenerated.
>
> **Display-item scheme (author decision, 2026-06-14):** a single **main table** = the
> prioritized repurposing shortlist (EGFR control block + Tier 1/Tier 2), assembled from the
> existing `Tab4` + `Tab5` `.tsv` (presentation reorganization, no re-analysis). Cohort/DE
> counts are reported **in-text**. All other tables move to **Supplementary**.

---

## Main table

## Table 1. Prioritized repurposing candidates for HNSCC.

Network-anchored prioritization of repurposing candidates, ordered by composite score
(`composite = 0.60·TargetPriority + 0.40·DrugViability`). The table is organized in two blocks.
The **EGFR control block** (top) lists the EGFR-axis agents that the unsupervised prioritization
recovers — including drugs with established approval in HNSCC (cetuximab, composite 0.73, #1) —
serving as the method's internal positive control. The **prioritized candidates block** lists the
leading non-EGFR candidates by tier: *hub-central* (anchor target is a network hub; e.g.,
metformin → NDUFS2, oxidative phosphorylation, composite 0.591; proteasome) and
*peripheral-differential* (target differentially abundant but not a hub; e.g., epigenetic
DNMT1/MAOA, antimetabolite IMPDH2). Columns: drug name, block/tier, functional module, anchor
hub/target, primary target subclass (where applicable), maximum clinical phase, approved in HNSCC
(yes/no), composite score, TargetPriority (TP), DrugViability (DV), number of supporting sources
(of 4), LOD-stable (yes/no), and weight-configuration robustness (n/6 configurations retaining the
candidate). Metformin's anchor (Complex I) is under-abundant in tumour; it is presented as a
metabolic-vulnerability / combination candidate, not as inhibition of an overexpressed target.
*Assembled from `results/tables/pub/main/Tab4_EGFR_validation.tsv` and
`Tab5_novel_candidates_by_module.tsv`.*

---

## Supplementary tables

## Table S1. Cohort and differential-abundance summary.

Summary of the proteomic dataset and differential-abundance analysis (also summarized in-text).
The DIA proteome (MaxQuant, 10 paired tumour/normal specimens) quantified 3,352 proteins, all
mapped to Entrez IDs (org.Hs.eg.db). Differential abundance was tested with a paired limma contrast
(tumour-vs-normal); thresholds |log₂FC| > 1 and FDR < 0.05 (Benjamini–Hochberg). 666 proteins
(19.9%) were differentially abundant — 329 over-abundant (49.4%) and 337 under-abundant (50.6%) in
tumour. HPV status: 6 HPV-positive patients (12 samples), 4 HPV-negative patients (8 samples).
Columns: parameter, value, note. *Source: `Tab1_resumen_DE.tsv`.*

## Table S2. Top differentially abundant proteins.

Most extreme differentially abundant proteins in tumour versus adjacent-normal tissue, ranked by
|log₂FC|. Columns: gene symbol, UniProt ID, log₂ fold change, FDR (limma / Benjamini–Hochberg),
direction (over-/under-abundant). *Source: `Tab2_top_proteinas.tsv`.*

## Table S3. Candidate provenance and multi-source support.

The leading repurposing candidates with their evidence provenance and clinical/regulatory context.
Columns: drug name, number of supporting sources (of 4), the sources (DGIdb / ChEMBL / Open Targets
/ L2S2), maximum clinical phase, HNSCC indication (yes/no), oncological indication (yes/no), ChEMBL
ID. *Source: `Tab3_top_candidatos.tsv`.*

## Table S4. External corroboration summary.

Cohort-level corroboration metrics for the two independent validation datasets: overlapping gene
count, Pearson correlation of log₂FC against the DIA proteome, and directional concordance for
CPTAC-HNSCC (proteome; r = 0.789, n = 636, 86.3%) and TCGA-HNSC (transcriptome; r = 0.601, n = 663,
76.2%), with per-anchor concordance counts. Columns: metric, value. *Source:
`Tab6_concordance_summary.tsv`.*

## Table S5. Extended candidate list by module.

Full list of non-EGFR anchor candidates across all modules (extended form of the Table 1
prioritized block). Columns: drug name, tier, functional module, anchor hub/target, maximum
clinical phase, composite score, TargetPriority (TP), DrugViability (DV), number of sources,
LOD-stable, weight-configuration robustness (n/6). *Source:
`supp/TabS1_extended_candidates_by_module.tsv`.*
