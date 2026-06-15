# Storyboard / Spine — HNSCC Drug Repurposing Manuscript

> **Step 0 of `WRITING_PLAN.md`.** Defines the central claim, the figure/table → message
> map, and the section outline with word budget. **Checkpoint before prose.**
> All numbers below must match `results/tables/pub/*` and `results/audit/*` (post-audit,
> commit `938d133`). No figure/table is regenerated.

---

## 1. Central claim

### One-sentence version
*Proteomics of paired HNSCC tumours yields a credible, mechanistically organized shortlist of
repurposable drugs for HNSCC — anchored by recovery of the clinically validated EGFR axis and
corroborated across two independent cohorts — offering actionable therapeutic candidates beyond
the current standard of care.*

> **Framing emphasis (per author, clinical/translational journal family):** the **hero is the
> translational finding** — an actionable, mechanism-organized set of repurposing candidates for
> HNSCC. The two-factor network prioritization is the **vehicle that produces and justifies the
> shortlist**, described honestly but kept subordinate; not pitched as a methods contribution.

### The argument, in five moves

**(1) The clinical context.** HNSCC has high morbidity and few systemic options beyond
platinum/EGFR/PD-1; new agents are slow and costly to develop. Drug repurposing — redeploying
approved drugs — is an attractive, faster route. A repurposing candidate list is arguably most
useful to a clinician when it is *credible* and *actionable* — which motivated us to prioritize
not only by how much a target changes, but by where it sits biologically and whether its drug is
deployable. *(Premise framed as motivation, not an asserted gap; any comparative claim against
existing pipelines stays hedged until supported by citation.)*

**(2) How we build the shortlist.** From 10 paired HNSCC tumour / adjacent-normal DIA proteomes
we prioritize candidates by *both* the biological position of the target in the tumour protein
network *and* the deployability of the drug that hits it (approval status, regulatory class,
transcriptomic reversal). Concretely this is a two-factor composite
(`0.60·TargetPriority + 0.40·DrugViability`) over a *network → module → druggable hub → drug*
chain — but the point for the reader is the *result*: a shortlist that weighs well-positioned
targets and deployable drugs together, rather than fold-change magnitude alone.

**(3) Why the shortlist is credible — EGFR recovery.** Without disease-specific tuning, the
analysis independently elevates the EGFR axis to the top (cetuximab #1; afatinib, gefitinib and
other EGFR agents nearby) — i.e., it *re-discovers the one molecular target with established
approved therapy in HNSCC*. This built-in recovery of known clinical truth is the evidence that
the **other** nominations are worth a clinician's attention.

**(4) The actionable candidates, organized by mechanism.** Beyond EGFR, the shortlist resolves
into clinically interpretable therapeutic axes: a **hub-central Tier 1** (OXPHOS/Complex I →
metformin; proteasome → carfilzomib) and a **peripheral-differential Tier 2** (epigenetic
DNMT1/MAOA, antimetabolite IMPDH2, and classical repurposing agents). The grouping is
data-driven, giving each candidate a mechanistic rationale a clinician/translational reader can
evaluate.

**(5) Independent corroboration.** The directionality of the prioritized targets replicates in
two independent cohorts — CPTAC-HNSCC proteome (r = 0.789, 86.3% concordance) and TCGA-HNSC
transcriptome (r = 0.601, 76.2%) — strengthening confidence in the shortlist. This
**corroborates target directionality; it is not a test of drug efficacy or of the rank order.**

### The contribution
The deliverable is an **actionable, mechanism-organized shortlist of repurposing candidates for
HNSCC**, made credible by EGFR recovery and cross-cohort corroboration. The network
prioritization is the means to that end, not the headline. Metformin is one Tier-1 candidate
among several, framed as a metabolic vulnerability and combination-therapy option (its anchor,
Complex I, is *under*-abundant), not as inhibition of an overexpressed target.

### What this paper is NOT
Not an efficacy claim for any drug; not a metformin paper; not a prognostic-biomarker study. The
prioritized targets are framed as therapeutic vulnerabilities, not survival markers. (The earlier
Kaplan–Meier survival analysis of the pillar genes was dropped from the manuscript: a non-
significant null that added little to the narrative.)

---

## 2. Figure → message map

| Figure | One-line message it carries | Anchored numbers |
|--------|------------------------------|------------------|
| **Fig1** | Study design: 10 paired DIA proteomes → DGIdb/ChEMBL/OpenTargets/L2S2 → PPI network → two-factor composite → 2-cohort corroboration. | 10 pairs; 4 drug DBs |
| **Fig2** | The tumour proteome is broadly remodelled with coherent pathway-level shifts. | 666 DE proteins (329↑/337↓, 19.9% of 3352); 11/50 Hallmarks FDR<0.05; OXPHOS NES≈−2.2 |
| **Fig3** | A reproducible, multi-source druggable space exists (not noise from one DB). | 458 multi-source candidates (≥2 DBs OR approved); clinical-phase & regulatory-class distributions; UpSet overlap |
| **Fig4** | Network structure organizes targets into modules with druggable hubs — the scaffold for prioritization. | Louvain modules; drugability tiers (approved / hubs-only / below-threshold); 99 druggable hubs |
| **Fig5** | The two-factor composite ranks candidates and decomposes *why* each ranks where it does. | composite = 0.60·TargetPriority + 0.40·DrugViability; tier = hub_central vs peripheral_diff; cetuximab 0.73 #1 |
| **Fig6** | Prioritized-target directionality replicates across an independent proteome and transcriptome. | CPTAC r=0.789/n=636/86.3%; TCGA r=0.601/n=663/76.2%; 14 anchors (CPTAC 11/12, TCGA 11/14 concordant) |
| **Fig S1** (robust) | The ranking is stable across 6 weight configs + LOD. | lod_stable panel |

## 3. Table → message map

**Display-item scheme (author decision 2026-06-14):** ONE main table = the prioritized shortlist
(EGFR control block + Tier1/Tier2, assembled from Tab4+Tab5). Cohort/DE counts in-text. Everything
else → Supplementary. Full legends in `tables.md`.

| Table | Placement | Message |
|-------|-----------|---------|
| **Table 1** (= Tab4+Tab5 merged) | **Main** | Prioritized repurposing shortlist: EGFR control block + hub-central/peripheral tiers (composite/TP/DV/phase/LOD/robustness). The deliverable. |
| **Table S1** (Tab1) | Suppl | Cohort + DE summary (3352 quantified, 666 DE, HPV+/− split) — also in-text. |
| **Table S2** (Tab2) | Suppl | Top DE proteins (most extreme up/down). |
| **Table S3** (Tab3) | Suppl | Candidate provenance / multi-source support (sources/phase/regulatory). |
| **Table S4** (Tab6) | Suppl | External corroboration summary (CPTAC + TCGA per-anchor). |
| **Table S5** (TabS1) | Suppl | Extended non-EGFR anchor list. |

---

## 4. Section outline + word budget (target body 3500–4500 w; abstract ≤250)

| Section | Budget | Spine content | Source files |
|---------|--------|---------------|--------------|
| **Abstract** | ~250 | Background → method (network-anchored 2-factor) → EGFR control → tiers → 2-cohort corroboration → contribution. Distilled last. | — |
| **Introduction** | 600–800 | HNSCC burden & limited systemic options → repurposing as a faster route to actionable therapy → motivation (hedged, not asserted as established gap): a clinically useful candidate list would ideally combine biological prioritization + drug deployability + independent support; any comparative claim vs existing pipelines kept tentative / citation-gated → objective: a credible, mechanism-organized shortlist of repurposing candidates for HNSCC, with EGFR as internal control + cross-cohort corroboration. (Method framed as means, not headline.) | thesis ARTICULO_2.0 L201–314 |
| **Methods** | 900–1100 | Proteomics & DE (limma, FDR<0.05, |log2FC|>1) → multi-source drug mapping (458, class A–D, manual overrides m4) → PPI network & Louvain → two-factor composite (M1 directional π bias) → robustness (LOD/weights) → external cohorts (CPTAC limma paired; TCGA DESeq2) → software. Language: "corroboration" not "validation" (M4). | METHODS.md + AUDIT_TEXT m4/M1 |
| **Results** | 1300–1600 | R1 Proteome remodelling (Fig2; DE counts in-text; B1 tissue-composition note). R2 Druggable space (Fig3). R3 Network modules & hubs (Fig4). R4 Two-factor ranking: EGFR control + Tier1/Tier2 (Fig5, **Table 1**). R5 External corroboration (Fig6; Table S4; M4 language). R6 Robustness (FigS). | pub tables + README notes + AUDIT B1/M4 |
| **Discussion** | 800–1000 | Lead with translational implications of the four therapeutic axes (EGFR/Proteasome/Epigenetic/OXPHOS) → metformin as metabolic vulnerability + combination candidate (B2, NOT inhibit-overexpression) → what credible prioritization may add over fold-change/DB-count lists (M6, brief, hedged — supports clinical actionability; any superiority claim citation-gated, not pitched as method) → limitations (n=10, tissue composition B1, imputation m2, π directional bias M1, corroboration ≠ efficacy) → perspectives (which candidates merit preclinical/trial follow-up). | thesis L701–1158 + AUDIT B1/B2/M1/M6 |
| **Title + keywords + highlights** | — | Reflect method + EGFR control + 2-tier shortlist; no single-drug headline. | — |

**Total body ≈ 4050–4750 w** → trim Results/Methods to land ≤4500.

---

## 5. Non-negotiable guardrails (carried from WRITING_PLAN)

- Numbers must match post-audit tables; if a number doesn't reconcile, **report, don't re-run**.
- Clinical/literature claims get inline PubMed-verified citations; log new criteria in
  `docs/CLINICAL_CRITERIA.md`.
- Metformin = metabolic-vulnerability + combination framing (B2); never "inhibit overexpression";
  never the sole headline.
- "Corroboration/replication of target directionality", never "validation" of efficacy (M4).
- Let the ranking speak; exclusions only by uniform criteria (confirmation-bias guard).

---

## 6. Open question for author (resolve at checkpoint)

- **Target journal family confirmed?** README lists Cancers / BMC Cancer / Front Oncol / J Transl
  Med / Cancer Cell Int. Draft stays journal-agnostic (per plan) — but confirming the family now
  helps set Intro length & Discussion depth. Default: keep agnostic, decide at Phase 10.
