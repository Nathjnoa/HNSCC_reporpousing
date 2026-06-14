# Introduction (draft, Step 4 · WRITING_PLAN.md)

> Condensed/translated from the thesis Introduction (`docs/thesis/ARTICULO_2.0_REVJC.md`,
> ≈L201–314). ~720 words. Clinical/literature claims carry placeholder citation keys
> `[@key]` and are flagged **⚠️ TO-VERIFY**; every clinical claim is verified against PubMed
> when references are compiled (Step 8; rule `06-clinical-evidence-rigor`). The repurposing
> gap is framed as motivation, not asserted (author guidance); method kept subordinate.

---

## Introduction

Head and neck squamous cell carcinoma (HNSCC) is among the leading causes of cancer morbidity and
mortality worldwide, and overall survival remains limited, particularly in locally advanced and
recurrent disease. ⚠️ TO-VERIFY `[@hnscc_epidemiology]` Despite advances across surgery,
radiotherapy, chemotherapy, and immunotherapy, outcomes are still poor, a reflection of how
biologically complex the disease is. ⚠️ TO-VERIFY `[@hnscc_outcomes]`

Over the past decade, management has shifted with the introduction of immune-checkpoint
inhibitors. The phase 3 KEYNOTE-048 trial established pembrolizumab (as monotherapy in
PD-L1-positive patients or combined with platinum/5-fluorouracil) as the first-line standard for
recurrent or metastatic disease, with patient selection guided by the combined positive score
(CPS). ⚠️ TO-VERIFY `[@keynote048]` More recently, KEYNOTE-689 reported that perioperative
pembrolizumab added to standard of care significantly improves event-free survival in
locally advanced, resectable disease, a setting where progress had been comparatively slow.
⚠️ TO-VERIFY `[@keynote689]` These advances highlight the need to find new molecular signatures
and therapeutic targets to improve systemic and perioperative treatment.

A central challenge in HNSCC is its molecular heterogeneity. The distinction between
human-papillomavirus (HPV)-associated tumors and those linked to tobacco and alcohol has defined
subgroups with substantial prognostic and biological differences, yet this stratification remains
insufficient to capture the functional complexity of the tumor. ⚠️ TO-VERIFY `[@hpv_subtypes]`
Proteomics helps here because it gives a direct picture of the tumor's functional state:
unlike genomic or transcriptomic profiles, the proteome reflects what cells are actually doing,
capturing both post-transcriptional and post-translational regulation. ⚠️ TO-VERIFY
`[@proteomics_phenotype]` Proteomic studies have repeatedly found changes in EGFR signaling,
energy metabolism, and the tumor microenvironment in HNSCC. ⚠️ TO-VERIFY
`[@hnscc_proteomics]`

Drug repurposing, the reuse of approved drugs whose safety is already known, is an efficient way
to speed up the search for new treatment options, and combining a tumor's molecular signature
with pharmacological-connectivity resources can point to drugs predicted to reverse the disease
state. ⚠️ TO-VERIFY `[@repurposing_rationale]` Connectivity-based approaches typically prioritize
candidates by how strongly a drug reverses the transcriptomic signature ⚠️ TO-VERIFY
`[@lamb2006; @subramanian2017; @evangelista2025]`, while target-centric selection often relies on the size of
the protein or expression change. Ranking by signature reversal or effect size alone, however,
does not capture where a target sits in the tumor's molecular network, a property that network
medicine links to biological importance ⚠️ TO-VERIFY `[@barabasi2011]`, and it leaves aside
whether the candidate drug can realistically be used in the clinic (approval status, regulatory
class) or whether the change in its target reproduces beyond a single cohort. A shortlist that is
both believable and usable for a translational reader needs an approach that brings these pieces
together: where the target sits in the network, how ready the drug is for the clinic, and whether
the signal reproduces across cohorts.

Our aim was to build and prioritize a credible, mechanistically organized shortlist of
drug-repurposing candidates for HNSCC by combining quantitative proteomics of paired tumors with
pharmacological evidence from several sources and the structure of the protein–protein
interaction network. Candidates were ranked not by the size of the protein change alone but by
where each target sits in the tumor network and how ready its drug is for clinical use. To check
that the approach is trustworthy, we asked whether it recovers the EGFR axis (the main target with
approved therapy in HNSCC) as a built-in positive control, and we tested whether the direction of
change in the prioritized targets agrees across two independent cohorts (a CPTAC proteome and a
TCGA transcriptome). The study was designed both to recover known targets and to bring forward
new, biologically interpretable vulnerabilities with translational potential.
