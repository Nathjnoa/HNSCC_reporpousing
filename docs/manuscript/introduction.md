# Introduction — draft (Step 4, WRITING_PLAN.md)

> Condensed/translated from the thesis Introduction (`docs/thesis/ARTICULO_2.0_REVJC.md`,
> ≈L201–314). ~720 words. Clinical/literature claims carry placeholder citation keys
> `[@key]` and are flagged **⚠️ TO-VERIFY** — every clinical claim is verified against PubMed
> when references are compiled (Step 8; rule `06-clinical-evidence-rigor`). The repurposing
> gap is framed as motivation, not asserted (author guidance); method kept subordinate.

---

## Introduction

Head and neck squamous cell carcinoma (HNSCC) is among the leading causes of cancer morbidity and
mortality worldwide, and overall survival remains limited — particularly in locally advanced and
recurrent disease. ⚠️ TO-VERIFY `[@hnscc_epidemiology]` Despite advances across surgery,
radiotherapy, chemotherapy, and immunotherapy, clinical outcomes are still suboptimal,
reflecting the biological complexity of the disease. ⚠️ TO-VERIFY `[@hnscc_outcomes]`

Over the past decade, management has shifted with the introduction of immune-checkpoint
inhibitors. The phase 3 KEYNOTE-048 trial established pembrolizumab — as monotherapy in
PD-L1-positive patients or combined with platinum/5-fluorouracil — as the first-line standard for
recurrent or metastatic disease, with patient selection guided by the combined positive score
(CPS). ⚠️ TO-VERIFY `[@keynote048]` More recently, KEYNOTE-689 reported that perioperative
pembrolizumab added to standard of care significantly improves event-free survival in
locally advanced, resectable disease — a setting where progress had been comparatively slow.
⚠️ TO-VERIFY `[@keynote689]` These milestones underscore the need to identify new molecular
signatures and therapeutic targets to optimize systemic and perioperative strategies.

A central challenge in HNSCC is its molecular heterogeneity. The distinction between
human-papillomavirus (HPV)-associated tumours and those linked to tobacco and alcohol has defined
subgroups with substantial prognostic and biological differences, yet this stratification remains
insufficient to capture the functional complexity of the tumour. ⚠️ TO-VERIFY `[@hpv_subtypes]`
In this context, proteomics provides a direct readout of the functional tumour phenotype:
unlike genomic or transcriptomic profiles, the proteome reflects cellular activity integrating
post-transcriptional and post-translational regulation. ⚠️ TO-VERIFY `[@proteomics_phenotype]`
Proteomic studies have implicated EGFR signalling, energy metabolism, epigenetic regulation, and
the tumour microenvironment as recurrently altered programmes in HNSCC. ⚠️ TO-VERIFY
`[@hnscc_proteomics]`

Drug repurposing — redeploying approved compounds with established safety profiles — is an
efficient route to accelerate the identification of new therapeutic options, and integrating
tumour molecular signatures with pharmacological-connectivity resources can nominate agents
predicted to reverse a disease state. ⚠️ TO-VERIFY `[@repurposing_rationale]` For such a
candidate list to be useful to a translational reader, however, it would ideally be both
*credible* — distinguishing targets that occupy a central position in the tumour's molecular
network from those that are merely strongly altered — and *actionable*, accounting for whether the
nominating drug is realistically deployable (approval status, regulatory class, evidence of
transcriptomic reversal) and whether its target directionality holds beyond a single cohort.
*(Any comparative claim regarding existing pipelines is kept tentative pending citation.)*

Here we address this by integrating quantitative proteomics of paired HNSCC tumours with
multi-source pharmacological evidence and protein–protein interaction network structure to derive
a credible, mechanistically organized shortlist of repurposing candidates. Candidates are
prioritized not by fold-change magnitude alone but by the biological position of their target in
the tumour network together with the deployability of the drug, and the resulting ranking is
evaluated against a built-in internal control — recovery of the EGFR axis, the principal target
with established approved therapy in HNSCC — and corroborated for target directionality in two
independent cohorts (a CPTAC proteome and a TCGA transcriptome). The aim is therefore not only to
recover validated targets but also to surface novel, mechanistically interpretable vulnerabilities
with immediate translational potential.
