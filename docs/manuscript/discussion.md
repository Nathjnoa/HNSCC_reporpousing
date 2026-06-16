# Discussion (draft, Step 5 · WRITING_PLAN.md)

> Full-prose English draft for author review. ~1,200 words. Follows the STORYBOARD §4
> five-move spine. Numbers verified against `results/tables/pub/*` and `results/audit/*`
> (post-audit, commit `938d133`). Clinical/literature claims carry PubMed-verified inline
> citations (`[@key]`; keys in `references.bib`, criteria logged in `docs/CLINICAL_CRITERIA.md`).
> Language: "corroboration / replication of target directionality", not "validation" of efficacy
> (M4); metformin framed as metabolic vulnerability, never inhibition of an over-expressed target
> (B2); associations stated as associations, not causation; hedged language for unproven claims.
> American English; lowercase p; no em dash.

---

## Discussion

Effective systemic options for HNSCC remain limited, and developing new agents de novo is slow and
costly; drug repurposing offers a faster, lower-risk route to deployable therapy
[@pushpakom2019; @tanoli2025]. Starting from the proteomes of paired HNSCC tumors, this study
delivers a mechanistically organized shortlist of repurposing candidates that resolves into four
interpretable therapeutic axes: EGFR (cetuximab and other anti-EGFR agents), epigenetic
(decitabine, azacitidine, tranylcypromine), proteasome (carfilzomib), and metabolic (metformin).
Two observations argue for taking the shortlist seriously. First, without disease-specific tuning,
the analysis re-discovered the EGFR axis, the only molecular target with established approved
therapy in HNSCC [@burtness2019; @carosi2026]. Second, the direction of change in the prioritized
targets was reproduced across two independent cohorts. The prioritization that produced this
shortlist weighs where a target sits in the tumor protein network together with how deployable its
drug is, rather than fold-change magnitude alone; it is the means by which the candidates are
ranked, not the claim of the work.

The unsupervised recovery of cetuximab at the top of the ranking (composite 0.73), accompanied by
other approved EGFR-directed agents, served as an internal positive control: a method that
independently elevates the established standard target is more credible when it also nominates less
obvious ones. The recovery is consistent with independent HNSCC analyses that converge on EGFR by
different routes, including single-cell nominations of anti-EGFR agents [@li2024] and proteogenomic
stratification in which EGFR-high tumors respond to cetuximab in patient-derived xenografts
[@wu2025]. Beyond this control, the result invites viewing EGFR not as an isolated receptor but
as a central hub within a broader signaling network. EGFR is frequently overexpressed in HNSCC, yet
single-agent anti-EGFR benefit is limited, a limitation that has been attributed to intrinsic and
acquired resistance, including compensatory receptor tyrosine kinase activation, downstream pathway
alterations, and epithelial-mesenchymal transition [@carosi2026]. Viewing EGFR as a network node
rather than a solitary target is consistent with the rationale for the combination strategies (dual
EGFR blockade, pan-HER inhibitors, or EGFR inhibition plus immunomodulation) currently being
explored to overcome resistance [@carosi2026; @bhatia2023]. The same network logic that recovers
EGFR also organizes the remaining candidates into interpretable mechanistic axes.

Two non-canonical axes are the most novel nominations. The epigenetic axis, anchored by DNMT1
(decitabine, azacitidine) and MAOA (tranylcypromine), implicates chromatin and DNA-methylation
regulators as an integral part of the HNSCC tumor state rather than only a secondary
epiphenomenon. This axis carries a plausible translational rationale: DNA methyltransferase
inhibition can increase tumor immunogenicity and remodel the immune microenvironment, which
provides a mechanistic basis for combining hypomethylating agents with immune-checkpoint inhibitors
[@bear2025]. A recent phase 2 window trial reported that decitabine followed by pembrolizumab
increased tumor-infiltrating lymphocytes and PD-L1 expression and reduced myeloid-derived
suppressor cells before standard therapy [@bear2025], illustrating the prime-then-checkpoint
strategy that our epigenetic nominations would support, an approach now being explored clinically.
The proteasome axis, anchored by the 20S subunits PSMA2 and PSMB3 (carfilzomib), points to a
possible dependence on protein homeostasis: a high biosynthetic burden and accumulation of
misfolded proteins could render HNSCC cells reliant on proteasomal clearance, a vulnerability that
next-generation proteasome inhibitors such as carfilzomib and oprozomib have been reported to
exploit in HNSCC preclinical models, in part through induction of the unfolded protein response
[@zang2012]. Proteostasis targeting is clinically established in hematologic malignancies but
comparatively under-explored in solid tumors such as HNSCC [@zang2012], which positions this axis
as a hypothesis worth dedicated preclinical testing rather than an established therapy. A feature
that distinguishes both axes from EGFR is that their anchors are differentially abundant nodes whose
proposed therapeutic logic would be reprogramming or stress-exploitation rather than the simple
inhibition of an overexpressed driver.

The metabolic axis requires the most careful framing. Unlike EGFR, which is over-abundant and
addressed by direct inhibition, Complex I subunits were robustly under-abundant in tumors
(NDUFS2 log2FC = -3.05), a direction reproduced in both the CPTAC proteome and the TCGA
transcriptome and consistent with a glycolytic, Warburg-like shift. Metformin rose to the top of
the metabolic module (composite 0.59) through the high network centrality of its anchor, NDUFS2,
independent of fold-change direction. Its nomination therefore rests on a metabolic-vulnerability
rationale rather than on target overexpression: residual oxidative-phosphorylation dependence may
remain a targetable feature of cancer cells, and Complex I inhibition has been reported to exert
immunometabolic effects on the tumor microenvironment [@pujaltemartin2024]. In HNSCC specifically,
metformin use has been associated with improved outcomes and may act in part through immune
mechanisms, including increased CD8+ infiltration and enhanced natural-killer-cell cytotoxicity
[@curry2018; @crist2022]. Because single-agent metformin trials have been less effective than
observational data suggested [@pujaltemartin2024], and early-phase HNSCC studies indicate
tolerability in combination with chemoradiation [@kemnade2023], we position metformin as a
combination candidate (with chemoradiation or immunotherapy) rather than as a single-agent
inhibitor of an over-expressed target.

More broadly, this study sits within, but is distinct from, prior repurposing efforts in HNSCC.
Experimental high-throughput screens nominate candidates directly from drug response but depend on
dedicated cell or patient-derived-cell resources [@gu2022], whereas existing in silico efforts have
typically operated on transcriptomic or genomic data from public repositories and converged on a
single hub target and a single drug [@kumar2025], or have applied network-medicine frameworks to
pan-cancer interactome modules rather than to a disease-specific proteome [@cheng2019; @barabasi2011].
Our contribution is the combination of three features: a paired tumor/normal DIA proteome as the
substrate, closer to the functional effector layer than mRNA; a network-centrality and
drug-deployability composite that resolves into several mechanistic axes rather than one target; and
directional corroboration across two external cohorts. By design the composite is dominated by
network centrality (Spearman rho approximately 0.78 with centrality alone): centrality carries the
biological signal, while the drug-level DrugViability term determines which agent to advance among
compounds that share a prioritized target. We therefore present the score as a transparent
aggregation of evidence rather than as a method whose superiority over fold-change or database-count
ranking is established, a comparison that would require dedicated benchmarking
[@pushpakom2019; @tanoli2025; @flanary2023]. Several prioritized agents are approved for
non-oncologic indications (valproate, acetazolamide, disulfiram) and could in principle offer
short-horizon, accessible repurposing options. Because single-agent repurposing has frequently
failed in HNSCC, as in a randomized trial of pantoprazole added to systemic therapy [@noronha2023],
and because the field is shifting toward perioperative, immunotherapy-based regimens in resectable
disease [@uppaluri2025], the prioritized agents are best positioned as combination candidates:
those that prime or sensitize the tumor microenvironment could complement checkpoint blockade, as
illustrated by repurposed fenofibrate in HPV+ models [@oneill2022], although our data speak to
target biology rather than to patient selection.

Several limitations bound these conclusions. The discovery cohort is small (10 paired specimens);
cross-cohort corroboration mitigates but does not replace adequately powered confirmation. Tumor
and adjacent-normal specimens differ in cellular composition, and muscle and respiratory-mucosa
signals in the down-regulated proteome likely reflect tissue context rather than tumor-intrinsic
biology; hallmark terms such as Myogenesis and Adipogenesis should be read accordingly, although
sensitivity analyses excluding tissue-composition markers left the prioritized targets and the
oxidative-phosphorylation axis essentially unchanged. Data-independent acquisition with MinProb
imputation can mildly bias proteins absent from tumor tissue toward apparent under-abundance, an
effect that appears to be largely genuine biology and affected only a small minority of
under-abundant proteins. The differential-abundance component of the score is directional,
rewarding inhibition of over-abundant targets, so under-abundant or potential tumor-suppressor
nodes are structurally down-weighted and reach prioritization mainly through centrality, as for the
metabolic axis. Crucially, the external cohorts corroborate the direction of change in the
prioritized targets; they do not test drug efficacy or the candidate rank order, and four anchors
(PDE6D, CA2, CDA, RPS11) were weakly concordant or discordant and should be treated as
lower-confidence nominations.

This is an in silico prioritization, so the candidates are hypotheses that require functional
validation, and the axes suggest concrete next experiments. The epigenetic axis is most directly
testable through a prime-then-checkpoint sequence, decitabine or azacitidine followed by PD-1
blockade, in immunocompetent HNSCC models with readouts of tumor-infiltrating lymphocytes and PD-L1
[@bear2025]. The proteasome axis warrants dose-response testing of carfilzomib in HNSCC cell lines
and patient-derived xenografts with an unfolded-protein-response readout, the under-explored axis
with the least precedent in solid tumors [@zang2012]. The metabolic axis is the closest to the
clinic and is best advanced as metformin combined with chemoradiation or checkpoint blockade rather
than as monotherapy, building on early-phase tolerability data [@kemnade2023]. Within these bounds,
the study offers a credible, mechanism-organized set of repurposing candidates and a prioritized
agenda for preclinical and combination-trial follow-up in HNSCC.
