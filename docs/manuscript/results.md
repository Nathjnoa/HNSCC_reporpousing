# Results — draft (Step 2, WRITING_PLAN.md)

> Full-prose English draft for author review. Numbers verified against
> `results/tables/pub/*`, `results/tables/network/*`, `results/tables/pathway_enrichment/*`
> and `results/audit/*` (post-audit, commit `938d133`). Figure/table callouts follow the
> display-item scheme (Fig 1–6, Fig S1–S3; Table 1 main; Tables S1–S5). Word count ≈ 1,450.
> Language: "corroboration / replication of target directionality", not "validation" (M4);
> let the ranking speak (confirmation-bias guard).

---

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
Open Targets, L2S2) yielded 458 multi-source candidate drugs — agents supported by at least two
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
**hub-central tier** comprised candidates whose anchor is a network hub: metformin, anchored to the
Complex I subunit NDUFS2 within the oxidative-phosphorylation module (composite 0.591), and
proteasome-directed agents anchored to PSMA2/PSMB3. Notably, Complex I subunits were robustly
*under*-abundant in tumour (NDUFS2 log₂FC = −3.05); metformin was prioritized through the high
network centrality of its anchor, independent of fold-change direction, and is framed here as a
metabolic-vulnerability and combination-therapy candidate rather than as inhibition of an
over-expressed target. The **peripheral-differential tier** comprised targets that are
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
limit-of-detection filter (Fig S2). Candidates flagged as LOD-stable — including the EGFR-axis
agents and the leading hub-central candidates — retained their positions across configurations
(Table 1, robustness columns), defining a robust core panel. Consistent with their nomination as
therapeutic vulnerabilities rather than prognostic markers, none of the four pillar genes (EGFR,
PSMB10, DNMT1, NDUFS3) was significantly associated with overall survival in TCGA-HNSC (all
p > 0.05; Fig S3).
