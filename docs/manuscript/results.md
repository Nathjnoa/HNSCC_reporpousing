# Results (draft, Step 2 · WRITING_PLAN.md)

> Full-prose English draft for author review. Numbers verified against
> `results/tables/pub/*`, `results/tables/network/*`, `results/tables/pathway_enrichment/*`
> and `results/audit/*` (post-audit, commit `938d133`). Figure/table callouts follow the
> display-item scheme (Fig 1–6, Fig S1; Table 1 main; Tables S1–S6). Word count ≈ 1,450.
> Language: "corroboration / replication of target directionality", not "validation" (M4);
> let the ranking speak (confirmation-bias guard).

---

## Results

### Differential protein abundance in the HNSCC proteome

Data-independent acquisition proteomics of ten paired tumor and adjacent-normal specimens
quantified 3,352 proteins, all mapped to Entrez gene identifiers. A paired differential-abundance
analysis (limma; |log₂FC| > 1, FDR < 0.05) identified 666 differentially abundant proteins (19.9%
of the quantified proteome), comprising 329 over-abundant and 337 under-abundant proteins in
tumor relative to adjacent-normal tissue (Fig 2A; Table S1). Gene-set enrichment analysis against
the 50 MSigDB Hallmark gene sets returned 11 significantly enriched programs (FDR < 0.05; Fig 2C).
The most strongly down-regulated signature was oxidative phosphorylation (NES = −2.21), followed by
myogenesis (NES = −1.99) and adipogenesis (NES = −1.83).

The most extreme under-abundant proteins were dominated by skeletal/cardiac muscle (e.g., MYH2,
MYL1/2/3, CASQ1, TCAP, ATP2A1) and airway/secretory markers (BPIFA1, BPIFB1). These
tissue-composition markers accounted for 12.8% of down-regulated proteins overall but 80% of the
20 most extreme. In a sensitivity analysis excluding these markers (Table S2), the hub
status of all prioritized targets was unchanged (e.g., NDUFS2 degree 51→51; EGFR 24→23), the
oxidative-phosphorylation signature was retained and even slightly stronger (NES −2.21→−2.32), and
only the myogenesis hallmark weakened appreciably (NES −2.00→−1.56), confirming that our
conclusions about the prioritized targets and the OXPHOS axis hold up regardless of tissue composition.

### Multi-source druggable candidate space

Mapping the differentially abundant proteins to four pharmacological resources (DGIdb, ChEMBL,
Open Targets, L2S2) yielded 458 multi-source candidate drugs: agents supported by at least two
databases, or already approved (regulatory class A/B, retained even with a single source), after
a curated exclusion list (Fig 3). The candidate set spanned the full clinical-development
spectrum, with a substantial fraction of approved agents (Fig 3A–B), and the cross-source UpSet
plot confirmed that the candidates were not an artifact of any single database (Fig 3C).

### Protein–protein interaction network and module structure

To prioritize by biological position rather than fold-change magnitude alone, we placed the
candidate targets in a high-confidence tumor protein–protein interaction network (STRING). Of the
differentially abundant proteins, 521 formed the network's largest connected component (2,714
interactions), which Louvain community detection partitioned into 13 functional modules (Fig 4A)
carrying 99 druggable hubs. Module-level Gene Ontology
enrichment assigned a coherent functional identity to each module (Fig 4B), and the per-module
hub counts (Fig 4C) defined the candidate anchors carried forward to prioritization. Modules were
classified by a data-driven drugability tier (containing an approved-drug target, druggable hubs
only, or below threshold). The most prominent druggable axes were anchored by their key hubs:
EGFR (module M2), proteasome subunits (M4, ubiquitin-dependent catabolism), Complex I of
oxidative phosphorylation (M8), and epigenetic regulators (M17, chromatin remodeling).

### Network-anchored two-factor prioritization

Each candidate was scored by a composite that combines how important its target is biologically
with how ready the drug is for clinical use (`composite = 0.60·TargetPriority + 0.40·DrugViability`),
with every drug anchored to its most credible central target through a curated ChEMBL/Open Targets
edge (Fig 5; Table 1). TargetPriority combines network centrality with directional differential
abundance, and DrugViability combines transcriptomic reversal (L2S2), regulatory class, and how
many sources support the drug.

Without any disease-specific tuning, the prioritization elevated the EGFR axis (module M2) to the
top of the ranking: cetuximab ranked first (composite 0.73), accompanied by other approved
EGFR-directed agents (nimotuzumab, panitumumab, gefitinib, afatinib) (Table 1, EGFR control
block). Recovering these clinically validated EGFR-directed agents without supervision served as
an internal positive control for the method.

Beyond this control, the shortlist split into two interpretable tiers (Table 1; Fig 5A). The
**hub-central tier** contained candidates whose anchor is a network hub: metformin, anchored to the
Complex I subunit NDUFS2 within the oxidative-phosphorylation module (composite 0.591), and
proteasome-directed agents anchored to PSMA2/PSMB3. Notably, Complex I subunits were strongly
*under*-abundant in tumor (NDUFS2 log₂FC = −3.05); metformin rose to the top through the high
network centrality of its anchor, independent of fold-change direction, and we present it as a
metabolic-vulnerability and combination-therapy candidate rather than as inhibition of an
over-expressed target. The **peripheral-differential tier** held targets that are
differentially abundant but not network hubs, including the epigenetic regulators DNMT1
(decitabine/azacitidine) and MAOA (tranylcypromine) and the antimetabolite target IMPDH2
(thioguanine), alongside classical repurposing agents. The extended candidate list across all
modules is provided in Table S5, and candidate provenance in Table S3.

### External corroboration in two independent cohorts

To see whether the findings generalize, we asked whether the direction of change in the
prioritized targets reproduced in independent data (Fig 6). At the global level, the DIA proteome correlated with the
CPTAC-HNSCC TMT proteome (Pearson r = 0.789 over 636 shared genes; 86.3% directional concordance;
Fig 6A) and with the TCGA-HNSC transcriptome (r = 0.601 over 663 shared genes; 76.2% concordance;
Fig 6B). For the 14 shortlist anchor targets specifically (Fig 6C; Table S4), directional
concordance was observed for 11 of 12 anchors with CPTAC data and 11 of 14 with TCGA, with 10/12
and 12/14 reaching FDR < 0.05, respectively. Concordance was strong for the principal candidates
(EGFR, NDUFS2, DNMT1, PSMA2, ALDH5A1, MAOA, MMP8, IMPDH2, ALDH2) and weaker or discordant for four
anchors (PDE6D, CA2, CDA, RPS11), which we treat as lower-confidence nominations. These cohorts
corroborate the direction of change in the prioritized targets; they are not a test of drug
efficacy or of the candidate rank order.

### Robustness of the ranking

Finally, the candidate ranking was stable across six independent weight configurations and a
limit-of-detection filter (Fig S1). Consistent with the network-anchored design, the composite
ranking was strongly correlated with target network centrality alone (Spearman ρ = 0.78),
confirming that network topology, rather than the drug-level tie-breakers, was the dominant driver
of prioritization. Candidates flagged as LOD-stable, including the EGFR-axis
agents and the leading hub-central candidates, retained their positions across configurations
(Table 1, robustness columns), defining a robust core panel.
