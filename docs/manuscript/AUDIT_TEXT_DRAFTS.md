# Borradores de texto derivados de la auditoría

> Párrafos sugeridos para integrar a Methods / Results / Limitations. **No** están insertados en el
> manuscrito; el autor decide su ubicación y redacción final. Las citas de B2 están verificadas en
> PubMed (DOIs incluidos). Ver `docs/AUDIT_PREMANUSCRIPT.md` para el sustento.

---

## B1 — Composición tisular del normal adyacente (Limitations + Results)

**Results (nota de interpretación):**
> The most strongly under-abundant proteins in tumour versus adjacent-normal tissue were dominated by
> skeletal/cardiac muscle (e.g., MYH2, MYL1/2/3, CASQ1, TCAP, ATP2A1) and airway/secretory markers
> (BPIFA1, BPIFB1), reflecting the muscle- and mucosa-rich composition of adjacent-normal
> tongue/larynx specimens. Tissue-composition markers accounted for 12.8% of down-regulated proteins
> overall but 80% of the 20 most extreme. To ensure that downstream conclusions were not driven by
> tissue composition, we repeated GSEA and network module/centrality analyses after excluding these
> markers: hub status of all prioritized targets was unchanged (e.g., NDUFS2 degree 51→51, EGFR
> 24→23), the oxidative-phosphorylation signature was retained and slightly strengthened
> (NES −2.21→−2.32), and only the Myogenesis hallmark was materially attenuated (NES −2.00→−1.56).

**Limitations:**
> Tumour and adjacent-normal specimens differ in cellular composition; muscle and respiratory-mucosa
> signals in the down-regulated proteome reflect tissue context rather than tumour-intrinsic biology,
> and hallmark terms such as Myogenesis/Adipogenesis should be interpreted accordingly. Sensitivity
> analyses (above) indicate that the prioritized targets and the OXPHOS axis are robust to this
> confound.

---

## B2 — Marco para dianas metabólicas down-reguladas (Results/Discussion sobre Metformina)

**Reencuadre obligatorio:** OXPHOS/Complejo I está **subexpresado** en el proteoma tumoral (NDUFS2
log2FC=−3.05; replicado en CPTAC y TCGA). No presentar como "sobreexpresado a inhibir".

> Unlike the EGFR axis, which is over-abundant and targeted by direct inhibition, Complex I subunits
> were robustly **under-abundant** in tumours (NDUFS2 log2FC = −3.05; replicated in CPTAC and
> TCGA-HNSC), consistent with a glycolytic (Warburg) metabolic shift. Metformin was prioritized here
> through the high network centrality of its anchor (NDUFS2), independent of fold-change direction,
> and its nomination rests on a **metabolic-vulnerability** rationale rather than target
> overexpression: residual OXPHOS dependence and the immunometabolic effects of Complex I inhibition,
> together with prior clinical precedent in HNSCC. Metformin is a mild Complex I inhibitor and OXPHOS
> constitutes a targetable vulnerability in cancer cells, although single-agent metformin trials have
> been less effective than observational data suggested (Pujalte-Martin et al., 2024, *Mol Oncol*,
> doi:10.1002/1878-0261.13583). In HNSCC specifically, metformin has been associated with improved
> outcomes and acts in part through tumour-microenvironment and immune mechanisms (NK-cell
> cytotoxicity, AMPK-independent), as shown in early-phase clinical trials (Curry et al., 2018,
> *Front Oncol*, doi:10.3389/fonc.2018.00436; Crist et al., 2022, *J Immunother Cancer*,
> doi:10.1136/jitc-2022-005632; Kemnade et al., 2023, *Oral Oncol*,
> doi:10.1016/j.oraloncology.2023.106536). We therefore position metformin as a combination-therapy
> candidate (chemo-radiation/immunotherapy), not as a single-agent inhibitor of an over-expressed
> target.

---

## M1 — Sesgo direccional del componente π (Methods)

> The differential-abundance component of TargetPriority (π_directional) is a signed, min–max-scaled
> statistic that assigns higher scores to over-abundant targets, encoding a preference for inhibiting
> up-regulated proteins. Consequently, under-abundant targets (including potential loss-of-function or
> tumour-suppressor nodes) are structurally down-weighted by this term. Network centrality is scored
> **independently of direction**; metabolic targets such as Complex I subunits, which are
> under-abundant, are therefore prioritized through centrality and external evidence rather than
> through π. This design choice favours small-molecule inhibition of over-expressed targets and should
> be considered when interpreting the ranking.

---

## M4 — Corroboración de dianas vs validación del método (Results, Fig6)

**Cambios de lenguaje:** sustituir "validation" por "corroboration/replication of target
directionality" donde aplique; no afirmar "14/14 validados".

> External proteomic (CPTAC-HNSCC) and transcriptomic (TCGA-HNSC) cohorts were used to **corroborate
> the directionality of the prioritized targets**, not to test drug efficacy or the rank order of
> candidates. Of the 14 anchor genes, directional concordance was observed for 11/12 with CPTAC data
> and 11/14 with TCGA, with 10/12 and 12/14 reaching FDR < 0.05, respectively. Concordance was strong
> for the principal candidates (EGFR, NDUFS2, DNMT1, PSMA2, ALDH5A1, MAOA, MMP8, IMPDH2, ALDH2) and
> weaker or discordant for four anchors (PDE6D, CA2, CDA, RPS11), which should be regarded as
> lower-confidence nominations. Per-gene log2FC and FDR for both cohorts are reported in Fig 6C /
> Table 6.

---

## M6 — Aporte del score sobre la centralidad (Methods/Limitations)

> The composite score is strongly correlated with network centrality alone (Spearman ρ ≈ 0.78) and
> with the differential-abundance term alone (ρ ≈ 0.78); network topology is the dominant driver, by
> design. The drug-level DrugViability term functions primarily as a clinical/pharmacological
> tie-breaker among agents that share a target rather than as the principal discriminator between
> targets.

---

## m2 — Sesgo de imputación MinProb (Methods + Limitations)

**Resultado empírico (verificado):** entre las 666 proteínas DE, las subexpresadas tienen mayor
fracción de valores ausentes que las sobreexpresadas (media 11.2% vs 6.8%; Wilcoxon p = 2.5e-7;
glm is_down ~ missingness β = 2.83, p = 1.5e-5; Spearman missingness-diferencial T−N vs log2FC =
−0.50). Solo 4/337 (1%) de las subexpresadas dependen de >50% de imputación.

> Quantification by data-independent acquisition with MinProb imputation can introduce a mild
> directional bias: proteins absent in tumour tissue are imputed at low values and therefore tend to
> appear under-abundant (under-abundant proteins had higher missingness than over-abundant proteins,
> 11.2% vs 6.8%; p = 2.5e-7). This largely reflects genuine biology — tissue-composition proteins
> (muscle, airway mucosa) are truly absent from tumour epithelium — and only 1% of under-abundant
> proteins depended on >50% imputed values. Sensitivity analyses excluding tissue-composition markers
> (see B1) confirmed that prioritized targets and pathway-level conclusions are robust to this effect.

---

## m3 — Transparencia de selección (Methods)

> Gene-set enrichment analysis tested all 50 MSigDB Hallmark gene sets; 11 reached FDR < 0.05 and are
> shown in Fig 2C. In the volcano plot (Fig 2A), gene labels were assigned to proteins with
> |log2FC| > 3.5 and FDR < 0.01 plus a fixed set of mechanistically relevant targets (e.g., EGFR,
> proteasome and Complex I subunits, DNMT1); labelling is illustrative and does not affect the
> analysis. The heatmap (Fig 2B) shows the 20 most up- and 20 most down-abundant proteins detected in
> all 20 samples (100% completeness, no imputation).

---

## m4 — Clasificación clínico-regulatoria (Methods/Supplementary)

> Drugs were assigned an ordinal regulatory class from approval status and indication: **A** =
> approved with HNSCC evidence; **B** = approved for another cancer; **C** = approved for a
> non-oncological indication; **D** = not approved/experimental, mapped to ordinal scores 1.00 / 0.75
> / 0.50 / 0.25. These weights encode rank order (regulatory maturity), not a pharmacological
> equivalence. Automated Open Targets indication matching was manually corrected for known mismatches:
> three agents were reassigned to class A (cetuximab, pembrolizumab, nivolumab — approved in HNSCC but
> missed by EFO matching), two from C to B (neratinib, lazertinib — approved in other cancers), and two
> from B to C (digoxin, digitoxin — cardiac glycosides with oncology trials but no oncological
> approval). All overrides are documented in `scripts/08_integrate_drug_targets.R` (blocks 4b–4d).

---

## m5 / m6 — Resueltos en la auditoría (no requieren texto)

- **m5** — Consistencia EGFR: el log2FC de EGFR (=2.23, FDR 6.8e-3) es único y consistente en datos,
  tablas y figuras; las menciones a "1.6/2.1" fueron confusión de exploración, no del material.
- **m6** — Claim de metilación (CDKN2A): marcado fuera de alcance en `docs/CLINICAL_CRITERIA.md`
  (script 19 removido); no debe entrar al manuscrito.
