# Manuscript — HNSCC Drug Repurposing (International Article)

Target: Q2/Q3 international journal (English, IMRAD format)
Candidates: Cancers, BMC Cancer, Frontiers in Oncology, Journal of Translational Medicine, Cancer Cell International

## Structure

| File | Section | Status |
|------|---------|--------|
| `manuscript.md` | Full IMRAD draft | In progress |
| `figures.md` | Figure legends | Pending |
| `tables.md` | Table legends | Pending |
| `supplementary.md` | Supplementary material | Pending |

## Central angle

**EGFR as internal validation + non-canonical vulnerabilities**
- Epigenetic: decitabine/azacitidine
- Proteostasis: carfilzomib
- Epigenetic modulator: tranylcypromine
- OXPHOS/Metformin: secondary finding (metabolic module in network)

## Figure plan

| Figure | Content |
|--------|---------|
| Fig1 | Workflow / study design |
| Fig2 | Volcano + heatmap + GSEA |
| Fig3 | Drug sources + phase + class |
| Fig4 | STRING PPI network + modules |
| Fig5 | Scores + LOD + sensitivity (EGFR vs non-EGFR) |
| Fig6 | External validation (CPTAC/TCGA concordance + survival) |

## Table plan

| Table | Content |
|-------|---------|
| T1 | Cohort / DE summary |
| T2 | Top prioritized candidates |
| T3 | EGFR LOD-stable (validation) |
| T4 | Non-EGFR LOD-stable (novelty) |

## Word count target: 3500–4500 words

---

## Análisis pendientes antes de submission

### Script 19 — Metilación de promotores TCGA-HNSC (PROPUESTO, no iniciado)

**Propósito:** Convertir la correlación de expresión DNMT1 en evidencia mecanística epigenética.
Eleva la narrativa de "DNMT1 sobreexpresado → hipótesis de inhibidores DNMT" a
"DNMT1 sobreexpresado → hipermetilación de supresores → los inhibidores DNMT revierten ese fenotipo".

**Datos:** Illumina HumanMethylation450K, TCGA-HNSC (~400 tumores + normales, GDC, público).
Mismo cohort que el SE ya cacheado en `data/intermediate/tcga/`.

**Tres partes del análisis:**
1. **DMPs tumor vs normal** — β-values en regiones promotoras (islas CpG ±1500 bp TSS).
   Mostrar hipermetilación diferencial en supresores conocidos (CDKN2A, CDH1, DAPK, RASSF1A).
2. **Correlación DNMT1↑ ↔ silenciamiento epigenético** — Para genes con promotor hipermetilado,
   correlacionar nivel de expresión RNA con β-value del promotor. Mostrar que la correlación
   es más fuerte en pacientes con DNMT1 alto.
3. **Methylation burden score ~ supervivencia** — Score agregado (media β en top DMPs) como
   predictor de OS en TCGA-HNSC. Une el hallazgo epigenético con relevancia clínica.

**Stack técnico (todo disponible en omics-R):**
- `TCGAbiolinks` — descarga array 450K desde GDC (~200 MB, mismo método que script 16)
- `minfi` o `ChAMP` — normalización de β-values
- `limma` — DMPs (differentially methylated positions)
- Output esperado: `results/figures/pub/main/OE4_Fig*`, `results/tables/pub/main/OE4_Tab*`

**Impacto estimado:** Diferencia entre target Cancers/BMC Cancer vs Molecular Cancer/EBioMedicine.
Revisor bio-pi: "mejor ratio impacto/costo de todos los análisis pendientes".

**Para iniciar en nueva sesión:** Crear `scripts/19_methylation_tcga.R` siguiendo el patrón de
`scripts/16_external_validation.R` (cache GDC, secciones numeradas, log, figuras pub/main/).
