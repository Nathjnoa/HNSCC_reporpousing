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

### Script 19 — Metilación de promotores TCGA-HNSC (IMPLEMENTADO 2026-06-10)

**Propósito:** Convierte la correlación de expresión DNMT1 en evidencia mecanística epigenética.
Narrativa: "DNMT1 sobreexpresado → hipermetilación de supresores → los inhibidores DNMT revierten ese fenotipo".

**Archivo:** `scripts/19_methylation_tcga.R`

**Datos:** Illumina HumanMethylation450K, TCGA-HNSC, GDC β-values armonizados (~200 MB primera descarga).
Cache: `data/intermediate/tcga/tcga_hnsc_meth450_se.rds`. Reutiliza SE RNA y clinical de script 16.

**Tres partes:**
1. **OE4_FigA** — Volcán de DMPs (limma sobre M-values): islas CpG en TSS1500/TSS200, tumor vs normal.
   Etiqueta CDKN2A, CDH1, DAPK1, RASSF1A (⚠️ claim literatura — ver CLINICAL_CRITERIA.md).
2. **OE4_FigB** — Scatter β-promotor vs RNA-seq, facetado por supresor, estratificado por DNMT1 (median-split).
   Spearman ρ por estrato en tabla `OE4_Tab2_meth_expr_corr.tsv`.
3. **OE4_FigC** — Methylation burden score (media β top 100 DMPs) ~ OS: KM + log-rank (median-split),
   idéntico a patrón OE3_FigB.

**Stack:** `TCGAbiolinks`, `limma`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`, `survival`, `survminer`.
Prerequisito: `BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")` en omics-R.

**Outputs:**
- `results/figures/pub/main/OE4_FigA_dmp_volcano.{pdf,png}`
- `results/figures/pub/main/OE4_FigB_dnmt1_meth_expr.{pdf,png}`
- `results/figures/pub/main/OE4_FigC_survival_burden.{pdf,png}`
- `results/tables/pub/main/OE4_Tab1_dmps_promoter.tsv`
- `results/tables/pub/main/OE4_Tab2_meth_expr_corr.tsv`
- `results/tables/pub/supp/OE4_TabS_survival_burden.tsv`

**Pendiente antes de submission:** Verificar claim de literatura (CDKN2A/CDH1/DAPK1/RASSF1A hipermetilados
en HNSCC) con PMIDs via @literature-scout y actualizar fila en CLINICAL_CRITERIA.md a ✅ VERIFICADO.
