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

| Figure | Content | Status |
|--------|---------|--------|
| Fig1 | Workflow / study design | Borrador (diagrama en español, pendiente EN) |
| Fig2 | (A) volcano · (B) Hallmarks GSEA · (C) top-40 DE heatmap | ✅ Finalizada — multipanel TIFF 600 DPI |
| Fig3 | Drug sources + phase + class | Pendiente revisión (hereda estilo nuevo) |
| Fig4 | STRING PPI network + modules | Pendiente revisión |
| Fig5 | Scores + LOD + sensitivity (EGFR vs non-EGFR) | Pendiente (paneles base en `results/figures/15_*`) |
| Fig6 | External validation (CPTAC/TCGA concordance + survival) | Pendiente revisión |

**Nota Fig2 (para figure legend):** los gene sets del panel B (GSEA) se ordenan por
el π-statistic = sign(log2FC) × |log2FC| × −log₁₀(FDR); incluir esta definición en
la leyenda. Layout multipanel: A (volcano) arriba-izq · B (GSEA) arriba-der · C
(heatmap, ancho completo) abajo. Estilo centralizado en `scripts/_fig_style.R`.

**Convención de naming:** paneles individuales `Fig2A_*`/`Fig2B_*`/`Fig2C_*` (insumos);
la figura publicable es `Fig2_multipanel.tif`. En el composite las letras de panel
las define el ensamblado (A=volcano, B=GSEA, C=heatmap) — el sufijo de archivo de los
paneles sueltos es solo organizativo.

## Table plan

| Table | Content |
|-------|---------|
| T1 | Cohort / DE summary |
| T2 | Top prioritized candidates |
| T3 | EGFR LOD-stable (validation) |
| T4 | Non-EGFR LOD-stable (novelty) |

## Word count target: 3500–4500 words
