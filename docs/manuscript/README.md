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
| Fig3 | (A) distribución de fase clínica · (B) clase regulatoria · (C) UpSet solapamiento de BD | ✅ Finalizada — multipanel TIFF 600 DPI; funnel completo → suppl. |
| Fig4 | (A) red coloreada por módulo Louvain (tier de drogabilidad) · (B) enriquecimiento GO por módulo · (C) hubs druggables por módulo | ✅ Finalizada — multipanel TIFF 600 DPI |
| Fig5 | Scores + LOD + sensitivity (EGFR vs non-EGFR) | Pendiente (paneles base en `results/figures/15_*`) |
| Fig6 | External validation (CPTAC/TCGA concordance + survival) | Pendiente revisión |

**Nota Fig2 (para figure legend):** los gene sets del panel B (GSEA) se ordenan por
el π-statistic = sign(log2FC) × |log2FC| × −log₁₀(FDR); incluir esta definición en
la leyenda. Layout multipanel: A (volcano) arriba-izq · B (GSEA) arriba-der · C
(heatmap, ancho completo) abajo. Estilo centralizado en `scripts/_fig_style.R`.

**Nota Fig3 (para figure legend):** "multi-source candidates" (n=458) = fármacos con
soporte en ≥2 bases de datos **O** ya aprobados (clase A/B, conservados aun con una
sola fuente), menos una lista de exclusión curada (`config`); definir así en la
leyenda (no solo "≥2 BD"). Los 3 paneles describen el mismo set de 458. Layout: A
(fase clínica) arriba-izq · B (clase) arriba-der · C (UpSet, ancho completo) abajo.
El **funnel completo** (3513→458→35→32) se movió a suplementario
(`FigS_selection_funnel`): con un solo filtro hasta multi-source no amerita embudo, y
las etapas top-ranked/LOD-stable pertenecen a la priorización (Fig5). Esto respeta la
dependencia analítica (la centralidad de red de Fig4 es un insumo del composite score).

**Nota Fig4 (para figure legend):** los módulos se colorean por **tier de drogabilidad
data-driven** (no por selección manual): *approved* = ≥1 fármaco aprobado dirigido a
un nodo del módulo (11 módulos, color propio); *hubs_only* = 0 aprobados pero ≥2 hubs
druggables (M10/mRNA processing, slate); *below threshold* = resto (M13/microtubule,
gris). **Importante para el texto:** las comunidades grises también contienen hubs
druggables (36/99 del total) y dianas priorizadas (p.ej. DNMT1→Decitabina en M17); el
color marca los ejes mecanísticos, NO exclusividad de dianas. Layout `stress` con pesos
intra-módulo (separa comunidades). A (red) ancho completo arriba · B (enriquecimiento
GO, 1 término/módulo) abajo-izq · C (hubs druggables × módulo) abajo-der.

**Convención de naming:** paneles individuales `Fig2A_*`/`Fig2B_*`/`Fig2C_*`,
`Fig3A_*`/`Fig3B_*`/`Fig3C_*` (insumos); las figuras publicables son
`Fig2_multipanel.tif` y `Fig3_multipanel.tif`. En el composite las letras de panel
las define el ensamblado — el sufijo de archivo de los paneles sueltos es solo
organizativo.

## Table plan

| Table | Content |
|-------|---------|
| T1 | Cohort / DE summary |
| T2 | Top prioritized candidates |
| T3 | EGFR LOD-stable (validation) |
| T4 | Non-EGFR LOD-stable (novelty) |

## Word count target: 3500–4500 words
