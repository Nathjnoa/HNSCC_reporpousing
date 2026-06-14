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

**EGFR as internal validation + priorización hub-céntrica de dos niveles (data-driven).**
La priorización (script 10 v3) ancla cada fármaco en la estructura de la red PPI
(red → módulo → hub drogable → fármaco) con un composite de dos factores
`TargetPriority × DrugViability`. El "titular" se lee del ranking, no se fuerza:

- **Validación:** el eje EGFR (módulo M2) encabeza y recupera fármacos ya aprobados en
  HNSCC (cetuximab, afatinib, gefitinib…) → coherencia biológica del método.
- **Tier 1 — hub-central (novedad):** OXPHOS/metformina (NDUFx), proteasoma/carfilzomib
  (PSMA2/PSMB3).
- **Tier 2 — periférico diferencial:** dianas en red diferencialmente abundantes pero
  no-hub: epigenético (DNMT1: decitabine/azacitidine; MAOA: tranylcypromine),
  antimetabolitos (IMPDH2: thioguanine), y repurposing clásico (valproato, acetazolamida,
  digoxina, disulfiram).

## Figure plan

| Figure | Content | Status |
|--------|---------|--------|
| Fig1 | Workflow / study design | Borrador (diagrama en español, pendiente EN) |
| Fig2 | (A) volcano · (B) Hallmarks GSEA · (C) top-40 DE heatmap | ✅ Finalizada — multipanel TIFF 600 DPI |
| Fig3 | (A) distribución de fase clínica · (B) clase regulatoria · (C) UpSet solapamiento de BD | ✅ Finalizada — multipanel TIFF 600 DPI; funnel completo → suppl. |
| Fig4 | (A) red coloreada por módulo Louvain (tier de drogabilidad) · (B) enriquecimiento GO por módulo · (C) hubs druggables por módulo | ✅ Finalizada — multipanel TIFF 600 DPI |
| Fig5 | (A) shortlist priorizado por módulo + descomposición composite (TP/DV), faceta por tier · (B) espacio bifactor TargetPriority × DrugViability | ✅ Finalizada — multipanel TIFF 600 DPI (`17b`+`17g`) |
| Fig6 | Validación externa en **dos cohortes** · **(A)** CPTAC-HNSCC proteoma TMT vs DIA (r=0.789, n=636, 86.3%) · **(B)** TCGA-HNSC RNA-seq vs DIA (r=0.601, n=663, 76.2%) · **(C)** 14 dianas-ancla del shortlist en CPTAC (12/14 con datos, 11 concord.) y TCGA (14/14, 11 concord.) + composite | ✅ `16b`+`16c`+`16`+`17i` — multipanel TIFF 600 DPI (A+B arriba / C ancho completo abajo) |
| FigS | Robustez: heatmap estabilidad de ranking × 6 configs de peso + LOD | ✅ `17h_figS_robustness` |
| FigS surviv. | KM OS para 4 genes-pilar (EGFR/PSMB10/DNMT1/NDUFS3); todos p>0.05 — dianas = vulnerabilidades terapéuticas, no biomarcadores pronósticos | ✅ `16` → `supp/FigS_survival_targets` |

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

**Nota Fig5 (para figure legend):** priorización hub-céntrica (scoring v3). Definir en
la leyenda:
- *Composite score* = 0.60·Target priority + 0.40·Drug viability.
- *Target priority* = 0.55·centralidad de red (grado/intermediación/eigenvector, min-max)
  + 0.45·abundancia diferencial direccional de la diana (π = sign(log2FC)·|log2FC|·−log₁₀FDR,
  min-max direccional; premia inhibir dianas UP).
- *Drug viability* = 0.40·reversión transcriptómica L2S2 + 0.40·clase regulatoria (A–D)
  + 0.20·nº de fuentes/4.
- Cada fármaco se ancla a su diana **creíble más central** (arista curada ChEMBL/OpenTargets,
  no DGIdb por su ruido); *tier* = hub_central (diana es hub de red) vs peripheral_diff.
- Líneas punteadas en B = isolíneas de composite constante.
- **No hay test de permutación** (se eliminó: para un score de agregación de evidencia la
  validez viene de recuperación de controles + robustez + validación externa, no de un p).

**Nota Fig6 (para figure legend):** validación externa en **dos cohortes
independientes**. **Panel A** — concordancia proteoma-vs-proteoma: log2FC del
proteoma DIA (tumor vs NAT, n=10 pares) vs CPTAC-HNSCC TMT (Huang et al. 2021;
limma pareado, 116 tumor / 66 normal, 66 pares); r=0.789, n=636 genes solapantes,
86.3% concordancia direccional; 4 genes-pilar etiquetados. **Panel B** —
concordancia proteoma-vs-transcriptoma: DIA vs TCGA-HNSC RNA-seq (STAR-Counts,
DESeq2, 520 tumor / 44 normal); r=0.601, n=663 genes solapantes, 76.2%
concordancia direccional; 4 genes-pilar etiquetados. **Panel C** — las 14
dianas-ancla del shortlist con log2FC en CPTAC y TCGA (símbolo relleno =
concordancia con el proteoma DIA; asteriscos = FDR<0.05/0.01/0.001) + barra de
composite descompuesta (Target priority vs Drug viability), ordenadas por
composite descendente. CPTAC: 12/14 con datos (11 concordantes, 10 FDR<0.05);
TCGA: 14/14 (11 concordantes, 11 FDR<0.05). Colores = pilar terapéutico
(EGFR/Proteasome/Epigenetic/OXPHOS); escala de grises en A/B = concordancia/
discordancia direccional de fondo.

**Convención de naming:** paneles individuales `Fig2A_*`/`Fig2B_*`/`Fig2C_*`,
`Fig3A_*`/`Fig3B_*`/`Fig3C_*` (insumos); las figuras publicables son
`Fig2_multipanel.tif` y `Fig3_multipanel.tif`. En el composite las letras de panel
las define el ensamblado — el sufijo de archivo de los paneles sueltos es solo
organizativo.

## Table plan

**Display-item scheme (decisión autor 2026-06-14):** una sola tabla **main** = el shortlist
priorizado (bloque-control EGFR + Tier1/Tier2, ensamblada de `Tab4`+`Tab5`); conteos cohorte/DE
**in-text**; el resto a **suplemento**. Leyendas en `tables.md`.

| Archivo fuente (results/tables/pub/) | Content | Placement en manuscrito |
|-------|---------|---------|
| `main/Tab4_EGFR_validation` + `main/Tab5_novel_candidates_by_module` | Shortlist priorizado (EGFR control + tiers) | **Table 1 (main)** — ensamblada |
| `main/Tab1_resumen_DE` | Cohort / DE summary | **Table S1** (+ in-text) |
| `main/Tab2_top_proteinas` | Top DE proteins | **Table S2** |
| `main/Tab3_top_candidatos` | Candidate provenance / multi-source support | **Table S3** |
| `main/Tab6_concordance_summary` | Validación externa TCGA/CPTAC (Fig6) | **Table S4** |
| `supp/TabS1_extended_candidates_by_module` | Lista extendida (todos los hubs-ancla no-EGFR) | **Table S5** |

## Word count target: 3500–4500 words
