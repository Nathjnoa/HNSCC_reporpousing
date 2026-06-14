# RUNBOOK — HNSCC Drug Repurposing Pipeline

Guía paso a paso para reproducir el análisis completo.

Pipeline principal: 01–03, 04–10, 15, 16b+16c+16 (validación externa), 17–18
Figuras de publicación: 17 + 17b (paneles), 17c (GSEA), 17d–17g (multipaneles Fig2–5),
17h (FigS robustez), 17i (Fig6), módulo de estilo `_fig_style.R`
Scripts archivados (fuera del manuscrito): scripts/archive/{11,12,13}
Auxiliar: supp/14_outputs_catalogue.R (regenera docs/OUTPUTS.md; METHODS.md se mantiene a mano)

Última actualización: 2026-06-13 — Validación externa en **DOS cohortes** independientes
(Fig6, 3 paneles): A = CPTAC-HNSCC proteoma TMT (`16b` fetch + `16c` limma DE), B = TCGA-HNSC
RNA-seq (`16`), C = dianas-ancla del shortlist en ambas cohortes. Antes: priorización
hub-céntrica v3 (script 10: composite `TargetPriority × DrugViability`, anclaje por arista
curada, tiers hub_central/peripheral); Fig5 (`17b`+`17g`) por módulo; robustez sin
permutación → FigS (`17h`); Tab4/Tab5 por tier/módulo. Scripts 11/12/13 archivados;
script 19 (metilación OE4) eliminado del paper.

---

## Requisitos previos

### Ambientes conda

| Ambiente | Scripts | Propósito |
| --- | --- | --- |
| `omics-R` | 01, 02, 03, 08, 09, 10, 15, 16c, 16, 17, 17b, 17c, 17d, 17e, 17f, 17g, 17h, 17i, 18 | Análisis proteómica, enriquecimiento, red, priorización, robustez, DE CPTAC, validación TCGA, figuras (paneles + multipaneles + suplementaria), tablas |
| `omics-py` | 04, 05, 06, 07, 16b | Consultas a bases de datos de fármacos (DGIdb, ChEMBL, OpenTargets, L2S2) + descarga proteómica CPTAC |

### Paquetes adicionales (instalar una sola vez)

```r
# En omics-R
BiocManager::install(c("org.Hs.eg.db", "clusterProfiler"))
install.packages(c("igraph", "ggraph", "tidygraph", "yaml", "httr2", "jsonlite",
                    "openxlsx", "fs", "patchwork", "ggrepel", "ComplexHeatmap"))
```

```bash
# En omics-py
conda activate omics-py
pip install requests pandas numpy
pip install cptac==1.5.14   # solo para 16b (descarga proteómica CPTAC-HNSCC)
```

> Validación externa (script 16) también requiere en `omics-R`:
> `BiocManager::install(c("TCGAbiolinks","DESeq2","SummarizedExperiment"))` y
> `install.packages(c("survival","survminer"))`.

### Datos requeridos

| Archivo | Ubicación | Descripción |
| --- | --- | --- |
| `proteomicacyc.csv` | `data/raw/` | Matriz de intensidades brutas |
| `metadata.csv` | `data/raw/` | Metadatos de muestras |
| `results_proteomica.tsv` | `data/raw/` | Resultados DE pre-computados (proteoDA/limma) |

### Configuración

Todos los parámetros del pipeline están centralizados en `config/analysis_params.yaml`.

---

## Ejecución del pipeline

### Directorio de trabajo

Todos los scripts esperan ejecutarse desde la raíz del proyecto:

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
```

### Nota sobre conda

Si `conda` no está en el PATH:

```bash
/home/jcarvajalv/anaconda3/bin/conda run -n omics-R Rscript scripts/01_parse_results_qc.R
# O activar primero:
source ~/.bashrc && conda activate omics-R
```

### Fase 1: Preparar datos (R)

```bash
conda activate omics-R

# 01 - Modelo limma HPV-ajustado + duplicateCorrelation (diseño pareado)
#      Exporta tablas DE
Rscript scripts/01_parse_results_qc.R
# Output: results/tables/de_limma/01_TVsS_*.tsv (666 sig: 329 up, 337 down)

# 02 - Mapear UniProt → Entrez/Symbol
Rscript scripts/02_id_mapping.R
# Output: data/intermediate/id_mapping/02_uniprot_to_ids.tsv
```

### Fase 1b: Enriquecimiento de vías (R)

```bash
conda activate omics-R

# 03 - GSEA (Hallmarks + GO BP + KEGG + Reactome). Semilla fija (set.seed 42).
#      NOTA v3: s_pathway eliminado del scoring (script 10 ya no usa GSEA).
#      GO/KEGG/Reactome quedan como contexto biológico; Fig4C usa enriquecimiento por módulo.
#      03_Hallmarks_GSEA.tsv -> panel GSEA de Fig2C (script 17c)
Rscript scripts/03_pathway_enrichment.R
# Output: results/tables/pathway_enrichment/03_*.tsv
#         results/figures/pathway_enrichment/03_*.pdf
```

### Fase 2: Consultar bases de datos de fármacos (Python)

Los scripts 04, 05, 06 son independientes y pueden ejecutarse en paralelo.

```bash
conda activate omics-py

# 04 - DGIdb GraphQL API v5
python scripts/04_query_dgidb.py
# Output: results/tables/drug_targets/04_dgidb_raw.tsv

# 05 - ChEMBL REST API (fármacos fase >= 3)
python scripts/05_query_chembl.py
# Output: results/tables/drug_targets/05_chembl_drugs.tsv

# 06 - Open Targets GraphQL API v4
python scripts/06_query_opentargets.py
# Output: results/tables/drug_targets/06_opentargets_gene_drugs.tsv

# 07 - L2S2 connectivity analysis (~2 min, requiere internet)
python scripts/07_l2s2_connectivity.py
# Output: results/tables/drug_targets/07_l2s2_results.tsv
```

### Fase 3: Integración y red PPI (R)

```bash
conda activate omics-R

# 08 - Integrar 4 fuentes en tabla maestra
Rscript scripts/08_integrate_drug_targets.R
# Output: results/tables/drug_targets/08_*.tsv

# 09 - Red STRING + Louvain modules + hubs por módulo
Rscript scripts/09_string_network.R
# Output: results/tables/network/09_*.tsv

# 10 - Priorización hub-céntrica v3: composite = 0.60·TargetPriority + 0.40·DrugViability
#      TP = 0.55·centralidad red + 0.45·pi direccional ; DV = 0.40·reversal L2S2 +
#      0.40·clase regulatoria + 0.20·breadth. Ancla = diana creíble (ChEMBL/OpenTargets)
#      más central; tier hub_central/peripheral_diff. Pesos en config: scoring_v3.
Rscript scripts/10_prioritization_scoring.R
# Output: results/tables/10_all_candidates_scored.tsv
#         results/tables/10_top20_candidates.tsv
#         results/tables/10_module_hub_candidates.tsv  ← resumen módulo→hub→fármaco
```

### Fase 4: Validación de robustez (R)

```bash
conda activate omics-R

# 15 - Robustez del composite v3: 6 configs de pesos + LOD (leave-one-database).
#      (Permutación eliminada: validez = recuperación de controles + robustez +
#       validación externa, no un p de permutación. Ver METHODS Fase 5.)
Rscript scripts/15_sensitivity_analysis.R
# Output: results/tables/15_sensitivity_ranks.tsv
#         results/tables/15_lod_stability.tsv
```

### Fase 4b: Evidencia clínica + reporte final (LEGACY — archivado)

> Esta capa (scripts 11, 12, 13) quedó **fuera del manuscrito actual** y se movió
> a `scripts/archive/`. Sus outputs (ensayos clínicos, drivers, Excel maestro
> `13_FINAL_drug_candidates.xlsx`, figuras de evidencia) no alimentan ninguna
> figura/tabla de publicación. Se conserva como referencia histórica. Para
> revivirla habría que regenerar sus inputs y reconectar rutas.

```bash
# scripts/archive/11_clinicaltrials_pubmed.py  — ClinicalTrials.gov + PubMed mining
# scripts/archive/12_cosmic_overlap.py         — overlap con drivers (COSMIC/NCG)
# scripts/archive/13_evidence_summary.R        — integración + reporte final
```

### Fase 5: Validación externa en dos cohortes (Python + R)

Fig6 valida las predicciones en **dos cohortes independientes**: CPTAC-HNSCC
(proteoma TMT, mismo nivel ómico que el descubrimiento) y TCGA-HNSC (RNA-seq,
validación cruzada multi-ómica). Orden: 16b → 16c → 16 → 17i.

```bash
# 16b - Descarga proteómica CPTAC-HNSCC (Huang et al. 2021) vía paquete `cptac`.
#       Solo descarga: el DE se hace en R (16c) con el MISMO método (limma).
conda activate omics-py
python scripts/16b_cptac_fetch.py
# Output: data/intermediate/cptac/{cptac_hnscc_proteomics.tsv, cptac_hnscc_samples.tsv}

conda activate omics-R

# 16c - DE proteómico CPTAC (limma pareado, duplicateCorrelation por patient_id).
#       Cohorte: 116 tumor / 66 normal (66 pares). Mismo método que el descubrimiento.
Rscript scripts/16c_cptac_de.R
# Output: data/intermediate/cptac/16c_cptac_hnscc_de.tsv

# 16 - Validación externa: ensambla Fig6 (3 paneles) + FigS supervivencia.
#      Lee el DE CPTAC (16c) y calcula/cachea el DE TCGA (DESeq2).
#      PRIMERA EJECUCION: descarga ~800 MB desde GDC (cache en data/intermediate/tcga/)
#      Dependencias: BiocManager::install(c("TCGAbiolinks","DESeq2","SummarizedExperiment"))
#                    install.packages(c("survival","survminer"))
#      Nota: S4Vectors enmascara dplyr::rename/count → script usa dplyr::select/count explícitos
Rscript scripts/16_external_validation.R
# Output (main):   results/figures/pub/main/Fig6A_cptac_concordance.{pdf,png}
#                  results/figures/pub/main/Fig6B_tcga_concordance.{pdf,png}
#                  results/figures/pub/main/Fig6C_targets_unified.{pdf,png}
# Output (supp):   results/figures/pub/supp/FigS_survival_targets.{pdf,png}
# Output (tables): results/tables/pub/main/Tab6_concordance_summary.tsv (CPTAC + TCGA)
#                  results/tables/pub/supp/TabS2_survival_genes.tsv
# Checkpoint en log: concordancia global + dianas-ancla con datos/concordantes/FDR<0.05

# 17i - Ensamble Fig6 multipanel (lee .rds cacheados por 16; NO re-ejecuta análisis)
#        Layout: A (CPTAC) + B (TCGA) arriba | C (dianas, ancho completo) abajo.
#        Correr DESPUÉS de 16_external_validation.R
Rscript scripts/17i_fig6_multipanel.R
# Output: results/figures/pub/main/Fig6_multipanel.{tif,pdf,png}  (TIFF 600 DPI LZW)
```

> **Panel A (Fig6A) — CPTAC proteoma vs proteoma:** concordancia del log2FC DIA vs
> CPTAC-HNSCC TMT (limma pareado, 116 tumor / 66 normal). r=0.789, n=636 genes,
> 86.3% concordancia direccional. Validación en el **mismo nivel ómico**, cohorte
> distinta.
>
> **Panel B (Fig6B) — TCGA proteoma vs RNA:** concordancia del log2FC DIA vs
> TCGA-HNSC RNA-seq (DESeq2, 520 tumor / 44 normal). r=0.601, n=663 genes, 76.2%
> concordancia direccional. Validación **cruzada multi-ómica**.
>
> **Panel C (Fig6C) — dianas-ancla en ambas cohortes:** las 14 anclas del shortlist
> (Fig5) con log2FC + FDR en CPTAC (12/14 con datos: 11 concordantes, 10 FDR<0.05)
> y TCGA (14/14: 11 concordantes, 11 FDR<0.05) + composite score. Cierra el lazo
> Fig5 → validación: las predicciones se sostienen en dos cohortes independientes.
>
> **FigS supervivencia (suplementario):** KM OS para EGFR/PSMB10/DNMT1/NDUFS3 en
> TCGA-HNSC (n=476 con datos OS), todos p>0.05 (p=0.098–0.422). Análisis
> exploratorio/contextual: las dianas son vulnerabilidades terapéuticas, **no
> biomarcadores pronósticos de OS**. Suplementario (no valida el método de repurposing).
>
> **Nota:** la antigua "Fase 5b" (metilación TCGA-HNSC, script 19 / OE4) fue
> **eliminada y excluida del manuscrito**. El script `19_methylation_tcga.R` y sus
> outputs `OE4_*` ya no existen (recuperables del historial git, commit `f3ff717`).

### Fase 6: Outputs de publicación (R)

> **Sistema de estilo centralizado.** Todos los scripts de figuras hacen
> `source("scripts/_fig_style.R")` — fuente única de paleta (Okabe-Ito,
> colorblind-safe), presets de tamaño/fuente, `theme_pub()` y funciones de
> guardado. Modificar el estilo en un solo lugar mantiene todas las figuras
> homogéneas. Sin títulos embebidos en los plots (van en los figure legends).
>
> **Convención de export:**
> - Paneles individuales → PNG (300 DPI, revisión) + PDF vectorial — `save_pub()` / `save_ch()`
> - Figura multipanel final → **TIFF 600 DPI LZW** (submission) + PDF — `save_tiff()`

```bash
conda activate omics-R

# 17 - Paneles Fig2/3/4 (PDF + PNG 300 DPI). Cachea objetos en results/figures/pub/.objects/
Rscript scripts/17_pub_figures.R

# 17b - Paneles Figura 5 (priorización hub-céntrica). Lee 10_all_candidates_scored.
Rscript scripts/17b_fig5_panels.R
# Output (cache): Fig5A_shortlist, Fig5B_tpdv_space (.rds + PDF/PNG)

# 17c - Panel GSEA Hallmarks (ordenado por π-statistic). CORRER DESPUÉS de 17.
Rscript scripts/17c_pub_figD_gsea.R

# 17d/17e/17f/17g - Ensamblado multipanel (leen objetos cacheados; NO re-grafican).
#   17d→Fig2 (volcano|GSEA·heatmap) · 17e→Fig3 (drug landscape) ·
#   17f→Fig4 (red·enrichment·hubs) · 17g→Fig5 (A shortlist por módulo | B espacio TP×DV)
Rscript scripts/17d_fig2_multipanel.R
Rscript scripts/17e_fig3_multipanel.R
Rscript scripts/17f_fig4_multipanel.R
Rscript scripts/17g_fig5_multipanel.R
# Output: results/figures/pub/main/Fig{2,3,4,5}_multipanel.{tif,pdf,png}  (TIFF 600 DPI LZW)

# 17h - Figura SUPLEMENTARIA de robustez (heatmap estabilidad × 6 configs + LOD).
Rscript scripts/17h_figS_robustness.R
# Output: results/figures/pub/supp/FigS_robustness.{tif,pdf,png}

# 17i - Figura 6 multipanel — CORRER DESPUÉS de 16_external_validation.R
Rscript scripts/17i_fig6_multipanel.R
# Output: results/figures/pub/main/Fig6_multipanel.{tif,pdf,png}  (TIFF 600 DPI LZW)

# 18 - Tablas de publicación (TSV)
Rscript scripts/18_pub_tables.R
# Output: results/tables/pub/main/{Tab1_resumen_DE, Tab2_top_proteinas,
#                                  Tab3_top_candidatos, Tab4_EGFR_validation,
#                                  Tab5_novel_candidates_by_module, Tab6_concordance_summary}.tsv
#         results/tables/pub/supp/TabS1_extended_candidates_by_module.tsv
```

---

## Dependencias entre scripts

```text
01 → 02 → 03 (pathway enrichment; panel GSEA para 17c)
       ↓
      04, 05, 06, 07 (paralelos)
               → 08 → 09 → 10 → 15 (robustez)
                          (10 anchor usa aristas curadas ChEMBL/OpenTargets de 08;
                           centralidad/módulos de 09)

Validación externa (Fig6, dos cohortes):
  16b (py: fetch CPTAC) → 16c (R: DE CPTAC limma) ─┐
                                                   ├→ 16 (paneles Fig6A/B/C + FigS) → 17i (Fig6 multipanel)
  09/10 (shortlist + DE TCGA DESeq2) ──────────────┘

Figuras (todas sourcean scripts/_fig_style.R; sin títulos embebidos):
  17  (paneles Fig2/3/4 + cachea objetos)   17b (paneles Fig5, lee 10)
  17c (GSEA, sobrescribe Fig2C)
  17d→Fig2 · 17e→Fig3 · 17f→Fig4 · 17g→Fig5   (ensamblan desde cache)
  17h→FigS_robustness (suplementaria, lee tablas 15)
  17i→Fig6 (ensambla A/B/C cacheados por 16)
  18  (tablas pub)

Opcional (archivado, no bloquea el pipeline): scripts/archive/{11,12,13}
```

---

## Outputs principales

| Archivo | Descripción |
| --- | --- |
| `results/figures/pub/main/` | Figuras de publicación: paneles (PDF + PNG 300 DPI) + multipaneles `Fig{2,3,4,5,6}_multipanel.tif` (TIFF 600 DPI) |
| `results/figures/pub/supp/FigS_robustness.{tif,pdf,png}` | Suplementaria: estabilidad ranking × 6 configs + LOD |
| `results/figures/pub/supp/FigS_survival_targets.{pdf,png}` | Suplementaria: KM OS de las 4 dianas-pilar en TCGA-HNSC (no-significativo) |
| `results/figures/pub/.objects/` | Caché de objetos de panel (`.rds`) — insumo de los `17d–17g`, regenerable, gitignored |
| `results/tables/pub/main/` | Tablas de publicación Tab1–Tab6 (TSV); Tab4=eje EGFR, Tab5=candidatos por módulo |
| `results/tables/pub/supp/TabS1_extended_candidates_by_module.tsv` | Lista extendida de hubs-ancla no-EGFR |
| `results/tables/10_module_hub_candidates.tsv` | Resumen módulo→hub→fármaco (priorización) |
| `results/tables/10_all_candidates_scored.tsv` | Todos los candidatos con scores (TP/DV/tier/módulo) |
| `results/tables/15_lod_stability.tsv` | Estabilidad leave-one-database (anotación de robustez) |
| `results/tables/network/09_network_giant.graphml` | Red PPI para Cytoscape |

---

## Troubleshooting

| Problema | Solución |
| --- | --- |
| `conda: command not found` | `source ~/.bashrc` o usar ruta completa: `~/anaconda3/bin/conda` |
| Script 07 lento | L2S2 requiere internet (l2s2.maayanlab.cloud). Cada ejecución consulta la API (~2 min). |
| Script 09 timeout STRING API | Reintentar. El script reintenta automáticamente con score=400 si score=700 no retorna aristas. |
| Script 06: error `Cannot query field 'knownDrugs'` | Verificar que se usa la versión actualizada del script (usa `drugAndClinicalCandidates`, no `knownDrugs`). |
| Script 17: columna no encontrada | Verificar que scripts 01, 08, 09, 15 se ejecutaron sin errores antes de 17. |
| Script 18: columna no encontrada | Verificar que scripts 01, 08, 10, 15 se ejecutaron sin errores antes de 18. |
| Script 17d–17g: `Faltan objetos de panel` | Correr `17` (Fig2/3/4) o `17b` (Fig5) y `17c` antes (generan la caché `.objects/*.rds`). |
| Heatmap colapsado/leyenda cortada en multipanel | Ajustar `width/height` del `grid.grabExpr()` y posición de leyendas en el ensamblador (el heatmap requiere ancho completo). |

---

## Checklist de reproducción

- [ ] Ambientes `omics-R` y `omics-py` con paquetes instalados
- [ ] Datos en `data/raw/` (3 archivos)
- [ ] `config/analysis_params.yaml` revisado
- [ ] Scripts 01–03, 04–10, 15 ejecutados en orden sin errores
- [ ] Validación externa: 16b (py, CPTAC fetch) → 16c (R, CPTAC DE) → 16 (R, TCGA + Fig6 paneles); requiere internet en primera ejecución
- [ ] Script 17i ejecutado después de 16 (ensambla Fig6_multipanel.tif)
- [ ] Figuras: 17 + 17b (paneles) → 17c (GSEA) → 17d/17e/17f/17g (multipaneles Fig2–5) → 17h (FigS robustez) → 16b/16c/16+17i (Fig6) y 18 (tablas)
- [ ] Scripts 11, 12 ejecutados si se necesita evidencia clínica/COSMIC (suplementario)
- [ ] Logs verificados en `logs/`
