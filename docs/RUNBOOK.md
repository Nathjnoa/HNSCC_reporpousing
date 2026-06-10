# RUNBOOK — HNSCC Drug Repurposing Pipeline

Guía paso a paso para reproducir el análisis completo (12 scripts: 01–02, 04–10, 15, 17–18).

Última actualización: 2026-04-18 — Pipeline reducido a 12 scripts; panel final = 32 candidatos LOD-stable

---

## Requisitos previos

### Ambientes conda

| Ambiente | Scripts | Propósito |
| --- | --- | --- |
| `omics-R` | 01, 02, 08, 09, 10, 15, 17, 18 | Análisis proteómica, red, scoring, sensibilidad, figuras, tablas |
| `omics-py` | 04, 05, 06, 07 | Consultas a bases de datos de fármacos |

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
```

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

# 10 - Scoring multi-criterio (5 dimensiones)
#      Pool top 35 para análisis LOD
Rscript scripts/10_prioritization_scoring.R
# Output: results/tables/10_all_candidates_scored.tsv
# Pesos: pi_stat=0.40, clinical=0.20, pathway=0.15, network=0.15, cmap=0.10
```

### Fase 4: Validación de robustez (R)

```bash
conda activate omics-R

# 15 - 6 configs de pesos + LOD (leave-one-database) + permutation test (n=1000)
#      LOD-stability es el criterio principal para definir el panel final
Rscript scripts/15_sensitivity_analysis.R
# Output: results/tables/15_sensitivity_ranks.tsv
#         results/tables/15_lod_stability.tsv      ← criterio de inclusión en panel final
```

### Fase 5: Outputs de publicación (R)

```bash
conda activate omics-R

# 17 - 7 figuras publication-ready (PDF + PNG 300 DPI)
Rscript scripts/17_pub_figures.R
# Output: results/figures/pub/main/{Sec0_FigB, Sec0_FigC, OE1_FigA, OE1_FigB,
#                                   OE2_FigA, OE2_FigB, OE2_FigC}.pdf/.png

# 18 - 6 tablas de publicación (TSV)
Rscript scripts/18_pub_tables.R
# Output: results/tables/pub/main/{Sec0_Tab1, Sec0_Tab2, OE1_Tab2,
#                                  OE2_Tab1, OE2_Tab2}.tsv
#         results/tables/pub/supp/OE2_TabS1_candidatos_extendidos_noEGFR.tsv
```

---

## Dependencias entre scripts

```text
01 → 02 → 04, 05, 06, 07 (paralelos)
               → 08 → 09 → 10 → 15 → 17
                                     18
```

---

## Outputs principales

| Archivo | Descripción |
| --- | --- |
| `results/figures/pub/main/` | 7 figuras de publicación (PDF + PNG) |
| `results/tables/pub/main/` | 5 tablas de publicación (TSV) |
| `results/tables/pub/supp/OE2_TabS1_*.tsv` | Tabla suplementaria |
| `results/tables/15_lod_stability.tsv` | Panel final: 32 candidatos LOD-stable |
| `results/tables/10_all_candidates_scored.tsv` | Todos los candidatos con scores |
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

---

## Checklist de reproducción

- [ ] Ambientes `omics-R` y `omics-py` con paquetes instalados
- [ ] Datos en `data/raw/` (3 archivos)
- [ ] `config/analysis_params.yaml` revisado
- [ ] Scripts 01–02, 04–10, 15 ejecutados en orden sin errores
- [ ] Scripts 17 y 18 ejecutados para outputs de publicación
- [ ] Logs verificados en `logs/`
