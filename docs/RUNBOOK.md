# RUNBOOK — HNSCC Drug Repurposing Pipeline

Guia paso a paso para reproducir el analisis completo (16 scripts: 01–15 + 17).

*Última actualización: 2026-04-05 — Revisión metodológica v3: pi-stat direccional, reclasificación DIGOXIN, exclusión miosinas cardíacas; panel final = 23 candidatos LOD-stable*

---

## Requisitos previos

### Ambientes conda

| Ambiente | Scripts | Proposito |
| -------- | ------- | --------- |
| `omics-R` | 01, 02, 03, 08, 09, 10, 13, 14, 15, 17 | Analisis proteomica, enriquecimiento, red, scoring, figuras pub |
| `omics-py` | 04, 05, 06, 07, 11, 12 | Consultas a bases de datos de farmacos (incluye L2S2) |

### Paquetes adicionales (instalar una sola vez)

```r
# En omics-R
BiocManager::install(c("ReactomePA"))
install.packages(c("igraph", "ggraph", "tidygraph", "yaml", "httr2", "jsonlite",
                    "openxlsx", "fs", "patchwork", "ggrepel"))
```

```bash
# En omics-py
conda activate omics-py
pip install chembl-webresource-client bioservices requests
```

### Datos requeridos

| Archivo | Ubicacion | Descripcion |
| ------- | --------- | ----------- |
| `proteomicacyc.csv` | `data/raw/` | Matriz de intensidades brutas |
| `metadata.csv` | `data/raw/` | Metadatos de muestras |
| `results_proteomica.tsv` | `data/raw/` | Resultados DE pre-computados (proteoDA/limma) |

### Configuracion

Todos los parametros del pipeline estan centralizados en `config/analysis_params.yaml`. Modificar ahi para analisis de sensibilidad.

---

## Ejecucion del pipeline

### Directorio de trabajo

Todos los scripts esperan ejecutarse desde la raiz del proyecto:

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
```

### Nota sobre conda

Si `conda` no está en el PATH, usar la ruta completa:

```bash
# En lugar de: conda activate omics-R
/home/jcarvajalv/anaconda3/bin/conda run -n omics-R Rscript scripts/01_parse_results_qc.R

# O activar primero:
source ~/.bashrc && conda activate omics-R
```

### Fase 1: Preparar datos (R)

```bash
conda activate omics-R

# 01 - Modelo limma HPV-ajustado + duplicateCorrelation (diseño pareado)
#      Exporta tablas DE + análisis de missingness
Rscript scripts/01_parse_results_qc.R
# Output: results/tables/de_limma/01_TVsS_*.tsv (666 sig: 329 up, 337 down)
#         results/tables/de_limma/01_missingness_analysis.tsv
#         results/figures/qc/01_*.pdf

# 02 - Mapear UniProt → Entrez/Symbol
Rscript scripts/02_id_mapping.R
# Output: results/tables/de_limma/02_TVsS_*_with_ids.tsv
```

### Fase 2: Enriquecimiento funcional (R)

```bash
# 03 - ORA (GO/KEGG/Reactome) + GSEA con ranking π-estadístico (Hallmarks, GO BP, KEGG, Reactome)
Rscript scripts/03_pathway_enrichment.R
# Output: results/tables/pathway_enrichment/03_*.tsv, results/figures/pathway_enrichment/03_*.pdf
```

### Fase 3: Consultar bases de datos de farmacos (Python + R)

Los scripts 04, 05, 06 son independientes y pueden ejecutarse en paralelo.

```bash
conda activate omics-py

# 04 - DGIdb GraphQL API v5
python scripts/04_query_dgidb.py
# Output: results/tables/drug_targets/04_dgidb_*.tsv

# 05 - ChEMBL REST API (farmacos fase >= 3)
python scripts/05_query_chembl.py
# Output: results/tables/drug_targets/05_chembl_drugs.tsv

# 06 - Open Targets GraphQL (filtra por score HNSCC >= 0.2, ver config/analysis_params.yaml)
python scripts/06_query_opentargets.py
# Output: results/tables/drug_targets/06_opentargets_*.tsv
```

```bash
conda activate omics-py

# 07 - L2S2 connectivity analysis (~2 min, requiere internet)
python scripts/07_l2s2_connectivity.py
# Output: results/tables/drug_targets/07_l2s2_results.tsv
#         results/tables/drug_targets/07_l2s2_top_reversors.tsv
```

```bash
conda activate omics-R

# 08 - Integrar 4 fuentes en tabla maestra
Rscript scripts/08_integrate_drug_targets.R
# Output: results/tables/drug_targets/08_*.tsv, results/figures/08_*.pdf
```

### Fase 4: Red PPI + Scoring (R)

```bash
conda activate omics-R

# 09 - Red STRING + Louvain modules + hubs por módulo (top 10% betweenness intra-módulo)
Rscript scripts/09_string_network.R
# Output: results/tables/network/09_*.tsv (incluye 09_modules.tsv), results/figures/09_*.pdf

# 10 - Scoring multi-criterio v2 (pi-stat + clinical + pathway + network + cmap)
#      Sin límite de diversidad por target; sin score_evidence en composite
Rscript scripts/10_prioritization_scoring.R
# Output: results/tables/10_top20_candidates.tsv (top 35), results/tables/10_all_candidates_scored.tsv
# Pesos v2: pi_stat=0.325, clinical=0.20, pathway=0.15, network=0.195, cmap=0.13
# Panel final = candidatos con lod_stable=TRUE en 15_lod_stability.tsv
```

### Fase 5: Validacion in silico (Python)

```bash
conda activate omics-py

# 11 - ClinicalTrials.gov + PubMed (requiere internet)
python scripts/11_clinicaltrials_pubmed.py
# Output: results/tables/evidence/11_clinical_evidence.tsv, results/figures/evidence/11_*.pdf

# 12 - Cruce con genes driver de cancer
python scripts/12_cosmic_overlap.py
# Output: results/tables/evidence/12_*.tsv, results/figures/evidence/12_*.pdf
```

### Fase 6: Evidencia final + Documentacion (R)

```bash
conda activate omics-R

# 13 - Ranking combinado + Excel final + figuras
Rscript scripts/13_evidence_summary.R
# Output: results/tables/13_FINAL_drug_candidates.xlsx, results/figures/final/13_*.pdf

# 14 - Generar METHODS.md y OUTPUTS.md
Rscript scripts/14_methods_summary.R
# Output: docs/METHODS.md, docs/OUTPUTS.md
```

### Fase 7: Validacion de robustez (R)

```bash
conda activate omics-R

# 15 - Análisis de sensibilidad: 6 configs de pesos + LOD (leave-one-database) + permutation test (n=1000)
#      LOD-stability es el criterio principal para definir el panel final de candidatos
Rscript scripts/15_sensitivity_analysis.R
# Output: results/tables/15_sensitivity_ranks.tsv
#         results/tables/15_lod_stability.tsv      ← criterio de inclusión en panel final
#         results/tables/15_permutation_test.tsv
```

### Fase 8: Figuras de publicacion (R)

```bash
conda activate omics-R

# 17 - Figuras calidad publicacion (PDF + PNG 300 DPI)
Rscript scripts/17_pub_figures.R
# Output: results/figures/pub/A1..F2_*.pdf/.png, results/figures/pub/FIG1..FIG5_*.pdf/.png
```

---

## Dependencias entre scripts

```text
01 → 02 → 03
        → 04, 05, 06 (paralelos)
        → 07 (L2S2)
              → 08 → 09 → 10 → 11, 12 (paralelos) → 13 → 14
                              ↓
                             15 (sensibilidad; requiere columna `bonus` de 10)
                              ↓
                             17 (figuras pub; requiere outputs de 01, 03, 09, 10, 13, 15)
```

Los scripts solo pueden ejecutarse si sus dependencias estan completas. Respetar el orden numerico garantiza que todos los inputs existan.

---

## Outputs principales

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/13_FINAL_drug_candidates.xlsx` | **Excel final** con Top 20, todos los candidatos, matriz de evidencia |
| `results/tables/10_top20_candidates.tsv` | Top 20 candidatos con scores desglosados |
| `results/tables/10_all_candidates_scored.tsv` | Candidatos con scoring completo (s_pi_stat, s_clinical, s_cmap, s_pathway, s_network) |
| `results/tables/network/09_network_giant.graphml` | Red PPI para visualizar en Cytoscape |
| `docs/METHODS.md` | Seccion de metodos para manuscrito |
| `docs/OUTPUTS.md` | Catalogo de todos los archivos generados |

---

## Troubleshooting

| Problema | Solucion |
| -------- | -------- |
| `conda: command not found` | `source ~/.bashrc` o usar ruta completa: `~/anaconda3/bin/conda` |
| Script 07 lento o sin resultados | L2S2 requiere conexión a internet (l2s2.maayanlab.cloud). No descarga datos localmente; cada ejecución consulta la API (~2 min). |
| Script 09 timeout en STRING API | Reintentar. STRING tiene rate limits; el script hace pausa entre batches. Si score >= 700 no retorna aristas, el script reintenta automáticamente con score = 400 y actualiza el umbral en consecuencia. |
| Scripts 11-12 sin resultados de internet | Verificar conexion a internet. Las APIs de ClinicalTrials.gov y PubMed requieren acceso HTTPS. |
| Script 11: advertencia de email NCBI | El script usa `NCBI_EMAIL` como identificador en PubMed E-utilities. Por defecto usa `jcarvajal@fucsalud.edu.co`. Para sobreescribir: `export NCBI_EMAIL="tu@email.com"` antes de ejecutar. |
| `STRINGdb` no disponible | Normal. Script 09 usa STRING REST API directamente con httr2+jsonlite. |
| COSMIC CGC vacio | Descarga manual requerida (login en cancer.sanger.ac.uk). Script 12 funciona sin el con NCG7 embebido. |
| Script 17: `avg_intensity_TVsS` no encontrada | Columna opcional. Script 17 detecta automáticamente `AveExpr` (columna estándar de salida limma) como alternativa para el MA plot. No requiere intervención. |
| Script 17: `drug_class_label` no encontrada | La columna exportada por script 10 se llama `drug_class` (valores A/B/C/D). Script 17 traduce internamente a etiquetas descriptivas; no existe columna `drug_class_label` en los archivos intermedios. |
| Script 15: advertencia `bonus` column not found | La columna `bonus` la exporta script 10. Si no está presente, script 15 usa 0 y continúa; verificar que script 10 se ejecutó sin errores y que `10_all_candidates_scored.tsv` está completo. |
| Script 17: `logFC_TVsS` no encontrada en nodos de red | El archivo `09_network_node_metrics.tsv` solo contiene métricas de red. Script 17 hace un `left_join` automático desde la tabla DE para añadir `logFC_TVsS` y `adj.P.Val_TVsS`. No requiere intervención. |

---

## Checklist de reproduccion

- [ ] Ambientes `omics-R` y `omics-py` activados y con paquetes instalados
- [ ] Datos en `data/raw/` (3 archivos)
- [ ] `config/analysis_params.yaml` revisado (parametros por defecto son los usados en el analisis)
- [ ] Scripts 01-15 ejecutados en orden sin errores
- [ ] Script 17 ejecutado para figuras de publicacion
- [ ] `results/tables/13_FINAL_drug_candidates.xlsx` generado
- [ ] Logs verificados en `logs/` (cada script genera log con timestamp)
