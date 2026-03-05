# RUNBOOK — HNSCC Drug Repurposing Pipeline

Guia paso a paso para reproducir el analisis completo (16 scripts: 01–15 + 17).

---

## Requisitos previos

### Ambientes conda

| Ambiente | Scripts | Proposito |
| -------- | ------- | --------- |
| `omics-R` | 01, 02, 03, 07, 08, 09, 10, 13, 14, 15, 17 | Analisis proteomica, enriquecimiento, red, scoring, figuras pub |
| `omics-py` | 04, 05, 06, 11, 12 | Consultas a bases de datos de farmacos |

### Paquetes adicionales (instalar una sola vez)

```r
# En omics-R
BiocManager::install(c("ReactomePA", "signatureSearch", "ExperimentHub"))
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

### Fase 1: Preparar datos (R)

```bash
conda activate omics-R

# 01 - Parsear TSV, QC plots, exportar tablas DE limpias
Rscript scripts/01_parse_results_qc.R
# Output: results/tables/de_limma/01_TVsS_*.tsv, results/figures/qc/01_*.pdf

# 02 - Mapear UniProt → Entrez/Symbol
Rscript scripts/02_id_mapping.R
# Output: results/tables/de_limma/02_TVsS_*_with_ids.tsv
```

### Fase 2: Enriquecimiento funcional (R)

```bash
# 03 - ORA (GO/KEGG/Reactome) + GSEA (Hallmarks, GO BP)
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

# 06 - Open Targets GraphQL
python scripts/06_query_opentargets.py
# Output: results/tables/drug_targets/06_opentargets_*.tsv
```

```bash
conda activate omics-R

# 07 - CMap2 connectivity analysis (requiere ~370MB de cache la primera vez)
Rscript scripts/07_cmap_connectivity.R
# Output: results/tables/drug_targets/07_cmap_*.tsv

# 08 - Integrar 4 fuentes en tabla maestra
Rscript scripts/08_integrate_drug_targets.R
# Output: results/tables/drug_targets/08_*.tsv, results/figures/08_*.pdf
```

### Fase 4: Red PPI + Scoring (R)

```bash
conda activate omics-R

# 09 - Red STRING + metricas de centralidad
Rscript scripts/09_string_network.R
# Output: results/tables/network/09_*.tsv, results/figures/09_*.pdf

# 10 - Scoring multi-criterio + Top 20 (aplica exclusiones de config)
Rscript scripts/10_prioritization_scoring.R
# Output: results/tables/10_top20_candidates.tsv, results/tables/10_all_candidates_scored.tsv
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

# 15 - Analisis de sensibilidad de pesos (6 configuraciones)
Rscript scripts/15_sensitivity_analysis.R
# Output: results/tables/15_sensitivity_ranks.tsv, results/figures/15_score_distribution.pdf
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
        → 07 (CMap)
              → 08 → 09 → 10 → 11, 12 (paralelos) → 13 → 14
```

Los scripts solo pueden ejecutarse si sus dependencias estan completas. Respetar el orden numerico garantiza que todos los inputs existan.

---

## Outputs principales

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/13_FINAL_drug_candidates.xlsx` | **Excel final** con Top 20, todos los candidatos, matriz de evidencia |
| `results/tables/10_top20_candidates.tsv` | Top 20 candidatos con scores desglosados |
| `results/tables/10_all_candidates_scored.tsv` | 177 candidatos con scoring completo |
| `results/tables/network/09_network_giant.graphml` | Red PPI para visualizar en Cytoscape |
| `docs/METHODS.md` | Seccion de metodos para manuscrito |
| `docs/OUTPUTS.md` | Catalogo de todos los archivos generados |

---

## Troubleshooting

| Problema | Solucion |
| -------- | -------- |
| `conda: command not found` | `source ~/.bashrc` o usar ruta completa: `~/anaconda3/bin/conda` |
| Script 07 lento la primera vez | CMap2 descarga ~370MB de ExperimentHub. Solo ocurre una vez (se cachea). |
| Script 09 timeout en STRING API | Reintentar. STRING tiene rate limits; el script hace pausa entre batches. |
| Scripts 11-12 sin resultados de internet | Verificar conexion a internet. Las APIs de ClinicalTrials.gov y PubMed requieren acceso HTTPS. |
| `STRINGdb` no disponible | Normal. Script 09 usa STRING REST API directamente con httr2+jsonlite. |
| COSMIC CGC vacio | Descarga manual requerida (login en cancer.sanger.ac.uk). Script 12 funciona sin el con NCG7 embebido. |

---

## Checklist de reproduccion

- [ ] Ambientes `omics-R` y `omics-py` activados y con paquetes instalados
- [ ] Datos en `data/raw/` (3 archivos)
- [ ] `config/analysis_params.yaml` revisado (parametros por defecto son los usados en el analisis)
- [ ] Scripts 01-15 ejecutados en orden sin errores
- [ ] Script 17 ejecutado para figuras de publicacion
- [ ] `results/tables/13_FINAL_drug_candidates.xlsx` generado
- [ ] Logs verificados en `logs/` (cada script genera log con timestamp)
