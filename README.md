# HNSCC Drug Repurposing — Reposicionamiento de Fármacos en HNSCC

## Descripción del proyecto

Análisis bioinformático de reposicionamiento de fármacos para el carcinoma escamocelular de cabeza y cuello (HNSCC), integrando datos proteómicos cuantitativos con bases de datos farmacológicas, análisis de conectividad y evidencia clínica.

**Tipo**: Proyecto R + Python
**Ambiente principal**: `omics-R` (análisis proteómico y vías), `omics-py` (consultas a bases de datos de fármacos)
**Datos**: Proteómica DIA, MaxQuant + proteoDA, Tumor vs Normal (n=10 pares pareados)

---

## Objetivos

**General**: Evaluar in silico el potencial de reposicionamiento de fármacos para HNSCC integrando datos proteómicos, análisis bioinformáticos y evidencia farmacológica.

**Específico 1**: Relacionar proteínas expresadas diferencialmente con fármacos reportados en bases de datos públicas, mediante herramientas de reposicionamiento y análisis de conectividad.

**Específico 2**: Priorizar fármacos candidatos integrando los resultados de conectividad con evidencia en rutas alteradas y análisis de redes.

**Específico 3**: Contrastar bibliográfica y clínicamente los fármacos priorizados mediante búsqueda en bases de datos especializadas.

---

## Datos

### Diseño experimental

| Variable | Detalle |
| -------- | ------- |
| Tipo de datos | Proteómica DIA (Data-Independent Acquisition) |
| Software | MaxQuant → proteoDA (wrapper de limma) |
| Comparación principal | Tumor vs Sano (TVsS) |
| Comparaciones adicionales | TpVsTn (Tumor HPV+ vs HPV-), TpVsSp (Tumor HPV+ vs Normal HPV+) |
| N muestras | 20 total = 10 pares paciente (M1-M10) |
| Variable adicional | VPH/HPV status (POSITIVE/NEGATIVE) |
| Distribución HPV | 4 pares HPV-, 6 pares HPV+ |

### Archivos de datos

| Archivo | Descripcion |
| ------- | ----------- |
| `data/raw/proteomicacyc.csv` | Matriz de intensidades brutas (peak name x muestras) |
| `data/raw/metadata.csv` | Metadatos de muestras (condition, vph, group) |
| `data/raw/results_proteomica.tsv` | Resultados DE pre-computados (output de `write_limma_tables()`) |

### Estructura del TSV de resultados

El archivo `results_proteomica.tsv` tiene 3 filas de encabezados fusionados (estilo Excel) + fila 4 con nombres de columnas reales:

- **Cols 1-3**: gene_symbol, uniprot_id, entry_name
- **Cols 4-23**: Intensidades log2 Cyclic Loess normalizadas (M1S, M1T, ..., M10S, M10T)
- **Cols 24-33**: Estadisticos TVsS (logFC, CI.L, CI.R, avg_intensity, t, B, P.Value, adj.P.Val, sig.PVal, sig.FDR)
- **Cols 34-43**: Estadisticos TpVsTn
- **Cols 44-53**: Estadisticos TpVsSp

**Codificacion sig.FDR**: `1` = upregulated significativo, `-1` = downregulated significativo, `0` = no significativo. Nota: proteoDA puede generar 2 valores borderline no enteros (0.73, 1.71) que se excluyen con filtro estricto.

---

## Metodologia del analisis proteomico (pre-computado)

El script `scripts/proteomicacyc.R` realizo el analisis en el laboratorio con **proteoDA**:

1. Filtrar muestras Pool
2. Convertir ceros a missing (NA)
3. Filtrar proteinas: `min_prop = 0.33` por grupo
4. Normalizar con **Cyclic Loess** (`cycloess`)
5. Modelo limma: `~0 + group` (4 grupos)
6. Contrastes:
   - **TVsS** = (TUMORAL_NEG + TUMORAL_POS) - (NORMAL_NEG + NORMAL_POS) <- principal
   - TpVsTn = TUMORAL_POS - TUMORAL_NEG
   - TpVsSp = TUMORAL_POS - NORMAL_POS
7. Criterios significancia: `adj.P.Val < 0.05`, `|logFC| > 1`, correccion BH

### Resultados TVsS

| Metrica | Valor |
| ------- | ----- |
| Total proteinas cuantificadas | 3,352 |
| Upregulated (sig.FDR = 1) | 248 |
| Downregulated (sig.FDR = -1) | 272 |
| **Total significativas** | **520** |
| logFC minimo (sig.) | -11.61 |
| logFC maximo (sig.) | +12.28 |

---

## Pipeline de analisis

**Nota**: El analisis proteomico DE (normalizacion, limma, contrastes) ya fue
realizado en el laboratorio con proteoDA. Los resultados estan en
`data/raw/results_proteomica.tsv`. Este pipeline parte de esos resultados.

### Estado de scripts

| Script | Fase | Estado | Descripcion |
| ------ | ---- | ------ | ----------- |
| `01_parse_results_qc.R` | 1 | COMPLETADO | Parsear TSV, QC plots, exportar tablas DE limpias |
| `02_id_mapping.R` | 1 | COMPLETADO | Agregar Entrez IDs a proteinas DE (org.Hs.eg.db). Input: `01_TVsS_significant.tsv`. Output: `02_TVsS_significant_with_ids.tsv` |
| `03_pathway_enrichment.R` | 2 | COMPLETADO | ORA: GO BP/MF/CC, KEGG, Reactome (clusterProfiler + ReactomePA). GSEA: Hallmarks MSigDB + GSEA GO BP. Input: `02_TVsS_significant_with_ids.tsv`. Output: tablas en `results/tables/pathway_enrichment/` |
| `04_query_dgidb.py` | 3 | COMPLETADO | Consultar DGIdb GraphQL API v5. 226/520 genes con interacciones. 2,252 farmacos unicos, 2,846 interacciones. Output: `results/tables/drug_targets/04_dgidb_raw.tsv` |
| `05_query_chembl.py` | 3 | COMPLETADO | Consultar ChEMBL REST API. 309/520 proteinas con target ChEMBL. 90 pares gen-farmaco fase>=3. 55 aprobados + 34 Fase III. Output: `results/tables/drug_targets/05_chembl_drugs.tsv` |
| `06_query_opentargets.py` | 3 | COMPLETADO | Consultar Open Targets GraphQL. 354/520 genes con evidencia HNSCC. 84 genes con drugs. 66 candidatos reposicionamiento novel (aprobados no-HNSCC). Top: Metformin (16 genes). Output: `results/tables/drug_targets/06_opentargets_gene_drugs.tsv` |
| `07_l2s2_connectivity.py` | 3 | COMPLETADO | L2S2 GraphQL API (l2s2.maayanlab.cloud, NAR 2025). Enriquecimiento bidireccional, filterFda=TRUE, pvalue<0.001. 1,044 drugs FDA-aprobados, 248 líneas celulares. Top: bortezomib (-1.0), calcitriol (-0.858), dasatinib (-0.855). Output: `results/tables/drug_targets/07_l2s2_results.tsv` |
| `07_cmap_connectivity_DEPRECATED.R` | — | DEPRECADO | Reemplazado por `07_l2s2_connectivity.py`. Conservado como referencia. |
| `08_integrate_drug_targets.R` | 3 | COMPLETADO | Unificar 4 fuentes (DGIdb+ChEMBL+OpenTargets+L2S2). Clasificar farmacos: A=aprobado HNSCC, B=otro cancer, C=no-cancer, D=experimental. 2,421 farmacos unicos; 187 candidatos multi-fuente. Output: `results/tables/drug_targets/08_drug_target_master_table.tsv` |
| `09_string_network.R` | 4 | COMPLETADO | Red PPI via STRING REST API v12 (score>=700). 403 nodos, 2,001 aristas. 41 hubs (top 10%), todos subunidades de Complejo I mitocondrial. 16 druggable hubs (Metformin). Export GraphML para Cytoscape. Output: `results/tables/network/09_network_node_metrics.tsv` |
| `10_prioritization_scoring.R` | 4 | COMPLETADO | Scoring multi-criterio (6 dimensiones) + filtro de exclusion de farmacos no-antitumorales. 177 candidatos evaluados (10 excluidos por nombre, 10 por target unico). Top 20 refinado. Output: `results/tables/10_top20_candidates.tsv` |
| `11_clinicaltrials_pubmed.py` | 5 | COMPLETADO | Para top 20: buscar en ClinicalTrials.gov API v2 (drug+HNSCC) y contar papers PubMed. 11/20 con trials HNSCC, 8 activos. Output: `results/tables/evidence/11_clinical_evidence.tsv` |
| `12_cosmic_overlap.py` | 5 | COMPLETADO | Cruzar targets DE con cancer driver databases (NCG7 + literatura HNSCC). 2 driver genes solapan (EGFR principal). Output: `results/tables/evidence/12_cosmic_overlap.tsv` |
| `13_evidence_summary.R` | 6 | COMPLETADO | Compilar evidencia final. Clasificar nivel 1-4 por candidato. Ranking combinado: 60% scoring + 40% evidencia. Output: `results/tables/13_FINAL_drug_candidates.xlsx` |
| `14_methods_summary.R` | 6 | COMPLETADO | Generar `docs/METHODS.md` con parametros exactos y `docs/OUTPUTS.md`. |
| `15_sensitivity_analysis.R` | 7 | COMPLETADO | Analisis de sensibilidad de pesos (6 configuraciones). 9 candidatos altamente estables (6/6 configs), top 5 ALTAMENTE ROBUSTO. Output: `results/tables/15_sensitivity_ranks.tsv` |
| `17_pub_figures.R` | 8 | COMPLETADO | Figuras de calidad publicacion (pub-figures double_col, Okabe-Ito, PDF+PNG 300 DPI). 14 figuras individuales + 5 multipanel. 5 figuras nuevas: MA plot, heatmap top DE, scatter logFC-vs-grado PPI, dot matrix scoring, bump chart. Output: `results/figures/pub/` |

### Scoring multi-criterio (script 10)

```text
Score = 0.20 x |log2FC| normalizado
      + 0.15 x significancia (-log10 adj.P)
      + 0.20 x fase clinica (aprobado=1, fase3=0.75, fase2=0.5, fase1=0.25)
      + 0.15 x L2S2 reversal score (normalizado, mas negativo = mejor)
      + 0.15 x relevancia de ruta (target en pathway enriquecido?)
      + 0.15 x centralidad en red (betweenness normalizado)
```

Pesos ajustables en `config/analysis_params.yaml`.

### Niveles de evidencia (script 13)

- **Nivel 1**: Aprobado para HNSCC
- **Nivel 2**: Aprobado otro cancer + trials HNSCC + gen COSMIC
- **Nivel 3**: Aprobado no-cancer + pathway enriquecido + hub + L2S2 reversor
- **Nivel 4**: Experimental + >=2 bases de datos de soporte

### Dependencias

```text
01 -> 02 -> 03
         -> 04, 05, 06 (paralelos)
         -> 07 (L2S2)
               -> 08 -> 09 -> 10 -> 11, 12 (paralelos) -> 13 -> 14
```

---

## Scripts implementados

### `scripts/01_parse_results_qc.R` (COMPLETADO - 2026-03-04)

**Proposito**: Parsear `results_proteomica.tsv`, validar estructura, generar figuras QC y exportar tablas DE limpias para analisis downstream.

**Input**: `data/raw/results_proteomica.tsv`, `data/raw/metadata.csv`

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/de_limma/01_TVsS_all_proteins.tsv` | 3,352 proteinas con estadisticos TVsS completos |
| `results/tables/de_limma/01_TVsS_significant.tsv` | 520 proteinas DE (up + down) |
| `results/tables/de_limma/01_TVsS_upregulated.tsv` | 248 proteinas upregulated |
| `results/tables/de_limma/01_TVsS_downregulated.tsv` | 272 proteinas downregulated |
| `results/figures/qc/01_boxplot_intensidades.pdf` | Distribucion log2 por muestra |
| `results/figures/qc/01_PCA_muestras.pdf` | PCA coloreado por condicion y VPH |
| `results/figures/qc/01_volcano_TVsS.pdf` | Volcano plot mejorado |
| `results/figures/qc/01_resumen_DE.pdf` | Conteos up/down |

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/01_parse_results_qc.R
```

---

## Paquetes adicionales requeridos

Instalar en `omics-R` antes de ejecutar scripts posteriores:

```r
# STRINGdb NO disponible (dependencia chron falla en compilacion)
# Script 09 usa STRING REST API v12 directamente con httr2 + jsonlite
BiocManager::install(c("ReactomePA"))
install.packages(c("igraph", "ggraph", "tidygraph", "yaml", "httr2", "jsonlite"))
```

Instalar en `omics-py`:

```bash
pip install chembl-webresource-client bioservices
```

---

## Controles de sesgo incorporados

1. **Multiple testing**: BH en cada paso de enriquecimiento y DE
2. **Filtro estricto sig.FDR**: solo valores exactos `1` y `-1` (excluye borderline)
3. **Multi-database**: minimo 2 bases de datos para candidatos finales
4. **Categorizacion repurposing**: A (HNSCC) / B (otro cancer) / C (no-cancer) / D (experimental)
5. **Publication bias**: conteo PubMed por candidato en Fase 5
6. **Parametros centralizados**: `config/analysis_params.yaml`

---

## Consideraciones biologicas importantes

- El diseno controla el **status HPV** (VPH): variable critica en HNSCC. El contraste TVsS promedia sobre HPV+/-, que es el enfoque correcto para el objetivo de reposicionamiento general.
- HPV+ HNSCC tiene mejor pronostico y diferente biologia que HPV-. Los contrastes TpVsTn y TpVsSp permiten explorar diferencias especificas.
- 1,112 proteinas sin NA en todas las muestras (base para PCA). Alta presencia de missing values tipica de proteomica DIA.

---

## Estructura del proyecto

```text
hnscc_drug_repurposing/
|-- data/
|   |-- raw/                    <- datos originales (NO modificar)
|   |   |-- proteomicacyc.csv
|   |   |-- metadata.csv
|   |   |-- results_proteomica.tsv
|   |   `-- cosmic/             <- para COSMIC Cancer Gene Census (descarga posterior)
|   |-- intermediate/
|   |   `-- id_mapping/
|   `-- processed/
|-- scripts/
|   |-- proteomicacyc.R         <- analisis original del laboratorio
|   `-- 01_parse_results_qc.R   <- primer script del pipeline de reposicionamiento
|-- config/
|   |-- analysis_params.yaml    <- todos los parametros centralizados
|   `-- sample_metadata.csv     <- metadatos de muestras
|-- results/
|   |-- figures/ {qc, de_limma, pathway_enrichment, drug_targets, network, final_candidates}
|   `-- tables/  {de_limma, pathway_enrichment, drug_targets, network, evidence}
|-- docs/
|   `-- RUNBOOK.md
|-- logs/
`-- env/
    `-- environment.yml
```

---

---

### `scripts/02_id_mapping.R` (COMPLETADO - 2026-03-04)

**Proposito**: Mapear UniProt IDs → Gene Symbol + Entrez ID usando `org.Hs.eg.db` para preparar inputs de enriquecimiento funcional.

**Inputs**: `results/tables/de_limma/01_TVsS_significant.tsv`, `01_TVsS_all_proteins.tsv`

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/de_limma/02_TVsS_significant_with_ids.tsv` | 520 proteinas DE con entrez_id y symbol_org |
| `results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv` | 3,352 proteinas con entrez_id (para GSEA ranked list) |
| `data/intermediate/id_mapping/02_uniprot_to_ids.tsv` | Tabla de referencia UniProt → Entrez/Symbol |

**Resultados del mapeo**:

| Metrica | Valor |
| ------- | ----- |
| UniProt IDs totales | 3,352 |
| Mapeados a Entrez | 3,263 (97.3%) |
| Sin mapeo | 89 (2.7%) |
| Significativas con Entrez | 518 / 520 (99.6%) |
| Sin mapeo en significativas | 2 (ALDH9A1 P49189, HNRNPAB Q99729) |

Nota: Las 2 proteinas significativas sin Entrez ID se retienen en la tabla (entrez_id = NA); scripts downstream las excluiran automaticamente donde sea necesario.

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/02_id_mapping.R
```

---

---

### `scripts/03_pathway_enrichment.R` (COMPLETADO - 2026-03-04)

**Proposito**: ORA y GSEA para caracterizar las vías biológicas alteradas en TVsS. Inputs: proteínas significativas con Entrez IDs (script 02) + lista rankeada completa para GSEA.

**Inputs**: `results/tables/de_limma/02_TVsS_significant_with_ids.tsv`, `02_TVsS_all_proteins_with_ids.tsv`, `config/analysis_params.yaml`

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/pathway_enrichment/03_GO_BP_ORA_full.tsv` | 77 terminos GO BP ORA |
| `results/tables/pathway_enrichment/03_GO_BP_ORA_simplified.tsv` | 19 terminos GO BP simplificados (cutoff=0.70) |
| `results/tables/pathway_enrichment/03_GO_MF_ORA.tsv` | 33 terminos GO MF ORA |
| `results/tables/pathway_enrichment/03_GO_CC_ORA.tsv` | 55 terminos GO CC ORA |
| `results/tables/pathway_enrichment/03_KEGG_ORA.tsv` | 16 rutas KEGG ORA |
| `results/tables/pathway_enrichment/03_Reactome_ORA.tsv` | 14 rutas Reactome ORA |
| `results/tables/pathway_enrichment/03_Hallmarks_GSEA.tsv` | 15 gene sets Hallmarks GSEA |
| `results/tables/pathway_enrichment/03_GO_BP_GSEA.tsv` | 455 terminos GO BP GSEA |
| `results/tables/pathway_enrichment/03_summary_top10_all.tsv` | Top 10 de cada analisis (resumen rapido) |
| `results/figures/pathway_enrichment/03_GO_BP_dotplot.pdf` | Dotplot GO BP simplificado |
| `results/figures/pathway_enrichment/03_GO_MF_dotplot.pdf` | Dotplot GO MF |
| `results/figures/pathway_enrichment/03_GO_CC_dotplot.pdf` | Dotplot GO CC |
| `results/figures/pathway_enrichment/03_KEGG_barplot.pdf` | Barplot KEGG |
| `results/figures/pathway_enrichment/03_Reactome_dotplot.pdf` | Dotplot Reactome |
| `results/figures/pathway_enrichment/03_GO_BP_cnetplot.pdf` | Cnetplot GO BP (top 8 terminos, coloreado por logFC) |
| `results/figures/pathway_enrichment/03_Hallmarks_GSEA_dotplot.pdf` | Dotplot Hallmarks GSEA (activados vs reprimidos) |
| `results/figures/pathway_enrichment/03_Hallmarks_GSEA_ridgeplot.pdf` | Ridgeplot Hallmarks GSEA |
| `results/figures/pathway_enrichment/03_GO_BP_GSEA_dotplot.pdf` | Dotplot GO BP GSEA |

**Resultados clave**:

| Analisis | N significativos | Universo |
| -------- | ---------------- | -------- |
| GO BP ORA (full) | 77 terminos | 518 sig / 3262 cuantif |
| GO BP ORA (simplif.) | 19 terminos | simplify cutoff=0.70 |
| GO MF ORA | 33 terminos | idem |
| GO CC ORA | 55 terminos | idem |
| KEGG ORA | 16 rutas | idem |
| Reactome ORA | 14 rutas | idem |
| Hallmarks GSEA | 15 gene sets | 3255 genes rankeados |
| GO BP GSEA | 455 terminos | idem |

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/03_pathway_enrichment.R
```

---

---

### `scripts/04_query_dgidb.py` (COMPLETADO - 2026-03-04)

**Proposito**: Consultar DGIdb GraphQL API v5 para identificar fármacos conocidos que interactúan con las proteínas DE. Nota: API v2 REST fue deprecada; se usa GraphQL v5 (`https://dgidb.org/api/graphql`).

**Inputs**: `results/tables/de_limma/02_TVsS_significant_with_ids.tsv`

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/drug_targets/04_dgidb_raw.tsv` | 2,846 interacciones gen-farmaco deduplicadas |
| `results/tables/drug_targets/04_dgidb_summary.tsv` | 2,252 farmacos unicos con n_targets, logFC, interaction_score |
| `results/tables/drug_targets/04_dgidb_no_interactions.txt` | 294 genes sin interacciones en DGIdb |

**Resultados clave**:

| Metrica | Valor |
| ------- | ----- |
| Genes con interacciones | 226 / 520 (43.5%) |
| Interacciones totales (dedup) | 2,846 |
| Farmacos unicos | 2,252 |
| Top farmacos por targets | Bortezomib (12), Carfilzomib (10), Ixazomib (10) — inhibidores de proteosoma |

**Ejecucion**:

```bash
conda activate omics-py
cd ~/bioinfo/projects/hnscc_drug_repurposing
python scripts/04_query_dgidb.py
```

---

---

### `scripts/05_query_chembl.py` (COMPLETADO - 2026-03-04)

**Nota tecnica**: Reescrito usando REST API directa (requests). La libreria `chembl_webresource_client` se cuelga por cache interno sin timeout adecuado.

**Resultados clave**:

| Metrica | Valor |
| ------- | ----- |
| Proteinas con target ChEMBL | 309 / 520 |
| Mecanismos totales | 174 |
| Pares gen-farmaco (dedup, fase>=3) | 90 |
| Farmacos aprobados (fase 4) | 55 |
| Farmacos Fase III | 34 |
| Genes con farmaco fase>=3 | 21 |

Farmacos destacados: Erlotinib, Cetuximab (EGFR), Lapatinib (ERBB2), Valproic acid (HDAC), Decitabine (DNMT).

---

### `scripts/06_query_opentargets.py` (COMPLETADO - 2026-03-04)

**Resultados clave**:

| Metrica | Valor |
| ------- | ----- |
| Genes DE con evidencia HNSCC (OT score) | 354 / 520 |
| Genes con Ensembl ID mapeado | 517 / 520 |
| Genes con knownDrugs | 84 |
| Pares gen-farmaco-indicacion (dedup) | 807 |
| Farmacos unicos | 125 |
| Farmacos aprobados | 67 |
| Candidatos reposicionamiento novel | 66 (aprobados, indicacion no-HNSCC) |
| Farmacos con indicacion HNSCC directa | 1 |

Top candidatos: Metformin (16 genes DE, AMPK/OXPHOS), Ocriplasmin (10 genes, proteasa), Collagenase Clostridium histolyticum (5 genes, ECM).

---

---

### `scripts/07_l2s2_connectivity.py` (COMPLETADO - 2026-03-08)

**Metodo**: L2S2 GraphQL API (l2s2.maayanlab.cloud, Evangelista et al., NAR 2025). Enriquecimiento bidireccional con Fisher's exact test. Firma query: top 150 proteinas UP + top 150 DOWN por logFC. Filtro: `filterFda=TRUE`, `pvalue < 0.001`. Paginacion: 500 nodos/request.

**Score**: `reversal_score = -log10(best_p) × mean_log2(OR) × log2(n_sigs+1)`, normalizado a [-1, 0]. Bonus ×1.5 para reversores bidireccionales.

**Resultados clave**:

| Metrica | Valor |
| ------- | ----- |
| Drugs FDA-aprobados evaluados | 1,044 |
| Lineas celulares cubiertas | 248 |
| Con reversión bidireccional consistente | 931 (89%) |
| Mejor reversor | bortezomib (score=-1.0) |
| Top 5 | bortezomib (-1.0), calcitriol (-0.858), dasatinib (-0.855), olaparib (-0.822), vorinostat (-0.766) |
| EGFR inhibidores en top 20 | gefitinib (#7), erlotinib (#14), afatinib (#17) |

**Outputs**:

- `results/tables/drug_targets/07_l2s2_results.tsv` — 1,044 drugs con reversal_score
- `results/tables/drug_targets/07_l2s2_top_reversors.tsv` — drugs con scaled_score < 0

> `scripts/07_cmap_connectivity_DEPRECATED.R` — Script CMap2 original, conservado como referencia. NO ejecutar.

---

Ultima actualizacion: 2026-03-08 | Pipeline: 14/14 COMPLETADO

---

### `scripts/08_integrate_drug_targets.R` (COMPLETADO - 2026-03-04)

**Proposito**: Unificar las 4 fuentes de farmacos (DGIdb, ChEMBL, Open Targets, L2S2) en una tabla maestra normalizada. Clasificar cada farmaco en categorias de reposicionamiento. Identificar candidatos multi-fuente.

**Inputs**: `04_dgidb_raw.tsv`, `05_chembl_drugs.tsv`, `06_opentargets_gene_drugs.tsv`, `07_l2s2_top_reversors.tsv`

**Clasificacion de farmacos**:

| Clase | Criterio | Ejemplo |
| ----- | -------- | ------- |
| A | Aprobado con indicacion directa HNSCC | Cedazuridine |
| B | Aprobado para otro cancer | Bortezomib, Lapatinib |
| C | Aprobado para indicacion no-oncologica | Metformin, Valproic acid |
| D | Experimental / en investigacion | Compuestos experimentales no-aprobados |

**Resultados**:

| Metrica | Valor |
| ------- | ----- |
| Farmacos unicos totales | 2,421 |
| Clase A (HNSCC directo) | 1 (Cedazuridine) |
| Clase B (otro cancer) | 51 |
| Clase C (no-oncologico, reposicionamiento novel) | 80 |
| Clase D (experimental) | 2,289 |
| Candidatos multi-fuente (>=2 bases de datos) | 187 |

**Nota tecnica**: Los booleans de Python ("True"/"False") no son parseables directamente por `as.logical()` en R. Solucion: `tolower(as.character(col)) == "true"`.

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/drug_targets/08_drug_target_master_table.tsv` | Tabla maestra de todos los farmacos con genes, fuentes, clase |
| `results/tables/drug_targets/08_drug_summary_per_drug.tsv` | 2,421 farmacos con n_sources, n_targets, clase, cmap_score |
| `results/tables/drug_targets/08_multi_source_candidates.tsv` | 187 candidatos en >=2 fuentes (prioridad alta) |
| `results/figures/drug_targets/08_class_distribution.pdf` | Distribucion de clases A-D |
| `results/figures/drug_targets/08_multisource_heatmap.pdf` | Heatmap de fuentes por candidato top |
| `results/figures/drug_targets/08_top_candidates_lollipop.pdf` | Lollipop de candidatos multi-fuente top |
| `results/figures/drug_targets/08_logfc_vs_nsources.pdf` | Scatter logFC vs numero de fuentes |

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/08_integrate_drug_targets.R
```

---

---

### `scripts/09_string_network.R` (COMPLETADO - 2026-03-04)

**Proposito**: Construir red PPI de las 518 proteinas DE con Entrez ID usando STRING REST API v12 (score combinado >=700 = alta confianza). Calcular metricas de centralidad. Identificar hubs druggables.

**Nota tecnica importante**: El paquete STRINGdb no esta disponible (dependencia `chron` falla en compilacion). Se uso STRING REST API v12 directamente con `httr2` + `jsonlite::fromJSON(resp_body_string())`. STRING devuelve `Content-Type: text/json` (no `application/json`), por lo que `resp_body_json()` falla y hay que usar `resp_body_string()` primero.

**Inputs**: `results/tables/de_limma/02_TVsS_significant_with_ids.tsv`, `results/tables/drug_targets/08_drug_summary_per_drug.tsv`

**Resultados de red**:

| Metrica | Valor |
| ------- | ----- |
| Genes consultados | 518 |
| Mapeados en STRING | 513 (99%) |
| Nodos en componente gigante | 403 |
| Aristas (score >= 700) | 2,001 |
| Hubs (top 10% grado) | 41 |
| Hubs druggables | 16 |
| Identidad biologica de los hubs | Subunidades de Complejo I mitocondrial |

**Hallazgo clave**: Los 41 hubs son exclusivamente subunidades de cadena respiratoria mitocondrial (NDUFA*, NDUFB*, NDUFS*, ATP5*, COX*, UQCR*), todos downregulados en tumor. Los 16 hubs druggables son targets conocidos de **Metformin** (inhibidor de Complejo I). Convergencia biologica notable.

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/network/09_network_node_metrics.tsv` | 513 nodos con degree, betweenness, eigenvector, is_hub |
| `results/tables/network/09_network_edges.tsv` | 2,001 aristas con combined_score |
| `results/tables/network/09_druggable_hubs.tsv` | 16 hubs con fármacos conocidos |
| `results/network/09_network_giant.graphml` | Red para importar en Cytoscape |
| `results/figures/network/09_hub_degree_distribution.pdf` | Distribucion de grado (ley de potencia) |
| `results/figures/network/09_hub_metrics_scatter.pdf` | Scatter degree vs betweenness |
| `results/figures/network/09_druggable_hubs_barplot.pdf` | Barplot de hubs druggables |
| `results/figures/network/09_network_full.pdf` | Visualizacion completa de red (ggraph) |
| `results/figures/network/09_network_hubs_only.pdf` | Subred de hubs (ggraph coloreado por logFC) |

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/09_string_network.R
```

---

---

### `scripts/10_prioritization_scoring.R` (COMPLETADO - 2026-03-04)

**Proposito**: Scoring multi-criterio para los candidatos multi-fuente. 6 dimensiones con pesos configurables en `analysis_params.yaml`. Filtro de exclusion de farmacos no-antitumorales (gabapentinoides, anticoagulantes, anti-amiloides, antidotos). Seleccion de Top 20 con diversidad biologica (max 3 por target primario).

**Formula de scoring**:

```text
composite = 0.20 x |logFC| normalizado (0-1)
          + 0.15 x significancia (-log10 adj.P normalizado)
          + 0.20 x fase clinica (aprobado=1, fase3=0.75, fase2=0.5, fase1=0.25, exp=0)
          + 0.15 x L2S2 reversal score (negativo invertido, normalizado)
          + 0.15 x pathway relevance (target en pathway enriquecido)
          + 0.15 x centralidad red (betweenness normalizado)
```

**Exclusiones**: 17 farmacos por nombre + 10 por target unico no-antitumoral (F2, APP, CACNA2D1, ADH1B). Configurables en `config/analysis_params.yaml` seccion `exclusions`.

**Resultados**:

| Metrica | Valor |
| ------- | ----- |
| Candidatos evaluados (post-exclusion) | 177 |
| Score maximo | 0.637 (Mavacamten) |
| Score minimo Top 20 | 0.494 |
| Clase A en Top 20 | 1 (Cedazuridine) |
| Clase B en Top 20 | 5 |
| Clase C en Top 20 | 13 |
| Clase D en Top 20 | 1 |

**Top 10 candidatos (refinados)**:

| Rank | Farmaco | Score | Target | Clase |
| ---- | ------- | ----- | ------ | ----- |
| 1 | Mavacamten | 0.637 | MYH7 | C |
| 2 | Erlotinib HCl | 0.607 | EGFR | C |
| 3 | Lapatinib | 0.607 | EGFR | C |
| 4 | Cetuximab | 0.607 | EGFR | C |
| 5 | Metformin HCl | 0.600 | NDUFA/NDUFB/NDUFS | C |
| 6 | Omecamtiv Mecarbil | 0.577 | MYH7 | B |
| 7 | Collagenase C.h. | 0.542 | COL/MMP8 | C |
| 8 | Cedazuridine | 0.541 | CDA | A |
| 9 | Ocriplasmin | 0.536 | COL/LAM | C |
| 10 | Danicamtiv | 0.527 | MYH7 | D |

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/10_top20_candidates.tsv` | Top 20 candidatos con scores desglosados |
| `results/tables/10_all_candidates_scored.tsv` | 187 candidatos con scoring completo |
| `results/figures/10_top20_heatmap.pdf` | Heatmap de componentes de score por candidato |
| `results/figures/10_top20_lollipop.pdf` | Lollipop de score final, coloreado por clase |
| `results/figures/10_score_radar.pdf` | Radar chart de dimensiones por Top 5 |

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/10_prioritization_scoring.R
```

---

---

### `scripts/11_clinicaltrials_pubmed.py` (COMPLETADO - 2026-03-04)

**Proposito**: Consultar ClinicalTrials.gov API v2 y PubMed E-utilities para cada candidato del Top 20. Medir evidencia clinica y bibliografica en HNSCC.

**Nota tecnica**: Nombres de farmacos con sales (e.g., "ERLOTINIB HYDROCHLORIDE") se simplifican automaticamente para PubMed (`simplify_drug_name()`) para evitar resultados artificialmente bajos.

**Resultados**:

| Metrica | Valor |
| ------- | ----- |
| Candidatos consultados | 20 |
| Con ensayos HNSCC | 11 |
| Con trials activos | 8 |
| Con papers HNSCC | 9 |
| Mejor evidence_score | 1.0 (Erlotinib, Cetuximab, Metformin) |

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/evidence/11_clinical_evidence.tsv` | 20 candidatos con trials, papers, evidence_score |
| `results/tables/evidence/11_trials_detail.tsv` | Detalle de todos los ensayos encontrados |
| `results/figures/evidence/11_evidence_bubble.pdf` | Bubble plot: trials vs papers, tamano=evidence_score |
| `results/figures/evidence/11_trials_phase_bar.pdf` | Barras de fases clinicas por candidato |

**Ejecucion**:

```bash
conda activate omics-py
cd ~/bioinfo/projects/hnscc_drug_repurposing
python scripts/11_clinicaltrials_pubmed.py
```

---

---

### `scripts/12_cosmic_overlap.py` (COMPLETADO - 2026-03-04)

**Proposito**: Cruzar proteinas DE con bases de datos de genes driver de cancer (NCG7, literatura HNSCC). Anotar Top 20 con informacion de drivers.

**Nota tecnica**: COSMIC CGC requiere descarga manual (login-gated); el script maneja su ausencia gracefully usando NCG7 embebido (~126 oncogenes) + 22 drivers HNSCC de literatura como fuentes alternativas.

**Resultados**:

| Metrica | Valor |
| ------- | ----- |
| Proteinas DE totales | 520 |
| Overlap con NCG7 | 2 genes |
| Overlap con literatura HNSCC | 1 gen |
| Driver principal | EGFR (logFC=4.33, adj.P=0.015) |
| Top 20 con target driver | 3/20 |

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/evidence/12_cosmic_overlap.tsv` | 2 genes DE que son drivers de cancer |
| `results/tables/evidence/12_cancer_drivers_annot.tsv` | Anotacion completa de drivers |
| `results/tables/evidence/12_top20_with_drivers.tsv` | Top 20 con columna targets_driver |
| `results/figures/evidence/12_driver_overlap_bar.pdf` | Barras de overlap por fuente |
| `results/figures/evidence/12_driver_logfc_plot.pdf` | logFC de genes driver |

**Ejecucion**:

```bash
conda activate omics-py
cd ~/bioinfo/projects/hnscc_drug_repurposing
python scripts/12_cosmic_overlap.py
```

---

---

### `scripts/13_evidence_summary.R` (COMPLETADO - 2026-03-04)

**Proposito**: Integrar scoring (script 10) + evidencia clinica (script 11) + drivers (script 12) en ranking final. Generar Excel multi-hoja y figuras finales.

**Ranking combinado**:

```text
combined_rank_score = 0.60 x composite_score + 0.40 x evidence_score
```

**Resultados**:

| Metrica | Valor |
| ------- | ----- |
| Top 5 ranking final | Erlotinib (0.764), Cetuximab (0.764), Metformin (0.760), Lapatinib (0.724), Doxycycline (0.638) |
| Niveles de evidencia | Nivel 1: 0, Nivel 2: 1 (Cedazuridine), Nivel 3: mayoria, Nivel 4: 1 |

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `results/tables/13_FINAL_drug_candidates.xlsx` | Excel 5 hojas: Top20, Todos, Matriz, Ensayos, Metodologia |
| `results/tables/13_evidence_matrix.tsv` | Matriz binaria de evidencia (10 dimensiones x 20 candidatos) |
| `results/figures/final/13_evidence_heatmap.pdf` | Heatmap binario de evidencia |
| `results/figures/final/13_final_ranking.pdf` | Barplot ranking combinado |
| `results/figures/final/13_multipanel_summary.pdf` | Multipanel resumen (4 paneles) |

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/13_evidence_summary.R
```

---

---

### `scripts/14_methods_summary.R` (COMPLETADO - 2026-03-04)

**Proposito**: Generar documentacion reproducible automatica. Lee parametros de `config/analysis_params.yaml` y versiones de paquetes para escribir seccion de metodos.

**Outputs**:

| Archivo | Descripcion |
| ------- | ----------- |
| `docs/METHODS.md` | Seccion de metodos lista para manuscrito |
| `docs/OUTPUTS.md` | Catalogo de todos los archivos generados por el pipeline |

**Ejecucion**:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/14_methods_summary.R
```

---

## Documentacion adicional

- `docs/RUNBOOK.md` — Orden de ejecucion y comandos
- `docs/PROGRESS_REPORT.md` — Explicacion didactica de la investigacion (teoria + resultados)
- `docs/METHODS.md` — Seccion de metodos para manuscrito (generado por script 14)
- `docs/OUTPUTS.md` — Catalogo de archivos generados (generado por script 14)
- `docs/FUTURE_WORK.md` — Proximos pasos planificados (TCGA survival, figuras publicacion, manuscrito)
