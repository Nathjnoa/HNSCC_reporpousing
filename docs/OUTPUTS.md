# Outputs Catalogue — HNSCC Drug Repurposing

*Última actualización: 2026-03-06 — verificado contra archivos en disco; figuras de scripts intermedios reemplazadas por versiones pub/ cuando aplica.*
*Organizado script por script, con rutas relativas desde la raíz del proyecto.*

Convención de la columna Tipo: `T` = tabla/archivo de datos, `F` = figura.

---

## Script 01 — `01_parse_results_qc.R`

Parseo de resultados proteómica (proteoDA/limma), QC y exportación de proteínas DE.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/de_limma/01_TVsS_all_proteins.tsv` | Todas las proteínas cuantificadas con estadísticos del contraste Tumor vs Sano (logFC, adj.P.Val, intensidades) |
| T | `results/tables/de_limma/01_TVsS_significant.tsv` | Proteínas DE significativas (sig.FDR = 1 o −1): up + down |
| T | `results/tables/de_limma/01_TVsS_upregulated.tsv` | Solo proteínas upreguladas (sig.FDR = 1) |
| T | `results/tables/de_limma/01_TVsS_downregulated.tsv` | Solo proteínas downreguladas (sig.FDR = −1) |
| F | `results/figures/qc/01_boxplot_intensidades.pdf` | Distribución log2 de intensidades por muestra, coloreada por condición y VPH |
| F | `results/figures/qc/01_resumen_DE.pdf` | Barplot de conteos de proteínas upreguladas vs downreguladas |

---

## Script 02 — `02_id_mapping.R`

Mapeo UniProt → Entrez ID + símbolo canónico via `org.Hs.eg.db`.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/de_limma/02_TVsS_significant_with_ids.tsv` | Proteínas significativas con columnas adicionales `entrez_id` y `symbol_org` |
| T | `results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv` | Todas las proteínas cuantificadas con IDs mapeados (input para GSEA en script 03) |
| T | `data/intermediate/id_mapping/02_uniprot_to_ids.tsv` | Tabla de referencia: UniProt ID → Entrez ID + Symbol (incluye isoformas limpias) |

---

## Script 03 — `03_pathway_enrichment.R`

Enriquecimiento funcional: ORA (GO/KEGG/Reactome) y GSEA (Hallmarks/GO BP).

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/pathway_enrichment/03_GO_BP_ORA_full.tsv` | ORA GO Biological Process completo (todos los términos significativos) |
| T | `results/tables/pathway_enrichment/03_GO_BP_ORA_simplified.tsv` | ORA GO BP simplificado (cutoff semántico 0.70, elimina redundancia) |
| T | `results/tables/pathway_enrichment/03_GO_MF_ORA.tsv` | ORA GO Molecular Function |
| T | `results/tables/pathway_enrichment/03_GO_CC_ORA.tsv` | ORA GO Cellular Component |
| T | `results/tables/pathway_enrichment/03_KEGG_ORA.tsv` | ORA KEGG Pathways |
| T | `results/tables/pathway_enrichment/03_Reactome_ORA.tsv` | ORA Reactome Pathways |
| T | `results/tables/pathway_enrichment/03_Hallmarks_GSEA.tsv` | GSEA MSigDB Hallmarks con NES y p.adjust |
| T | `results/tables/pathway_enrichment/03_GO_BP_GSEA.tsv` | GSEA GO Biological Process |
| T | `results/tables/pathway_enrichment/03_summary_top10_all.tsv` | Top 10 resultados de cada análisis agrupados (vista rápida integrada) |
| F | `results/figures/pathway_enrichment/03_GO_BP_dotplot.pdf` | Dotplot top 20 términos GO BP ORA (versión simplificada) |
| F | `results/figures/pathway_enrichment/03_GO_MF_dotplot.pdf` | Dotplot top 15 términos GO MF ORA |
| F | `results/figures/pathway_enrichment/03_GO_CC_dotplot.pdf` | Dotplot top 15 términos GO CC ORA |
| F | `results/figures/pathway_enrichment/03_KEGG_barplot.pdf` | Barplot top 15 rutas KEGG ORA |
| F | `results/figures/pathway_enrichment/03_Reactome_dotplot.pdf` | Dotplot top 15 rutas Reactome ORA |
| F | `results/figures/pathway_enrichment/03_GO_BP_emapplot.pdf` | Red de similitud semántica entre términos GO BP (top 30) |
| F | `results/figures/pathway_enrichment/03_GO_BP_cnetplot.pdf` | Red gen-concepto GO BP: top 8 términos con genes coloreados por logFC |
| F | `results/figures/pathway_enrichment/03_GO_BP_GSEA_dotplot.pdf` | Dotplot GSEA GO BP separado por dirección de enriquecimiento |

---

## Script 04 — `04_query_dgidb.py`

Consulta a DGIdb API v4 para interacciones gen-fármaco de los genes DE.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/drug_targets/04_dgidb_raw.tsv` | Interacciones gen-fármaco crudas descargadas de DGIdb (todas las fuentes) |
| T | `results/tables/drug_targets/04_dgidb_summary.tsv` | Resumen deduplicado: un registro por par gen-fármaco con interaction_score |
| T | `results/tables/drug_targets/04_dgidb_no_interactions.txt` | Lista de genes DE sin interacciones en DGIdb |

---

## Script 05 — `05_query_chembl.py`

Consulta a ChEMBL API por fármacos en fase clínica >= 3 para genes DE.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/drug_targets/05_chembl_drugs.tsv` | Fármacos ChEMBL con max_phase >= 3 para genes DE, con SMILES y actividad |
| T | `results/tables/drug_targets/05_chembl_summary.tsv` | Resumen agregado por fármaco con fase máxima y genes diana |

---

## Script 06 — `06_query_opentargets.py`

Consulta a Open Targets Platform (GraphQL) para evidencia HNSCC de genes DE.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/drug_targets/06_opentargets_gene_drugs.tsv` | Interacciones gen-fármaco con indicaciones, aprobación y evidencia HNSCC (is_hnscc_indication) |
| T | `results/tables/drug_targets/06_opentargets_hnscc_scores.tsv` | Scores de asociación Open Targets para genes DE respecto a HNSCC (overall_score) |

---

## Script 07 — `07_cmap_connectivity.R`

Análisis de conectividad CMap2/LINCS con `signatureSearch` para identificar compuestos que revierten la firma tumoral.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/drug_targets/07_cmap_results.tsv` | Todos los compuestos CMap2 evaluados, ordenados por scaled_score (más negativo = mejor reversor) |
| T | `results/tables/drug_targets/07_cmap_top_reversors.tsv` | Top 200 compuestos reversores (bottom 5% del score de conectividad) |

---

## Script 08 — `08_integrate_drug_targets.R`

Integración de 4 fuentes (DGIdb, ChEMBL, Open Targets, CMap2) en tabla maestra; clasificación A/B/C/D.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/drug_targets/08_drug_target_master_table.tsv` | Tabla larga maestra: todas las interacciones de las 4 fuentes con metadatos de genes y fármacos |
| T | `results/tables/drug_targets/08_drug_summary_per_drug.tsv` | Resumen por fármaco: n_sources, fuentes, max_phase, is_approved, cmap_score, genes diana, clase |
| T | `results/tables/drug_targets/08_multi_source_candidates.tsv` | Candidatos con soporte en >= 2 fuentes (o clase A/B por criterio clínico) |
| F | `results/figures/08_drug_class_barplot.pdf` | Barplot horizontal de fármacos por clase (A: HNSCC, B: Phase III+, C: Repurposing, D: Experimental) |
| F | `results/figures/08_cmap_vs_targets.pdf` | Scatter CMap scaled_score vs número de genes diana (solo candidatos con score CMap) |

---

## Script 09 — `09_string_network.R`

Red PPI con STRING API v12 + igraph; métricas de centralidad y hubs druggables.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/network/09_network_edges.tsv` | Aristas de la red PPI (gene_A, gene_B, score STRING) filtradas por score >= 700 |
| T | `results/tables/network/09_network_node_metrics.tsv` | Métricas por nodo: degree, betweenness normalizado, eigenvector, clustering_coeff, logFC, is_hub |
| T | `results/tables/network/09_druggable_hubs.tsv` | Hubs (top 10% por grado) con candidatos farmacológicos asociados desde script 08 |
| T | `results/tables/network/09_network_giant.graphml` | Componente mayor de la red en formato GraphML para importar en Cytoscape |
| F | `results/figures/09_network_degree_dist.pdf` | Histograma de distribución de grados con umbral de hub marcado |
| F | `results/figures/09_top25_hubs_degree.pdf` | Barplot horizontal top 25 hubs por grado, coloreado por dirección de expresión |
| F | `results/figures/09_hub_druggable_scatter.pdf` | Scatter grado vs betweenness de hubs; etiquetados los druggables |
| F | `results/figures/09_network_full.pdf` | Red completa del componente mayor (ggraph layout FR): nodos por dirección y tamaño por grado |
| F | `results/figures/09_network_hubs_only.pdf` | Subred de hubs únicamente: diamantes = hubs con candidato farmacológico |

---

## Script 10 — `10_prioritization_scoring.R`

Scoring multi-criterio (6 dimensiones) y selección del Top 20 candidatos con diversidad de targets.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/10_all_candidates_scored.tsv` | Todos los candidatos con s_logfc, s_sig, s_clinical, s_cmap, s_pathway, s_network, composite_score, final_score |
| T | `results/tables/10_top20_candidates.tsv` | Top 20 candidatos finales seleccionados con criterio de diversidad (max. 3 por target primario) |
| F | `results/figures/10_score_vs_targets.pdf` | Scatter score compuesto vs número de genes diana DE; etiquetados los Top 20 |

---

## Script 11 — `11_clinicaltrials_pubmed.py`

Consulta a ClinicalTrials.gov API v2 y PubMed para evidencia clínica de los Top 20.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/evidence/11_clinical_evidence.tsv` | Por candidato: n_ensayos HNSCC, n_ensayos activos, n_publicaciones PubMed HNSCC, evidence_score |
| T | `results/tables/evidence/11_trials_detail.tsv` | Detalle de cada ensayo clínico encontrado: título, fase, estado, NCT ID |
| F | `results/figures/evidence/11_evidence_bubble.pdf` | Bubble plot: candidatos por evidencia clínica y bibliográfica |

---

## Script 12 — `12_cosmic_overlap.py`

Solapamiento de genes diana DE con Cancer Gene Census (COSMIC) y NCG7.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/evidence/12_cancer_drivers_annot.tsv` | Todos los genes DE con anotación de driver (COSMIC tier, NCG7, hallmark) |
| T | `results/tables/evidence/12_cosmic_overlap.tsv` | Solo genes DE que solapan con el Cancer Gene Census de COSMIC |
| T | `results/tables/evidence/12_top20_with_drivers.tsv` | Top 20 candidatos con columnas `n_target_drivers` y `has_driver_target` |
| F | `results/figures/evidence/12_driver_overlap_bar.pdf` | Barplot de solapamiento genes DE con genes driver por categoría (COSMIC/NCG7) |
| F | `results/figures/evidence/12_driver_logfc_plot.pdf` | Scatter logFC de genes driver DE (cuando hay solapamiento) |

---

## Script 13 — `13_evidence_summary.R`

Integración final de toda la evidencia y generación del reporte maestro.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/13_FINAL_drug_candidates.xlsx` | Archivo maestro con 5 hojas: Top20_Final, Todos_Candidatos, Evidencia_Matriz, Ensayos_Clinicos, Metodologia |
| T | `results/tables/13_evidence_matrix.tsv` | Matriz binaria 20 candidatos × 10 dimensiones de evidencia (0/1 por dimensión) |
| F | `results/figures/final/13_evidence_radar.pdf` | Barplot grouped del perfil de scoring de los Top 5 candidatos por dimensión |

---

## Script 14 — `14_methods_summary.R`

Generación automática de documentación del pipeline (no produce resultados experimentales).

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `docs/METHODS.md` | Sección de métodos para manuscrito, generada automáticamente desde los resultados |
| T | `docs/OUTPUTS.md` | Catálogo de outputs (este archivo) |

---

## Script 15 — `15_sensitivity_analysis.R`

Análisis de sensibilidad de los pesos del scoring en 6 configuraciones distintas.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| T | `results/tables/15_sensitivity_ranks.tsv` | Rankings de todos los candidatos en 6 configuraciones de pesos, con n_configs_top20 |
| F | `results/figures/15_score_distribution.pdf` | Boxplot + jitter de la distribución de scores finales entre configuraciones (candidatos robustos) |

---

## Script 17 — `17_pub_figures.R`

Figuras de calidad publicación. Exporta PDF + PNG 300 DPI en `results/figures/pub/`.

Secciones: A = Proteómica QC/DE · B = Enriquecimiento · C = Bases de datos · D = Red PPI + Scoring · E = Evidencia · F = Sensibilidad · FIG = Figuras multipanel para manuscrito.

| Tipo | Archivo | Descripción |
| --- | --- | --- |
| F | `results/figures/pub/A1_volcano.pdf/.png` | Volcano plot publicación con paleta Okabe-Ito y etiquetas ggrepel |
| F | `results/figures/pub/A2_MA_plot.pdf/.png` | MA plot: logFC vs intensidad media por proteína |
| F | `results/figures/pub/A3_PCA.pdf/.png` | PCA biplot con líneas de pares por condición |
| F | `results/figures/pub/A4_heatmap_topDE.pdf/.png` | Heatmap top 30 proteínas DE (15 up + 15 down) × 20 muestras |
| F | `results/figures/pub/A_multipanel_QC_DE.pdf/.png` | Multipanel: Volcano + MA + PCA |
| F | `results/figures/pub/B1_hallmarks_barplot.pdf/.png` | Hallmarks GSEA en barplot horizontal ordenado por NES |
| F | `results/figures/pub/C1_drug_sources_bar.pdf/.png` | Número de candidatos detectados por cada base de datos (DGIdb/ChEMBL/OT/CMap) |
| F | `results/figures/pub/C2_drug_phase_dist.pdf/.png` | Distribución de fases clínicas de todos los candidatos |
| F | `results/figures/pub/D1_degree_vs_logFC.pdf/.png` | Scatter grado PPI vs logFC de genes DE |
| F | `results/figures/pub/D2_scoring_components_dot.pdf/.png` | Dot matrix: 6 componentes de score × Top 20 candidatos |
| F | `results/figures/pub/D3_top20_lollipop.pdf/.png` | Lollipop de ranking final mejorado para publicación |
| F | `results/figures/pub/E1_evidence_heatmap.pdf/.png` | Heatmap binario de evidencia (ComplexHeatmap) |
| F | `results/figures/pub/F1_bump_chart.pdf/.png` | Bump chart de estabilidad de ranks entre configuraciones |
| F | `results/figures/pub/F2_stability_bar.pdf/.png` | Barplot de robustez mejorado |
| F | `results/figures/pub/FIG1_QC_DE.pdf/.png` | Figura 1 manuscrito: Volcano + MA plot |
| F | `results/figures/pub/FIG2_pathways.pdf/.png` | Figura 2 manuscrito: Hallmarks GSEA |
| F | `results/figures/pub/FIG3_scoring.pdf/.png` | Figura 3 manuscrito: Lollipop Top 20 + Dot matrix de componentes |
| F | `results/figures/pub/FIG4_network.pdf/.png` | Figura 4 manuscrito: Scatter logFC vs grado PPI |
| F | `results/figures/pub/FIG5_sensitivity.pdf/.png` | Figura 5 manuscrito: Bump chart + Barplot de robustez |

---

## Archivo clave del proyecto

`results/tables/13_FINAL_drug_candidates.xlsx` — Resultado maestro con 5 hojas:

1. **Top20_Final**: Top 20 candidatos con toda la evidencia integrada y ranking final
2. **Todos_Candidatos**: Todos los candidatos evaluados con sus scores (~187 fármacos)
3. **Evidencia_Matriz**: Matriz binaria 10 dimensiones × 20 candidatos
4. **Ensayos_Clinicos**: Detalle de ensayos ClinicalTrials.gov por candidato
5. **Metodologia**: Resumen del pipeline con conteos de resultados por etapa

---

*Para regenerar todos los outputs, seguir `docs/RUNBOOK.md` en orden.*
