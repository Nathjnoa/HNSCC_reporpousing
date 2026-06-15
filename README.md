# HNSCC Drug Repurposing

Drug repurposing analysis for head and neck squamous cell carcinoma (HNSCC)
integrating quantitative DIA proteomics with four pharmacological databases,
PPI network analysis, multi-criteria scoring, and LOD stability validation.

**Estado**: Análisis completado (2026-03-08) · En preparación para artículo internacional (IMRAD, inglés)
**Ambientes**: `omics-R` (R scripts), `omics-py` (Python scripts)
**Datos**: Proteómica DIA · MaxQuant · 10 pares tumor/normal · 6 HPV+, 4 HPV−

Top candidatos LOD-stable: Erlotinib/Cetuximab (EGFR), Metformina (Complejo I OXPHOS).
Workspace del artículo: `docs/manuscript/` — Material de tesis archivado: `docs/thesis/`

---

## Estructura del proyecto

```text
hnscc_drug_repurposing/
├── data/
│   ├── raw/                 # Inmutables — no modificar
│   │   ├── proteomicacyc.csv
│   │   ├── metadata.csv
│   │   └── results_proteomica.tsv
│   ├── intermediate/id_mapping/
│   └── processed/
├── scripts/                 # Pipeline principal (01–03, 04–10, 15, 16b/16c/16, 17–18) + supp/ + archive/
├── config/
│   └── analysis_params.yaml # Todos los parámetros centralizados
├── results/
│   ├── figures/pub/{main,supp}/    # Figuras de publicación
│   ├── figures/pathway_enrichment/ # Enriquecimiento GO/KEGG/Reactome/Hallmarks
│   └── tables/pub/{main,supp}/     # Tablas de publicación
├── tests/                   # Tests de scripts Python
└── docs/
    ├── RUNBOOK.md
    ├── METHODS.md
    ├── manuscript/          # Workspace artículo internacional (IMRAD, en construcción)
    └── thesis/              # Tesis FUCS archivada (solo referencia)
```

---

## Pipeline

| Script | Fase | Ambiente | Descripción |
| --- | --- | --- | --- |
| `01_parse_results_qc.R` | 1 | omics-R | Parsear limma TVsS, exportar tablas DE (666 sig.) |
| `02_id_mapping.R` | 1 | omics-R | UniProt → Entrez/Symbol (org.Hs.eg.db) |
| `03_pathway_enrichment.R` | 1b | omics-R | GSEA (GO/KEGG/Reactome/Hallmarks); Hallmarks → panel GSEA de Fig2C |
| `04_query_dgidb.py` | 2 | omics-py | DGIdb GraphQL API v5 |
| `05_query_chembl.py` | 2 | omics-py | ChEMBL REST API, fase ≥ 3 |
| `06_query_opentargets.py` | 2 | omics-py | Open Targets GraphQL API v4 |
| `07_l2s2_connectivity.py` | 2 | omics-py | L2S2 reversión transcriptómica (1 044 drugs FDA) |
| `08_integrate_drug_targets.R` | 3 | omics-R | Unificar 4 fuentes · clasificar A/B/C/D |
| `09_string_network.R` | 4 | omics-R | Red PPI STRING v12 · módulos Louvain · hubs |
| `10_prioritization_scoring.R` | 4 | omics-R | Priorización hub-céntrica v3: composite = 0.60·TargetPriority + 0.40·DrugViability · anclaje por arista curada · tiers |
| `15_sensitivity_analysis.R` | 5 | omics-R | Robustez: LOD (leave-one-database) + 6 configs de pesos |
| `16b_cptac_fetch.py` | 5 | omics-py | Descarga proteómica CPTAC-HNSCC (paquete `cptac`) |
| `16c_cptac_de.R` | 5 | omics-R | DE proteómico CPTAC (limma pareado, mismo método que descubrimiento) |
| `16_external_validation.R` | 5 | omics-R | Validación 2 cohortes: concordancia CPTAC (Fig6A) + TCGA RNA-seq (Fig6B) + dianas (Fig6C) + KM supervivencia (supp) |
| `17 / 17b / 17c` | 6 | omics-R | Paneles Fig2/3/4 · Fig5 · GSEA Hallmarks (PDF + PNG 300 DPI, cache `.objects/`) |
| `17d–17i` | 6 | omics-R | Ensamblado multipanel Fig2/3/4/5/6 (TIFF 600 DPI LZW) |
| `17h_figS_robustness.R` | 6 | omics-R | Figura suplementaria de robustez (estabilidad × configs + LOD) |
| `18_pub_tables.R` | 6 | omics-R | Tablas de publicación Tab1–Tab6 (TSV) |
| `archive/{11,12,13}` | — | — | Evidencia clínica/COSMIC/reporte (archivado, fuera del manuscrito) |

Cadena de dependencias y comandos de ejecución: ver `docs/RUNBOOK.md`.

---

## Outputs de publicación

**Figuras multipanel** (`results/figures/pub/main/`, TIFF 600 DPI + PDF + PNG):

| Archivo | Contenido |
| --- | --- |
| `Fig1_workflow.*` | Workflow / diseño del estudio (esquema horizontal de 5 carriles; embudo 3,352→666→458) |
| `Fig2_multipanel.*` | (A) volcano · (B) Hallmarks GSEA · (C) heatmap top-40 DE |
| `Fig3_multipanel.*` | (A) fase clínica · (B) clase regulatoria · (C) UpSet solapamiento de BD |
| `Fig4_multipanel.*` | (A) red PPI por módulo Louvain · (B) enriquecimiento GO por módulo · (C) hubs druggables |
| `Fig5_multipanel.*` | (A) shortlist priorizado por módulo (TP/DV, faceta por tier) · (B) espacio TargetPriority × DrugViability |
| `Fig6_multipanel.*` | Validación 2 cohortes: (A) CPTAC proteoma · (B) TCGA RNA-seq · (C) dianas-ancla en ambas |

**Figuras suplementarias** (`results/figures/pub/supp/`):

| Archivo | Contenido |
| --- | --- |
| `FigS_robustness.*` | Estabilidad del ranking × 6 configs de peso + LOD |
| `FigS_selection_funnel.*` | Funnel de selección (candidatos → panel LOD-estable) |
| `FigS_survival_targets.*` | KM OS de 4 genes-pilar en TCGA-HNSC (no-significativo) |

**Tablas** (`results/tables/pub/`):

| Archivo | Contenido |
| --- | --- |
| `main/Tab1_resumen_DE.tsv` | Resumen estadístico del análisis DE |
| `main/Tab2_top_proteinas.tsv` | Top proteínas DE (up + down) |
| `main/Tab3_top_candidatos.tsv` | Top candidatos priorizados (composite) |
| `main/Tab4_EGFR_validation.tsv` | Eje EGFR (validación interna) |
| `main/Tab5_novel_candidates_by_module.tsv` | Candidatos novel por módulo |
| `main/Tab6_concordance_summary.tsv` | Concordancia externa (CPTAC + TCGA) |
| `supp/TabS1_extended_candidates_by_module.tsv` | Lista extendida de hubs-ancla por módulo |
| `supp/TabS2_survival_genes.tsv` | Genes para análisis de supervivencia |
| `supp/TabS3_candidatos_100_3fuentes.tsv` | Candidatos con ≥3 fuentes |
| `supp/TableS_exclusions.tsv` | Criterios de exclusión |

---

## Priorización hub-céntrica v3 (script 10)

```text
composite = 0.60 × TargetPriority + 0.40 × DrugViability

TargetPriority = 0.55 × centralidad en red + 0.45 × π direccional
DrugViability  = 0.40 × reversal L2S2 + 0.40 × clase regulatoria + 0.20 × breadth
```

Ancla = diana creíble (arista curada ChEMBL/OpenTargets) más central del módulo;
tiers `hub_central` / `peripheral_diff`. Pesos configurables en
`config/analysis_params.yaml` (sección `scoring_v3`).
Panel robusto = candidatos con `lod_stable = TRUE` en `results/tables/15_lod_stability.tsv`.
