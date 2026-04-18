# HNSCC Drug Repurposing

Drug repurposing analysis for head and neck squamous cell carcinoma (HNSCC)
integrating quantitative DIA proteomics with four pharmacological databases,
PPI network analysis, multi-criteria scoring, and LOD stability validation.

**Ambientes**: `omics-R` (R scripts), `omics-py` (Python scripts)
**Datos**: Proteómica DIA · MaxQuant · 10 pares tumor/normal · 6 HPV+, 4 HPV−

Top candidatos LOD-stable: Erlotinib/Cetuximab (EGFR), Metformina (Complejo I OXPHOS).

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
├── scripts/                 # 12 scripts numerados (ver pipeline)
├── config/
│   └── analysis_params.yaml # Todos los parámetros centralizados
├── results/
│   ├── figures/pub/main/    # 7 figuras de publicación
│   └── tables/pub/          # 5 tablas main + 1 supp
└── docs/
    ├── RUNBOOK.md
    ├── METHODS.md
    └── HNSCC_DrugRepurposing_Figuras_reviewed.docx
```

---

## Pipeline

| Script | Fase | Ambiente | Descripción |
| --- | --- | --- | --- |
| `01_parse_results_qc.R` | 1 | omics-R | Parsear limma TVsS, exportar tablas DE (666 sig.) |
| `02_id_mapping.R` | 1 | omics-R | UniProt → Entrez/Symbol (org.Hs.eg.db) |
| `04_query_dgidb.py` | 2 | omics-py | DGIdb GraphQL API v5 |
| `05_query_chembl.py` | 2 | omics-py | ChEMBL REST API, fase ≥ 3 |
| `06_query_opentargets.py` | 2 | omics-py | Open Targets GraphQL API v4 |
| `07_l2s2_connectivity.py` | 2 | omics-py | L2S2 reversión transcriptómica (1 044 drugs FDA) |
| `08_integrate_drug_targets.R` | 3 | omics-R | Unificar 4 fuentes · clasificar A/B/C/D |
| `09_string_network.R` | 4 | omics-R | Red PPI STRING v12 · módulos Louvain · hubs |
| `10_prioritization_scoring.R` | 4 | omics-R | Scoring multi-criterio 5D · pool top 35 |
| `15_sensitivity_analysis.R` | 5 | omics-R | LOD + 6 configs pesos + permutation test |
| `17_pub_figures.R` | 6 | omics-R | 7 figuras de publicación (PDF + PNG 300 DPI) |
| `18_pub_tables.R` | 6 | omics-R | 6 tablas de publicación (TSV) |

Cadena de dependencias y comandos de ejecución: ver `docs/RUNBOOK.md`.

---

## Outputs de publicación

**Figuras** (`results/figures/pub/main/`):

| Archivo | Contenido |
| --- | --- |
| `Sec0_FigB_volcano.png` | Volcano plot DE proteómica |
| `Sec0_FigC_heatmap_topDE.png` | Heatmap top 40 proteínas DE |
| `OE1_FigA_drug_sources_bar.png` | Candidatos por número de fuentes |
| `OE1_FigB_drug_phase_dist.png` | Distribución de fases clínicas |
| `OE2_FigA_ppi_network.png` | Red PPI con módulos y hubs |
| `OE2_FigB_module_barplot.png` | Candidatos por módulo de red |
| `OE2_FigC_class_distribution.png` | Distribución de clases A/B/C/D |

**Tablas** (`results/tables/pub/`):

| Archivo | Contenido |
| --- | --- |
| `main/Sec0_Tab1_resumen_DE.tsv` | Resumen estadístico del análisis DE |
| `main/Sec0_Tab2_top_proteinas.tsv` | Top 40 proteínas DE (20 up + 20 down) |
| `main/OE1_Tab2_top_candidatos.tsv` | Top 30 candidatos por número de fuentes |
| `main/OE2_Tab1_EGFR_LOD_stable.tsv` | Candidatos EGFR LOD-stable |
| `main/OE2_Tab2_noEGFR_LOD_stable.tsv` | Candidatos no-EGFR LOD-stable |
| `supp/OE2_TabS1_candidatos_extendidos_noEGFR.tsv` | Extendidos no-EGFR (robustos a pesos) |

---

## Scoring multi-criterio (script 10)

```text
composite = 0.325 × π-estadístico (logFC × -log10 adj.P)
          + 0.200 × fase clínica
          + 0.195 × centralidad en red (betweenness)
          + 0.150 × relevancia de ruta
          + 0.130 × L2S2 reversal score
```

Pesos configurables en `config/analysis_params.yaml`.
Panel final = candidatos con `lod_stable = TRUE` en `results/tables/15_lod_stability.tsv` (32 candidatos).
