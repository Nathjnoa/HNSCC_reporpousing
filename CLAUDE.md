# hnscc_drug_repurposing

Pipeline R+Python de drug repurposing para cáncer de cabeza y cuello (HNSCC) basado
en proteómica cuantitativa DIA (10 pares tumor/normal, MaxQuant). Integra múltiples
bases de datos (DGIdb, ChEMBL, OpenTargets, L2S2) y genera ranking de candidatos.

Top candidatos: Erlotinib/Cetuximab (EGFR), Metformina (Complejo I OXPHOS).

## Ambientes

| Fase | Ambiente |
|------|----------|
| Scripts R (.R) | `conda activate omics-R` |
| Scripts Python (.py) | `conda activate omics-py` |

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
```

## Ejecución (orden de scripts)

```bash
# Fase 1: Datos
conda activate omics-R
Rscript scripts/01_parse_results_qc.R        # QC proteómica MaxQuant
Rscript scripts/02_id_mapping.R              # UniProt → Ensembl → Entrez
Rscript scripts/03_pathway_enrichment.R      # GO/KEGG/Reactome/Hallmarks

# Fase 2: Consulta de drogas
conda activate omics-py
python scripts/04_query_dgidb.py             # DGIdb
python scripts/05_query_chembl.py            # ChEMBL IC50/Ki
python scripts/06_query_opentargets.py       # OpenTargets scores
python scripts/07_l2s2_connectivity.py       # L2S2 conectividad

# Fase 3: Integración y red
conda activate omics-R
Rscript scripts/08_integrate_drug_targets.R  # Integrar fuentes
Rscript scripts/09_string_network.R          # Red STRING
Rscript scripts/10_prioritization_scoring.R  # Scoring multi-criterio

# Fase 4: Evidencia clínica
conda activate omics-py
python scripts/11_clinicaltrials_pubmed.py   # ClinicalTrials + PubMed
python scripts/12_cosmic_overlap.py          # COSMIC cancer genes

# Fase 5: Resumen y figuras
conda activate omics-R
Rscript scripts/13_evidence_summary.R        # Tabla resumen
Rscript scripts/14_methods_summary.R         # Métodos
Rscript scripts/15_sensitivity_analysis.R    # Análisis de sensibilidad
Rscript scripts/17_pub_figures.R             # Figuras publicación
```

## Rutas clave
- `data/raw/` — outputs MaxQuant (no modificar)
- `results/tables/` — tablas intermedias y finales
- `results/figures/` — figuras
- `docs/RUNBOOK.md` — guía completa con troubleshooting
- `docs/METHODS.md` — métodos para manuscrito

## Notas
- `scripts/07_cmap_connectivity_DEPRECATED.R` — no usar (reemplazado por L2S2)
- `scripts/proteomicacyc.R` — script exploración inicial, no parte del pipeline
- Análisis completado: 2026-03-08
