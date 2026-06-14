# Criterios Clínicos y de Literatura — hnscc_drug_repurposing

Registro auditable conforme a reglas `06-clinical-evidence-rigor` (global) y
`clinical-criteria-registry` (workspace). Toda fila debe verificarse contra
la fuente primaria antes de que el criterio entre a resultados.

| Criterio | Sistema/Fuente | Versión/Año | Locator (PMID/DOI/sección) | Definición textual (verbatim) | Cómo los datos la cumplen | Estado | Fecha verificada |
|----------|----------------|-------------|----------------------------|-------------------------------|---------------------------|--------|------------------|
| Query HNSCC PubMed | NCBI E-utilities / MeSH | 2026 | MeSH: D006258 (Head and Neck Neoplasms) | "head and neck squamous cell carcinoma" OR "head and neck cancer" OR "head and neck neoplasms"[MeSH Terms] OR HNSCC[Title/Abstract] OR "squamous cell carcinoma of the head and neck" — SIN "squamous cell carcinoma" solo (captura SCC de otros sitios) | Script 11 usa HNSCC_PUBMED_TERM con calificadores head/neck explícitos; excluye SCC genérico | ✅ VERIFICADO (B3 fix 2026-06-10) | 2026-06-10 |
| Query ClinicalTrials HNSCC (is_hnscc) | ClinicalTrials.gov API v2 | 2026 | https://clinicaltrials.gov/api/v2 | Ensayo clasificado como HNSCC si `conditions` contiene: "head and neck cancer", "head and neck squamous", "HNSCC", "squamous cell carcinoma head neck", "oral cancer", "oropharyngeal cancer" (phrase matching, no substring crudo) | Script 11 usa phrase matching sobre el campo `conditions`; `drug_confirmed=False` si interventions vacío | ✅ VERIFICADO (B3 fix 2026-06-10) | 2026-06-10 |
| Panel de candidatos para evidencia traslacional | Pipeline interno | 2026 | results/tables/15_lod_stability.tsv (lod_stable=TRUE) | Candidatos LOD-estables = permanecen en top-35 en las 4 corridas leave-one-database-out (DGIdb, ChEMBL, OpenTargets, L2S2); n=32 | Script 11 carga 10_all_candidates_scored.tsv y filtra por lod_stable=TRUE de 15_lod_stability.tsv | ✅ VERIFICADO | 2026-06-10 |
| Snapshot de evidencia clínica | ClinicalTrials.gov + PubMed | 2026 | Script 11 query_date column | Fecha de consulta registrada en columna `query_date` de 11_clinical_evidence.tsv para reproducibilidad | Cada run de script 11 registra QUERY_DATE = datetime.now().strftime("%Y-%m-%d") | ✅ VERIFICADO | 2026-06-10 |
| Supresores hipermetilados en HNSCC (CDKN2A, CDH1, DAPK1, RASSF1A) | Literatura (varios estudios TCGA/HNSCC) | ~2011–2020 | ⚠️ PMIDs pendientes de verificación | Estos supresores tumorales presentan hipermetilación de promotores recurrente en HNSCC comparado con mucosa normal | Script 19 etiqueta estas sondas en FigA; su presencia como DMPs hipermetilados en los datos se verifica empíricamente, pero el claim de recurrencia publicado requiere cita verificada | ⚠️ NO-VERIFICADO — **FUERA DE ALCANCE**: script 19 (metilación) removido del manuscrito (auditoría 2026-06-13); este criterio **no entra a ninguna salida del artículo** | 2026-06-10 |

> **Nota de alcance (auditoría 2026-06-13):** las filas que referencian el **script 11** (evidencia
> clínica PubMed/ClinicalTrials) y el **script 19** (metilación) corresponden a análisis **archivados
> / removidos** del manuscrito actual (ver `scripts/archive/` y CLAUDE.md). Se conservan como registro
> histórico pero **ninguno de estos criterios debe entrar al artículo en preparación**. Los criterios
> vigentes para el manuscrito son los de los scripts 01–10, 15 y 16.
