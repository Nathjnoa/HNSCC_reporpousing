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

## Afirmaciones de literatura — Discussion (Step 5, manuscrito)

Cada afirmación clínica/de literatura de la Discusión, con su fuente primaria verificada vía
PubMed (2026-06-15) y su citekey en `docs/manuscript/references.bib`. Estas son afirmaciones de
**soporte de literatura** (no sistemas multi-criterio); se registran para trazabilidad por la
regla `06-clinical-evidence-rigor`.

| Afirmación (Discussion) | Fuente | Año | Locator (PMID / DOI) | Citekey | Estado | Fecha verificada |
|-------------------------|--------|-----|----------------------|---------|--------|------------------|
| Beneficio anti-EGFR de agente único limitado por resistencia intrínseca/adquirida (RTK compensatorias, vías downstream, EMT); racional de combinación en HNSCC | Carosi et al., Crit Rev Oncol Hematol | 2026 | PMID 41690534 / 10.1016/j.critrevonc.2026.105207 | `carosi2026` | ✅ VERIFICADO | 2026-06-15 |
| Combinaciones EGFR + inmunoterapia y panorama ICI en HNSCC | Bhatia & Burtness, Drugs | 2023 | PMID 36645621 / 10.1007/s40265-023-01835-2 | `bhatia2023` | ✅ VERIFICADO | 2026-06-15 |
| Inhibición de DNMT (decitabina) + pembrolizumab aumenta TILs/PD-L1 y reduce MDSCs (sensibilización a ICI; ensayo fase 2 ventana) | Bear et al., J Immunother Cancer | 2025 | PMID 40021215 / 10.1136/jitc-2024-010294 | `bear2025` | ✅ VERIFICADO | 2026-06-15 |
| Inhibidores de proteasoma de nueva generación (carfilzomib/oprozomib) activos en modelos preclínicos de HNSCC vía UPR/ATF4 | Zang et al., Autophagy | 2012 | PMID 22995770 / 10.4161/auto.22185 | `zang2012` | ✅ VERIFICADO | 2026-06-15 |
| OXPHOS/Complejo I como vulnerabilidad metabólica targeteable; metformina inhibidor leve de Complejo I con efectos inmunometabólicos; ensayos de agente único poco efectivos | Pujalte-Martin et al., Mol Oncol | 2024 | PMID 38214418 / 10.1002/1878-0261.13583 | `pujaltemartin2024` | ✅ VERIFICADO | 2026-06-15 |
| Metformina en HNSCC: apoptosis tumoral + aumento de infiltrado CD8+ (ensayo clínico HPV+/−) | Curry et al., Front Oncol | 2018 | PMID 30364350 / 10.3389/fonc.2018.00436 | `curry2018` | ✅ VERIFICADO | 2026-06-15 |
| Metformina aumenta citotoxicidad de células NK en HNSCC vía inhibición de CXCL1 | Crist et al., J Immunother Cancer | 2022 | PMID 36328378 / 10.1136/jitc-2022-005632 | `crist2022` | ✅ VERIFICADO | 2026-06-15 |
| Metformina tolerable como quimio-radiosensibilizador en HNSCC (fase I/II) | Kemnade et al., Oral Oncol | 2023 | PMID 37562095 / 10.1016/j.oraloncology.2023.106536 | `kemnade2023` | ✅ VERIFICADO | 2026-06-15 |
| HNSCC resecable vira hacia inmunoterapia perioperatoria (pembrolizumab neoadyuvante/adyuvante, KEYNOTE-689) | Uppaluri et al., N Engl J Med | 2025 | DOI 10.1056/NEJMoa2415434 | `uppaluri2025` | ✅ VERIFICADO (ya en bib) | 2026-06-15 |
| Principios de network medicine (efecto terapéutico vía módulos funcionales) | Barabási et al., Nat Rev Genet | 2011 | DOI 10.1038/nrg2918 | `barabasi2011` | ✅ VERIFICADO (ya en bib) | 2026-06-15 |
| Reposicionamiento de fármacos: progreso y marco general | Pushpakom et al. 2019 / Tanoli et al. 2025 | 2019/2025 | 10.1038/nrd.2018.168 · 10.1038/s41573-025-01164-x | `pushpakom2019` · `tanoli2025` | ✅ VERIFICADO (ya en bib) | 2026-06-15 |
| Repurposing experimental de alto rendimiento en HNSCC (screen de PDC, 2248 compuestos; fedratinib, mitoxantrona) — contraste experimental al enfoque in silico | Gu et al., Sci Transl Med | 2022 | PMID 36070368 / 10.1126/scitranslmed.abo5987 | `gu2022` | ✅ VERIFICADO | 2026-06-15 |
| Repurposing in silico previo single-target/single-drug en HNSCC (transcriptoma → hub COL1A1 → docking → bleomicina) — vecino metodológico más cercano | Kumar et al., OMICS | 2025 | PMID 40662975 / 10.1177/15578100251359275 | `kumar2025` | ✅ VERIFICADO | 2026-06-15 |
| Framework network-medicine de repurposing in silico pan-cáncer (GPSnet; interactoma PPI + módulos de TCGA) — antecedente conceptual de la lógica de centralidad | Cheng et al., Nat Commun | 2019 | PMID 31375661 / 10.1038/s41467-019-10744-6 | `cheng2019` | ✅ VERIFICADO | 2026-06-15 |
| Revisión de métodos computacionales de predicción de terapia combinada en cáncer | Flanary et al., JCO Precis Oncol | 2023 | PMID 37824797 / 10.1200/PO.23.00261 | `flanary2023` | ✅ VERIFICADO | 2026-06-15 |
| Repurposing de agente único negativo en HNSCC avanzado (pantoprazol + terapia sistémica, fase I/II RCT, sin mejora de PFS) — respalda framing de combinación | Noronha et al., Med Oncol | 2023 | PMID 38129716 / 10.1007/s12032-023-02234-z | `noronha2023` | ✅ VERIFICADO | 2026-06-15 |
| Fenofibrato reposicionado reprograma el microambiente tumor-inmune en HNSCC HPV+ (prime/sensibilizar + cisplatino/checkpoint) | O'Neill et al., Cancers | 2022 | PMID 35053444 / 10.3390/cancers14020282 | `oneill2022` | ✅ VERIFICADO | 2026-06-15 |
| Estratificación proteogenómica de OSCC; tumores EGFR-high responden a cetuximab en PDX — corroboración independiente del control positivo EGFR | Wu et al., Cancer Lett | 2025 | PMID 39842500 / 10.1016/j.canlet.2025.217482 | `wu2025` | ✅ VERIFICADO | 2026-06-15 |
| scRNA-seq de LSCC; nominación bioinformática de anti-EGFR (erlotinib) contra genes de CSC — recuperación independiente del eje EGFR | Li et al., Genomics Proteomics Bioinformatics | 2024 | PMID 39107908 / 10.1093/gpbjnl/qzae056 | `li2024` | ✅ VERIFICADO | 2026-06-15 |
