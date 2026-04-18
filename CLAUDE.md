# hnscc_drug_repurposing

Pipeline R+Python de drug repurposing para HNSCC basado en proteómica DIA (10 pares tumor/normal, MaxQuant).
Integra DGIdb, ChEMBL, OpenTargets, L2S2. **Análisis completado: 2026-03-08.**

Top candidatos: Erlotinib/Cetuximab (EGFR), Metformina (Complejo I OXPHOS).

## Ambientes

```bash
conda activate omics-R   # scripts R
conda activate omics-py  # scripts Python
cd ~/bioinfo/projects/hnscc_drug_repurposing
```

## Ejecución

Ver `docs/RUNBOOK.md` para el orden completo y troubleshooting.

Resumen de fases:

1. R: QC proteómica, ID mapping, pathway enrichment (scripts 01-03)
2. Python: consulta DGIdb, ChEMBL, OpenTargets, L2S2 (scripts 04-07)
3. R: integración, red STRING, scoring (scripts 08-10)
4. Python: ClinicalTrials + PubMed, COSMIC (scripts 11-12)
5. R: resumen, figuras publicación (scripts 13-15, 17)

## Notas

- `scripts/07_cmap_connectivity_DEPRECATED.R` — no usar (reemplazado por L2S2)
- `scripts/proteomicacyc.R` — exploración inicial, no parte del pipeline
- `data/raw/` — outputs MaxQuant (no modificar)
