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

1. R: QC proteómica, ID mapping (scripts 01–02)
2. Python: consulta DGIdb, ChEMBL, OpenTargets, L2S2 (scripts 04–07)
3. R: integración, red STRING, scoring (scripts 08–10)
4. R: sensibilidad LOD, figuras publicación, tablas (scripts 15, 17, 18)

## Notas

- `data/raw/` — outputs MaxQuant (no modificar)
- Panel final = candidatos con `lod_stable = TRUE` en `results/tables/15_lod_stability.tsv`
- Pesos de scoring configurables en `config/analysis_params.yaml`
