# hnscc_drug_repurposing

Pipeline R+Python de drug repurposing para HNSCC basado en proteómica DIA (10 pares tumor/normal, MaxQuant).
Integra DGIdb, ChEMBL, OpenTargets, L2S2. Análisis completado 2026-03-08 · **En preparación para artículo internacional**.

Top candidatos: Erlotinib/Cetuximab (EGFR), Metformina (Complejo I OXPHOS).
Workspace artículo: `docs/manuscript/` · Tesis archivada: `docs/thesis/`

## Ambientes

```bash
conda activate omics-R   # scripts R
conda activate omics-py  # scripts Python
cd ~/bioinfo/projects/hnscc_drug_repurposing
```

## Ejecución

Ver `docs/RUNBOOK.md` para el orden completo y troubleshooting.

Resumen de fases (pipeline principal):

1. R: QC proteómica, ID mapping, enriquecimiento de vías (scripts 01–03)
2. Python: consulta DGIdb, ChEMBL, OpenTargets, L2S2 (scripts 04–07)
3. R: integración, red STRING, scoring (scripts 08–10)
4. R: sensibilidad LOD, figuras publicación, tablas (scripts 15, 17, 18)
5. Python (suplementario): evidencia clínica y COSMIC (scripts 11, 12)

## Notas

- `data/raw/` — outputs MaxQuant (no modificar)
- Panel final = candidatos con `lod_stable = TRUE` en `results/tables/15_lod_stability.tsv`
- Pesos de scoring configurables en `config/analysis_params.yaml`
