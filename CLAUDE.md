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

1. R: QC proteómica, ID mapping, enriquecimiento de vías GSEA (scripts 01–03)
2. Python: consulta DGIdb, ChEMBL, OpenTargets, L2S2 (scripts 04–07)
3. R: integración, red STRING, priorización hub-céntrica v3 (scripts 08–10)
4. R: robustez (LOD + configs de pesos) (script 15)
5. Validación externa en dos cohortes: CPTAC proteoma (16b py + 16c R) + TCGA RNA-seq (16 R)
6. R: figuras de publicación (17 + 17b/17c paneles → 17d–17i multipaneles) y tablas (18)

> Scripts `11`/`12`/`13` (evidencia clínica, COSMIC, reporte) **archivados** en
> `scripts/archive/` — fuera del manuscrito actual. Script 19 (metilación OE4) eliminado.

## Notas

- `data/raw/` — outputs MaxQuant (no modificar)
- Panel robusto = candidatos con `lod_stable = TRUE` en `results/tables/15_lod_stability.tsv`
- Priorización v3: composite = 0.60·TargetPriority + 0.40·DrugViability (pesos en `config/analysis_params.yaml`, sección `scoring_v3`)
