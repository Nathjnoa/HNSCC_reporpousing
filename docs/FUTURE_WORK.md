# Trabajo Futuro — HNSCC Drug Repurposing

Este documento registra los próximos pasos planificados para fortalecer el análisis hacia publicación.
El pipeline computacional (scripts 01-15) está completo. Lo que sigue son análisis de validación independiente.

---

## Nivel 2 — Tesis final / preprint

### Script 16: Validación pronóstica en TCGA-HNSC

**Archivo**: `scripts/16_tcga_survival.R` (pendiente)
**Ambiente**: `omics-R`
**Paquetes requeridos**: `TCGAbiolinks`, `survival`, `survminer`, `SummarizedExperiment`

**Objetivo**: Verificar de forma independiente que los targets de los candidatos top tienen impacto pronóstico en HNSCC, usando la cohorte TCGA-HNSC (~500 pacientes con RNA-seq + datos clínicos).

**Plan de análisis**:

1. Descargar TCGA-HNSC (GDC portal) via TCGAbiolinks:
   - RNA-seq HTSeq-FPKM normalizado
   - Datos clínicos: OS, OS.time, HPV status, estadio, sitio anatómico

2. Estratificar pacientes por expresión alta/baja de cada target (mediana como punto de corte):
   - `EGFR` → Erlotinib/Cetuximab/Lapatinib
   - `NDUFA13` (representativo de Complejo I) → Metformin
   - `MMP8` → Doxycycline/Collagenase/Marimastat
   - `MYH7` → Mavacamten (¿tiene impacto pronóstico?)

3. Para cada target:
   - Kaplan-Meier (log-rank test)
   - Cox univariado (HR, 95% CI, p-valor)
   - Cox multivariado ajustando por: HPV status, estadio (I-II vs III-IV), sitio (orofaringe vs otros)

4. Pregunta clave: ¿targets de candidatos con alta evidencia clínica (Metformin, Doxycycline) también tienen impacto pronóstico independiente?

**Outputs esperados**:
- `results/figures/tcga/16_km_EGFR.pdf` — Kaplan-Meier EGFR
- `results/figures/tcga/16_km_NDUFA13.pdf` — Kaplan-Meier Complejo I
- `results/figures/tcga/16_km_MMP8.pdf` — Kaplan-Meier MMP8
- `results/tables/16_cox_univariado.tsv` — resultados Cox univariado
- `results/tables/16_cox_multivariado.tsv` — HR ajustados

**Nota**: Si NDUFA13 (downregulado en tumor, logFC=-6) tiene HR > 1 (alta expresión = mejor pronóstico), refuerza directamente la lógica de reactivar el Complejo I con Metformin.

---

### Script 17: Figuras de calidad publicación ✅ COMPLETADO (2026-03-04)

**Archivo**: `scripts/17_pub_figures.R` (COMPLETADO)
**Ambiente**: `omics-R`
**Skill usado**: `/pub-figures` (preset `double_col`, 180×120 mm, Okabe-Ito)

**19 PDFs + 19 PNGs en `results/figures/pub/`**:

| Código | Figura | Nueva? |
| ------ | ------ | ------ |
| A1 | Volcano plot mejorado (labels top 20) | Mejorada |
| A2 | MA plot (logFC vs intensidad media) | **NUEVA** |
| A3 | PCA biplot con líneas de pares | Mejorada |
| A4 | Heatmap top 30 DE × 20 muestras | **NUEVA** |
| B1 | Hallmarks GSEA barplot horizontal | Mejorada |
| C1 | Drug candidates por fuente DB | Mejorada |
| C2 | Distribución fases clínicas | **NUEVA** |
| D1 | Scatter logFC vs grado PPI | **NUEVA** |
| D2 | Dot matrix 6 componentes × Top 20 | **NUEVA** |
| D3 | Top 20 lollipop (por clase farmacológica) | Mejorada |
| E1 | Evidence heatmap (ComplexHeatmap, categorías) | Mejorada |
| F1 | Bump chart estabilidad de ranks | **NUEVA** |
| F2 | Stability bar por candidato | Mejorada |
| FIG1-5 | Multipanel listos para manuscrito | — |

**Pendiente (Fig 6 manuscrito)**: Kaplan-Meier → ver Script 16 abajo.

---

## Nivel 3 — Manuscrito formal

### Sección de limitaciones (para Discussion)

Puntos a incluir en la discusión del manuscrito:

1. **Proteómica vs transcriptómica**: L2S2 (LINCS L1000) usa firmas transcriptómicas de líneas celulares. La correlación proteoma-transcriptoma es ~0.6 en promedio; los resultados de reversión transcriptómica son sugestivos pero no concluyentes a nivel proteómico.

2. **Tamaño de muestra**: n=10 pares pareados es estadísticamente robusto para detectar logFC>1 con FDR<5%, pero la generalización requiere validación en cohortes más grandes.

3. **Ausencia de validación experimental**: El análisis es in silico. El siguiente paso lógico es validar los candidatos emergentes (Metformin, Doxycycline) en líneas celulares HNSCC (Cal-27, FaDu, SCC-25) con ensayos de viabilidad (MTT), migración, e invasión.

4. **Un solo centro**: Los datos proteómicos provienen de una sola cohorte. La validación en TCGA (script 16) proporciona evidencia independiente a nivel transcriptómico.

5. **Scoring subjetivo**: Los pesos del scoring compuesto son razonables pero arbitrarios. El análisis de sensibilidad (script 15) muestra que el top 5 es robusto, pero los candidatos en posiciones 10-20 son más dependientes de los pesos.

### Comparación con literatura existente

**Tarea pendiente** — buscar y sintetizar papers de drug repurposing en HNSCC:

Búsqueda sugerida en PubMed:
```
("head and neck" OR "HNSCC") AND ("drug repurposing" OR "drug repositioning") AND ("proteomic" OR "transcriptomic") AND (2015:2026[dp])
```

Preguntas a responder:
- ¿Erlotinib/Cetuximab aparecen en estudios similares? (esperado: sí)
- ¿Metformin ha sido reportado en HNSCC computacional? (buscar evidencia)
- ¿Doxycycline como candidato de reposicionamiento en HNSCC? (evaluar novedad)
- ¿Hay candidatos en nuestro top 20 completamente no reportados? (potencial novedad)

Papers clave a revisar:
- Estudios de drug repurposing con datos TCGA-HNSC (transcriptómica)
- Metformina en HNSCC: epidemiología y mecanismo (varios grupos han publicado)
- Inhibidores de MMP en cáncer de cabeza y cuello (ECM como target terapéutico)

---

## Cronograma sugerido

| Etapa | Tarea | Tiempo estimado |
| ----- | ----- | --------------- |
| Tesis | Script 16 (TCGA survival) | 1-2 semanas |
| Tesis | Script 17 (pub figures) | 1 semana |
| Tesis | Redacción capítulo resultados | 2-3 semanas |
| Publicación | Sección limitaciones | 2-3 días |
| Publicación | Búsqueda y síntesis literatura | 1 semana |
| Publicación | Revisión por co-autores | variable |

---

Creado: 2026-03-04 | Actualizar cuando se complete cada etapa
