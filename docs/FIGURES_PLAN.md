# Plan de Figuras y Tablas — HNSCC Drug Repurposing
**Destino primario:** Tesis de grado  
**Destino final:** Artículo científico  
**Última actualización:** 2026-04-11  
**Estado:** En construcción — Sección 0 y OE1 definidas. OE2 y OE3 pendientes.

---

## Convención de nomenclatura y directorios

Todos los archivos finales (tesis/artículo) viven en `results/figures/pub/`.
Los archivos de exploración intermedia permanecen en sus carpetas originales (`qc/`, `pathway_enrichment/`, `figures/`) — no se tocan.

```
results/figures/pub/
├── main/          # Figuras principales (van en el cuerpo del texto)
│   ├── Sec0_FigA_[nombre].pdf
│   ├── Sec0_FigB_[nombre].pdf
│   ├── OE1_FigA_[nombre].pdf
│   └── ...
└── supp/          # Figuras y tablas suplementarias
    ├── Sec0_S1_[nombre].pdf
    ├── OE1_S1_[nombre].pdf
    └── ...
```

**Archivos actuales en `pub/`** (A1–F2) se renombran siguiendo este esquema al momento de ensamblar la tesis.  
Mientras tanto el plan usa el nombre actual como referencia.

---

## Objetivo General

> Evaluar *in silico* el potencial de reposicionamiento de fármacos existentes para el tratamiento
> del cáncer escamocelular de cabeza y cuello, integrando datos proteómicos, análisis bioinformáticos
> y evidencia farmacológica que permita identificar candidatos con relevancia terapéutica.

---

## Sección 0 — Caracterización del perfil proteómico diferencial de HNSCC

*Base metodológica. En la tesis va como primera sección de Resultados, antes de los objetivos
específicos. Responde: ¿Cuál es el estado del proteoma en HNSCC comparado con tejido normal?
¿Qué procesos biológicos están alterados?*

---

### Figuras principales (`pub/main/Sec0_Fig*.pdf`)

| ID propuesto | Archivo actual | Descripción | Notas / Mejoras |
|---|---|---|---|
| `Sec0_FigA` | `pub/A3_PCA.pdf` | PCA de las 20 muestras. Separación tumor/normal en PC1; variabilidad HPV en PC2. | Añadir forma de punto para HPV+/−. Verificar etiquetas de ejes con varianza explicada. |
| `Sec0_FigB` | `pub/A1_volcano.pdf` | Volcano plot TVsS: 666 proteínas DE (329 up / 337 down). Umbral: adj.P < 0.05, \|log2FC\| > 1. | Etiquetar proteínas clave (EGFR, TOP2A, PSMB5 u otras del proteosoma). |
| `Sec0_FigC` | `pub/A4_heatmap_topDE.pdf` | Heatmap top proteínas DE. Patrón expresión tumor vs. normal con anotación HPV. | Confirmar barra de anotación HPV+/−. Verificar orden de clusters. |
| `Sec0_FigD` | `pub/B1_hallmarks_barplot.pdf` | GSEA Hallmarks MSigDB. Identifica los procesos biológicos activados/reprimidos y justifica la estrategia terapéutica (ej.: proliferación, OXPHOS, EMT). | Ordenar por NES. Paleta: rojo = activado, azul = reprimido. Agregar línea de corte FDR. |

> **Nota:** B1 (Hallmarks) va en Sec0 porque el enriquecimiento hace parte de la caracterización
> biológica previa, no de la búsqueda de fármacos. Es el puente entre "qué está alterado" y
> "qué fármacos buscar".

---

### Figuras suplementarias (`pub/supp/Sec0_S*.pdf`)

| ID propuesto | Archivo actual | Descripción |
|---|---|---|
| `Sec0_S1` | `pub/A2_MA_plot.pdf` | MA plot (intensidad media vs. log2FC). Confirma ausencia de sesgo por abundancia. |
| `Sec0_S2` | `qc/01_boxplot_intensidades.pdf` | Boxplot de intensidades normalizadas. Valida éxito de normalización Cyclic Loess. |
| `Sec0_S3` | `qc/01_missingness_vs_direction.pdf` | Datos faltantes vs. dirección del cambio. Evidencia de missing-not-at-random. |
| `Sec0_S4` | `pathway_enrichment/03_GO_BP_dotplot.pdf` | ORA GO Biological Process (19 términos simplificados). |
| `Sec0_S5` | `pathway_enrichment/03_KEGG_barplot.pdf` | ORA KEGG (16 rutas). |
| `Sec0_S6` | `pathway_enrichment/03_Reactome_dotplot.pdf` | ORA Reactome (14 rutas). |

---

### Tablas — Sección 0

| ID | Tipo | Contenido | Fuente de datos | Estado |
|---|---|---|---|---|
| `Sec0_Tab1` | **Principal** | Resumen del análisis DE: n proteínas cuantificadas, n significativas, n up, n down, umbrales usados, coeficiente de variación inter-muestra. | Script 01 | Existe en logs/texto — **falta formalizar como tabla** |
| `Sec0_Tab2` | **Principal** | Top 20–30 proteínas más sobreexpresadas y 20–30 más reprimidas (Gene Symbol, log2FC, adj.P, función anotada). | `results/tables/de_limma/01_TVsS_*.tsv` | Existe — **falta seleccionar y formatear** |
| `Sec0_TabS1` | Suplementaria | Tabla completa de las 666 proteínas DE con todos los estadísticos. | `results/tables/de_limma/` | Existe — dar formato final |
| `Sec0_TabS2` | Suplementaria | Tabla de enriquecimiento Hallmarks con NES, p-valor, FDR y genes líderes (leading edge). | `results/tables/` (script 03) | Existe — revisar si está guardada |

---

### Figuras a eliminar (Sección 0)

| Archivo | Motivo |
|---|---|
| `qc/01_resumen_DE.pdf` | Redundante — los números están en Sec0_Tab1 y en el volcano. |
| `pathway_enrichment/03_GO_MF_dotplot.pdf` | GO Molecular Function — menos informativa para historia de drug repurposing. |
| `pathway_enrichment/03_GO_CC_dotplot.pdf` | GO Cellular Component — descriptiva, no suma al argumento terapéutico. |
| `pathway_enrichment/03_GO_BP_emapplot.pdf` | Demasiado complejo; no agrega sobre el dotplot. |
| `pathway_enrichment/03_GO_BP_cnetplot.pdf` | Idem. |
| `pathway_enrichment/03_GO_BP_GSEA_dotplot.pdf` | Cubierta por B1 (Hallmarks). |
| `pathway_enrichment/03_KEGG_GSEA_dotplot.pdf` | Idem. |
| `pathway_enrichment/03_Reactome_GSEA_dotplot.pdf` | Idem. |
| `pathway_enrichment/03_Hallmarks_GSEA_ridgeplot.pdf` | El barplot (B1) comunica mejor para manuscrito. |

---

## OE1 — Relacionar proteínas DE con fármacos reportados en bases de datos públicas

> Relacionar proteínas expresadas diferencialmente en cáncer escamocelular de cabeza y cuello con
> fármacos reportados en bases de datos públicas, utilizando herramientas bioinformáticas de
> reposicionamiento y análisis de conectividad.

*Responde: ¿Qué fármacos existentes tienen como blanco las proteínas DE en HNSCC?
¿Cuántos candidatos encontró cada base de datos? ¿Cuántos tienen respaldo clínico?*

---

### Figuras principales (`pub/main/OE1_Fig*.pdf`)

| ID propuesto | Archivo actual | Descripción | Notas / Mejoras |
|---|---|---|---|
| `OE1_FigA` | `pub/C1_drug_sources_bar.pdf` | Candidatos identificados por cada fuente (DGIdb, ChEMBL, OpenTargets, L2S2). Muestra la estrategia multi-fuente. | Agregar encima de cada barra el n° de genes consultados. Unificar paleta con OE1_FigB. |
| `OE1_FigB` | `pub/C2_drug_phase_dist.pdf` | Distribución de fases clínicas. Muestra que la mayoría de candidatos son fármacos aprobados (Fase IV). | Separar o colorear por fuente de datos si es posible. |
| `OE1_FigC` | **CREAR** | UpSet plot o diagrama de overlap entre las 4 fuentes de datos. Muestra qué candidatos tienen respaldo de múltiples fuentes (evidencia convergente). Candidatos en múltiples fuentes → mayor confianza. | Usar `UpSetR` en R. Datos disponibles en `results/tables/` (script 08). **Esta figura falta y es importante para OE1.** |

---

### Figuras suplementarias (`pub/supp/OE1_S*.pdf`)

| ID propuesto | Archivo actual | Descripción |
|---|---|---|
| `OE1_S1` | `figures/08_drug_sources_support.pdf` | Número de fuentes que respaldan cada candidato (versión alternativa a C1). Útil como validación de OE1_FigC. |

---

### Tablas — OE1

| ID | Tipo | Contenido | Fuente de datos | Estado |
|---|---|---|---|---|
| `OE1_Tab1` | **Principal** | Resumen de consulta a bases de datos: n genes DE consultados, n genes con hit, n fármacos únicos, n fármacos Fase ≥ III — una fila por base de datos (DGIdb, ChEMBL, OpenTargets, L2S2). | Scripts 04–07 | Datos existen — **falta construir tabla resumen** |
| `OE1_Tab2` | **Principal** | Top candidatos por número de fuentes que los respaldan: fármaco, mecanismo, blanco(s), fases clínicas, n° fuentes. Top 15–20. | Script 08 (`results/tables/`) | Probablemente existe — **revisar y formatear** |
| `OE1_TabS1` | Suplementaria | Tabla completa de pares gen-fármaco de DGIdb (2,846 pares, 226 genes, 2,252 fármacos). | Script 04 | Existe como TSV — dar formato |
| `OE1_TabS2` | Suplementaria | Tabla de fármacos ChEMBL Fase ≥ 3 (90 pares, 21 genes). | Script 05 | Existe como TSV |
| `OE1_TabS3` | Suplementaria | Tabla de candidatos Open Targets con association score HNSCC (354 genes DE con evidencia). | Script 06 | Existe como TSV |

---

### Figuras a eliminar (OE1)

| Archivo | Motivo |
|---|---|
| `figures/08_cmap_vs_targets.pdf` | Análisis CMap deprecado (reemplazado por L2S2). No aplica. |
| `figures/08_drug_class_barplot.pdf` | Cubierto por OE1_FigA (C1) y OE1_Tab1. |
| `figures/08_top_candidates_lollipop.pdf` | Versión preliminar del lollipop — reemplazado por D3 en OE2. |

---

## Pendiente — OE2 y OE3

*(Se planificará en la siguiente iteración)*

- **OE2:** Priorizar candidatos por conectividad, rutas y redes → figuras D1, D2, D3 + red STRING
- **OE3:** Contrastar clínica y bibliográficamente → figura E1, bump chart, stability bar

---

## Figuras confirmadas para eliminar (global)

Las siguientes se pueden borrar del directorio sin pérdida de información:

```
pathway_enrichment/03_GO_MF_dotplot.pdf
pathway_enrichment/03_GO_CC_dotplot.pdf
pathway_enrichment/03_GO_BP_emapplot.pdf
pathway_enrichment/03_GO_BP_cnetplot.pdf
pathway_enrichment/03_GO_BP_GSEA_dotplot.pdf
pathway_enrichment/03_KEGG_GSEA_dotplot.pdf
pathway_enrichment/03_Reactome_GSEA_dotplot.pdf
pathway_enrichment/03_Hallmarks_GSEA_ridgeplot.pdf
qc/01_resumen_DE.pdf
figures/08_cmap_vs_targets.pdf
figures/08_drug_class_barplot.pdf
figures/08_top_candidates_lollipop.pdf
```
