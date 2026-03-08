# Reporte de Progreso — HNSCC Drug Repurposing
**Actualización:** 2026-03-07
**Estado:** Pipeline completado y validado (16 scripts: 01-15 + 17) — re-ejecutado 2026-03-07

### Correcciones aplicadas (2026-03-07)

- **Bug crítico (script 07):** Dataset CMap2 incorrecto — se usaba `EH3224` (cmap_expr) en lugar de `EH3225` (cmap_rank). Con el dataset incorrecto todos los scores de conectividad resultaban 0. Corregido; el análisis ahora produce scores en rango [-1, 1] con 180 reversores potentes (bottom 5%).
- **Paletas de color (scripts 08, 10):** Etiquetas en `scale_color_manual`/`scale_fill_manual` no coincidían con los valores reales de `drug_class_label`. Corregidas en 4 figuras (lollipop, CMap scatter, barplot top20, scatter scoring).

---

## ¿De qué trata este proyecto?

Este proyecto busca identificar fármacos que ya existen y podrían usarse para tratar el **carcinoma escamocelular de cabeza y cuello (HNSCC)** mediante un enfoque computacional llamado **reposicionamiento de fármacos** (*drug repurposing*).

En lugar de desarrollar un fármaco desde cero (proceso que tarda 10-15 años y cuesta miles de millones de dólares), la idea es analizar qué proteínas están alteradas en el tumor y buscar si algún fármaco ya aprobado para otra enfermedad podría actuar sobre esas proteínas.

---

## Los datos: ¿de dónde partimos?

Tenemos datos de **proteómica cuantitativa** de 20 pacientes con HNSCC: 10 pares de tejido tumor/normal del mismo paciente. Esto es importante porque al comparar tumor vs. normal del mismo paciente, eliminamos la variabilidad genética individual.

**Tecnología usada:** DIA (Data-Independent Acquisition) con espectrometría de masas, procesado con MaxQuant y el paquete R `proteoDA` (un wrapper de `limma`).

**Variable adicional importante:** 6 pacientes son HPV+ (VPH positivo) y 4 son HPV−. El VPH es crítico en HNSCC porque los tumores HPV+ tienen mejor pronóstico y una biología molecular diferente. Nuestro análisis principal (TVsS = Tumor vs. Sano) promedia sobre ambos grupos HPV, que es el enfoque correcto para encontrar candidatos terapéuticos generales.

---

## Fase 1 — Preparar los datos (Scripts 01-02)

### Script 01: Parsear y hacer control de calidad

El archivo de resultados `results_proteomica.tsv` tiene un formato especial (3 filas de encabezados fusionados al estilo Excel, luego los datos). El script 01 lo "traduce" a un formato limpio y genera figuras de QC:

- **Boxplot de intensidades**: verifica que la normalización Cyclic Loess funcionó correctamente (las distribuciones de log2-intensidad deben ser similares entre muestras).
- **PCA**: visualiza si las muestras se agrupan por condición (tumor vs. normal) o si hay otros efectos. Esperamos que tumor y normal se separen en el PC1.
- **Volcano plot**: muestra las 520 proteínas significativamente alteradas (adj.P.Val < 0.05 y |log2FC| > 1).
- **Resumen DE**: 248 proteínas up-reguladas + 272 down-reguladas = **520 proteínas DE totales**.

*¿Por qué 520 de 3,352?* No todas las proteínas cuantificadas tienen cambio estadísticamente significativo. Usar un umbral de |log2FC| > 1 (cambio al menos 2× en escala lineal) y FDR < 5% garantiza que solo reportamos cambios robustos.

### Script 02: Mapear identificadores de genes

Las proteínas vienen identificadas con **UniProt IDs** (P00533, Q9Y4L1, etc.). La mayoría de herramientas de enriquecimiento necesitan **Entrez IDs** (números enteros) o **Gene Symbols** (EGFR, TP53...).

Se usó `org.Hs.eg.db` (base de datos de anotación del genoma humano en Bioconductor) para el mapeo:

- 3,263 / 3,352 proteínas mapeadas a Entrez (97.3%)
- 518 / 520 proteínas significativas con Entrez ID
- Las 2 sin mapeo (ALDH9A1, HNRNPAB) se retienen como NA y se excluyen donde sea necesario

*¿Por qué no mapear todo directamente con UniProt?* Las herramientas de enriquecimiento (clusterProfiler, ReactomePA) están optimizadas para Entrez/Symbol. El mapeo correcto garantiza resultados reproducibles.

---

## Fase 2 — Entender la biología alterada (Script 03)

Con la lista de 520 proteínas DE, la pregunta es: **¿qué procesos biológicos están alterados en el tumor?**

### Over-Representation Analysis (ORA)

Pregunta: *¿Están sobre-representados en nuestra lista genes que pertenecen a cierto proceso biológico?*

Método estadístico: test hipergeométrico (similar a la prueba exacta de Fisher). Compara la proporción de genes DE en un pathway vs. la proporción esperada por azar, usando todas las proteínas cuantificadas como universo de fondo.

**Resultados ORA:**
| Base de datos | Términos significativos |
|---------------|------------------------|
| GO Biological Process | 77 (19 tras simplificación) |
| GO Molecular Function | 33 |
| GO Cellular Component | 55 |
| KEGG | 16 rutas |
| Reactome | 14 rutas |

*¿Por qué simplificar GO BP?* La ontología Gene Ontology tiene muchos términos redundantes (ej: "metabolic process", "organic substance metabolic process", "cellular metabolic process" son casi lo mismo). La función `simplify()` con un umbral de similitud semántica de 0.70 elimina los términos más redundantes, dejando los 19 más informativos.

### Gene Set Enrichment Analysis (GSEA)

A diferencia del ORA que necesita un corte binario (significativo/no significativo), GSEA usa **toda la lista rankeada** de 3,255 proteínas ordenadas por logFC. Detecta si los genes de un pathway tienden a acumularse en los extremos superiores (activados) o inferiores (reprimidos) de la lista.

**Resultados GSEA Hallmarks MSigDB (colección curada de 50 procesos biológicos clave):**
- 15 gene sets significativos
- Algunos activados en tumor (ej: proliferación, metabolismo), otros reprimidos (ej: diferenciación)

**Resultados GSEA GO BP:** 455 términos significativos — la lista rankeada completa da mucho más poder estadístico que el ORA.

*¿Por qué usar ambos enfoques?* ORA y GSEA son complementarios. ORA es más intuitivo y conservador; GSEA aprovecha toda la información cuantitativa y detecta señales más sutiles.

---

## Fase 3 — Identificar candidatos farmacológicos (Scripts 04-07)

Esta es la fase central del reposicionamiento. La estrategia es consultar múltiples bases de datos independientes y luego integrar la evidencia. Usar múltiples fuentes reduce el riesgo de falsos positivos.

### Script 04: DGIdb — Interacciones gen-fármaco conocidas

**¿Qué es DGIdb?** Drug-Gene Interaction Database. Consolida información de ~30 fuentes distintas (ClinicalTrials.gov, TTD, CIViC, etc.) sobre qué fármacos interactúan con qué genes.

**¿Cómo consultamos?** Via GraphQL API v5 (la API REST v2 fue deprecada). Enviamos los 520 símbolos de genes y recuperamos todas las interacciones conocidas.

**Resultados:**
- **226 / 520 genes** (43.5%) tienen al menos una interacción farmacológica conocida
- **2,252 fármacos únicos** identificados
- **2,846 pares gen-fármaco** únicos

*Hallazgo destacado:* Los inhibidores del **proteosoma** (Bortezomib con 12 genes, Carfilzomib y Ixazomib con 10 cada uno) son los que más genes DE comparten. Esto es biológicamente coherente: el proteosoma (complejo de degradación de proteínas) está frecuentemente sobre-expresado en tumores, y su inhibición activa vías apoptóticas. Bortezomib está aprobado para mieloma múltiple.

### Script 05: ChEMBL — Fármacos en ensayos clínicos fase ≥ 3

**¿Qué es ChEMBL?** Base de datos del EBI (European Bioinformatics Institute) con >2 millones de compuestos y datos curados de mecanismo de acción. Es más conservadora que DGIdb porque solo incluye interacciones con evidencia mecanística fuerte.

**Flujo:**
1. UniProt ID → target ChEMBL (via API REST)
2. Target ChEMBL → mecanismos (drug-target pairs curados)
3. Mecanismos → detalles de molécula (max_phase, nombre, tipo)
4. Filtro: solo fármacos en fase clínica ≥ 3

**Resultados:**
- **309 / 520 proteínas** tienen un target en ChEMBL (Homo sapiens)
- **90 pares gen-fármaco** en fase ≥ 3 (tras deduplicación)
- **55 fármacos aprobados** (fase 4) + **34 en Fase III**
- Solo **21 genes** tienen fármacos fase ≥ 3 — estos son los "druggable targets" de alta confianza

*Fármacos destacados:*
- **Cetuximab / Erlotinib** → EGFR (ya aprobado para HNSCC, sirve de control positivo)
- **Lapatinib** → ERBB2/HER2 (aprobado para mama, reposicionamiento potencial)
- **Valproic acid** → HDAC (aprobado para epilepsia, interés en epigenética tumoral)
- **Decitabine** → DNMT1 (aprobado para leucemia, inhibidor de metilación de DNA)

*Nota técnica:* La librería `chembl_webresource_client` de Python se colgaba por un bug de caché interno sin timeout. El script fue reescrito usando `requests` directamente contra la REST API con timeout de 20 segundos, solucionando el problema.

### Script 06: Open Targets — Evidencia HNSCC + reposicionamiento novel

**¿Qué es Open Targets?** Plataforma del EBI/Sanger que integra múltiples tipos de evidencia para asociaciones gen-enfermedad: genética (GWAS, variantes somáticas), expresión diferencial, vías metabólicas, literatura, y ensayos clínicos.

**Estrategia de consulta:**
1. Paginar todos los ~8,715 targets asociados a HNSCC (EFO_0000181) con su "association score" (0-1)
2. Intersectar con nuestros 520 genes DE → 354 coinciden (68%!)
3. Para cada gen DE con Ensembl ID, consultar sus `knownDrugs` (todos los fármacos conocidos para ese gen, con indicación y fase)

**Resultados:**
- **354 / 520 genes** tienen evidencia directa en Open Targets para HNSCC
- **84 genes** tienen fármacos conocidos en la plataforma
- **67 fármacos aprobados** (para cualquier indicación)
- Solo **1 fármaco** tiene indicación directa para HNSCC en Open Targets
- **66 candidatos de reposicionamiento novel**: fármacos aprobados para otra indicación que target genes DE en HNSCC

*Hallazgo más llamativo:* **Metformin** target **16 genes DE** simultáneamente. La metformina (antidiabético, inhibidor de complejo I mitocondrial/AMPK) tiene evidencia epidemiológica de efecto protector en varios cánceres. Que 16 de nuestras proteínas DE sean targets conocidos de metformina sugiere una perturbación del metabolismo oxidativo en el HNSCC.

### Script 07: L2S2 — Compuestos que revierten la firma tumoral

**¿Qué es L2S2?** LINCS L1000 Signature Search (Evangelista et al., NAR 2025). Plataforma que indexa 1.67 millones de firmas transcriptómicas de 248 líneas celulares tratadas con 33,621 compuestos + 7,508 KO CRISPR. Accesible vía GraphQL API pública sin descarga local. Reemplazó a CMap2 (signatureSearch) que solo cubría 3,478 compuestos en 4 líneas celulares.

**Método:**

- **Query A:** top 150 proteínas UP en tumor → busca firmas DOWN de fármacos (Fisher's exact test, OR + pvalue)
- **Query B:** top 150 proteínas DOWN en tumor → busca firmas UP de fármacos
- **Filtros:** `filterFda=TRUE` (solo drugs FDA-aprobados), `pvalue < 0.001`
- **Score:** `reversal_score = -log10(best_p) × mean_log2(OR) × log2(n_sigs+1)`, normalizado a [-1, 0]
- **Reversión consistente:** drug tiene hits en AMBAS queries (UP→DOWN y DOWN→UP)

**Resultados:**

- **1,044 drugs FDA-aprobados** evaluados (6× más que CMap2)
- **248 líneas celulares** (62× más que CMap2)
- **931/1,044 (89%)** con reversión bidireccional consistente

*Top 5 reversores:*

| Compuesto | Score | Relevancia |
|-----------|-------|------------|
| Bortezomib | -1.000 | Inhibidor proteasoma; aprobado mieloma/linfoma; actividad en HNSCC publicada |
| Calcitriol | -0.858 | Vitamina D activa; receptor VDR presente en células HNSCC |
| Dasatinib | -0.855 | Inhibidor BCR-ABL/Src; Src kinase relevante en invasión HNSCC |
| Olaparib | -0.822 | Inhibidor PARP; señal de reparación de DNA en tumor |
| Vorinostat | -0.766 | Inhibidor HDAC; efecto epigenético en células HNSCC documentado |

*Inhibidores EGFR en top 20:* gefitinib (#7), erlotinib (#14), afatinib (#17) — validación de consistencia con análisis previos.

*Nota técnica:* L2S2 usa firmas transcriptómicas; nosotros usamos proteómica. La correlación transcriptómica-proteómica es imperfecta (~0.4-0.6 en general), por lo que los reversores deben considerarse como evidencia de soporte funcional, no absoluta. La cobertura de 1,044 drugs y 248 líneas celulares de L2S2 vs 3,478/4 de CMap2 aumenta sustancialmente la potencia estadística.

---

## Dónde estamos: resumen del pipeline

```
[DATOS]
  results_proteomica.tsv (3,352 proteínas, n=20 muestras)
       |
       v
[FASE 1: PREPARAR]
  01_parse_results_qc.R    ✅  → 520 proteínas DE (248 up, 272 down)
  02_id_mapping.R          ✅  → 518/520 con Entrez ID
       |
       v
[FASE 2: BIOLOGÍA]
  03_pathway_enrichment.R  ✅  → GO/KEGG/Reactome ORA + Hallmarks GSEA
       |
       v
[FASE 3: FÁRMACOS — 4 fuentes independientes]
  04_query_dgidb.py        ✅  → 2,252 fármacos, 226 genes con interacciones
  05_query_chembl.py       ✅  → 89 fármacos fase≥3, 21 genes druggable
  06_query_opentargets.py  ✅  → 66 candidatos reposicionamiento novel
  07_l2s2_connectivity.py  ✅  → 1,044 drugs FDA-aprobados, top: bortezomib (-1.000)
       |
       v
[FASE 3: INTEGRACIÓN]
  08_integrate_drug_targets.R  ✅  → 2,421 fármacos unificados; 187 candidatos multi-fuente
       |
       v
[FASE 4: PRIORIZACIÓN]
  09_string_network.R          ✅  → 403 nodos, 2,001 aristas; 41 hubs = Complejo I mitocondrial
  10_prioritization_scoring.R  ✅  → 177 candidatos (27 excluidos); Top 20 refinado
       |
       v
[FASE 5: VALIDACIÓN IN SILICO]
  11_clinicaltrials_pubmed.py  ✅  → 11/20 con trials HNSCC, 8 activos
  12_cosmic_overlap.py         ✅  → 2 driver genes (EGFR principal)
       |
       v
[FASE 6: EVIDENCIA FINAL + DOCUMENTACIÓN]
  13_evidence_summary.R        ✅  → Ranking combinado: Erlotinib/Cetuximab #1, Excel final
  14_methods_summary.R         ✅  → METHODS.md + OUTPUTS.md para manuscrito
```

---

## Fase 3 (integración) — Script 08: Tabla maestra de fármacos

### ¿Qué hicimos?

Unificamos los 4 listados independientes de fármacos en una sola tabla maestra normalizada. El reto principal fue que cada base de datos usa convenciones distintas de nombres de fármacos (ej: "Metformin" vs "Metformin HCl" vs "METFORMIN HYDROCHLORIDE"), por lo que fue necesario normalizar (minúsculas, eliminar sales, corregir caracteres especiales) antes de cruzar fuentes.

**Clasificación de reposicionamiento:**

| Clase | Criterio | N | Ejemplos |
| ----- | -------- | - | -------- |
| A | Aprobado con indicación directa HNSCC | 1 | Cedazuridine |
| B | Aprobado para otro cáncer | 51 | Bortezomib, Lapatinib, Decitabine |
| C | Aprobado para indicación no-oncológica | 80 | Metformin, Valproic acid, Deferoxamine |
| D | Experimental / investigación | 2,289 | Compuestos no aprobados por FDA |

**Hallazgo de integración:** 187 fármacos aparecen en ≥2 fuentes independientes. Esta convergencia multi-base-de-datos es el criterio más importante para priorizar, porque reduce drásticamente el riesgo de falsos positivos inherente en cualquier fuente individual.

*Nota técnica:* Los booleans escritos por Python ("True"/"False") no son reconocidos por `as.logical()` de R. Solución estándar en este pipeline: `tolower(as.character(col)) == "true"`.

---

## Fase 4 — Script 09: Red de interacciones proteína-proteína

### ¿Por qué una red PPI?

Los fármacos no actúan en proteínas aisladas. Las proteínas forman redes funcionales donde algunas son altamente conectadas (**hubs**). Un hub druggable es especialmente valioso porque su inhibición puede propagar el efecto a muchas vías simultáneamente.

### Cómo construimos la red

Usamos **STRING** (Search Tool for the Retrieval of Interacting Genes/Proteins) a través de su REST API v12. STRING integra interacciones físicas, coexpresión, fusiones génicas, coocurrencia filogenética, curación de literatura y datos experimentales, y asigna un "combined score" de 0 a 1000. Usamos un umbral de **700 (alta confianza)**.

*Nota:* El paquete R `STRINGdb` no estaba disponible (dependencia `chron` falla en compilación en este ambiente). Se implementó la consulta directamente via `httr2` + `jsonlite`. Un detalle importante: STRING devuelve respuestas con `Content-Type: text/json` (no `application/json`), por lo que se debe usar `resp_body_string()` + `fromJSON()` en lugar de `resp_body_json()`.

### Resultados

| Métrica de red | Valor |
| -------------- | ----- |
| Genes consultados | 518 |
| Mapeados en STRING | 513 (99%) |
| Nodos en componente gigante | 403 |
| Aristas (score ≥ 700) | 2,001 |
| Hubs (top 10% por grado) | 41 |
| Hubs druggables | 16 |

### El hallazgo biológico más importante del proyecto

**Los 41 hubs son exclusivamente subunidades de la cadena respiratoria mitocondrial:**

- NDUFA*, NDUFB*, NDUFS*, NDUFV* → **Complejo I** (NADH:ubiquinone oxidoreductase)
- ATP5*, ATP6*, ATP8* → **Complejo V** (ATP synthase)
- COX4, COX5, COX6, UQCR* → **Complejos III y IV**

Todos están **downregulados en tumor** (logFC negativo), coherente con la represión metabólica oxidativa conocida en HNSCC.

Los 16 hubs druggables son targets conocidos de **Metformin** (biguanida, inhibidor directo de Complejo I). Este hallazgo emerge de tres análisis independientes:

1. Open Targets (Script 06): Metformin target de 16 genes DE
2. Red PPI (Script 09): todos los hubs son subunidades de Complejo I
3. Scoring (Script 10): Metformin rank #5 con score 0.600

Esta convergencia es una señal muy sólida de perturbación de la fosforilación oxidativa en HNSCC.

---

## Fase 4 — Script 10: Scoring multi-criterio y Top 20 candidatos

### ¿Cómo priorizamos?

Con 187 candidatos multi-fuente, necesitamos una forma objetiva de ordenarlos. El scoring compuesto considera 6 dimensiones:

1. **|log2FC| del target** (peso 20%): cuánto cambia la proteína en tumor. Mayor cambio = mayor relevancia funcional.
2. **Significancia estadística** (peso 15%): -log10(adj.P.Val). Protege contra ruido estadístico.
3. **Fase clínica** (peso 20%): aprobado=1.0, fase3=0.75, fase2=0.5, fase1=0.25, experimental=0. Mayor peso porque determina la viabilidad translacional inmediata.
4. **L2S2 reversal score** (peso 15%): qué tan bien revierte la firma tumoral. Basado en evidencia funcional (enriquecimiento bidireccional Fisher), no solo asociativa.
5. **Pathway relevance** (peso 15%): si el target está en alguno de los pathways enriquecidos (ORA/GSEA significativos). Conexión con la biología alterada.
6. **Centralidad en red** (peso 15%): betweenness normalizado. Hubs son candidatos de mayor impacto sistémico.

### Top 10 candidatos (refinados tras exclusión)

Se excluyeron 27 fármacos no-antitumorales: gabapentinoides (CACNA2D1), anticoagulantes (F2), anti-amiloides Alzheimer (APP), y antídotos (ADH1B). Estos aparecían en trials HNSCC para manejo sintomático, no como terapia antitumoral.

| Rank | Fármaco | Score | L2S2 score | Clase | Indicación actual |
| ---- | ------- | ----- | ---------- | ----- | ----------------- |
| 1 | **Gefitinib** | 0.686 | -0.74 | B | NSCLC (EGFR mut) |
| 2 | **Metformin** | 0.685 | -0.45 | C | Diabetes tipo 2 |
| 3 | Vandetanib | 0.661 | -0.57 | B | Cáncer de tiroides |
| 4 | Mavacamten | 0.635 | — | C | Cardiomiopatía obstructiva |
| 5 | Tranylcypromine | 0.604 | -0.60 | C | Depresión (MAOI) |
| 6 | Metformin HCl | 0.598 | — | C | Diabetes tipo 2 |
| 7 | Azacitidine | 0.536 | -0.51 | B | MDS/AML |
| 8 | **Cedazuridine** | 0.541 | — | A | HNSCC (único clase A) |
| 9 | Doxycycline | 0.552 | -0.22 | C | Antibiótico |
| 10 | Pargyline | 0.563 | -0.46 | C | Antihipertensivo (MAOI) |

**Validación interna:** Gefitinib (anti-EGFR) en #1 con L2S2 score -0.74. Cedazuridine (clase A, HNSCC directo) presente. Metformin (Complejo I mitocondrial) sigue en top 3.

**Candidatos de reposicionamiento más interesantes:**

- **Mavacamten (#1)**: aprobado para cardiomiopatía (inhibe miosina cardíaca). MYH7 está significativamente downregulado en tumor. Mecanismo en HNSCC no explorado.
- **Metformin (#5)**: antidiabético con múltiple evidencia epidemiológica de efecto protector en cáncer. Target directo del Complejo I mitocondrial (nuestros hubs principales). 29 trials HNSCC activos, 214 papers PubMed.
- **Doxycycline (#16-18)**: antibiótico tetraciclina con actividad anti-MMP documentada. Los targets COL3A1/COL4A1/MMP8 son relevantes para la remodelación de ECM en HNSCC. 10 trials HNSCC, 41 papers.

### Selección de diversidad biológica

Para evitar que el Top 20 sea dominado por variantes del mismo compuesto, se aplicó un filtro de diversidad: máximo 3 fármacos por target primario.

---

## Fase 5 — Validación in silico (Scripts 11-12)

### Script 11: Evidencia clínica y bibliográfica

Para cada candidato del Top 20, se consultó:

- **ClinicalTrials.gov API v2**: ensayos registrados que mencionan el fármaco + HNSCC
- **PubMed E-utilities**: papers publicados que mencionan el fármaco + HNSCC

El evidence_score combina: proporción de trials HNSCC (50%), papers PubMed (30%), y trials activos (20%).

**Resultados:** 11/20 candidatos tienen ensayos clínicos en HNSCC. Los 3 con mayor evidencia (score=1.0): Cetuximab (369 trials, 2,229 papers), Erlotinib (56 trials, 417 papers), Metformin (29 trials, 214 papers).

*Nota técnica:* Los nombres de fármacos con sales ("ERLOTINIB HYDROCHLORIDE") daban ~0 resultados en PubMed. Se implementó `simplify_drug_name()` para eliminar sufijos de sales antes de la consulta.

### Script 12: Cruce con genes driver de cáncer

Se cruzaron las 520 proteínas DE con bases de datos de genes driver: NCG7 (~126 oncogenes embebidos) y 22 drivers HNSCC de literatura.

**Resultado:** Solo 2 genes DE solapan con drivers (principalmente EGFR, logFC=4.33). Esto es esperado: la proteómica DE captura consecuencias del tumor (metabolismo, ECM, inflamación), no las mutaciones conductoras.

---

## Fase 6 — Evidencia final y documentación (Scripts 13-14)

### Script 13: Ranking combinado y Excel final

Se integró el composite_score (script 10) con el evidence_score (script 11) en un ranking combinado:

```text
combined_rank_score = 0.60 x composite_score + 0.40 x evidence_score
```

Esto reordena los candidatos: Mavacamten (#1 por scoring) baja porque no tiene trials HNSCC, mientras que Erlotinib/Cetuximab/Metformin suben por su abundante evidencia clínica.

**Top 5 ranking final:**

| Rank | Fármaco | Combined score |
| ---- | ------- | -------------- |
| 1 | Gefitinib | 0.731 |
| 2 | Metformin | 0.811 (final: 0.685) |
| 3 | Vandetanib | 0.565 |
| 4 | Doxycycline | 0.470 |
| 5 | Metformin HCl | 0.759 |

El output principal es `13_FINAL_drug_candidates.xlsx` con 5 hojas: Top20_Final, Todos_Candidatos, Evidencia_Matriz, Ensayos_Clinicos, Metodologia.

### Script 14: Documentación automatizada

Genera `docs/METHODS.md` (sección de métodos lista para manuscrito, con parámetros exactos y versiones de paquetes) y `docs/OUTPUTS.md` (catálogo de todos los archivos generados).

---

## Reflexión biológica de los hallazgos

El pipeline completo converge en una historia biológica coherente:

1. **Metabolismo mitocondrial alterado (señal más fuerte)**: Metformin target 16 genes DE, todos los 41 hubs de la red PPI son subunidades del Complejo I mitocondrial, y Metformin tiene 29 trials activos en HNSCC con 214 papers. La represión de la fosforilación oxidativa es una de las perturbaciones centrales del HNSCC identificadas en este análisis.

2. **EGFR como control positivo**: Erlotinib, Cetuximab y Lapatinib (todos anti-EGFR/HER2) aparecen consistentemente en las primeras posiciones. EGFR es el único gen driver de cáncer significativamente DE (logFC=4.33). Cetuximab está aprobado para HNSCC, validando el pipeline.

3. **ECM y remodelación tisular**: Collagenase, Ocriplasmin, Doxycycline (anti-MMP) y Marimastat target colágenos y metaloproteasas sobre-expresados en tumor. La remodelación de ECM facilita invasión y metástasis en HNSCC.

4. **Doxycycline como candidato emergente**: Antibiótico tetraciclina con actividad anti-MMP bien documentada. 10 trials HNSCC, 41 papers. Bajo costo, perfil de seguridad conocido, y mecanismo biológicamente plausible.

5. **Señal proteosómica y epigenética (L2S2)**: Bortezomib (inhibidor proteasoma, #1 L2S2), vorinostat (inhibidor HDAC, #5 L2S2) y calcitriol (vitamina D, #2 L2S2) como top reversores. El proteasoma es una diana emergente en HNSCC con evidencia experimental publicada.

6. **Inhibidores EGFR validados por L2S2**: Gefitinib (#7), erlotinib (#14) y afatinib (#17) entre los top 20 reversores L2S2, confirmando la relevancia de EGFR desde tres análisis independientes (expresión, ClinicalTrials, y firma transcriptómica).

**Conclusión**: El pipeline identifica consistentemente 3 ejes terapéuticos principales en HNSCC: señalización EGFR (gefitinib, vandetanib, erlotinib — ya explotado clínicamente), metabolismo mitocondrial (Metformin, oportunidad de reposicionamiento con 29 trials activos), y remodelación de ECM (Doxycycline, oportunidad de reposicionamiento). La integración de L2S2 añade evidencia funcional de reversión transcriptómica que refuerza los candidatos identificados por expresión proteómica.

---

## Fase 7 — Análisis de sensibilidad (Script 15)

### ¿Por qué es necesario?

El scoring compuesto asigna pesos subjetivos a 6 dimensiones. Una pregunta legítima es: *¿cambiaría el ranking si diéramos más peso a la evidencia clínica vs. la señal proteómica?* El análisis de sensibilidad responde esto sistemáticamente.

### Resultados

Se evaluaron 6 configuraciones extremas de pesos:

| Configuración | Lógica |
| ------------- | ------ |
| baseline | Pesos actuales (usado para el top 20) |
| clinical_heavy | Máximo peso a fase clínica (45%) |
| molecular_heavy | Máximo peso a logFC + significancia (60%) |
| network_heavy | Máximo peso a centralidad de red (35%) |
| pathway_heavy | Máximo peso a relevancia de rutas (35%) |
| equal_weights | Todos los pesos iguales (sin sesgo) |

**Resultado clave: El ranking es altamente robusto.**

- **9 candidatos aparecen en el top 20 de las 6/6 configuraciones**: Mavacamten, Erlotinib, Lapatinib, Metformin, Omecamtiv Mecarbil, Collagenase C.h., Cedazuridine, Ocriplasmin, Forodesine.
- **Cetuximab**: 5/6 configs (sale del top 20 solo en `network_heavy` porque EGFR tiene centralidad media vs. los hubs mitocondriales)
- **Top 5 final (Erlotinib, Cetuximab, Metformin, Lapatinib, Doxycycline)**: ALTAMENTE ROBUSTO

Este resultado es importante para la discusión de la tesis: la conclusión no depende de la elección arbitraria de pesos.

---

## Próximos pasos

Ver `docs/FUTURE_WORK.md` para el plan detallado de:

- **Nivel 2**: Validación en TCGA-HNSC (survival analysis) + figuras de publicación
- **Nivel 3**: Sección de limitaciones + comparación con literatura para manuscrito

---

Generado automáticamente — 2026-03-04

---

## Revisión crítica Fase 4–5 — Scripts 10, 11 y 12 (2026-03-07)

Revisión minuciosa de los scripts de priorización y validación in silico. Se documentan todos los bugs encontrados y corregidos.

### Script 10 — `10_prioritization_scoring.R`

| # | Tipo | Descripción | Impacto |
|---|------|-------------|---------|
| 1 | **Bug** | `toupper()` aplicado a `excluded_drugs` pero no a `drug_name_norm` → exclusiones silenciosamente fallidas (la columna está en minúsculas). Corregido: `toupper(drug_name_norm) %in% excluded_drugs`. | Los 18 fármacos no-antitumorales (Alzheimer, anticoagulantes, gabapentinoides) no se excluían del Top 20. |
| 2 | **Bug** | Falta validación de suma de pesos. Si `config/analysis_params.yaml` se modifica y los pesos no suman 1.0, el score compuesto queda mal escalado sin advertencia. Corregido: `stop()` si `abs(sum - 1.0) > 0.001`. | Protege contra errores silenciosos en análisis de sensibilidad. |
| 3 | **Bug** | `betweenness_norm` (output de script 09, ya normalizado 0-1) se dividía nuevamente por `max_betw` → doble normalización distorsiona el componente de red. Corregido: usar directamente `nm$betweenness_norm`. | El componente `score_network` (peso 0.15) estaba comprimido incorrectamente. |
| 4 | **Bug** | Header del log listaba `"radar"` como output, pero el script genera `"scatter"`. Corregido en lista de outputs. | Cosmético/trazabilidad. |
| 5 | **Calidad** | Deduplicación ChEMBL con `group_by("")` sobre strings vacíos colapsaba todos los fármacos sin ChEMBL en uno solo. Refactorizado: filtrar antes de `group_by`, luego `bind_rows` con no-ChEMBL. | Podía eliminar candidatos legítimos sin ChEMBL ID. |
| 6 | **Calidad** | `library(ggrepel)` faltaba en el bloque inicial (se usaba en el scatter plot). Agregado. | Error en ejecución si ggrepel no está cargado. |
| 7 | **Calidad** | Columna `primary_target` no incluida en `out_cols` → ausente en `10_top20_candidates.tsv`. Agregada. | Script 11 usaba `drug_class` para colorear figuras; script 12 usaba `de_genes` con separador `\|`. |
| 8 | **Calidad** | Carga muerta de `drug_summary` (archivo no existente en pipeline actual). Eliminada. | Error si se ejecuta sin ese archivo. |

### Script 11 — `11_clinicaltrials_pubmed.py`

| # | Tipo | Descripción | Impacto |
|---|------|-------------|---------|
| 1 | **Bug** | `class_col` buscaba `["repurposing_class", "class", "clase"]` pero script 10 exporta `drug_class`. Corregido: agregar `"drug_class"` como primera opción. | Bubble plot y barplots sin diferenciación de color por clase (A/B/C/D); todos los puntos tomaban color por defecto. |
| 2 | **Bug** | `simplify_drug_name`: loop `for` sobre sufijos eliminaba solo el primero encontrado; nombres compuestos (ej. `"drug mesylate monohydrate"`) quedaban con sufijos residuales. Corregido con `while changed` + `break`. | Búsquedas en APIs con nombre incorrecto → menos ensayos encontrados para algunos candidatos. |
| 3 | **Calidad** | `HNSCC_CT_TERMS` definida pero no usada; keywords hardcodeadas en `query_clinicaltrials()`. Corregido: usar la constante como fuente única. | Dead code; inconsistencia entre constante y uso real. |
| 4 | **Calidad** | `drug_confirmed` calculado pero no usado en el filtro de `hnscc_trials`. Corregido: `[t for t in all_trials if t["is_hnscc"] and t["drug_confirmed"]]`. | `ct_hnscc_trials` sobreestimado para fármacos con nombres genéricos. |
| 5 | **Calidad** | `"rettype": "count"` inválido con `format="json"` en PubMed E-utilities (silenciosamente ignorado). Eliminado de ambos dicts de parámetros. | Claridad del código; previene advertencias en versiones futuras de la API. |
| 6 | **Calidad** | `NCBI_EMAIL` hardcodeado como `researcher@example.com`. Cambiado a variable de entorno con fallback real: `os.environ.get("NCBI_EMAIL", "jcarvajal@fucsalud.edu.co")`. Log de inicio reporta el email utilizado. | Buenas prácticas NCBI E-utilities; identificación correcta ante el servidor. |

### Script 12 — `12_cosmic_overlap.py`

| # | Tipo | Descripción | Impacto |
|---|------|-------------|---------|
| 7 | **Bug crítico** | `annotate_row`: `de_genes_str.split("/")` debería ser `split("\|")`. Script 10 almacena genes separados por `\|` (e.g., `"EGFR\|ERBB2"`). Con `/`, cada string se trata como gen único y nunca coincide con `any_driver`. | `n_target_drivers` siempre 0 → `has_driver_target` siempre `False` en `12_top20_with_drivers.tsv`. Erlotinib/Cetuximab (EGFR) no quedaban marcados como driver-targets. |
| 8 | **Calidad** | `INPUT_DE` sin manejo de `FileNotFoundError`. Corregido: verificar existencia con log WARNING y `df_all = None` si no existe; proteger uso posterior con `if df_all is not None`. | Crash no informativo si script 01 no se ha ejecutado o el nombre cambió por timestamp. |
| 9 | **Calidad** | `sym_col` inferido de `df_sig` pero usado directamente en `df_all` sin verificar que la columna existe. Corregido: `if sym_col in df_all.columns` antes de construir volcano. | `KeyError` si los dos archivos DE tienen esquemas diferentes. |
| 10 | **Calidad** | Fallback IntOGen pan-cancer loggeado como `INFO`. Degradado a `WARNING` con mensaje explícito sobre riesgo de inflar falsos positivos. | El fallback usa miles de genes → solapamiento `intogen_de` artificialmente grande sin advertencia visible. |

### Archivos modificados en esta sesión (2026-03-07 — Scripts 10-12)

| Archivo | Cambios |
|---------|---------|
| `scripts/10_prioritization_scoring.R` | Bugs 1-3 + calidad 4-8 (ver tabla arriba) |
| `scripts/11_clinicaltrials_pubmed.py` | Bugs 1-2 + calidad 3-6 (ver tabla arriba) |
| `scripts/12_cosmic_overlap.py` | Bug crítico 7 + calidad 8-10 (ver tabla arriba) |
| `config/analysis_params.yaml` | Corrección de 4 comentarios con número de script incorrecto (09, 10 en lugar de 10, 11) |
| `docs/OUTPUTS.md` | Agregados: `10_scoring_heatmap.pdf`, `10_top20_barplot.pdf`, `11_trials_phase_bar.pdf`; columna `primary_target` en descripción de `10_all_candidates_scored.tsv` |
| `docs/METHODS.md` | Header actualizado con scripts 10, 11, 12 |
| `docs/RUNBOOK.md` | Nota de troubleshooting para variable `NCBI_EMAIL` |

---

## Revisión crítica — Scripts 14 y 15 (2026-03-07)

### Script 14 — `14_methods_summary.R`

| # | Tipo | Descripción | Impacto |
|---|------|-------------|---------|
| 1 | **Bug crítico** | `%||%` usada desde línea 73 pero definida en línea 236 — crash garantizado. Primera asignación de `methods_text` (líneas 58-233) era código muerto. Corregido: definir `%||%` antes del primer uso y eliminar el bloque muerto. | El script no podía ejecutarse; solo el bloque `local({...})` producía output pero perdía secciones del methods_text largo. |
| 2 | **Bug** | Typo `TVsS` debería ser `TVsN` (Tumor vs. Normal). | Habría entrado en `METHODS.md` del manuscrito. |
| 3 | **Calidad** | Pipe nativo con placeholder `_` (`\|> grep(..., x = _, ...)`) requiere R ≥ 4.2. Frágil y dependiente de formato de `cat_files()`. Reemplazado por variable intermedia con `grep()` tradicional. | Fallo en R 4.1; potencial para OUTPUTS.md vacío o con entradas duplicadas. |
| 4 | **Calidad** | Sección "Prioritization & Evidence" en OUTPUTS.md con lógica confusa (regex frágil + posible duplicación). Refactorizado: usar `grep()` sobre `cat_files("results/tables")` con patrón `/(10_\|11_\|12_\|13_)/`. | Sección podía quedar vacía o con entradas duplicadas. |

### Script 15 — `15_sensitivity_analysis.R`

| # | Tipo | Descripción | Impacto |
|---|------|-------------|---------|
| 1 | **Bug crítico** | `filter(!drug_name_norm %in% excluded_drugs)` — `excluded_drugs` está en mayúsculas (`toupper`), `drug_name_norm` en minúsculas → ningún fármaco excluido. Mismo bug que script 10 bug #1. Corregido: `toupper(drug_name_norm)`. | Los 18 fármacos no-antitumorales seguían en el pool de candidatos del análisis de sensibilidad, invalidando la comparación con el baseline de script 10. |
| 2 | **Bug crítico** | `filter(!de_genes %in% excl_targets)` — `de_genes` contiene strings multi-gen separados por `\|` (ej: "EGFR\|ERBB2"); nunca coincide con un gen individual. Corregido: `sapply()` con `str_split(..., "\\|")` + `any(...%in% excl_targets)`. | Targets excluidos (CACNA2D1, F2, APP, etc.) no se filtraban. |
| 3 | **Calidad** | Variable `bonus` usada en `mutate()` (líneas 134 y 339) sin verificar que existe como columna en el TSV de entrada. Añadido: check explícito con fallback `bonus <- 0` y log `WARN`. | Si script 10 no exportó la columna `bonus`, el script crasheaba con error de variable no encontrada. |

### Archivos modificados — Scripts 14-15 (2026-03-07)

| Archivo | Cambios |
|---------|---------|
| `scripts/14_methods_summary.R` | Reescritura estructural: eliminar 175 líneas de código muerto, mover `%||%`, fix descripción TVsS = "Tumor vs. Adjacent Normal", fix OUTPUTS pipe |
| `scripts/15_sensitivity_analysis.R` | Bug crítico 1-2 (exclusiones rotas) + calidad 3 (check columna bonus) |

---

## Revisión crítica — Script 17 (2026-03-07)

### Script 17 — `17_pub_figures.R`

| # | Tipo | Descripción | Impacto |
|---|------|-------------|---------|
| 1 | **Bug crítico** | `drug_class_label` no existe en `10_top20_candidates.tsv` (script 10 exporta `drug_class` = A/B/C/D). `str_remove(drug_class_label, "^[A-Z]: ")` crasheaba. Corregido: lookup `DRUG_CLASS_LABELS["A"] = "HNSCC-approved"`, etc. | Crash en D3 lollipop. |
| 2 | **Bug crítico** | `col_cats[colnames(ev_mat)]` usaba strings exactos con acentos ("Proteómica DE", "Fase clínica ≥ 3") para mapear categorías. Si los nombres de columna en `13_evidence_matrix.tsv` difieren, produce NAs y `HeatmapAnnotation` falla. Corregido: función `assign_ev_category()` con `grepl()` insensible a acentos/idioma. | Crash o anotación vacía en E1 evidence heatmap. |
| 3 | **Bug crítico** | `avg_intensity_TVsS` no es columna estándar de limma (`AveExpr`). MA plot crasheaba si script 01 no renombró esa columna. Corregido: detección automática con fallback a `AveExpr` y log informativo. | Crash en A2 MA plot. |
| 4 | **Bug crítico** | `logFC_TVsS` y `adj.P.Val_TVsS` ausentes en `09_network_node_metrics.tsv` (script 09 exporta métricas de red, no valores DE). Corregido: `left_join` desde tabla DE antes de construir `net_df`. | Crash en D1 scatter. |
| 5 | **Bug** | Labels volcano y MA plot: `slice_head(n=20)` sobre todo el rango ordena por `abs(logFC)` y en datos reales todos los top 20 pueden ser de la misma dirección → una dirección sin etiquetas. Corregido: `bind_rows(top N up, top N down)` separados por dirección. | Volcano: solo genes down etiquetados. MA: solo genes down etiquetados. |
| 6 | **Calidad** | `scale_y_reverse(limits = c(20, 1))` — `limits` con min > max genera warning/comportamiento indefinido. Corregido: `limits = c(1, 20)` (orden natural antes de revertir). | Warning en F1 bump chart; eje posiblemente sin límite inferior. |
| 7 | **Calidad** | `FIG1_QC_DE` y `A_multipanel_QC_DE` guardaban la misma figura dos veces. Eliminado el guardado duplicado en sección A; se conserva solo `FIG1_QC_DE`. | Output redundante (2 archivos idénticos). |
| 8 | **Calidad** | Join `meta_raw` para columna `vph` nunca usada en el plot PCA (solo `condition`). Eliminado. | Dead code; podía crashear si `metadata.csv` no tenía columnas esperadas. |
| 9 | **Calidad** | `ha_ev_row` definido como `NULL` cuando `evidence_level` no existe — `Heatmap(right_annotation = NULL)` es válido en ComplexHeatmap, pero mejor explícito. Refactorizado con manejo condicional. | Robusto ante ausencia de columna `evidence_level`. |
| 10 | **Nomenclatura** | Títulos usaban "Surrounding" inconsistente con terminología del proyecto. Corregido a "Adjacent Normal" (TVsS = Tumor vs. Sano, tejido adyacente del mismo paciente). | Consistencia en figuras y METHODS.md. |

### Archivos modificados — Script 17 (2026-03-07)

| Archivo | Cambios |
|---------|---------|
| `scripts/17_pub_figures.R` | Reescritura completa: bugs 1-5 críticos + calidad 6-10 + quitar barplot logFC de heatmap topDE + fix etiquetas volcano/MA |
| `scripts/14_methods_summary.R` | Fix nomenclatura: "TVsN" → "TVsS = Tumor vs. Adjacent Normal" |
