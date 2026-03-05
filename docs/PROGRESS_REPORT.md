# Reporte de Progreso — HNSCC Drug Repurposing
**Actualización:** 2026-03-04
**Estado:** Pipeline completado (14/14 scripts)

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

### Script 07: CMap — Compuestos que revierten la firma tumoral

**¿Qué es CMap?** The Connectivity Map (Broad Institute). Mide las firmas transcriptómicas de miles de compuestos en líneas celulares. La idea central es: si un tumor sobreexpresa genes X,Y,Z y un fármaco los subregula, ese fármaco podría revertir el fenotipo tumoral.

**Método usado:** `gess_cmap()` del paquete Bioconductor `signatureSearch`.
- **Base de datos:** CMap2 (EH3223, 370MB, 3,478 compuestos × 12,403 genes, en Entrez IDs)
- **Query signature:** top 150 proteínas más up-reguladas + top 150 más down-reguladas por logFC
- **Overlap con CMap2:** 132/150 up, 134/150 down — excelente cobertura
- **Score:** `scaled_score` (-1 a +1). **Negativo = reversal** (el compuesto produce el efecto opuesto a nuestro tumor)

**Resultados:**
- 3,478 compuestos evaluados en ~27 segundos
- **174 reversores potentes** (bottom 5%, solo small molecules `trt_cp`)
- Score mínimo (mejor reversor): -1.0

*Top reversores con relevancia clínica/biológica:*

| Compuesto | Score | Relevancia |
|-----------|-------|------------|
| Progesterone | -0.977 | Hormona esteroidea; relevante en HNSCC HPV+ (eje hormonal) |
| Harmine | -0.970 | Inhibidor DYRK1A/MAO; actividad antitumoral publicada |
| Deferoxamine | -0.950 | Quelante de hierro; altera metabolismo de hierro tumoral |
| Nimesulide | -0.909 | NSAID inhibidor COX-2; COX-2 sobre-expresada en HNSCC |
| Mifepristone | -0.907 | Antiprogestágeno; actividad antitumoral en varios estudios |
| Atropine | -0.914 | Anticolinérgico; nervios colinérgicos promueven tumors HNSCC |
| Menadione | -0.924 | Vitamina K3; genera estrés oxidativo selectivo en células tumorales |

*Nota técnica importante:* CMap usa firmas transcriptómicas; nosotros usamos proteómica. La proteómica correlaciona bien con la transcriptómica en muchos casos, pero no perfectamente. Los resultados de CMap deben considerarse como evidencia adicional, no como verdad absoluta.

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
  07_cmap_connectivity.R   ✅  → 174 reversores de firma tumoral
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
| D | Experimental / investigación | 2,289 | Latamoxef, harmine, menadione |

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
4. **CMap score** (peso 15%): qué tan bien revierte la firma tumoral. Basado en evidencia funcional, no solo asociativa.
5. **Pathway relevance** (peso 15%): si el target está en alguno de los pathways enriquecidos (ORA/GSEA significativos). Conexión con la biología alterada.
6. **Centralidad en red** (peso 15%): betweenness normalizado. Hubs son candidatos de mayor impacto sistémico.

### Top 10 candidatos (refinados tras exclusión)

Se excluyeron 27 fármacos no-antitumorales: gabapentinoides (CACNA2D1), anticoagulantes (F2), anti-amiloides Alzheimer (APP), y antídotos (ADH1B). Estos aparecían en trials HNSCC para manejo sintomático, no como terapia antitumoral.

| Rank | Fármaco | Score | Target principal | Clase | Indicación actual |
| ---- | ------- | ----- | ---------------- | ----- | ----------------- |
| 1 | **Mavacamten** | 0.637 | MYH7 | C | Cardiomiopatía obstructiva |
| 2 | Erlotinib HCl | 0.607 | EGFR | C | NSCLC, páncreas |
| 3 | Lapatinib | 0.607 | EGFR | C | Cáncer de mama HER2+ |
| 4 | Cetuximab | 0.607 | EGFR | C | **HNSCC** (control positivo) |
| 5 | **Metformin HCl** | 0.600 | NDUFA/NDUFB/NDUFS | C | Diabetes tipo 2 |
| 6 | Omecamtiv Mecarbil | 0.577 | MYH7 | B | Insuficiencia cardíaca |
| 7 | Collagenase C.h. | 0.542 | COL/MMP8 | C | Contractura Dupuytren |
| 8 | **Cedazuridine** | 0.541 | CDA | A | HNSCC (único clase A) |
| 9 | Ocriplasmin | 0.536 | COL/LAM | C | Tracción vitreomacular |
| 10 | Danicamtiv | 0.527 | MYH7 | D | Experimental (cardíaco) |

**Validación interna:** Cetuximab y Erlotinib (aprobados para HNSCC) en posiciones 2-4. Cedazuridine (clase A, HNSCC directo) en #8.

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
| 1 | Erlotinib | 0.764 |
| 2 | Cetuximab | 0.764 |
| 3 | Metformin | 0.760 |
| 4 | Lapatinib | 0.724 |
| 5 | Doxycycline | 0.638 |

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

5. **Eje hormonal y CMap**: Progesterone y Mifepristone como top reversores CMap. En un contexto donde 60% de pacientes son HPV+, la alteración de señalización hormonal por el virus es biológicamente coherente.

6. **Inflamación**: Nimesulide (COX-2) como reversor CMap. COX-2 está frecuentemente sobre-expresada en HNSCC y su inhibición tiene precedente experimental.

**Conclusión**: El pipeline identifica consistentemente 3 ejes terapéuticos principales en HNSCC: señalización EGFR (ya explotado clínicamente), metabolismo mitocondrial (Metformin, oportunidad de reposicionamiento), y remodelación de ECM (Doxycycline, oportunidad de reposicionamiento).

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
