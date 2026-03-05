# Metodología — HNSCC Drug Repurposing Pipeline

*Guía académica y didáctica del pipeline computacional*
*Última actualización: 2026-03-06*

---

## Introducción y Lógica del Proyecto

### ¿Por qué reposicionamiento de fármacos?

El desarrollo de un fármaco nuevo desde cero toma en promedio 12–17 años y costos superiores a los mil millones de dólares. El **reposicionamiento de fármacos** (*drug repurposing*) es una estrategia alternativa que busca nuevos usos terapéuticos para compuestos ya aprobados o en desarrollo clínico avanzado. Al partir de moléculas con perfiles de seguridad conocidos, se reduce significativamente el tiempo y el riesgo del desarrollo clínico.

### ¿Por qué HNSCC?

El **carcinoma de células escamosas de cabeza y cuello** (HNSCC) es el sexto cáncer más frecuente a nivel mundial. Pese a avances en inmunoterapia, la mayoría de los pacientes con enfermedad avanzada o recurrente tienen opciones limitadas y pronóstico desfavorable. Esta situación clínica justifica la búsqueda activa de candidatos farmacológicos alternativos mediante estrategias computacionales.

### Lógica general del pipeline

El pipeline parte de datos proteómicos de tumores HNSCC comparados con tejido normal adyacente. Mediante análisis diferencial se identifican las proteínas alteradas en el tumor. Estas proteínas son la "firma molecular" que guía todo lo demás: el enriquecimiento funcional responde *¿en qué procesos biológicos están involucradas?*, las consultas a bases de datos responden *¿hay fármacos que actúen sobre ellas?*, la red de interacciones responden *¿cuáles son las más centrales e influyentes?*, y finalmente el scoring integra toda esta información para priorizar candidatos con el mayor respaldo multimodal.

```text
Proteómica tumoral
       ↓
   Análisis DE (limma)
       ↓
  Enriquecimiento funcional ← ¿Qué procesos están alterados?
       ↓
  Consulta a BD de fármacos ← ¿Hay fármacos para esos genes?
       ↓
  Red PPI (STRING)          ← ¿Qué genes son más influyentes en la red?
       ↓
  Scoring multi-criterio    ← ¿Qué candidatos tienen mayor soporte?
       ↓
  Validación clínica        ← ¿Hay ensayos o publicaciones que los respalden?
       ↓
  Ranking final integrado
```

---

## Datos de Entrada

### Proteómica DIA (Script 01)

El estudio utiliza datos de **proteómica cuantitativa** adquiridos por el método **DIA** (*Data-Independent Acquisition*). A diferencia del método tradicional DDA (*Data-Dependent Acquisition*), en DIA todos los péptidos presentes en la muestra son fragmentados sistemáticamente, lo que produce cuantificaciones más reproducibles y con menor número de datos faltantes entre muestras. Esto es especialmente importante en estudios comparativos como este.

Los datos provienen de **20 muestras** de tejido: 10 pares tumor/normal adyacente de pacientes con HNSCC. El diseño pareado es metodológicamente sólido porque permite controlar la variabilidad inter-paciente: cada paciente actúa como su propio control, de modo que las diferencias observadas se atribuyen más directamente al estado tumoral que a características individuales del paciente.

El conjunto incluye pacientes **VPH-positivos** y **VPH-negativos**, lo que refleja la heterogeneidad biológica real del HNSCC. El análisis primario promedia sobre el estado viral (contraste Tumor vs. Normal), lo que busca identificar alteraciones proteicas compartidas independientemente del mecanismo de carcinogénesis.

Las intensidades de proteínas fueron preprocesadas con **MaxQuant** (software estándar para análisis de espectros MS) antes de ingresar al pipeline de R.

---

## Fase 1: Control de Calidad y Análisis Diferencial

### Script 01 — `01_parse_results_qc.R`

**¿Qué hace?** Carga los resultados de expresión diferencial pre-calculados, realiza QC visual de las intensidades y exporta tablas clasificadas por dirección de cambio.

**¿Por qué QC?** Antes de cualquier análisis funcional es indispensable verificar que los datos tienen la distribución esperada y que no hay muestras con anomalías técnicas. El boxplot de intensidades por muestra permite detectar muestras con distribuciones atípicas que podrían indicar fallas en la preparación de la muestra o en la adquisición.

El **PCA** (Análisis de Componentes Principales) proyecta las 20 muestras en un espacio reducido de varianza máxima. Si el análisis funciona correctamente, se espera que las muestras se agrupen por condición (tumor vs. normal), lo que validaría que la señal biológica dominante en los datos es la diferencia tumoral y no el efecto batch ni otras variables técnicas.

**¿Por qué limma/proteoDA?** El análisis diferencial se realizó con **limma** (*Linear Models for Microarray and RNA-seq Data*), adaptado para proteómica mediante el wrapper **proteoDA**. limma utiliza modelos lineales con estimación empírica de Bayes del error (*empirical Bayes shrinkage*), lo que estabiliza las estimaciones de varianza proteína a proteína especialmente en experimentos con pocos replicados (como los n=10 pares aquí).

**Criterios de significancia:** |log2FC| > 1 (al menos el doble de cambio) y FDR ajustado < 0.05 (corrección Benjamini-Hochberg). El umbral de |log2FC| > 1 busca un balance entre sensibilidad y relevancia biológica: cambios menores son estadísticamente detectables pero pueden carecer de impacto funcional. La corrección FDR controla la tasa de falsos positivos esperada en pruebas múltiples (una por cada proteína cuantificada).

---

## Fase 2: Mapeo de Identificadores

### Script 02 — `02_id_mapping.R`

**¿Qué hace?** Convierte los identificadores de proteínas de formato UniProt a Entrez ID y símbolo génico canónico (HGNC).

**¿Por qué es necesario?** Los datos proteómicos usan **UniProt IDs** (sistema de identificación de proteínas). Sin embargo, las herramientas de enriquecimiento funcional, las bases de datos de interacciones génicas y las redes PPI usan predominantemente **Entrez IDs** (sistema NCBI) o **símbolos génicos** (como *TP53*, *EGFR*). Sin este paso de conversión, no sería posible conectar los resultados proteómicos con los recursos de anotación funcional.

Se utiliza el paquete **`org.Hs.eg.db`** (base de datos de anotación del genoma humano en R/Bioconductor) a través de la función `bitr()` de clusterProfiler. Este paquete contiene los mapeos entre todos los sistemas de identificación principales para *Homo sapiens* y se actualiza con cada release de Bioconductor.

Un aspecto técnico importante: algunas proteínas UniProt corresponden a isoformas (identificadores con sufijo `-1`, `-2`, etc.). El script limpia estos sufijos antes del mapeo para evitar pérdidas innecesarias de anotación.

---

## Fase 3: Enriquecimiento Funcional

### Script 03 — `03_pathway_enrichment.R`

**¿Qué hace?** Realiza análisis de enriquecimiento funcional para interpretar qué procesos biológicos y rutas de señalización están sobre-representados entre las proteínas DE.

**¿Por qué enriquecimiento funcional?** Una lista de cientos de proteínas alteradas es difícil de interpretar gene por gene. El análisis de enriquecimiento permite elevar el nivel de interpretación de *genes individuales* a *conceptos biológicos* (procesos, funciones, rutas). Esto facilita entender la fisiopatología subyacente y también guía la identificación de dianas farmacológicas con relevancia en vías conocidas del cáncer.

#### ORA — Over-Representation Analysis

El **ORA** (análisis de sobre-representación) evalúa si la proporción de genes DE que pertenecen a un conjunto génico dado es mayor de lo esperado por azar, usando la prueba exacta de Fisher (o hipergeométrica). Para esto se necesita:

1. **Lista de genes de interés**: proteínas DE significativas
2. **Universo de fondo**: todas las proteínas cuantificadas (no solo las DE)

El universo de fondo es metodológicamente crítico: si se usara el genoma humano completo como fondo, se sobreestimaría el enriquecimiento de categorías de proteínas que tienden a ser detectables por proteómica (proteínas abundantes, de expresión alta). Al usar solo las proteínas cuantificadas como fondo, la comparación es justa.

Se aplica ORA a tres ontologías de Gene Ontology (GO):

- **BP** (*Biological Process*): procesos celulares como "respuesta al daño en el ADN" o "ciclo celular"
- **MF** (*Molecular Function*): actividades moleculares como "actividad quinasa"
- **CC** (*Cellular Component*): compartimentos celulares como "núcleo" o "membrana plasmática"

Y a dos bases de rutas:

- **KEGG**: rutas metabólicas y de señalización con representación visual de circuitos moleculares
- **Reactome**: base de datos curada manualmente con alta cobertura de procesos humanos

Los términos GO BP se simplifican semánticamente con `simplify()` (umbral 0.7) para eliminar términos redundantes: GO es una ontología jerárquica donde términos padre e hijo suelen aparecer enriquecidos simultáneamente. La simplificación retiene el término más específico o representativo.

#### GSEA — Gene Set Enrichment Analysis

A diferencia del ORA (que trabaja con un umbral binario de significancia), el **GSEA** usa el *ranking continuo* de todas las proteínas cuantificadas, ordenadas por su estadístico de cambio (logFC × -log10 p.adj). Esto tiene dos ventajas:

1. No depende de un umbral arbitrario de significancia
2. Detecta enriquecimientos en genes con cambios modestos pero coordinados en el mismo set génico

El GSEA calcula un **Normalized Enrichment Score (NES)** positivo si el set tiende a aparecer en el extremo superior del ranking (genes up-regulados) o negativo si tiende al extremo inferior (genes down-regulados).

Se aplica GSEA a:

- **MSigDB Hallmark gene sets**: 50 conjuntos curados que representan estados biológicos bien definidos en cáncer (e.g., *Hallmark_E2F_Targets*, *Hallmark_Hypoxia*). Su nivel de curación los hace ideales para interpretación biológica de alto nivel.
- **GO Biological Process**: para mayor granularidad mecanística

Parámetros clave: `minGSSize=10, maxGSSize=500` limitan el análisis a sets que tengan representación cuantificable (ni muy pequeños para ser ruidosos, ni tan grandes que sean inespecíficos). Se usan 1000 permutaciones para estimar la distribución nula.

---

## Fase 4: Consulta a Bases de Datos de Fármacos

La identificación de candidatos farmacológicos se apoya en **cuatro fuentes complementarias**. Ninguna base de datos cubre el espacio farmacológico completo; la integración multi-fuente maximiza la cobertura y permite identificar candidatos con soporte independiente en múltiples sistemas.

### Script 04 — `04_query_dgidb.py`

**DGIdb** (*Drug Gene Interaction Database*) agrega información de interacciones gen-fármaco de más de 30 fuentes primarias, incluyendo PharmGKB, TTD, ChEMBL, DrugBank, entre otras. Para cada gen DE, se consulta qué fármacos conocidos lo tienen como diana, junto con el tipo de interacción (inhibidor, activador, modulador) y la fuente de la evidencia.

**¿Por qué DGIdb?** Es el repositorio de interacciones gen-fármaco más comprehensivo y de acceso libre. Su cobertura es amplia pero la calidad de la evidencia es heterogénea (mezcla fuentes computacionales y experimentales), por lo que se complementa con fuentes de mayor curación clínica.

Se usa la **GraphQL API v5** para consultas estructuradas que retornan directamente los campos necesarios (gene, drug, interaction_score, interaction_types, sources), minimizando el procesamiento posterior.

### Script 05 — `05_query_chembl.py`

**ChEMBL** es la base de datos de moléculas bioactivas de EMBL-EBI, con más de 2.4 millones de compuestos. Se consulta específicamente por fármacos en **fase clínica >= 3** que tengan como diana proteínas DE del estudio.

**¿Por qué filtrar por fase >= 3?** La fase clínica es un proxy del nivel de evidencia de seguridad y actividad de un compuesto. Los fármacos en fase 3 o aprobados han superado pruebas extensas de toxicidad y eficacia inicial, lo que los hace candidatos más viables para reposicionamiento (menor riesgo clínico).

De ChEMBL se extrae además el **SMILES** del compuesto (representación estructural química), útil para análisis de similitud química posterior.

### Script 06 — `06_query_opentargets.py`

**Open Targets Platform** es una asociación entre EMBL-EBI, Sanger Institute, GSK y otras instituciones, diseñada específicamente para integrar evidencia de asociación gen-enfermedad y gen-fármaco. Proporciona:

1. **Scores de asociación gen-HNSCC**: cuantifican cuánta evidencia existe de que un gen está involucrado en HNSCC, integrando genómica, datos de expresión, ensayos clínicos y literatura.
2. **Interacciones gen-fármaco con indicación clínica**: incluye si el fármaco está aprobado para HNSCC específicamente (`is_hnscc_indication`).

**¿Por qué Open Targets?** Es la única fuente que conecta directamente la información de target con la indicación oncológica de interés (HNSCC). Un fármaco que ya está aprobado para HNSCC en Open Targets recibe la clasificación de mayor prioridad en este pipeline.

Se usa la **GraphQL API** de Open Targets (más eficiente que REST para este tipo de consultas relacionales).

### Script 07 — `07_cmap_connectivity.R`

**CMap** (*Connectivity Map*) es un proyecto del Broad Institute que mide los cambios transcriptómicos en líneas celulares tratadas con miles de compuestos. La idea central es que si la firma transcriptómica de un fármaco es **opuesta** a la firma del tumor, ese fármaco podría revertir el estado tumoral.

Se usa la función `gess_cmap()` del paquete **signatureSearch** con las **top 150 proteínas up-reguladas y 150 down-reguladas** (por |logFC|) como firma de consulta. El resultado es un **scaled_score**: valores muy negativos indican que el fármaco tiende a inducir la expresión opuesta a la firma tumoral (reversores potenciales).

**¿Por qué CMap?** A diferencia de las otras tres fuentes (que buscan fármacos que actúen directamente sobre los genes DE), CMap identifica compuestos que revierten la *firma global* del tumor, independientemente del mecanismo. Esto permite descubrir candidatos no obvios.

**Limitación importante**: CMap2 (EH3223) usa datos transcriptómicos de líneas celulares, no de tejido tumoral primario. La correlación entre firma proteómica tumoral primaria y perfil transcriptómico en líneas celulares es imperfecta, por lo que el score CMap se usa como una dimensión adicional de evidencia y no como único criterio.

---

## Fase 5: Integración de Fuentes de Fármacos

### Script 08 — `08_integrate_drug_targets.R`

**¿Qué hace?** Une los resultados de las cuatro fuentes en una tabla maestra unificada y asigna una clasificación de prioridad clínica a cada fármaco.

**¿Por qué integrar en lugar de tomar solo la "mejor" fuente?** Cada base de datos cubre una fracción distinta del espacio farmacológico y usa criterios de evidencia distintos. Un fármaco respaldado por múltiples fuentes independientes tiene mayor validez que uno presente en una sola base. La integración permite también combinar tipos de evidencia complementarios: DGIdb aporta amplitud, ChEMBL aporta fase clínica, Open Targets aporta relevancia en HNSCC, y CMap aporta reversión transcriptómica.

**Clasificación A/B/C/D:**

| Clase | Criterio | Interpretación |
| ----- | -------- | -------------- |
| A | Aprobado para HNSCC | Mayor prioridad: evidencia directa en la indicación |
| B | Aprobado para otro cáncer, fase III+ | Alta prioridad: evidencia oncológica sólida |
| C | Aprobado en indicación no oncológica | Media: evidencia de seguridad, mecanismo a explorar |
| D | Solo experimental o CMap | Menor: hipótesis computacional, sin validación clínica |

Esta clasificación es útil porque un reposicionamiento de clase A requiere mucho menos trabajo de validación que uno de clase D.

---

## Fase 6: Red de Interacciones Proteína-Proteína

### Script 09 — `09_string_network.R`

**¿Qué hace?** Construye una red de interacciones proteína-proteína (PPI) entre los genes DE usando STRING y calcula métricas de centralidad en la red.

**¿Por qué una red PPI?** Las proteínas no actúan de forma aislada sino en redes funcionales. Una proteína que interacciona con muchas otras (alto grado de conexión) es candidata a ser un nodo crítico de la red, cuya inhibición o activación puede tener efectos cascada sobre múltiples vías. Estos nodos se denominan **hubs** y son dianas farmacológicas de especial interés porque la modulación de un hub puede afectar múltiples procesos disfuncionales simultáneamente.

**STRING v12** (*Search Tool for the Retrieval of Interacting Genes/Proteins*) agrega evidencia de interacciones de múltiples fuentes: co-expresión, fusiones génicas, vecindad genómica, co-mención en literatura, bases de datos experimentales y predicciones computacionales. Se usa un umbral de **combined_score >= 700** (escala 0-1000), que corresponde al nivel de **alta confianza** según los propios umbrales de STRING. Esto filtra interacciones poco respaldadas y retiene solo las más fiables.

**Métricas de centralidad calculadas con igraph:**

- **Degree (grado)**: número de proteínas con las que interacciona directamente. Mide la conectividad local.
- **Betweenness**: número de caminos más cortos entre pares de nodos que pasan por este nodo. Mide la influencia como "puente" en la red; nodos con alta betweenness controlan el flujo de información.
- **Eigenvector centrality**: similar al algoritmo PageRank; un nodo es más central si sus vecinos también son altamente conectados. Mide la influencia global considerando la calidad de las conexiones.

**Definición de hub**: top 10% de proteínas por grado. Los hubs druggables (hubs para los que existe un candidato farmacológico) son de especial interés ya que representan dianas accesibles farmacológicamente que además son centrales en la red tumoral.

---

## Fase 7: Scoring Multi-Criterio

### Script 10 — `10_prioritization_scoring.R`

**¿Qué hace?** Integra seis dimensiones de evidencia en un score compuesto para cada candidato fármaco y selecciona el Top 20.

**¿Por qué un scoring multi-criterio?** No existe una métrica única que capture todos los aspectos relevantes de un candidato para reposicionamiento. El scoring permite combinar evidencias de naturaleza distinta (cambio de expresión, relevancia clínica, conectividad en red, etc.) en una escala común, y asignar pesos que reflejen la importancia relativa de cada dimensión.

**Las seis dimensiones y su justificación:**

| Dimensión | Peso | ¿Qué captura? | ¿Por qué este peso? |
| --------- | ---- | ------------- | ------------------- |
| `s_logfc` — magnitud del cambio | 0.20 | Cuánto está alterada la diana en el tumor | Mayor cambio sugiere mayor relevancia de la proteína en el fenotipo tumoral |
| `s_sig` — significancia estadística | 0.15 | Confianza estadística del cambio | Complementa logFC; penaliza cambios con alta variabilidad |
| `s_clinical` — fase clínica | 0.20 | Avance del fármaco en desarrollo clínico | Fármacos más avanzados tienen mayor probabilidad de éxito futuro |
| `s_cmap` — score CMap | 0.15 | Potencial de reversión transcriptómica | Evidencia mecanística de que el fármaco puede revertir el estado tumoral |
| `s_pathway` — relevancia en rutas | 0.15 | Si el gen diana está en rutas enriquecidas en el tumor | Favorece dianas con contexto biológico relevante en HNSCC |
| `s_network` — centralidad en PPI | 0.15 | Si el gen es hub en la red tumoral | Favorece dianas cuya modulación puede tener impacto sistémico en la red |

Cada dimensión se normaliza a escala 0-1 antes de ponderar, para que dimensiones con distintas unidades sean comparables.

**Filtro de diversidad de targets**: para evitar que el Top 20 esté dominado por múltiples fármacos que actúan sobre el mismo gen (e.g., 5 inhibidores de EGFR), se aplica un **máximo de 3 candidatos por gen diana primario**. Esto asegura diversidad de mecanismos en el ranking final.

---

## Fase 8: Validación de Evidencia Clínica

### Script 11 — `11_clinicaltrials_pubmed.py`

**¿Qué hace?** Consulta ClinicalTrials.gov y PubMed para cada uno de los Top 20 candidatos, buscando ensayos clínicos y publicaciones científicas que los relacionen con HNSCC.

**¿Por qué validar con evidencia clínica?** El scoring anterior es puramente computacional y se basa en datos proteómicos y de bases de datos. La evidencia clínica externa (ensayos registrados, publicaciones) aporta validación independiente de que un candidato es farmacológicamente activo en un contexto relevante. Un fármaco con alto score computacional *y* ensayos clínicos activos en HNSCC tiene mucho mayor prioridad que uno con solo score computacional.

**ClinicalTrials.gov API v2**: se consulta por el nombre del fármaco en combinación con HNSCC o head and neck cancer. Se registra el número de ensayos, su fase (I/II/III) y si están activos o completados. La API v2 de ClinicalTrials.gov (lanzada en 2023) ofrece mejor estructuración de datos que la v1 y permite búsquedas más precisas.

**PubMed E-utilities**: se consultan publicaciones usando el nombre simplificado del fármaco (sin formas de sal como "-HCl" o "-sodium", que varían entre estudios) combinado con términos MeSH de HNSCC. Esto captura literatura preclínica y clínica relevante.

### Script 12 — `12_cosmic_overlap.py`

**¿Qué hace?** Cruza los genes diana DE con el **Cancer Gene Census** de COSMIC y la base de datos **NCG7** para identificar cuáles son genes driver conocidos del cáncer.

**¿Por qué el solapamiento con genes driver?** Un gen driver es aquel cuya mutación o alteración de expresión contribuye activamente a la progresión tumoral (a diferencia de los genes "passenger" que se alteran sin consecuencia funcional). Si la diana de un fármaco candidato es además un driver conocido del cáncer (o específicamente de HNSCC), la justificación biológica del reposicionamiento es mucho más sólida.

**COSMIC CGC** (*Cancer Gene Census*) es el catálogo más curado de genes driver en cáncer humano, mantenido por el Sanger Institute. **NCG7** (*Network of Cancer Genes*) es una base de datos complementaria que incluye criterios adicionales de clasificación como driver. El pipeline usa NCG7 de forma embebida y descarga el CGC si está disponible localmente (requiere login en cancer.sanger.ac.uk).

---

## Fase 9: Integración Final y Ranking

### Script 13 — `13_evidence_summary.R`

**¿Qué hace?** Combina el score multi-criterio con la evidencia clínica en un **score final integrado** y genera el reporte maestro.

**Fórmula del score final:**

```text
score_final = 0.60 × score_multicriterio + 0.40 × score_evidencia_clínica
```

El peso mayor al score multi-criterio (0.60) refleja que la base proteómica y bioinformática es el fundamento del análisis. El score de evidencia clínica (0.40) actúa como modificador: eleva candidatos con respaldo externo fuerte y penaliza los que carecen de evidencia publicada.

**Niveles de evidencia (1–4):**

| Nivel | Criterio |
| ----- | -------- |
| 1 | Aprobado específicamente para HNSCC |
| 2 | Aprobado para otro cáncer de cabeza y cuello, o ensayo fase III+ activo en HNSCC |
| 3 | Ensayo fase I/II en HNSCC, o aprobado en otra indicación oncológica |
| 4 | Solo evidencia preclínica o computacional en HNSCC |

Esta clasificación permite comunicar de forma clara a clínicos o revisores cuánta evidencia clínica tiene cada candidato, más allá del score numérico.

---

## Fase 10: Análisis de Sensibilidad

### Script 15 — `15_sensitivity_analysis.R`

**¿Qué hace?** Evalúa la **robustez** del ranking final repitiendo el scoring con 6 configuraciones distintas de pesos.

**¿Por qué un análisis de sensibilidad?** Los pesos asignados al scoring multi-criterio son decisiones metodológicas razonadas pero no absolutas. Distintos expertos podrían asignar pesos distintos. El análisis de sensibilidad responde la pregunta: *¿cambia significativamente el ranking si modificamos los pesos?* Si un candidato aparece en el Top 20 en casi todas las configuraciones de pesos probadas, es un candidato **robusto** cuya priorización no depende de decisiones metodológicas específicas. Si solo aparece en una configuración, es más **sensible** a los parámetros y debe interpretarse con más cautela.

Las 6 configuraciones varían sistemáticamente los pesos para favorecer distintos aspectos (e.g., una configuración prioriza más la fase clínica, otra la centralidad en red, etc.). El estadístico de robustez es simplemente el número de configuraciones en que cada candidato aparece en el Top 20.

---

## Fase 11: Figuras de Publicación

### Script 17 — `17_pub_figures.R`

**¿Qué hace?** Genera versiones de alta calidad de las figuras clave del análisis, listas para incluir en manuscritos o presentaciones.

**¿Por qué un script separado para figuras de publicación?** Los scripts de análisis priorizan la generación rápida de figuras exploratorias para verificar resultados. Las figuras de publicación requieren mayor control sobre paletas de color, tamaño de fuentes, proporciones, y eliminación de elementos visuales redundantes. Separar estas responsabilidades mantiene los scripts de análisis simples y reproducibles, y concentra las decisiones estéticas en un único lugar.

Se utilizan paletas **colorblind-safe** (Okabe-Ito, okabe8) para garantizar que las figuras sean accesibles para personas con deficiencia en la percepción del color, requisito estándar en journals científicos modernos. Todas las figuras se exportan en **PDF** (vectorial, para edición posterior) y **PNG a 300 DPI** (resolución mínima estándar para publicación impresa).

---

## Software y Reproducibilidad

### Lenguajes y entornos

| Ambiente | Lenguaje | Scripts |
| -------- | -------- | ------- |
| `omics-R` (conda) | R 4.4.3 | 01, 02, 03, 07, 08, 09, 10, 13, 14, 15, 17 |
| `omics-py` (conda) | Python 3 | 04, 05, 06, 11, 12 |

### Paquetes R principales

| Paquete | Versión | Uso |
| ------- | ------- | --- |
| limma | 3.62.2 | Análisis diferencial |
| clusterProfiler | 4.14.6 | ORA y GSEA |
| ReactomePA | 1.50.0 | Enriquecimiento Reactome |
| signatureSearch | 1.20.0 | Análisis CMap2 |
| igraph | 2.2.1 | Red PPI y métricas |
| ggraph / tidygraph | 2.2.2 | Visualización de redes |
| ComplexHeatmap | 2.22.0 | Heatmaps anotados |
| openxlsx | 4.2.8.1 | Exportación Excel |

### Bases de datos consultadas

| Base de datos | Versión/Fecha | Acceso |
| ------------- | ------------- | ------ |
| DGIdb | v5 | GraphQL API |
| ChEMBL | v33 | REST API |
| Open Targets Platform | Mar 2026 | GraphQL API |
| CMap2 | EH3223 | ExperimentHub (R) |
| STRING | v12 | REST API |
| ClinicalTrials.gov | API v2 | REST API |
| NCBI PubMed | - | E-utilities |
| MSigDB Hallmarks | v2024.1 | msigdbr (R) |
| NCG7 | v7.1 | Embebido en script 12 |

### Parámetros del análisis

Todos los parámetros (umbrales, pesos, filtros) están centralizados en `config/analysis_params.yaml`. Modificar ese archivo es la forma correcta de reproducir variantes del análisis sin editar los scripts directamente.

### Reproducción completa

Para reproducir el análisis desde cero:

```bash
# Ver docs/RUNBOOK.md para instrucciones paso a paso
cd ~/bioinfo/projects/hnscc_drug_repurposing
# Seguir las fases 1–8 del RUNBOOK en orden
```

El orden de ejecución de los scripts respeta sus dependencias (ver sección "Dependencias entre scripts" del RUNBOOK). Cada script genera un log con timestamp en `logs/`.

---

*Para el catálogo completo de archivos generados por cada script, ver `docs/OUTPUTS.md`.*
*Para instrucciones de ejecución paso a paso, ver `docs/RUNBOOK.md`.*
