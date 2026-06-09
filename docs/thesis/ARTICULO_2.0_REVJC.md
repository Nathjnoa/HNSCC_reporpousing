FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

I. TÍTULO

Reposicionamiento in silico de fármacos para el tratamiento del cáncer

escamocelular de cabeza y cuello.

II. AUTORES

Jonathan  Carvajal  Veloza  Instructor  Asistente,  Investigador  grupo  de  Ciencias
Básicas en Salud CBS-FUCS. Universidad FUCS

Luz  Dary  Gutierrez  Castañeda  Instituto  de  investigación,  Vicerrectoría  de

investigaciones. Grupo de Ciencias Básicas-FUCS, Universidad FUCS

Alvaro Eduardo Granados Calixto  Docente especialización Cirugía de Cabeza y
Cuello. Universidad FUCS.

Juan Manuel Marquez Duque Estudiante especialización Cirugía Cabeza y Cuello.
FUCS
Universidad

Jaime Rafael Montero Arrieta Estudiante de la Especialización en Cirugía general.
Universidad FUCS.

Ana Maria Martin Pinilla Estudiante de Pregrado Medicina. Universidad FUCS .

III. DEPARTAMENTO

- Instituto de Investigaciones. Instituto de ciencias Básicas Universidad FUCS.
- Departamento de Cirugía de Cabeza y Cuello Universidad FUCS. - Sociedad

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

de Cirugía de Bogotá, Hospital de San José

IV. DIRECCIÓN – CONTACTO

Jonathan Carvajal Veloza

Email:  jcarvajal@fucsalud.edu.co

Instituto de Ciencias Básicas. Universidad FUCS
Teléfono: 3016837292

V. CARACTERÍSTICAS

Número de palabras máximas: 4.500 entre resumen y cuerpo.
Número de figuras 07
Número de tablas 07
Número de referencias 21

VI. RESUMEN

Introducción:

El carcinoma escamocelular de cabeza y cuello (CECC) constituye una neoplasia
altamente  heterogénea  con
insatisfactorios,
especialmente  en  enfermedad  avanzada  o  recurrente.  La  identificación  de
nuevas  estrategias  terapéuticas mediante  enfoques  de  medicina  de  precisión
es una necesidad urgente.

resultados  clínicos  aún

Comentado [1]: Quitar referencias, los resumenes no
llevan referencias

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

Objetivo:
Evaluar  in  silico  el  potencial  de  reposicionamiento  de  fármacos  mediante  la
integración  de  datos  proteómicos,  análisis  de  conectividad  farmacológica  y

modelos de priorización multicriterio en CECC.

Métodos:

in  silico  utilizando  datos
Se  realizó  un  estudio  observacional  analítico
proteómicos  cuantitativos  derivados  de  muestras  pareadas  tumor-tejido

normal.  Se  identificaron  proteínas  diferencialmente  expresadas  (FDR  <  0,05;
|log2FC| > 1) y se construyeron firmas moleculares tumorales. Estas firmas fueron

interrogadas  en  plataformas  de  conectividad  farmacológica  (L2S2).  Se
realizaron análisis de enriquecimiento funcional (FGSEA), mapeo fármaco–diana
(DrugBank,  ChEMBL,  OpenTargets)  y  análisis  de  redes.  Los  candidatos  se
priorizaron  mediante  un  modelo  multicriterio  integrando  conectividad,  soporte
mecanístico y factibilidad clínica.

Resultados:
Se identificaron  666  proteínas  diferencialmente  expresadas,  evidenciando  una
profunda reprogramación proteómica tumoral, consistente con estudios previos

en  CECC.  El  análisis  de  redes  reveló  módulos
funcionales  altamente
interconectados,  destacando  señalización  EGFR,  proteostasis,  metabolismo  y

remodelación  de  matriz  extracelular.  Se
identificaron  552  compuestos
candidatos, de los cuales una proporción significativa correspondió a fármacos
aprobados.  El  eje  EGFR/ErbB  emergió  como  el  principal  nodo  terapéutico,
validando la aproximación en concordancia con su papel conocido en CECC. De
manera  relevante,  se
incluyendo
inhibidores  de  DNMT  (decitabina,  azacitidina),  inhibidores  del  proteasoma
(carfilzomib)  y  moduladores  metabólicos,  sugiriendo  vulnerabilidades

identificaron  candidatos  no  canónicos,

terapéuticas previamente subexploradas.

Conclusiones:
El  reposicionamiento farmacológico  guiado  por  proteómica  permite  identificar

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

de  manera  robusta  blancos  terapéuticos  en  CECC,  integrando  coherencia
biológica  con  aplicabilidad  clínica.  Este  enfoque  representa  una  herramienta
estratégica  para  acelerar  la  implementación  de  medicina  de  precisión  en

oncología de cabeza y cuello.

Palabras clave:  cáncer  de  cabeza  y  cuello,  proteómica,  reposicionamiento  de

fármacos, EGFR, bioinformática, medicina de precisión

VII. PALABRAS CLAVE

Cáncer  de cabeza  y cuello,  proteómica,  reposicionamiento  de  fármacos,  EGFR,

bioinformática, medicina de precisión

VIII. INTRODUCCIÓN

El carcinoma escamocelular de cabeza y cuello (CECC) representa una de las

principales causas de morbilidad y mortalidad por cáncer a nivel mundial, con
una incidencia en aumento y una supervivencia global que permanece limitada,

especialmente  en  estadios  localmente  avanzados  y  enfermedad  recurrente
(1,3,19). A pesar de los avances en cirugía oncológica, radioterapia, quimioterapia
e  inmunoterapia,  los  resultados  clínicos  siguen  siendo  subóptimos  (5,6,20),

reflejando la complejidad biológica de esta enfermedad.

En  la  última  década,  el  manejo  del  CECC  ha  experimentado  un  cambio  de
paradigma  con  la  llegada  de  los  inhibidores  de  puntos  de  control  inmunitario
(ICIs)(23). El ensayo clínico fase 3 KEYNOTE-048 estableció a pembrolizumab, ya
sea  como  monoterapia  en  pacientes  con  expresión  positiva  del  ligando  1  de

muerte programada (PD-L1) o en combinación con quimioterapia (platino y 5-
FU),  como  el  nuevo  estándar  de  primera  línea  para  enfermedad  recurrente  o

metastásica  (4).  La  selección  de  estos  pacientes  depende  crucialmente  de  la
Puntuación Positiva Combinada (CPS), un biomarcador que integra la expresión

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

de PD-L1 tanto en células tumorales como en células inmunes asociadas (7). No
obstante, en el escenario del CECC localmente avanzado y resecable —campo
de acción primordial del cirujano de cabeza y cuello—, los avances han sido más

lentos (10,11). Recientemente, el estudio KEYNOTE-689 ha demostrado que el uso
de pembrolizumab perioperatorio (neoadyuvante y adyuvante) en combinación

con el estándar de cuidado mejora significativamente la supervivencia libre de
eventos (SLE) (12). Estos hitos clínicos validan la necesidad urgente de explorar

nuevas firmas moleculares y dianas terapéuticas que permitan optimizar estas
estrategias de tratamiento sistémico y perioperatorio (11,15).

Uno  de  los  principales  desafíos  en  el  manejo  del  CECC  radica  en  su  marcada
heterogeneidad  molecular.  La  distinción  entre  tumores  asociados  al  virus  del
papiloma  humano  (VPH)  y  aquellos  relacionados  con  carcinógenos  clásicos
como  tabaco  y  alcohol  ha  permitido  reconocer  subgrupos  con  diferencias
pronósticas y biológicas sustanciales (6,21). Sin embargo, esta estratificación es

aún insuficiente para capturar la complejidad funcional del tumor.

En  este  contexto,
la  proteómica  ha  emergido  como  una  herramienta
fundamental para caracterizar el fenotipo funcional tumoral. A diferencia de la

genómica  o  transcriptómica,  la  proteómica  refleja  directamente  la  actividad
integrando  múltiples  niveles  de  regulación  postranscripcional  y
celular,

postraduccional (7,8). Estudios recientes han demostrado que el CECC presenta
rutas  críticas  como  señalización  EGFR,
alteraciones  significativas  en
metabolismo energético, regulación epigenética y microambiente tumoral (9–

13,16–18).

Paralelamente,  el  reposicionamiento  farmacológico  se  ha  consolidado  como
una  estrategia  eficiente  para  acelerar  la  identificación  de  nuevas  opciones
terapéuticas,  aprovechando  compuestos  ya  aprobados  con  perfiles  de
seguridad conocidos (7,8). La integración de firmas moleculares tumorales con

plataformas  de  conectividad  farmacológica  permite

identificar  fármacos

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

capaces  de  revertir  estados  patológicos,  ofreciendo  una  aproximación
sistemática a la medicina de precisión.

El  presente  estudio  propone  un  enfoque  integrativo  que  combina  datos
proteómicos,  análisis  bioinformáticos  avanzados  y  evidencia  farmacológica
para identificar candidatos terapéuticos en CECC. Esta estrategia busca no solo

validar blancos conocidos, sino también descubrir nuevas vulnerabilidades con
potencial traslacional inmediato.

IX. MÉTODOS

Se realizó un estudio observacional analítico in silico basado en la integración de
datos proteómicos y farmacológicos de fuentes secundarias (7,8). Se analizaron
datos proteómicos cuantitativos obtenidos de 10 pares de muestras tumor-tejido
normal,  con  un  total  de  3352  proteínas  cuantificadas.  Se  realizó  análisis
diferencial  utilizando  corrección  por  Benjamini-Hochberg,  considerando

significativo FDR < 0,05 y |log2FC| > 1 (7–9).

Se  utilizó  FGSEA  para  identificar  rutas  biológicas  alteradas,  empleando  bases

Hallmark  y  Reactome,  ampliamente  utilizadas  en  análisis  funcional  de  datos
ómicos (8,17). Las proteínas diferencialmente expresadas se cruzaron con cuatro

bases de datos farmacológicas: DGIdb v5, ChEMBL 34 (fármacos en fase clínica
≥
III),  Open  Targets  Platform  y  L2S2  (reversión  de  firma  molecular).  Se
seleccionaron  compuestos  con  puntaje  de  conectividad  inversa  ≤ −90  en L2S2,
indicativo  de  reversión  del  fenotipo  tumoral.  Se  priorizaron  candidatos
independientes.
respaldados

fuentes

menos

dos

por

al

Adicionalmente,  con  las  proteínas  diferencialmente  expresadas,  se  construyó
una red de interacción proteína-proteína con STRING v12 (umbral de confianza ≥

700) para evaluar la posición topológica de las dianas moleculares(13,17).

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

Para integrar información heterogénea proveniente de múltiples fuentes en una
métrica  única,  se  construyó  un  puntaje  compuesto  de  cinco  dimensiones  con
pesos  diferenciados  según  su  relación  directa  con  la  biología  tumoral  e

información en bases de datos (Tabla X). El π-estadístico proteómico recibió el
peso  más  alto  (0.40)  por  ser  la  información  base  del  estudio  (las  proteínas

obtenidas  de
tumorales),  este  estadístico  captura
simultáneamente  la  magnitud,  la  significancia  estadística  y  la  dirección  del

las  muestras

cambio  de  expresión  en  tejido  tumoral  de  HNSCC,  por  lo  cual  es  el  más
importante para la priorización. La fase clínica FDA (0.20) refleja la seguridad y
viabilidad  traslacional  del  compuesto  con
independencia  del  contexto
oncológico específico. La relevancia de vías enriquecidas (0.15) y la centralidad
en  la  red  de  interacción  proteína-proteína  (0.15)  aportan  contexto  biológico
sistémico: el primero evalúa si el blanco actúa dentro de procesos desregulados

en  el  tumor,  y  el  segundo  si  ocupa  una  posición  topológicamente  estratégica
susceptible de amplificar el efecto terapéutico. La conectividad transcriptómica

derivada  de  L2S2  (base  de  datos,  https://l2s2.maayanlab.cloud/)  (0.10)
contribuye  evidencia  de  reversión  farmacológica  de  la  firma  molecular,  pero

recibió un peso menor dado que opera sobre transcriptómica de líneas celulares,
una  capa  de  datos  más  distal  a  la  proteómica  tisular  del  estudio.  Los  pesos
fueron fijados a priori antes de ordenar los candidatos y su robustez fue evaluada
mediante análisis de sensibilidad de pesos y análisis leave-one-dimension-out
(LOD).

Tabla  X.  Dimensiones  y  pesos  del  puntaje  compuesto  de  priorización  de
candidatos farmacológicos.

Dimensión

Peso  Descripción

π-estadístico
proteómico

0.40  sign(log₂FC) × |log₂FC| × −log₁₀(adj.P); dirección,
la  expresión

magnitud  y  significancia  de
diferencial

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

Fase clínica FDA

0.20  Aprobado=1;  Fase  III=0.75;  Fase  II=0.50;  Fase

I=0.25; preclínico=0.00

Relevancia de vías

0.15

Presencia  del  blanco  primario  en  una  vía
significativamente enriquecida en HNSCC

Centralidad de red PPI  0.15  Centralidad  de  intermediación  (betweenness)

del  blanco  primario  en
normalizada al intervalo [0, 1]

la

red  STRING,

Conectividad

0.10

Puntaje  de  reversión  de  firma  génica  (más

transcriptómica
(L2S2)

negativo = mayor reversión); derivado de L2S2

X. RESULTADOS

El  análisis  proteómico  cuantitativo  incluyó  3.352  proteínas  identificadas  en  10

pares  de  muestras  tumor-tejido  normal.  Tras  aplicar  criterios  de  significancia
(FDR  <  0,05;  |log2FC|  >  1),  se  identificaron  666  proteínas  diferencialmente
expresadas, correspondientes al 19,9% del proteoma cuantificado.

De  estas,  329  proteínas  (49,4%)  se  encontraron  sobreexpresadas  en  tejido
tumoral, mientras que 337 (50,6%) se encontraban subexpresadas, evidenciando

una reprogramación proteómica global del CECC (Tabla 1).
El volcano plot evidenció una clara separación entre proteínas diferencialmente

expresadas,  destacando  moléculas  con  relevancia  terapéutica  como  EGFR,
DNMT1 y componentes del proteasoma (Figura 1). Asimismo, el mapa de calor de

las 40 proteínas más diferencialmente expresadas permitió discriminar de forma

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

consistente entre muestras tumorales y normales, evidenciando la robustez del
perfil molecular (Figura 2).

las  proteínas  más  sobreexpresadas  se

Entre
identificaron  componentes
relacionados con matriz extracelular y señalización tumoral, incluyendo LAMB3,

LAMA3, SERPINE1 y PTK7, mientras que las proteínas subexpresadas se asociaron
predominantemente a funciones musculares y metabólicas, como MYH2, MYL2,

ATP2A1 y CASQ1 (Tabla 2).

2.  Identificación  de  candidatos  farmacológicos  mediante  integración  multi-
fuente

La integración  de  bases  de  datos  farmacológicas  (DGIdb,  OpenTargets,  L2S2  y

ChEMBL) permitió identificar un conjunto de fármacos candidatos asociados a
proteínas diferencialmente expresadas.

El  análisis  comparativo  entre  fuentes  mostró  una  cobertura  complementaria,

con  variabilidad  en  el  número  de  compuestos  identificados  por cada  base  de
datos (Figura 3). Al restringir el análisis a compuestos presentes en al menos dos
fuentes,  se  identificaron  458  candidatos  farmacológicos,  lo  que  incrementa  la
confiabilidad de los hallazgos.
De estos, el 75% correspondieron a fármacos aprobados (Fase IV), lo que resalta
el alto potencial de reposicionamiento clínico inmediato (Figura 4).
El ranking por número de fuentes de respaldo identificó múltiples fármacos con
evidencia consistente, destacando:

● Inhibidores EGFR: afatinib, lapatinib, gefitinib
● Agentes epigenéticos: azacitidina, decitabina
● Fármacos reposicionables no oncológicos: disulfiram, acetazolamida

3. Arquitectura de red y módulos funcionales en CECC

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

El análisis de interacción proteína-proteína reveló una red altamente organizada
con múltiples nodos hub, sugiriendo la existencia de ejes biológicos críticos en la

tumorogénesis del CECC (Figura 5).

A partir de esta red, se identificaron cinco módulos funcionales principales que
concentran la mayoría de blancos terapéuticos:

1. Señalización EGFR
2. Fosforilación oxidativa (OXPHOS)
3. Proteasoma y degradación proteica
4. Respuesta inmune y matriz extracelular
5. Metabolismo de pequeñas moléculas

El módulo de respuesta inmune/matriz extracelular (n=64) y el de EGFR (n=57)

concentraron  el  mayor  número  de  fármacos  candidatos,  lo  que  sugiere  su
relevancia como nodos terapéuticos prioritarios (Figura 6).

Esta notable concentración de candidatos en el módulo de 'Respuesta inmune y
matriz extracelular' (n=64) es biológicamente coherente con el éxito clínico de
los  inhibidores  de  PD-1  observados  en  la  práctica  oncológica  actual  (16,17).  La
identificación  de  proteínas  relacionadas  con  la  remodelación  de  la  matriz
sugiere  que  el  microambiente  tumoral  no  solo  actúa  como  una  barrera  física,
sino como un modulador activo que facilita la evasión inmune en el CECC (17,18).
Asimismo, la estabilidad en el ranking de fármacos como la decitabina (puntaje

0.715)  adquiere  una  relevancia  traslacional  inmediata  ante
la  evidencia
emergente  que  propone  el  uso  de  agentes  epigenéticos  para  'sensibilizar'  el

tumor (19,20). Estos fármacos podrían potenciar la inmunogenicidad y mejorar
la respuesta a los ICIs, abriendo una vía estratégica para superar la resistencia
intrínseca  observada  en  una  proporción  importante  de  pacientes  tratados
actualmente con pembrolizumab (20,21).

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

4. Clasificación clínico-regulatoria de candidatos

El  cruce  entre  proteínas  hub  y  bases  farmacológicas  permitió  identificar  552
fármacos únicos, los cuales fueron clasificados según su estado regulatorio:

● Clase A (aprobados en CECC): 11 (2,0%)
● Clase B (aprobados en otros cánceres): 31 (5,6%)
● Clase C (aprobados no oncológicos): 62 (11,2%)
● Clase D (experimentales): 448 (81,2%)

5. Priorización multicriterio de candidatos farmacológicos
Se  implementó  un  modelo  de  priorización  integrando  señal  proteómica,  fase

clínica,  relevancia  de  vías,  centralidad  en  red  y  conectividad  transcriptómica
(Tabla 4).

5.1 Dominancia del eje EGFR como validación del modelo

El análisis identificó una clara predominancia de fármacos dirigidos al eje
EGFR/ErbB, los cuales se mantuvieron consistentemente en el top del ranking bajo
diferentes condiciones de análisis (Tabla 5).

Entre los candidatos más robustos se incluyen:
● Gefitinib (puntaje 0,777)
● Erlotinib (0,772)
● Neratinib (0,753)
● Dacomitinib (0,739)

Además, la presencia de cetuximab y afatinib, ambos con indicación en CECC,
valida la coherencia biológica del pipeline.

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

5.2 Identificación de candidatos no-EGFR: principal aporte del estudio
De particular relevancia, el análisis permitió identificar fármacos no dirigidos a
EGFR que permanecieron estables en el ranking (Tabla 6):

● Decitabina (0,715) – inhibidor DNMT
● Azacitidina (0,669) – epigenético
● Carfilzomib (0,649) – inhibidor del proteasoma
● Tranylcypromine (0,641) – modulador epigenético

Estos  hallazgos
particularmente en regulación epigenética y proteostasis.

vulnerabilidades

sugieren

terapéuticas  alternativas,

5.3 Candidatos exploratorios y reposicionamiento clásico

El análisis extendido identificó fármacos robustos en análisis de sensibilidad pero

dependientes de parámetros específicos (Tabla 7), incluyendo:

● Disulfiram
● Ácido valproico
● Digoxina
● Acetazolamida

Estos  compuestos  representan  oportunidades  de  reposicionamiento  clásico,
especialmente relevantes por su disponibilidad clínica.

XI. DISCUSIÓN

El carcinoma de células escamosas de cabeza y cuello (HNSCC) es un término

que agrupa diversas neoplasias malignas que se originan en la cavidad oral, la
faringe,  la  hipofaringe,  la  laringe,  la  cavidad  nasal  y  las  glándulas  salivales
(20,35,36).  En  este  contexto,  las  características  mutacionales  y  los  perfiles  de

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

expresión  génica  asociados  a  estas  neoplasias  presentan  una  marcada
heterogeneidad. Por lo tanto, la identificación de biomoléculas que puedan ser
utilizadas como potenciales blancos farmacológicos resulta de gran relevancia

y  utilidad  clínica  (9–11,16).  En  el  presente  trabajo,  el  análisis  de  las  proteínas
diferencialmente  expresadas  en  muestras  tumorales  de  pacientes,  en

comparación  con  muestras  de  tejido  sano  provenientes  del  mismo  individuo,
permitió identificar un conjunto de proteínas con potencial para ser propuestas

como blancos moleculares.

Las  vías  de  señalización  en  las  que  se  encuentran  los  principales  blancos

moleculares  hacen  parte  de  señalización  para  proliferación  (EGFR)  y  de
degradación de proteínas. La identificación de módulos funcionales altamente
interconectados, particularmente aquellos relacionados con señalización EGFR,
remodelación  de  matriz  extracelular,  metabolismo  energético  y  proteostasis,
refuerza  la  noción  de  que  la  biología  del  CECC  está  gobernada  por  redes

dinámicas más que por nodos aislados (17,18,35,36).

La  consistencia  del  eje  EGFR/ErbB  como  principal  nodo  terapéutico  constituye
una  validación  interna  crítica  del  enfoque  metodológico.  La  identificación

reiterada de inhibidores de EGFR en los primeros lugares del ranking, junto con la
recuperación de fármacos ya aprobados en CECC como cetuximab, sugiere que

el pipeline no solo es capaz de capturar dependencias biológicas conocidas, sino
también  de  jerarquizarlas  adecuadamente.  Sin  embargo,  más  allá  de  esta
validación, estos resultados invitan a reconsiderar el papel de EGFR no como un
blanco  aislado,  sino  como  un  nodo  central  dentro  de  una  red  de  señalización
más  amplia,  potencialmente  susceptible  a  estrategias  combinatorias  que

aborden mecanismos de resistencia intrínseca y adquirida (5,6,14,35).

Esta  validación  adquiere  una  relevancia  crítica  ante  el  cambio  de  paradigma
terapéutico  en  el  CECC  recurrente  o  metastásico  (R/M).  Históricamente,  el

régimen EXTREME (cetuximab, platino y 5-FU) fue el estándar de cuidado (24). Sin
embargo, los resultados del ensayo KEYNOTE-048 han desplazado este enfoque,

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

demostrando  que  pembrolizumab  (monoterapia  o  con  quimioterapia)  ofrece
una supervivencia global superior (25,27). Es imperativo notar que el beneficio de
pembrolizumab  es  dependiente  del  biomarcador  PD-L1, medido  a  través  de  la

Puntuación Positiva Combinada (CPS) (22,25). Nuestros hallazgos, que destacan
el  módulo  de  “respuesta  inmune  y  matriz  extracelular”,  se  alinean  con  esta

transición hacia la inmunoterapia, sugiriendo que las proteínas identificadas en
este estudio podrían complementar la selección de pacientes más allá del CPS

(15,17).  La  priorización  de  inhibidores  de  DNMT  cobra  un  sentido  traslacional
mayor al considerar la resistencia a inhibidores de puntos de control inmunitario
(ICIs). Estudios sugieren que la reprogramación epigenética puede aumentar la
tumoral  (30,31).  La  combinación  de  decitabina  con
inmunogenicidad
inmunoterapia  busca  no  solo  inhibir  la  proliferación,  sino  “sensibilizar”  el
respuesta  en  pacientes  que
microambiente

tumoral,  potenciando

la

actualmente no se benefician de los ICIs de manera aislada (15,30).

El  aporte  más  significativo  del  estudio  radica  en

la

identificación  de

vulnerabilidades  no  canónicas,  particularmente  aquellas  relacionadas  con
regulación epigenética y homeostasis proteica. La priorización de inhibidores de

DNA  metiltransferasa  como  decitabina  y  azacitidina  sugiere  que  la  disrupción
epigenética no es simplemente un fenómeno secundario, sino un componente

estructural  del  estado  tumoral  en  CECC.  Este  hallazgo  es  consistente  con  la
creciente  evidencia  de  que
la
adaptación tumoral, la evasión inmune y la resistencia terapéutica (17,30,31). En
este sentido, el reposicionamiento de agentes epigenéticos podría representar
una  estrategia  para  reprogramar  estados  celulares  malignos,  más  que

la  plasticidad  epigenética  contribuye  a

simplemente inhibir vías proliferativas.

La  identificación  de  inhibidores  del  proteasoma  como  carfilzomib  sugiere  una
dependencia  del  CECC  de  mecanismos  de  proteostasis.  La  acumulación  de

proteínas mal plegadas y la necesidad de mantener un equilibrio proteico en un
entorno  de  alta  demanda  biosintética  podrían  constituir  un  talón  de  Aquiles

explotable terapéuticamente. Este eje ha sido relativamente poco explorado en

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

CECC  en  comparación  con  otras  neoplasias  hematológicas,  lo  que  posiciona
estos hallazgos como una oportunidad particularmente innovadora (17,18).

Otro  aspecto  destacable  es  la  identificación  de  fármacos  aprobados  para
indicaciones no oncológicas dentro de los candidatos priorizados. La presencia
de compuestos como disulfiram, ácido valproico o acetazolamida sugiere que

la  biología  tumoral  comparte  nodos  funcionales  con  procesos  fisiológicos
aparentemente no relacionados. Esta convergencia abre la puerta a estrategias

de  reposicionamiento  clásico  con  un  alto  potencial  de  traslación  clínica,
especialmente en contextos donde el acceso a terapias innovadoras es limitado.

Desde una perspectiva de salud global, este punto adquiere especial relevancia

(7,8,33).

Un  aspecto  relevante  que  refuerza  la  plausibilidad  traslacional  de  estos
hallazgos es la evidencia emergente en ensayos clínicos en curso. La consulta de
registros  como  ClinicalTrials.gov  evidencia  que  varios  de
los  fármacos
priorizados no solo han sido evaluados en otros contextos oncológicos, sino que
actualmente  se  encuentran  en  investigación  en  cáncer  escamocelular  de
la
cabeza  y  cuello  (32).  Por  ejemplo,  estudios  fase

I  están  evaluando

combinación  de  decitabina  con  terapias  estándar  en  enfermedad  resecable
HPV-negativa, con el objetivo de potenciar la sensibilidad tumoral a tratamientos

convencionales. De igual forma, ensayos recientes exploran la combinación de
decitabina  con  inmunoterapia,  como  nivolumab,  en  enfermedad  recurrente  o
metastásica,  buscando  aumentar  la  inmunogenicidad  tumoral  y  mejorar  la
respuesta  clínica  (26,32).  Estos  enfoques  reflejan  una  tendencia  actual  en
oncología  hacia  estrategias  combinatorias  basadas  en  la  modulación  del
microambiente  tumoral  y  la  reprogramación  epigenética.  En  paralelo,  aunque
fármacos como disulfiram aún no cuentan con ensayos específicos en CECC, sí
han sido evaluados en estudios clínicos en otros tumores sólidos y sarcomas, lo

que respalda su seguridad y su potencial de reposicionamiento (33). En conjunto,
estos  hallazgos  sugieren  que  varios  de  los  candidatos  identificados  en  este

estudio no se encuentran en una etapa puramente teórica, sino que ya forman

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

parte del ecosistema actual de investigación clínica en oncología, lo que facilita
su eventual traslación al manejo del cáncer de cabeza y cuello.

El enfoque  basado  en  redes  utilizado  en este  estudio  también  permite  superar
una  limitación  conceptual  importante  en  oncología:  la  tendencia  a  priorizar
blancos  moleculares  individuales.  La  evidencia  generada  sugiere  que  la

efectividad  terapéutica  podría  depender  más  de  la  modulación  de  módulos
funcionales  completos  que  de  la  inhibición  de  proteínas  específicas  (13,17,35).

Esto  tiene  implicaciones  directas  en  el  diseño  de  terapias  combinadas
racionales, particularmente en tumores complejos como el CECC.

informar

relevantes  para

la  práctica  quirúrgica  oncológica.

Desde  una  perspectiva  traslacional,  los  hallazgos  de  este  estudio  tienen
implicaciones
La
identificación  de  subgrupos  moleculares  con  vulnerabilidades  específicas
la  selección  de  terapias  neoadyuvantes  o  adyuvantes,
podría
optimizando  resultados  quirúrgicos  y  reduciendo  recurrencias.  Asimismo,  el
reposicionamiento  de  fármacos  con  perfiles  de  seguridad  conocidos  podría
facilitar la implementación de estrategias terapéuticas en tiempos clínicamente
relevantes  (6,20,35).  Para  el  cirujano  de  cabeza  y  cuello,  el  escenario

perioperatorio  está  evolucionando  rápidamente.  Estudios  recientes  sobre
inmunoterapia  perioperatoria  con  pembrolizumab  han  demostrado  beneficios

en  supervivencia  libre  de  eventos  y  respuesta  patológica  en  CECC  resecable,
validando  el  papel  creciente  de  las  estrategias  neoadyuvantes  basadas  en
inmunoterapia (28,29). Este avance respalda nuestra propuesta de utilizar firmas
proteómicas para informar la selección de terapias neoadyuvantes, optimizando
la viabilidad de la resección y reduciendo la probabilidad de recurrencia local y

a distancia (17,28).

Además  de  los  marcadores  moleculares,  factores  clínicos  como  la  carga
tumoral  desempeñan  un  papel  pronóstico  esencial.  Evidencia  reciente  indica

que  el  número  de  lesiones  metastásicas  (BNML),  el  diámetro  de  las  lesiones
(BSLD)  y  el  valor  de  captación  estándar  máximo  (SUVmax)  en  PET/CT  están

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

inversamente  correlacionados  con  la  eficacia  de  los  ICIs  (23).  Por  tanto,  la
integración  de  la  proteómica  funcional  con  la  evaluación  imagenológica
detallada pretratamiento es fundamental para predecir qué pacientes con alta

carga  tumoral  podrían  requerir  terapias  citotóxicas  concomitantes  con
inmunoterapia para lograr el control de la enfermedad (22,23).

No  obstante,  estos  resultados  deben  interpretarse  en  el  contexto  de  las
limitaciones inherentes a un estudio in silico. Aunque la integración de múltiples

capas de evidencia proteómica, farmacológica y de redes refuerza la robustez
de  los  hallazgos,  la  validación  funcional  en  modelos  preclínicos  y  clínicos  es

indispensable  (9,16).  Adicionalmente,  la  naturaleza  estática  de  los  datos
proteómicos limita la captura de dinámicas temporales y de interacción con el
microambiente tumoral (17,18).

En conjunto, este estudio propone un cambio de paradigma en la identificación
de  blancos terapéuticos en  CECC,  pasando  de  un  enfoque centrado  en  genes
individuales a uno basado en estados funcionales integrados. La convergencia
entre  coherencia  biológica,  robustez  computacional  y  viabilidad  clínica
posiciona  al  reposicionamiento  farmacológico  guiado  por  proteómica  como

una estrategia prometedora  para acelerar la implementación  de medicina  de
precisión en oncología de cabeza y cuello (7,8,15,35,36).

En  este  sentido,  la  integración  de  evidencia  molecular  con  datos  clínicos
emergentes  posiciona  a  los  candidatos  identificados  no  solo  como  hallazgos

teóricos,  sino  como  oportunidades  reales  para  el  desarrollo  de  estrategias
terapéuticas en cáncer escamocelular de cabeza y cuello.

Conclusión

El presente estudio demuestra que el análisis integrativo basado en proteómica
permite  caracterizar  de  manera
la  biología  del  carcinoma
escamocelular  de  cabeza  y  cuello,  identificando  no  solo  las  dependencias

funcional

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

oncogénicas  clásicamente  descritas  (5,6,20),  sino  también  vulnerabilidades
terapéuticas  previamente  subexploradas  (7,8).  La  consistencia  del  eje  EGFR
como nodo central valida la solidez del enfoque, mientras que la identificación

de  rutas  emergentes,  particularmente  aquellas  relacionadas  con  regulación
epigenética  y  proteostasis,  amplía  de  manera  significativa  el  panorama

terapéutico potencial en esta enfermedad (17,18).

Nuestros resultados posicionan al reposicionamiento farmacológico guiado por

proteómica no  solo  como  una  alternativa  teórica,  sino  como  una  herramienta
práctica de medicina de precisión quirúrgica (22,29). La identificación de dianas

epigenéticas  y  de  proteostasis  ofrece  una  oportunidad  estratégica  para  el
tratamiento  neoadyuvante,  buscando  convertir  escenarios  quirúrgicos  de  alto
riesgo en intervenciones con mayores probabilidades de éxito curativo (12,30). En
conclusión,  la  alineación  entre  las  firmas  moleculares  identificadas  y  los
estándares  establecidos  por  los  ensayos  KEYNOTE-048  y  KEYNOTE-689  valida
este  flujo  de  trabajo  bioinformático  como  un  puente  sólido  hacia  la  práctica

clínica traslacional en el carcinoma escamocelular de cabeza y cuello (29).

Desde  una  perspectiva  traslacional,  estos  hallazgos  adquieren  relevancia  al

evidenciar  que  una  proporción  considerable  de  los  candidatos  identificados
corresponde a fármacos ya aprobados, lo que abre la posibilidad de estrategias

de  reposicionamiento  con  aplicabilidad  clínica  en  el  corto  plazo  (7,8).  Este
aspecto  resulta  especialmente  pertinente  en  el  contexto  del  CECC,  donde  las
opciones  terapéuticas  siguen  siendo  limitadas  y  los  desenlaces  clínicos

continúan siendo subóptimos en estadios avanzados (1,3,6).

Más allá de la identificación de blancos individuales, este estudio propone una
visión  de  la  enfermedad  basada  en  redes  funcionales,  sugiriendo  que  la
intervención terapéutica efectiva podría depender de la modulación coordinada
de  múltiples  procesos  biológicos  interconectados  (13,17).  En  este  sentido,  el

reposicionamiento  farmacológico  guiado  por  proteómica  se  posiciona  como
una  herramienta  estratégica  para  el  desarrollo  de  terapias  combinadas

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

racionales  y  para  la  implementación  de  enfoques  de  medicina  de  precisión
(7,8,15).

Finalmente, aunque estos resultados requieren validación experimental y clínica
(9,16),  la  coherencia  biológica  y  la  plausibilidad  farmacológica  observadas
respaldan su relevancia. En conjunto, este trabajo aporta un marco conceptual

y metodológico sólido para la identificación de oportunidades terapéuticas en
CECC, con potencial impacto en la optimización de estrategias multimodales y

en la mejora de los resultados oncológicos en esta población (6,15,20).

XII. AGRADECIMIENTOS

Expresamos nuestro reconocimiento a los grupos de investigación cuyos datos
proteómicos  contribuyeron  al  desarrollo  de  este
trabajo.  Asimismo,
agradecemos el apoyo académico y formativo recibido durante el proceso de

elaboración del presente estudio, el cual fue fundamental para su ejecución.

Finalmente,  agradecemos  a  nuestros  colegas  y  mentores  por  sus  valiosos

aportes, discusiones y revisión crítica del manuscrito, que permitieron fortalecer
la calidad científica de este trabajo

XIV. DECLARACIÓN DE FINANCIACIÓN DEL PROYECTO

Declaramos que los investigadores no tenemos conflictos de interés

Declaramos  que  el  archivo  no  requiere  financiación  por  parte  de  ninguna
identidad

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

XV. REFERENCIAS BIBLIOGRÁFICAS

1.  Bray  F,  Laversanne  M,  Sung  H,  Ferlay  J,  Siegel  RL,  Soerjomataram  I,  et  al.
Global  cancer  statistics  2022:  GLOBOCAN  estimates  of  incidence  and

mortality  worldwide  for  36  cancers  in  185  countries.  CA  Cancer  J  Clin.
2024;74(3):229-63. doi:10.3322/caac.21834.

2.  Liu  Y,  Zhang  N,  Wen  Y,  Wen  J.  Head  and  neck  cancer:  pathogenesis  and

targeted therapy. MedComm. 2024;5(9):e702. doi:10.1002/mco2.702.

3.  Siegel RL, Miller KD, Wagle NS, Jemal A. Cancer statistics, 2023. CA Cancer J

Clin. 2023;73(1):17-48. doi:10.3322/caac.21763.

4.  Budach V, Tinhofer I. Novel prognostic clinical factors and biomarkers for
outcome prediction in head and neck cancer: a systematic review. Lancet
Oncol. 2019;20(6):e313-e326.

5.  Chow LQM, Longo DL. Head and neck cancer. N Engl J Med. 2020;382(1):60-

72.

6.  Mody MD, Rocco JW, Yom SS, Haddad RI, Saba NF. Head and neck cancer.

Lancet. 2021;398(10318):2289-99.

7.  Russell C, Rahman A, Mohammed AR. Application of genomics, proteomics
and metabolomics in drug discovery, development and clinic. Ther Deliv.
2013;4(3):395-413. doi:10.4155/tde.13.4.

8.  Gheybi E, Hosseinzadeh P, Tayebi-Khorrami V, Rostami M, Soukhtanloo M.
Proteomics in decoding cancer: A review. Clin Chim Acta. 2025;574:120302.

doi:10.1016/j.cca.2025.120302.

9.  Babu  N,  Mohan  S,  Nanjappa  V,  Chavan  S,  Advani  J,  Khan  AA,  et  al.
Identification  of  potential  biomarkers  of  head  and  neck  squamous  cell

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

carcinoma  using  iTRAQ-based  quantitative  proteomic  approach.  Data
Brief. 2018;19:1124-30. doi:10.1016/j.dib.2018.05.100.

10.  Sudha  R,  Kawachi  N,  Du  P,  Nieves  E,  Belbin  TJ,  Negassa  A,  et  al.  Global
proteomic  analysis  distinguishes  biologic  differences  in  head  and  neck

squamous carcinoma. Lab Invest. 2007;87(8):755-66.

11.  Marimuthu A, Chavan S, Sathe G, Sahasrabuddhe NA, Srikanth SM, Renuse
S,  et  al.  Identification  of  head  and  neck  squamous  cell  carcinoma

biomarker  candidates  through  proteomic  analysis  of  cancer  cell
secretome. Biochim Biophys Acta. 2013;1834(11):2308-16.

12.  Rivera  C,  Oliveira  AK,  Costa  RAP,  De  Rossi  T,  Paes  Leme  AF.  Prognostic
biomarkers  in  oral  squamous  cell  carcinoma:  A  systematic  review.  Oral

Oncol. 2017;72:38-47.

13.  Li  H,  Wheeler  S,  Park  Y,  Ju  Z,  Thomas  SM,  Fichera  M,  et  al.  Proteomic
characterization of head and neck cancer patient-derived xenografts. Mol

Cancer Res. 2015;14(3):278-86.

14.  Cívico-Ortega  JL,  González-Ruiz  I,  Ramos-García  P,  Cruz-Granados  D,
Samayoa-Descamps
and
clinicopathological significance of EGFR expression in oral squamous cell
Int  J  Mol  Sci.
carcinoma:  systematic  review  and  meta-analysis.

V,  González-Moles  MÁ.

Prognostic

2023;24(15):11888.

15.  Ruffin AT, Li H, Vujanovic L, Zandberg DP, Ferris RL, Bruno TC. Improving head
immunomodulation  of  the  tumour

and  neck  cancer  therapies  by
microenvironment. Nat Rev Cancer. 2023;23:173-88.

16.  Harris TM, Du P, Kawachi N, Belbin TJ, Wang Y, Schlecht NF, et al. Proteomic
analysis  of  oral  cavity  squamous  cell  carcinoma  specimens  identifies

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

patient  outcome-associated  proteins.  Arch
2015;139(4):494-507.

Pathol

Lab  Med.

17.  Huang  C,  Chen  L,  Savage  SR,  Vargas  Eguez  R,  Dou  Y,  Li  Y,  et  al.
Proteogenomic  insights  into  the  biology  and  treatment  of  HPV-negative

head and neck squamous cell carcinoma. Cancer Cell. 2021;39(3):361-379.

18.  Kaneko T, Zeng PYF, Liu X, Abdo R, Barrett JW, Zhang Q, et al. Proteome and
phosphoproteome  signatures  of  recurrence  for  HPV+  head  and  neck

squamous cell carcinoma. Commun Med (Lond). 2022;2:95.

19.  Sung H, Ferlay J, Siegel RL, Laversanne M, Soerjomataram I, Jemal A, et al.
Global  cancer  statistics  2020:  GLOBOCAN  estimates  of  incidence  and
mortality  worldwide  for  36  cancers  in  185  countries.  CA  Cancer  J  Clin.

2021;71(3):209-49. doi:10.3322/caac.21660.

20. Marur S, Forastiere AA. Head and neck squamous cell carcinoma: update
on epidemiology, diagnosis, and treatment. Mayo Clin Proc. 2016;91(3):386-

96. doi:10.1016/j.mayocp.2015.12.017.

21.  Zhang  Y,  Fakhry  C,  D’Souza  G.  Projected  association  of  human
papillomavirus  vaccination  with  oropharynx  cancer  incidence  in  the  US,
2022;8(3):426-33.
2020-2045.

Oncol.

JAMA

doi:10.1001/jamaoncol.2021.6553.

22. Haddad RI, Seiwert TY, Chow LQM, Gupta S, Weiss J, Gluck I, et al. Influence
of  tumor  mutational  burden,  inflammatory  gene  expression  profile,  and

PD-L1  expression  on  response  to  pembrolizumab  in  head  and  neck
Immunother  Cancer.  2022;10:e003026.
squamous  cell  carcinoma.  J

doi:10.1136/jitc-2021-003026.

23. Matoba T, Minohara K, Kawakita D, Takano G, Oguri K, Murashima A, et al.
Impact of tumor burden on survival in patients with recurrent or metastatic

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

head and neck cancer treated with immune checkpoint inhibitors. Sci Rep.
2022;12:14319. doi:10.1038/s41598-022-18611-z.

24. Vermorken  JB,  Mesia  R,  Rivera  F,  Remenar  E,  Kawecki  A,  Rottey  S,  et  al.
Platinum-based chemotherapy plus cetuximab in head and neck cancer.

N Engl J Med. 2008;359(11):1116-27. doi:10.1056/NEJMoa0802656.

25. Burtness B, Harrington KJ, Greil R, Soulières D, Tahara M, de Castro G Jr, et
al.  Pembrolizumab  alone  or  with  chemotherapy  versus  cetuximab  with

chemotherapy for recurrent or metastatic head and neck squamous cell
carcinoma  (KEYNOTE-048):  overall  survival  results  from  a  randomised,

open-label,
doi:10.1016/S0140-6736(19)32591-7.

phase

3

study.

Lancet.

2019;394(10212):1915-28.

26. Ferris RL, Blumenschein G Jr, Fayette J, Guigay J, Colevas AD, Licitra L, et al.
Nivolumab for recurrent squamous-cell carcinoma of the head and neck.
N Engl J Med. 2016;375(19):1856-67. doi:10.1056/NEJMoa1602252.

27. Cohen  EEW,  Soulières  D,  Le  Tourneau  C,  Dinis  J,  Licitra  L,  Ahn  MJ,  et  al.
for
Pembrolizumab  versus  methotrexate,  docetaxel,  or  cetuximab
recurrent  or  metastatic  head-and-neck  squamous  cell  carcinoma
(KEYNOTE-040):  a  randomised,  open-label,  phase  3  study.  Lancet.

2019;393(10167):156-67. doi:10.1016/S0140-6736(18)31999-8.

28. Uppaluri R, Campbell KM, Egloff AM, Zolkind P, Skidmore ZL, Nussenbaum B,
et  al.  Neoadjuvant  and  adjuvant  pembrolizumab  in  resectable  locally

advanced,  human  papillomavirus-unrelated  head  and  neck  cancer:  a
trial.  Clin  Cancer  Res.  2020;26(19):5140-52.
multicenter  phase

II

doi:10.1158/1078-0432.CCR-20-1695.

29. Wise-Draper TM, Old MO, Worden FP, O'Neill A, Cohen EEW, Dunlap N, et al.
trial  of  neoadjuvant
investigator-initiated
Phase
pembrolizumab  and  adjuvant  concurrent  radiation  and  pembrolizumab

II  multi-site

FUNDACIÓN UNIVERSITARIA DE CIENCIAS  DE LA

SALUD

VERSIÓN 01

FORMULACIÓN Y EJECUCIÓN DE PROYECTOS DE

INVESTIGACIÓN

CÓDIGO: F-PI-FEP-09

GUIA DE ELABORACION DE UN ARTÍCULO DE

INVESTIGACIÓN

FECHA 14-02-2018

with or without cisplatin in resectable stage III-IV head and neck squamous
cell carcinoma. J Clin Oncol. 2022;40(2):138-49. doi:10.1200/JCO.21.01256.

30. Topper MJ, Vaz M, Chiappinelli KB, DeStefano Shields CE, Niknafs N, Yen RC,
et al. Epigenetic therapy ties MYC depletion to reversing immune evasion
2017;171(6):1284-1300.e21.
and

treating

cancer.

lung

Cell.

doi:10.1016/j.cell.2017.10.022.

31.  Wrangle J, Wang W, Koch A, Easwaran H, Mohammad HP, Vendetti F, et al.
Alterations  of  immune  response  of  non-small  cell  lung  cancer  with
azacytidine. Oncotarget. 2013;4(11):2067-79. doi:10.18632/oncotarget.1542.

32. ClinicalTrials.gov  [Internet].  Bethesda  (MD):  National  Library  of  Medicine

(US). Available from: https://clinicaltrials.gov/

33. Skrott  Z,  Mistrik  M,  Andersen  KK,  Friis  S,  Majera  D,  Gursky  J,  et  al.  Alcohol-
abuse  drug  disulfiram  targets  cancer  via  p97  segregase  adaptor  NPL4.
Nature. 2017;552(7684):194-9. doi:10.1038/nature25016.

34. González-Moles  MÁ,  Scully  C,  Ruiz-Ávila  I,  Plaza-Campillo  JJ.  The  cancer
stem cell hypothesis applied to oral carcinoma. Oral Oncol. 2013;49(8):738-

46. doi:10.1016/j.oraloncology.2013.05.011.

35. Leemans CR, Snijders PJF, Brakenhoff RH. The molecular landscape of head
and neck cancer. Nat Rev Cancer. 2018;18(5):269-82. doi:10.1038/nrc.2018.11.

36. Johnson DE, Burtness B, Leemans CR, Lui VWY, Bauman JE, Grandis JR. Head
and  neck  squamous  cell  carcinoma.  Nat  Rev  Dis  Primers.  2020;6(1):92.

doi:10.1038/s41572-020-00224-3.

