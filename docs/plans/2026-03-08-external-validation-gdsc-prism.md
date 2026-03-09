# HNSCC Drug Repurposing — External Cross-Validation Plan (GDSC2 + PRISM)

> **Para Claude (nueva sesión):** Ejecutar este plan con el approach **Subagent-Driven**:
> 1. Invocar `superpowers:subagent-driven-development` ANTES de cualquier acción
> 2. Despachar un subagente por tarea
> 3. Revisar output del subagente (spec compliance + code quality) antes de continuar
> 4. Usar `TodoWrite` para trackear progreso
>
> **Directorio de trabajo:** `~/bioinfo/projects/hnscc_drug_repurposing`
> **Ambiente Python:** `omics-py` | **Ambiente R:** `omics-R`
> **Plan completo:** este archivo — leerlo completo antes de empezar
> **Orden de ejecución:** Task 1 → Task 2 → Task 3 → Task 4

**Goal:** Validate the Top 20 drug candidates from the repurposing pipeline against
publicly available drug sensitivity data (GDSC2 and PRISM), computing
HNSCC-preferential sensitivity z-scores as orthogonal in silico validation.
Additionally, investigate the Mavacamten/M11 muscle module as a potential
proteomics artifact.

**Context:**
The pipeline produces candidate drugs ranked by a composite score based on
proteomics DE, pathway enrichment, network centrality, clinical phase, and
literature evidence. A key gap is lack of external validation — we have not
checked whether our top candidates actually show preferential anti-cancer
activity in HNSCC cell lines. GDSC2 and PRISM provide IC50 and viability
data across hundreds of cancer cell lines including HNSCC lines.

**Key inputs:**
- `results/tables/10_top20_candidates.tsv` — Top 20 candidates with `drug_name_norm` column
- `results/tables/10_all_candidates_scored.tsv` — All scored candidates

**Key outputs (to generate):**
- `results/tables/validation/16_gdsc_sensitivity.tsv`
- `results/tables/validation/16_prism_sensitivity.tsv`
- `results/tables/validation/16_validation_summary.tsv`
- `results/figures/validation/16_gdsc_zscore_dotplot.pdf`
- `results/figures/validation/16_prism_zscore_dotplot.pdf`
- `results/figures/validation/16_validation_concordance.pdf`
- `results/tables/validation/16_muscle_module_assessment.tsv`

---

## Dependency Map

```
PHASE I:  Task 1 (GDSC2 validation script)
          Task 2 (PRISM validation script)   ← independientes, paralelos
               │
PHASE II: Task 3 (muscle module artifact assessment)  ← independiente
               │
PHASE III: Task 4 (update docs + RUNBOOK + METHODS)
```

Tasks 1 y 2 pueden despacharse en paralelo. Task 3 también es independiente.
Task 4 requiere que 1, 2, y 3 hayan terminado.

---

## PHASE I — External Validation Scripts

### Task 1: Script `16_gdsc_validation.py` (GDSC2)

**Ambiente:** `omics-py`
**Archivo a crear:** `scripts/16_gdsc_validation.py`

**Descripción:** Descargar datos de sensibilidad farmacológica de GDSC2,
mapear nuestros Top 20 candidatos a los nombres de fármacos de GDSC2, y
calcular z-scores de sensibilidad preferencial en HNSCC.

**Implementación detallada:**

**Paso 1: Descarga de datos GDSC2**

Descargar el archivo de fitted dose-response de GDSC2. Usar requests con
fallback a múltiples URLs. Guardar en `data/raw/gdsc2_dose_response.csv`
(solo si no existe — cache local para evitar re-descargar).

URLs a intentar en orden:
```
https://cog.sanger.ac.uk/cancerrxgene/GDSC_data_BULK/drug_screening_matrices/GDSC2_fitted_dose_response_27Oct23.xlsx
```

Si descarga falla, el script debe continuar con un aviso y producir outputs
vacíos (no romper el pipeline). Logging explícito de si los datos vienen
de cache o de descarga fresca.

Columnas relevantes del archivo GDSC2:
- `DRUG_NAME`: nombre del fármaco
- `DRUG_ID`: ID numérico GDSC
- `COSMIC_ID`: ID de línea celular
- `CELL_LINE_NAME`: nombre de la línea
- `TCGA_DESC`: tipo de cáncer (usar "HNSC" para HNSCC)
- `LN_IC50`: log natural del IC50 (en µM) — menor = más sensible
- `AUC`: área bajo la curva — menor = más sensible

**Paso 2: Identificar líneas celulares HNSCC**

Filtrar por `TCGA_DESC == "HNSC"`. Reportar en log cuántas líneas HNSCC
hay en el dataset y sus nombres.

**Paso 3: Normalización de nombres de fármacos**

Función `normalize_drug_name(name)`:
```python
import re

SALT_SUFFIXES = [
    "hydrochloride", "sodium", "sulfate", "citrate", "mesylate",
    "acetate", "phosphate", "tartrate", "maleate", "fumarate",
    "tosylate", "ditosylate", "monohydrate", "dihydrate", "hemihydrate",
]

def normalize_drug_name(name):
    if not isinstance(name, str):
        return ""
    name = name.lower().strip()
    # Eliminar sufijos de sal
    for suffix in SALT_SUFFIXES:
        name = re.sub(r'\s+' + suffix + r'(\s+|$)', ' ', name)
    # Eliminar espacios extra y puntuación no alfanumérica
    name = re.sub(r'[^a-z0-9\s-]', '', name)
    name = re.sub(r'\s+', ' ', name).strip()
    return name
```

Aplicar a `DRUG_NAME` de GDSC2 y a `drug_name_norm` de nuestro Top 20
(que ya está normalizado). Hacer join por nombre normalizado. Para casos
sin match exacto, intentar matching parcial (substring de longitud ≥ 6
que esté contenido en ambos nombres).

Reportar en log: N candidatos del Top 20 encontrados en GDSC2, N no encontrados.

**Paso 4: Calcular z-scores de sensibilidad HNSCC**

Para cada fármaco encontrado en GDSC2:

```python
# Para cada drug encontrado:
ic50_hnscc = df[(df['drug_norm'] == drug) & (df['TCGA_DESC'] == 'HNSC')]['LN_IC50']
ic50_others = df[(df['drug_norm'] == drug) & (df['TCGA_DESC'] != 'HNSC')]['LN_IC50']

# z-score: negativo = HNSCC más sensible que promedio (BUENO para repurposing)
z_hnscc = (ic50_hnscc.mean() - ic50_others.mean()) / ic50_others.std()

# También calcular: ¿es HNSCC significativamente más sensible? (Mann-Whitney U)
from scipy import stats
stat, pval = stats.mannwhitneyu(ic50_hnscc, ic50_others, alternative='less')
```

Calcular también para AUC con la misma lógica.

**Paso 5: Exportar**

`results/tables/validation/16_gdsc_sensitivity.tsv` con columnas:
```
drug_name_norm | gdsc_drug_name | n_hnscc_lines | n_other_lines |
mean_ln_ic50_hnscc | mean_ln_ic50_others | z_hnscc_ic50 |
mean_auc_hnscc | mean_auc_others | z_hnscc_auc |
pval_mannwhitney | fdr_bh | validated_hnscc
```

`validated_hnscc = True` si z_hnscc_ic50 < -0.5 Y pval < 0.1.

**Paso 6: Figura**

Dotplot horizontal: eje X = z_hnscc_ic50, eje Y = drug_name_norm ordenado
por z-score. Colorear por: verde oscuro (z < -1, validado), verde claro
(-1 ≤ z < -0.5), gris (z ≥ -0.5), rojo (z > 0.5, resistente).
Línea vertical en z=0. Puntos con error bars (IC95 bootstrap si N≥3).

Guardar como `results/figures/validation/16_gdsc_zscore_dotplot.pdf`.

**Logging requerido:**
```
GDSC2 validation date: YYYY-MM-DD
Datos: [cache/fresh] desde [URL/path]
Líneas HNSCC: N (ej: FaDu, Cal-27, SCC-25, ...)
Candidatos encontrados en GDSC2: N/20
Candidatos validados (z<-0.5, p<0.1): N
Top validados: DRUG1 (z=X), DRUG2 (z=Y), ...
```

---

### Task 2: Script `16_prism_validation.py` (PRISM)

**Ambiente:** `omics-py`
**Archivo a crear:** `scripts/16_prism_validation.py`

**Descripción:** Descargar datos del PRISM Repurposing Secondary Screen de
DepMap y calcular sensibilidad preferencial HNSCC para nuestros candidatos.

**Paso 1: Descarga de datos PRISM**

El PRISM secondary screen está disponible en DepMap. Intentar descargar de:

```
https://depmap.org/portal/api/download/file?releasename=PRISM+Repurposing+19Q4&filename=secondary-screen-dose-response-curve-parameters.csv
```

Si falla, intentar la URL del archivo de viabilidad media:
```
https://depmap.org/portal/api/download/file?releasename=PRISM+Repurposing+19Q4&filename=secondary-screen-replicate-collapsed-logfold-change.csv
```

Cache en `data/raw/prism_secondary_screen.csv`. Si descarga falla, producir
outputs vacíos con aviso.

Columnas relevantes (dose-response parameters):
- `compound_name`: nombre del fármaco
- `broad_id`: ID Broad
- `depmap_id`: ID de línea celular (ej: ACH-000001)
- `auc`: área bajo la curva de viabilidad (menor = más letal)
- `ec50`: concentración de efecto 50%

Para el archivo de log-fold-change (alternativo):
- Columnas = compound broad_ids
- Filas = cell line depmap_ids
- Valores = log2(viabilidad) — más negativo = más letal

**Paso 2: Metadata de líneas celulares**

Descargar metadata para identificar líneas HNSCC:
```
https://depmap.org/portal/api/download/file?releasename=DepMap+Public+24Q4&filename=Model.csv
```

Filtrar por `OncotreeLineage == "Head and Neck"` o `OncotreePrimaryDisease`
que contenga "Head and Neck Squamous". Cache en `data/raw/depmap_model_metadata.csv`.

**Paso 3: Normalización y matching**

Usar la misma función `normalize_drug_name()` del Task 1. Aplicar a
`compound_name` de PRISM y a `drug_name_norm` de Top 20. Join por nombre
normalizado.

**Paso 4: Calcular z-scores**

Para cada fármaco en PRISM encontrado en nuestro Top 20:
```python
auc_hnscc = df_drug[df_drug['is_hnscc']]['auc']
auc_others = df_drug[~df_drug['is_hnscc']]['auc']

# AUC en PRISM: menor = más letal = MEJOR para repurposing
z_hnscc = (auc_hnscc.mean() - auc_others.mean()) / auc_others.std()
# Negativo = HNSCC más sensible que promedio
```

Mann-Whitney U test `alternative='less'` (HNSCC AUC < others AUC).

**Paso 5: Exportar**

`results/tables/validation/16_prism_sensitivity.tsv` con columnas:
```
drug_name_norm | prism_compound_name | n_hnscc_lines | n_other_lines |
mean_auc_hnscc | mean_auc_others | z_hnscc_auc |
pval_mannwhitney | fdr_bh | validated_hnscc
```

`validated_hnscc = True` si z < -0.5 Y pval < 0.1.

**Paso 6: Figura de concordancia GDSC2 vs PRISM**

Scatter plot: eje X = z_hnscc GDSC2, eje Y = z_hnscc PRISM.
Solo candidatos presentes en ambas bases. Anotar con nombre del fármaco.
Correlación de Spearman en el plot. Cuadrante Q3 (x<0, y<0) = validado
en ambas bases (colorear en verde oscuro).

Guardar como `results/figures/validation/16_validation_concordance.pdf`.

**Paso 7: Tabla resumen integrada**

`results/tables/validation/16_validation_summary.tsv` con columnas:
```
drug_name_norm | rank_pipeline | composite_score |
z_gdsc | gdsc_validated | z_prism | prism_validated |
validated_both | validation_status
```

`validation_status`: "CONFIRMED" (ambas), "PARTIAL" (una), "NOT_FOUND" (ninguna DB),
"DISCORDANT" (z positivo en alguna).

---

## PHASE II — Artifact Assessment

### Task 3: Muscle Module Artifact Assessment

**Ambiente:** `omics-R`
**Archivo a crear:** `scripts/16b_muscle_module_assessment.R`

**Contexto:** El módulo M11 del análisis Louvain fue anotado como
"muscle_contraction". Mavacamten (#3 en el ranking) es un inhibidor de
miosina cardíaca que conecta con este módulo. En proteómica tumoral, las
muestras quirúrgicas de HNSCC frecuentemente incluyen tejido muscular
peritumoral, lo que puede contaminar la señal proteómica con proteínas
sarcoméricas. Si M11 es un artefacto de contaminación, Mavacamten podría
ser un falso positivo.

**Implementación:**

```r
# Input: results/tables/network/09_modules.tsv
#        results/tables/de_limma/01_TVsS_all_proteins.tsv
```

**Paso 1: Extraer proteínas del módulo M11**

```r
modules <- read_tsv("results/tables/network/09_modules.tsv")
m11 <- modules %>% filter(grepl("muscle|M11", module_name))
```

**Paso 2: Caracterizar el módulo**

Para cada proteína en M11:
- Obtener su logFC y adj.P.Val del DE analysis (join con `01_TVsS_all_proteins.tsv`)
- Consultar si es una proteína sarcomérica/muscular específica vs. ubiqua
- Usar una lista de marcadores musculares conocidos:
  ```r
  MUSCLE_MARKERS <- c("MYH7", "MYH6", "MYH3", "MYH1", "MYH2", "MYH4",
                       "MYL1", "MYL2", "MYL3", "MYL4", "TNNT2", "TNNI3",
                       "TPM1", "TPM2", "ACTA1", "TTN", "MYBPC3", "MYBPC2")
  ```
- Calcular: % de proteínas M11 que son marcadores musculares

**Paso 3: Análisis de missingness**

¿Las proteínas de M11 tienen más datos faltantes que el promedio?
(Cargar `01_missingness_analysis.tsv`). Proteínas detectadas en pocas
muestras tienen más probabilidad de ser contaminación que expresión real.

**Paso 4: Distribución de logFC**

¿Las proteínas de M11 son upregulated o downregulated consistentemente
en tumoral vs normal? Si son consistentemente downregulated (como
esperaríamos si normal tissue tiene más músculo que tumor), refuerza
la hipótesis de contaminación. Si están upregulated, podría ser señal real
(aunque improbable para músculo cardíaco en HNSCC).

**Paso 5: Exportar y concluir**

`results/tables/validation/16_muscle_module_assessment.tsv` con columnas:
```
gene_symbol | module_id | is_muscle_marker | logFC | adj.P.Val |
n_samples_detected | direction | contamination_flag
```

Generar un reporte de texto en el log:
```
=== MUSCLE MODULE ARTIFACT ASSESSMENT ===
Module: M11 (muscle_contraction)
Proteins in module: N
Muscle marker proteins: N (X%)
Mean logFC in M11: X (positive = upregulated in tumor)
Mean missingness rate M11: X% vs pipeline mean: Y%
CONCLUSION: [LIKELY ARTIFACT | LIKELY REAL | INCONCLUSIVE]
```

`LIKELY ARTIFACT` si: >50% son muscle markers Y logFC medio es negativo
(menos en tumor) Y missingness > promedio pipeline.

`LIKELY REAL` si: <20% muscle markers O logFC es positivo.

`INCONCLUSIVE` en otro caso.

Guardar la conclusión también en
`results/tables/validation/16_muscle_module_assessment.tsv`.

---

## PHASE III — Documentation Update

### Task 4: Update docs

**Archivos a modificar:**
- `docs/RUNBOOK.md`: Agregar Fase 8b (validación externa) entre Fase 7 y Fase 8 (figuras)
- `docs/METHODS.md`: Agregar sección "External Cross-Validation" después de Sensitivity Analysis
- `docs/OUTPUTS.md`: Si existe sección de tablas, añadir los nuevos outputs

**RUNBOOK.md — añadir:**
```markdown
### Fase 8b: Validación externa GDSC2 + PRISM (Python)

```bash
conda activate omics-py

# 16 - Cross-validación contra sensibilidad farmacológica en cell lines HNSCC
python scripts/16_gdsc_validation.py
# Output: results/tables/validation/16_gdsc_sensitivity.tsv
#         results/figures/validation/16_gdsc_zscore_dotplot.pdf

python scripts/16_prism_validation.py
# Output: results/tables/validation/16_prism_sensitivity.tsv
#         results/tables/validation/16_validation_summary.tsv
#         results/figures/validation/16_validation_concordance.pdf
```

```bash
conda activate omics-R

# 16b - Evaluación de artefacto de contaminación muscular en módulo M11
Rscript scripts/16b_muscle_module_assessment.R
# Output: results/tables/validation/16_muscle_module_assessment.tsv
```
```

**METHODS.md — agregar sección:**
```markdown
## Phase 6b: External Cross-Validation (Scripts 16, 16b)

To provide orthogonal validation of the top candidates, we queried two
publicly available drug sensitivity databases:

**GDSC2** (Genomics of Drug Sensitivity in Cancer, v2): Fitted dose-response
data (LN_IC50 and AUC) from the Sanger Institute for HNSCC cell lines
(TCGA_DESC = "HNSC"). For each candidate found in GDSC2, a HNSCC-preferential
sensitivity z-score was computed: z = (mean_IC50_HNSCC − mean_IC50_others) /
sd_IC50_others. Negative z indicates preferential HNSCC sensitivity.
Statistical significance tested with Mann-Whitney U (one-sided, alternative='less').

**PRISM Repurposing Secondary Screen** (Broad Institute, 19Q4): Log-fold-change
viability data for 1,448 FDA-approved compounds across 578 cancer cell lines.
HNSCC lines identified via DepMap Model metadata (OncotreeLineage = "Head and Neck").
Same z-score framework applied using AUC values.

Drug name matching between our pipeline and each database used normalized
names (lowercase, salt suffixes removed) with substring fallback for
partial matches.

Candidates validated in both databases (z < −0.5, p < 0.1) were designated
"CONFIRMED". Candidates not found in either database were designated "NOT_FOUND"
and their validation status remains unknown.

**Muscle module artifact assessment (Script 16b):** The Louvain module M11
(muscle_contraction, containing MYH7/MYL2/MYL3) was assessed for contamination
with peritumoral muscle tissue — a known artifact in surgical HNSCC proteomics.
Criteria evaluated: fraction of muscle-specific markers, mean logFC direction,
and per-protein missingness rates.
```

---

## Re-ejecución tras implementación

Ejecutar en este orden:

```bash
# Paralelo: 16 y 16b son independientes
conda activate omics-py
python scripts/16_gdsc_validation.py 2>&1 | tee logs/16_gdsc_$(date +%Y%m%d_%H%M%S).log
python scripts/16_prism_validation.py 2>&1 | tee logs/16_prism_$(date +%Y%m%d_%H%M%S).log

conda activate omics-R
Rscript scripts/16b_muscle_module_assessment.R 2>&1 | tee logs/16b_muscle_$(date +%Y%m%d_%H%M%S).log
```

Revisar:
- `results/tables/validation/16_validation_summary.tsv` — status por candidato
- Log del 16b — conclusión sobre M11 (LIKELY ARTIFACT / LIKELY REAL / INCONCLUSIVE)
- Si M11 = LIKELY ARTIFACT: considerar excluir Mavacamten del Top 20 y documentarlo

---

## Notas para el implementador

1. **Si GDSC2 o PRISM no descargan** (timeout, URL caída): el script debe
   producir un TSV vacío con las columnas correctas y loggear "DATA_NOT_AVAILABLE".
   No romper el pipeline.

2. **Drug name matching**: Es la parte más propensa a fallar. Si un candidato
   no matchea, intentar:
   - Coincidencia exacta normalizada
   - Substring match (len ≥ 6)
   - Reportar todos los no-matches en el log para revisión manual

3. **Mavacamten en GDSC/PRISM**: Es un fármaco muy nuevo (aprobado 2022).
   Es probable que NO esté en estas bases de datos (que son de 2019-2023).
   El script debe manejar esto graciosamente ("NOT_FOUND").

4. **Interpretación de resultados**: El análisis es exploratorio. Cell lines
   ≠ primary tumors. Un z-score negativo valida la hipótesis pero no garantiza
   eficacia clínica. Un z positivo no invalida el candidato necesariamente.
   Documentar esta limitación en el log.

5. **Outputs vacíos no son error**: Si ningún candidato matchea en GDSC2,
   producir igualmente la tabla con columnas correctas y N=0 rows, no romper.
