# HNSCC Scoring Fixes — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Corregir tres problemas encontrados en la revisión crítica del pipeline que afectan la coherencia y defensibilidad de la lista de candidatos antes de la validación bibliográfica manual.

**Architecture:** Tres cambios independientes en scripts R del pipeline. El más importante es la corrección del π-statistic en script 10, que pasa de magnitud absoluta a dirección explícita (upregulated > downregulated). Los otros dos son correcciones de clasificación en script 08 y de exclusión en config. Cada cambio termina en re-ejecución del script afectado y verificación del output. Se termina con una re-ejecución completa de scripts 10 → 15 para que el panel LOD-estable refleje los cambios.

**Tech Stack:** R (dplyr, yaml), `conda activate omics-R`, proyecto en `~/bioinfo/projects/hnscc_drug_repurposing/`

---

## Contexto del problema

### Problema central: el π-statistic usa abs() — scoring "desarticulado"

El pipeline tiene **dos hipótesis de repurposing** que en teoría deben complementarse:

1. **Hipótesis de inversión de firma** (L2S2/CMap): busca fármacos que down-regulen lo que el tumor tiene up-regulado, y vice versa. Solo reversores (score < 0) son considerados. Peso en composite: 0.13.

2. **Hipótesis de magnitud DE** (π-statistic): asume que proteínas con mayor cambio diferencial en tumor vs. normal son los targets más relevantes. Actualmente usa `abs(π)`. Peso: 0.325.

**El problema:** `abs(π)` hace que inhibir un gen *fuertemente downregulado* (como MYH7, logFC ≈ −5) dé el mismo puntaje que inhibir uno *fuertemente upregulado* (como EGFR). Esto es inconsistente con la hipótesis de "inverter la firma": idealmente, un fármaco que target genes *upregulados* en el tumor (y que la dimensión L2S2 también confirma como reversor) debería puntuar más alto.

Con la corrección propuesta (π signed, min-max normalizado):
- Un fármaco apuntando a genes **up** → s_pi_stat cercano a 1.0
- Un fármaco apuntando a genes **down** (con mecanismo válido, como metformin vs OXPHOS) → s_pi_stat ~0.3–0.5 (sigue en la lista, pero por debajo)
- Falsos positivos por contaminación tisular (miosinas cardíacas) → s_pi_stat cercano a 0.0

**Efecto esperado:** EGFR inhibitors suben en ranking de s_pi_stat. Metformin mantiene posición alta por pathway y network. Fármacos cardíacos ya no LOD-stables igual — se confirma su caída.

---

## Task 1: Corregir directionality del π-statistic en script 10

**Files:**
- Modify: `scripts/10_prioritization_scoring.R:200-207`
- Verify: `results/tables/10_top20_candidates.tsv` (columna s_pi_stat)
- Log: `logs/10_prioritization_scoring_*.log`

### Step 1: Leer el bloque actual (verificar líneas antes de editar)

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
sed -n '195,210p' scripts/10_prioritization_scoring.R
```

Expected: ver `max_pi <- max(abs(...))` y `abs(gene_de$pi_stat[...])` en `score_pistat_fn`.

### Step 2: Reemplazar el bloque de normalización (líneas 200–207)

Reemplazar exactamente este bloque:

```r
max_pi <- max(abs(gene_de$pi_stat), na.rm = TRUE)

score_pistat_fn <- function(gene_str) {
  genes   <- get_target_genes(gene_str)
  pi_vals <- abs(gene_de$pi_stat[gene_de$symbol_org %in% genes])
  if (length(pi_vals) == 0) return(0)
  mean(pi_vals, na.rm = TRUE) / max_pi
}
```

Por este bloque (signed, min-max normalizado):

```r
# Normalización direccional: 0 = target más downregulado del dataset,
# 1 = target más upregulado. Mantiene plausibilidad para mecanismos
# metabólicos (metformin/OXPHOS) pero prioriza inhibición de genes UP.
pi_min <- min(gene_de$pi_stat, na.rm = TRUE)
pi_max <- max(gene_de$pi_stat, na.rm = TRUE)
pi_range <- pi_max - pi_min

score_pistat_fn <- function(gene_str) {
  genes   <- get_target_genes(gene_str)
  pi_vals <- gene_de$pi_stat[gene_de$symbol_org %in% genes]  # signed
  if (length(pi_vals) == 0) return(0)
  mean_pi <- mean(pi_vals, na.rm = TRUE)
  (mean_pi - pi_min) / pi_range  # 0–1, mayor = target más upregulado
}
```

También actualizar el mensaje de log en línea ~86 para reflejar el cambio:
- Cambiar `"pi_stat=%.3f"` a `"pi_stat=%.3f (directional, signed min-max)"` en el sprintf.

### Step 3: Verificar la lógica manualmente antes de correr

Valores de referencia del dataset actual (del log de script 03):
- `pi_stat` range: −20.91 a +12.43
- Con la normalización: `pi_range = 12.43 − (−20.91) = 33.34`
- EGFR (pi ≈ +5): `(5 − (−20.91)) / 33.34 ≈ 0.78` → correcto (alto)
- NDUF target (pi ≈ −8): `(−8 − (−20.91)) / 33.34 ≈ 0.39` → correcto (medio)
- MYH7 (pi ≈ −14): `(−14 − (−20.91)) / 33.34 ≈ 0.21` → correcto (bajo)

### Step 4: Ejecutar script 10

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/10_prioritization_scoring.R 2>&1 | tail -40
```

### Step 5: Verificar cambio en ranking de s_pi_stat

```bash
# Comparar s_pi_stat de EGFR drugs vs antes (era ~0.23, ahora debe ser ~0.70-0.80)
awk -F'\t' 'NR>1 {print $1, $16}' results/tables/10_top20_candidates.tsv | \
  grep -i "gefitinib\|lapatinib\|afatinib\|erlotinib\|metformin\|mavacamten" | column -t
```

**Criterio de aceptación:**
- `s_pi_stat` de GEFITINIB/LAPATINIB/AFATINIB > 0.70
- `s_pi_stat` de METFORMIN entre 0.30 y 0.50 (baja un poco pero mantiene alta posición por otras dimensiones)
- `s_pi_stat` de MAVACAMTEN < 0.30

### Step 6: Commit

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
git add scripts/10_prioritization_scoring.R
git commit -m "fix: pi_stat directional (signed min-max) — up-regulated targets score higher

Previously abs(pi_stat) gave equal weight to strongly downregulated targets
(e.g. cardiac miosins MYH7/MYL2/MYL3) as to upregulated oncogenes. This
inflated scores of tissue-composition artifacts.

New normalization: 0 = most downregulated target in dataset, 1 = most
upregulated. EGFR inhibitors now score ~0.78 vs 0.23 before. Metabolic
repurposing candidates (metformin/OXPHOS) score ~0.39, still meaningful.

Ref: critical review session 2026-04-04"
```

---

## Task 2: Reclasificar DIGOXIN (y DIGITOXIN) de Clase B → C en script 08

**Contexto:** DIGOXIN tiene `has_cancer_indication=TRUE` en OpenTargets porque hay clinical trials en Kaposi's/NSCLC/pancreático. Sin embargo, la FDA **no aprueba** digoxin para ningún cáncer. La clasificación correcta es Clase C (aprobado para indicación no-oncológica).

**Impacto en el panel:** DIGOXIN está en el panel LOD-estable (n=26). Cambiar su clase de B→C no modifica su score compuesto (solo añadiría +0.01 de bonus, irrelevante). El cambio afecta únicamente la etiqueta en la tabla del manuscrito.

**Files:**
- Modify: `scripts/08_integrate_drug_targets.R` — agregar bloque override después de `4c`
- Verify: `results/tables/drug_targets/08_drug_summary_per_drug.tsv`

### Step 1: Leer contexto del bloque de overrides actuales

```bash
sed -n '398,440p' scripts/08_integrate_drug_targets.R
```

Expected: ver el bloque `CANCER_APPROVED_OVERRIDE` (neratinib, lazertinib) terminando en `~línea 437`.

### Step 2: Agregar bloque NOT_CANCER_APPROVED_OVERRIDE después del bloque 4c

Insertar después de la línea que cierra el `if (n_b_overridden > 0) { ... }` del bloque 4c (aproximadamente después de la línea 437):

```r
# ---------------------------------------------------------------------------
# 4d. Override C: fármacos aprobados para indicaciones NO oncológicas cuya
#     clasificación como B es un artefacto de Open Targets (CT ≠ FDA approval).
# ---------------------------------------------------------------------------
NOT_CANCER_APPROVED_OVERRIDE <- c(
  "DIGOXIN",   # FDA: insuficiencia cardíaca/arritmia. OT muestra CTs en K. sarcoma/NSCLC, no aprobación oncológica.
  "DIGITOXIN"  # Igual — glucósido cardíaco sin aprobación oncológica FDA.
)

n_c_overridden <- sum(
  toupper(drug_summary$drug_name_norm) %in% NOT_CANCER_APPROVED_OVERRIDE &
    drug_summary$drug_class == "B"
)
if (n_c_overridden > 0) {
  cat(sprintf("\n  Corrigiendo %d drogas de Clase B → C (clinical trial ≠ FDA cancer approval):\n",
              n_c_overridden))
  drug_summary <- drug_summary %>%
    mutate(
      has_cancer_indication = ifelse(
        toupper(drug_name_norm) %in% NOT_CANCER_APPROVED_OVERRIDE,
        FALSE, has_cancer_indication
      ),
      drug_class = ifelse(
        toupper(drug_name_norm) %in% NOT_CANCER_APPROVED_OVERRIDE & drug_class == "B",
        "C", drug_class
      ),
      drug_class_label = case_when(
        drug_class == "A" ~ "A: Approved + HNSCC evidence",
        drug_class == "B" ~ "B: Approved other cancer",
        drug_class == "C" ~ "C: Approved non-oncology",
        drug_class == "D" ~ "D: Not approved / Experimental"
      )
    )
  cat(paste(
    drug_summary$drug_name_norm[toupper(drug_summary$drug_name_norm) %in%
                                   NOT_CANCER_APPROVED_OVERRIDE],
    collapse = ", "
  ), "\n")
}
```

### Step 3: Ejecutar script 08

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/08_integrate_drug_targets.R 2>&1 | grep -E "Clase|DIGOXIN|DIGITOXIN|override|Corrigiendo"
```

**Criterio de aceptación:** log muestra `Corrigiendo 2 drogas de Clase B → C: DIGOXIN, DIGITOXIN`

### Step 4: Verificar clasificación

```bash
grep -i "^digoxin\|^digitoxin" \
  results/tables/drug_targets/08_drug_summary_per_drug.tsv | \
  awk -F'\t' '{print $1, $15, $16}'
```

Expected: `DIGOXIN   C   C: Approved non-oncology`

### Step 5: Commit

```bash
git add scripts/08_integrate_drug_targets.R
git commit -m "fix: reclassify DIGOXIN/DIGITOXIN Class B→C (OT clinical trials ≠ FDA cancer approval)

OpenTargets shows cardiac glycosides in cancer clinical trials (Kaposi,
NSCLC, pancreatic) which triggered has_cancer_indication=TRUE and Class B.
FDA has no approved oncology indication for these drugs. Correct class is C.

Does not change composite score (no scoring dimension uses drug_class directly).
Affects manuscript table label only."
```

---

## Task 3: Excluir fármacos de miosina cardíaca del análisis

**Contexto:** MAVACAMTEN, OMECAMTIV MECARBIL, DANICAMTIV y RELDESEMTIV/TIRASEMTIV aparecen en el top 35 porque sus targets (MYH7, MYL2, MYL3 — miosinas cardíacas/esqueléticas) están fuertemente downregulados. Esta señal es un artefacto de composición tisular: el tejido "normal" adyacente en HNSCC incluye fibras musculares esqueléticas que no están presentes en el tumor. Ninguno de estos fármacos tiene mecanismo anticancerígeno conocido.

> **Nota:** Estos fármacos ya NO son LOD-estables y no aparecen en el panel final de 26. La exclusión es profiláctica: elimina ruido del top 35 y del análisis de sensibilidad.

**Files:**
- Modify: `config/analysis_params.yaml` — agregar entradas a `exclusions.drugs`

### Step 1: Ver estado actual de la lista de exclusiones

```bash
grep -A30 "^exclusions:" config/analysis_params.yaml
```

### Step 2: Agregar entradas de miosina cardíaca a `exclusions.drugs`

Agregar después de la última entrada de `exclusions.drugs` (actualmente termina con `"OCRIPLASMIN"` y `"TECHNETIUM..."`):

```yaml
    - "MAVACAMTEN"            # Inhibidor miosina cardíaca MYH7 — target muscular (artefacto tisular HNSCC)
    - "OMECAMTIV MECARBIL"    # Activador miosina cardíaca MYH7 — target muscular (artefacto tisular HNSCC)
    - "DANICAMTIV"            # Activador miosina esquelética MYH7 — target muscular (artefacto tisular HNSCC)
    - "RELDESEMTIV"           # Activador troponina esquelética TNNT3 — target muscular (artefacto tisular HNSCC)
    - "TIRASEMTIV"            # Activador troponina esquelética TNNT3 — target muscular (artefacto tisular HNSCC)
    - "LEVOSIMENDAN"          # Sensibilizador calcio TNNI3/TNC — target muscular (artefacto tisular HNSCC)
```

### Step 3: Verificar que el YAML es válido

```bash
conda activate omics-R
Rscript -e "yaml::read_yaml('config/analysis_params.yaml'); cat('YAML OK\n')"
```

Expected: `YAML OK`

### Step 4: Commit (antes de re-ejecutar scripts)

```bash
git add config/analysis_params.yaml
git commit -m "config: exclude cardiac myosin drugs from candidates

MAVACAMTEN, OMECAMTIV MECARBIL, DANICAMTIV, RELDESEMTIV, TIRASEMTIV,
LEVOSIMENDAN target skeletal/cardiac myosin proteins (MYH7, MYL2, MYL3,
TNNI3, TNNT3) that are downregulated in HNSCC vs. normal due to tissue
composition (adjacent normal contains muscle fibers, tumor does not).

These drugs are not LOD-stable and have no anticancer mechanism.
Exclusion removes tissue-composition artifacts from top-35 and sensitivity
analysis without affecting the 26-drug LOD-stable panel."
```

---

## Task 4: Re-ejecución completa del pipeline desde script 08

Después de los tres cambios anteriores, es necesario re-ejecutar scripts 08 → 10 → 15 en orden para que los resultados finales reflejen todas las correcciones.

**Files:**
- Outputs actualizados: `results/tables/drug_targets/08_*.tsv`, `results/tables/10_*.tsv`, `results/tables/15_lod_stability.tsv`

### Step 1: Re-ejecutar script 08

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing
Rscript scripts/08_integrate_drug_targets.R 2>&1 | tail -20
```

Verificar en output:
- `Corrigiendo 2 drogas de Clase B → C: DIGOXIN, DIGITOXIN`
- `Candidatos multi-fuente (>=2): ~398` (6 menos que antes por las exclusiones de miosina)

### Step 2: Re-ejecutar script 09 (red STRING)

Script 09 no usa la clasificación A/B/C/D, pero sí cruza candidatos con hubs. Re-ejecutar para que el cruce de hubs druggables use la tabla actualizada.

```bash
Rscript scripts/09_string_network.R 2>&1 | tail -10
```

### Step 3: Re-ejecutar script 10 (scoring)

```bash
Rscript scripts/10_prioritization_scoring.R 2>&1 | tail -30
```

Verificar en el desglose de scores:
- GEFITINIB: `s_pi_stat > 0.70`
- METFORMIN: `s_pi_stat` entre 0.30 y 0.55
- MAVACAMTEN/OMECAMTIV/DANICAMTIV: **no deben aparecer** (excluidos por config)

### Step 4: Re-ejecutar script 15 (sensitivity + LOD stability)

```bash
Rscript scripts/15_sensitivity_analysis.R 2>&1 | tail -20
```

Verificar:
- El número de LOD-stable candidates debe ser similar al actual (26 ± 2)
- DIGOXIN debe seguir en el panel pero ahora como Clase C
- Los fármacos cardíacos de miosina **no deben aparecer** en ningún ranking

```bash
awk -F'\t' 'NR==1 || $NF=="TRUE"' results/tables/15_lod_stability.tsv | \
  cut -f1-3 | column -t
```

### Step 5: Commit final

```bash
git add results/tables/
git commit -m "results: re-run 08→09→10→15 after scoring and classification fixes

- pi_stat now directional (signed, min-max normalized)
- DIGOXIN/DIGITOXIN reclassified to Class C
- Cardiac myosin drugs excluded from candidates
LOD-stable panel updated accordingly."
```

---

## Task 5: Actualizar METHODS.md con las decisiones metodológicas

**Files:**
- Modify: `docs/METHODS.md` — sección de Análisis estadístico y Scoring

### Párrafos a agregar/modificar:

**En sección de DE (script 01):**
> "El diseño incluyó bloqueo intra-paciente mediante `duplicateCorrelation` (limma). La correlación de consenso intra-paciente fue de 0.016, lo que indica que el efecto del bloqueo fue mínimo en este contexto de proteómica DIA. Esto es consistente con la alta heterogeneidad biológica inter-tumoral y puede reflejar la normalización por mediana aplicada en los datos MaxQuant, que reduce parcialmente la variabilidad inter-individuo."

**En sección de Scoring (script 10):**
> "El π-statistic se calculó como π = signo(logFC) × |logFC| × −log₁₀(adj.P.Val + ε) para cada proteína significativa. Para el score de prioritización, se usó la media del π-statistic *con signo* de los genes diana del fármaco, normalizada al rango global del dataset (min-max). Esta normalización asigna puntuaciones más altas a fármacos que target proteínas *upreguladas* en el tumor y más bajas a los que target proteínas *downreguladas*. La hipótesis complementaria — que el fármaco *revierta* la firma tumoral — es capturada de forma independiente por la dimensión de conectividad L2S2, que solo considera fármacos cuya firma transcriptómica en líneas celulares es opuesta a la firma proteómica tumoral."

**En sección de Limitaciones:**
> "La dimensión de relevancia de pathways (score_pathway) mostró saturación: la mayoría de los candidatos alcanzaron el valor máximo (≈1.0) porque casi todos sus genes diana pertenecen a las mismas rutas enriquecidas dominantes (OXPHOS, señalización EGFR). Esta dimensión confirma relevancia biológica pero tiene bajo poder discriminativo en el ranking."

### Step 1: Abrir y editar METHODS.md

Localizar las secciones relevantes y añadir los párrafos arriba en su contexto correcto.

### Step 2: Commit

```bash
git add docs/METHODS.md
git commit -m "docs: update METHODS.md with pi_stat directionality rationale and limitations

- Clarify paired design block (consensus correlation 0.016 in DIA proteomics)
- Document pi_stat signed min-max normalization and dual-hypothesis architecture
- Document s_pathway saturation as acknowledged limitation"
```

---

## Verificación Final

El pipeline está listo para búsqueda bibliográfica manual cuando:

```bash
# 1. Miosinas cardíacas ausentes del top 35
grep -i "mavacamten\|omecamtiv\|danicamtiv\|reldesemtiv\|tirasemtiv\|levosimendan" \
  results/tables/10_top20_candidates.tsv
# Expected: sin resultados

# 2. DIGOXIN como Clase C en panel LOD-stable
awk -F'\t' 'NR==1 || /DIGOXIN/' results/tables/15_lod_stability.tsv | cut -f1,3 | column -t
# Expected: DIGOXIN  TRUE  con clase C en tabla de scoring

# 3. s_pi_stat de EGFR inhibitors > 0.70
awk -F'\t' 'NR>1 {print $1, $16}' results/tables/10_top20_candidates.tsv | \
  grep -i "gefitinib\|lapatinib\|afatinib" | column -t
# Expected: todos > 0.70

# 4. Panel LOD-stable de 25-27 candidatos
awk -F'\t' '$NF=="TRUE"' results/tables/15_lod_stability.tsv | wc -l
# Expected: entre 23 y 28
```
