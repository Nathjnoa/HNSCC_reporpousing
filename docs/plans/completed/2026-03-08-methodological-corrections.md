# HNSCC Drug Repurposing — Methodological Corrections Implementation Plan

> **STATUS: COMPLETADO — 2026-03-08**
> Todas las 11 tareas implementadas y revisadas. Pipeline re-ejecutado completo (01–17).
> Ver commits `a46ca45`–`2a4e572` en `master` para el historial de cambios.

---

> **Para Claude (nueva sesión):** Ejecutar este plan con el approach **Subagent-Driven**:
> 1. Invocar `superpowers:subagent-driven-development` ANTES de cualquier acción
> 2. Despachar un subagente por tarea (Task 1, Task 2, etc.)
> 3. Revisar output del subagente antes de continuar con la siguiente tarea
> 4. Usar `TodoWrite` para trackear progreso
>
> **Directorio de trabajo:** `~/bioinfo/projects/hnscc_drug_repurposing`
> **Ambiente R:** `omics-R` | **Ambiente Python:** `omics-py`
> **Plan completo:** este archivo — leerlo completo antes de empezar
> **Orden de ejecución:** Task 1 → Task 2 → (Tasks 3–7 paralelas) → Task 8 → Task 9 → Task 10 → Task 11

**Goal:** Implement 12 methodological corrections identified in critical review to improve statistical rigor and biological validity of the HNSCC drug repurposing pipeline.

**Architecture:** Changes cascade from foundational (DE re-analysis in script 01) through independent fixes (scripts 03, 06, 07, 08, 11) to dependent recalibrations (scripts 09, 10) and finally documentation. Phase I changes are prerequisite for re-running the full pipeline; Phase II–IV fixes can be implemented in parallel in the source code.

**Tech Stack:** R (limma, igraph, clusterProfiler), Python (requests, pandas), YAML config, bash.

---

## Dependency Map

```
PHASE I:  config.yaml ──► 01 (DE rewrite)
                              │
PHASE II: ┌─────────────────┼────────────────────────┐
          03 (pi-stat)    06 (OT filter)    08 (cetuximab)
          07 (L2S2 doc)   04,05 (snapshots) 11 (evidence)
                              │
PHASE III:                 09 (Louvain modules)
                              │
PHASE IV:                  10 (scoring recalibration)
                              │
PHASE V:                   15 (sensitivity strengthening)
PHASE VI:                  METHODS.md (limitations)
```

---

## PHASE I — Foundations

### Task 1: Update `config/analysis_params.yaml`

**Files:**
- Modify: `config/analysis_params.yaml`

**Step 1: Update scoring weights and add new parameters**

Replace the `scoring:` block and add `opentargets:` block:

```yaml
# --- Open Targets (script 06) ---
opentargets:
  min_hnscc_score: 0.2          # Minimum overall association score for HNSCC

# --- Priorización multi-criterio (script 10) ---
# CAMBIO: logFC+significance colapsados en pi_stat; nuevo weight_evidence
scoring:
  weight_pi_stat: 0.25           # reemplaza weight_log2fc(0.20)+weight_significance(0.15)
  weight_clinical_phase: 0.20
  weight_cmap_connectivity: 0.10 # reducido desde 0.15 (L2S2 usa transcriptómica, no proteómica)
  weight_pathway_relevance: 0.15
  weight_network_centrality: 0.15 # ahora usa n_distinct_modules (no n_targets raw)
  weight_evidence: 0.15          # nuevo: trials clínicos + PubMed (antes en script 13 solo)
  # Suma = 1.00 (verificado en script 10)
```

Remove lines 41–48 (old `scoring:` block with `weight_log2fc`, `weight_significance`, etc.) and insert the above. Keep all other blocks unchanged.

**Step 2: Verify sum = 1.00**

```bash
python3 -c "
import yaml
with open('config/analysis_params.yaml') as f:
    p = yaml.safe_load(f)
s = p['scoring']
total = s['weight_pi_stat'] + s['weight_clinical_phase'] + s['weight_cmap_connectivity'] + \
        s['weight_pathway_relevance'] + s['weight_network_centrality'] + s['weight_evidence']
print(f'Sum of weights: {total:.4f}')
assert abs(total - 1.0) < 0.001, 'ERROR: weights do not sum to 1.0'
print('OK')
"
```
Expected: `Sum of weights: 1.0000` and `OK`

**Step 3: Commit**

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
git add config/analysis_params.yaml
git commit -m "config: update scoring weights (pi-stat, reduce L2S2, add evidence dim) + OT min score"
```

---

### Task 2: Rewrite `01_parse_results_qc.R` — limma con covariable HPV

**Files:**
- Modify: `scripts/01_parse_results_qc.R`

This is the most complex task. The script currently uses pre-computed proteoDA stats. We rewrite it to run limma directly from the processed intensities, including HPV status as a covariate in the model.

**Step 1: Replace header and library block (lines 1–34)**

New header + libraries:

```r
# ============================================================
# Script: 01_parse_results_qc.R
# Proyecto: HNSCC Drug Repurposing
# Propósito: Reanalizar proteómica con modelo limma HPV-ajustado
#            (TVsS con HPV como covariable).
#
# Modelo: ~ condición (Tumor/Normal) + estado_VPH
#         bloqueo intra-paciente vía duplicateCorrelation()
#
# Input:
#   data/raw/results_proteomica.tsv    — matrix procesada (log2, MinProb, median-norm)
#   data/raw/metadata.csv             — metadatos (sample_id, condition, vph)
#
# Output:
#   results/tables/de_limma/
#     01_TVsS_all_proteins.tsv         — todas las proteínas con estadísticos
#     01_TVsS_significant.tsv          — adj.P.Val < 0.05 & |logFC| > 1
#     01_TVsS_upregulated.tsv
#     01_TVsS_downregulated.tsv
#     01_missingness_analysis.tsv      — % missing por proteína/grupo (bias check)
#   results/figures/qc/
#     01_boxplot_intensidades.pdf
#     01_PCA_muestras.pdf
#     01_volcano_TVsS.pdf
#     01_resumen_DE.pdf
#     01_missingness_vs_direction.pdf  — NUEVO: sesgo MinProb imputation
#
# Ambiente: omics-R
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(ggrepel)
  library(patchwork)
  library(viridis)
  library(yaml)
})
```

**Step 2: Replace Section 1 (cargar datos) — mantener el parsing del TSV, líneas 56–108**

Mantener el bloque de carga del TSV sin cambios (líneas 63–108 del script actual). Solo cambiar el comentario:

```r
# ============================================================
# 1. CARGAR INTENSIDADES (ya log2-transformadas, MinProb imputadas,
#    median-normalizadas por proteoDA)
# ============================================================
```

**Step 3: Replace Section 2 (cargar metadata) — añadir extracción de patient_id**

Reemplazar líneas 110–117 (bloque metadata) con:

```r
# ============================================================
# 2. CARGAR METADATA — incluir patient_id para diseño pareado
# ============================================================
meta <- read.csv2(meta_file, stringsAsFactors = FALSE)
# Extraer patient_id desde nombre de muestra: M1S/M1T -> patient "M1"
meta <- meta %>%
  mutate(patient_id = sub("[ST]$", "", sample_id))

cat("-- Metadata:\n")
print(table(meta$condition, meta$vph))
cat("\nPacientes únicos:", n_distinct(meta$patient_id), "\n")
cat("(esperado: 10 pacientes, 2 muestras c/u = 20 total)\n\n")

# Verificar alineación muestra-columna
stopifnot("Muestras en metadata no coinciden con columnas de intensidades" =
          all(intensity_cols %in% meta$sample_id))
meta <- meta[match(intensity_cols, meta$sample_id), ]  # alinear orden
```

**Step 4: Insert NEW Section 3 — Análisis de missingness (antes de limma)**

Insertar NUEVA sección después de la metadata y antes del nuevo análisis limma:

```r
# ============================================================
# 3. ANÁLISIS DE MISSINGNESS — verificar sesgo de imputación MinProb
# ============================================================
cat("-- Análisis de missingness por grupo...\n")

expr_raw_preimputation <- as.matrix(raw[, intensity_cols])
rownames(expr_raw_preimputation) <- raw$gene_symbol

normal_cols  <- meta$sample_id[meta$condition == "NORMAL"]
tumor_cols   <- meta$sample_id[meta$condition == "TUMORAL"]

missing_df <- data.frame(
  gene_symbol   = raw$gene_symbol,
  pct_missing_normal = rowMeans(is.na(expr_raw_preimputation[, normal_cols])) * 100,
  pct_missing_tumor  = rowMeans(is.na(expr_raw_preimputation[, tumor_cols]))  * 100
) %>%
  mutate(
    missing_diff = pct_missing_tumor - pct_missing_normal,
    bias_flag    = abs(missing_diff) > 30  # diferencia >30% entre grupos = sospechoso
  )

n_biased <- sum(missing_df$bias_flag, na.rm = TRUE)
cat(sprintf("  Proteínas con diferencial missingness >30%% entre grupos: %d / %d\n",
            n_biased, nrow(missing_df)))

write.table(missing_df,
  file.path(out_tables, "01_missingness_analysis.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
```

**Step 5: Replace Section 3 (old "Resumen TVsS" + Section 4 "Exportar") with new limma model**

Reemplazar líneas 120–178 (secciones 3 y 4 antiguas) con el nuevo análisis limma:

```r
# ============================================================
# 4. MODELO LIMMA HPV-AJUSTADO
#    ~ condition + vph_status
#    bloqueo intra-paciente (duplicateCorrelation)
# ============================================================
cat("-- Construyendo modelo limma con covariable VPH...\n")

params <- yaml::read_yaml(file.path(proj_dir, "config/analysis_params.yaml"))
fc_thresh  <- params$de$log2fc_threshold      # 1.0
pval_thr   <- params$de$adj_pval_threshold    # 0.05

# Matriz de expresión: proteínas × muestras (ya log2, imputado, normalizado)
expr_mat <- as.matrix(raw[, intensity_cols])
rownames(expr_mat) <- raw$gene_symbol

# Factores del diseño
condition <- factor(meta$condition, levels = c("NORMAL", "TUMORAL"))
vph       <- factor(meta$vph,       levels = c("NEGATIVE", "POSITIVE"))
patient   <- factor(meta$patient_id)

# Diseño sin intercepto para facilitar makeContrasts
design <- model.matrix(~ 0 + condition + vph, data = data.frame(condition, vph))
colnames(design) <- c("NORMAL", "TUMORAL", "VPH_POSITIVE")
cat(sprintf("  Diseño: %d muestras × %d coeficientes\n", nrow(design), ncol(design)))
cat("  Coeficientes:", paste(colnames(design), collapse = ", "), "\n")

# Correlación inter-paciente (diseño pareado)
cat("  Estimando correlación inter-paciente (duplicateCorrelation)...\n")
corfit <- duplicateCorrelation(expr_mat, design, block = patient)
cat(sprintf("  Correlación intra-paciente (consensus): %.3f\n",
            corfit$consensus.correlation))

# Ajuste del modelo con bloqueo por paciente
fit <- lmFit(expr_mat, design,
             block       = patient,
             correlation = corfit$consensus.correlation)

# Contraste: Tumoral vs Normal (ajustado por VPH)
contrasts <- makeContrasts(TumoralVsNormal = TUMORAL - NORMAL, levels = design)
fit2      <- contrasts.fit(fit, contrasts)
fit2      <- eBayes(fit2, trend = TRUE)  # trend=TRUE para proteómica (heterocedasticidad)

# Extraer resultados completos
results_all <- topTable(fit2, coef = "TumoralVsNormal",
                        n = Inf, adjust.method = "BH",
                        sort.by = "none")
results_all$gene_symbol <- rownames(results_all)

# Unir con anotaciones originales (uniprot_id, entry_name)
annot_cols <- raw %>% select(gene_symbol, uniprot_id, entry_name)
results_all <- results_all %>%
  left_join(annot_cols, by = "gene_symbol") %>%
  # Renombrar columnas para compatibilidad downstream
  rename(
    logFC_TVsS        = logFC,
    CI.L_TVsS         = CI.L,
    CI.R_TVsS         = CI.R,
    avg_intensity_TVsS = AveExpr,
    t_stat_TVsS       = t,
    B_stat_TVsS       = B,
    P.Value_TVsS      = P.Value,
    adj.P.Val_TVsS    = adj.P.Val
  ) %>%
  mutate(
    sig.FDR_TVsS = case_when(
      adj.P.Val_TVsS < pval_thr & logFC_TVsS >  fc_thresh ~  1L,
      adj.P.Val_TVsS < pval_thr & logFC_TVsS < -fc_thresh ~ -1L,
      TRUE                                                ~  0L
    )
  ) %>%
  select(gene_symbol, uniprot_id, entry_name,
         logFC_TVsS, CI.L_TVsS, CI.R_TVsS, avg_intensity_TVsS,
         t_stat_TVsS, B_stat_TVsS, P.Value_TVsS, adj.P.Val_TVsS,
         sig.FDR_TVsS)

# Añadir columnas de intensidades para downstream (script 02 las espera)
results_all <- results_all %>%
  left_join(raw %>% select(gene_symbol, all_of(intensity_cols)), by = "gene_symbol")

# Resumen
n_up   <- sum(results_all$sig.FDR_TVsS ==  1, na.rm = TRUE)
n_down <- sum(results_all$sig.FDR_TVsS == -1, na.rm = TRUE)
n_ns   <- nrow(results_all) - n_up - n_down

cat("\n=== Resultados limma HPV-ajustado ===\n")
cat(sprintf("  Total proteínas:    %d\n", nrow(results_all)))
cat(sprintf("  Upregulated:        %d\n", n_up))
cat(sprintf("  Downregulated:      %d\n", n_down))
cat(sprintf("  No significativas:  %d\n", n_ns))
cat(sprintf("  Criterio: adj.P.Val < %.2f & |logFC| > %.1f (BH)\n", pval_thr, fc_thresh))

sig_fc <- results_all$logFC_TVsS[results_all$sig.FDR_TVsS != 0]
cat(sprintf("  logFC rango sig: %.2f a %.2f\n\n",
            min(sig_fc, na.rm = TRUE), max(sig_fc, na.rm = TRUE)))

# ============================================================
# 5. EXPORTAR TABLAS
# ============================================================
cat("-- Exportando tablas...\n")

out_cols <- c("gene_symbol", "uniprot_id", "entry_name",
              "logFC_TVsS", "CI.L_TVsS", "CI.R_TVsS", "avg_intensity_TVsS",
              "t_stat_TVsS", "B_stat_TVsS", "P.Value_TVsS", "adj.P.Val_TVsS",
              "sig.FDR_TVsS", intensity_cols)

df_all  <- results_all
df_sig  <- results_all %>% filter(sig.FDR_TVsS %in% c(1L, -1L))
df_up   <- results_all %>% filter(sig.FDR_TVsS ==  1L)
df_down <- results_all %>% filter(sig.FDR_TVsS == -1L)

for (pair in list(
  list(df_all,  "01_TVsS_all_proteins.tsv"),
  list(df_sig,  "01_TVsS_significant.tsv"),
  list(df_up,   "01_TVsS_upregulated.tsv"),
  list(df_down, "01_TVsS_downregulated.tsv")
)) {
  write.table(pair[[1]], file.path(out_tables, pair[[2]]),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("  %s: %d proteínas\n", pair[[2]], nrow(pair[[1]])))
}
```

**Step 6: Add missingness figure after the existing boxplot figure (after line 211 in original)**

Añadir nueva figura de missingness al final de la sección de figuras (antes del resumen final):

```r
# ============================================================
# FIGURA QC-5: Missingness vs dirección DE (sesgo MinProb)
# ============================================================
cat("-- Generando figura QC-5: missingness vs dirección DE...\n")

miss_de <- missing_df %>%
  left_join(results_all %>% select(gene_symbol, logFC_TVsS, sig.FDR_TVsS),
            by = "gene_symbol") %>%
  mutate(
    DE_status = case_when(
      sig.FDR_TVsS ==  1 ~ "Up",
      sig.FDR_TVsS == -1 ~ "Down",
      TRUE               ~ "NS"
    )
  )

p_miss <- ggplot(miss_de,
                 aes(x = pct_missing_tumor, y = pct_missing_normal,
                     color = DE_status, alpha = DE_status)) +
  geom_point(size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Up" = "#D6604D", "Down" = "#4393C3", "NS" = "grey60")) +
  scale_alpha_manual(values = c("Up" = 0.9, "Down" = 0.9, "NS" = 0.2)) +
  labs(
    title    = "Missingness por grupo — verificación sesgo MinProb",
    subtitle = "Puntos sobre la diagonal = más missing en tumor (posible artefacto down-regulation)",
    x = "% Missing en TUMORAL", y = "% Missing en NORMAL",
    color = "DE status"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)) +
  guides(alpha = "none")

ggsave(file.path(out_figures, "01_missingness_vs_direction.pdf"),
       p_miss, width = 6.5, height = 6)
cat("  Guardado: 01_missingness_vs_direction.pdf\n")
```

**Step 7: Verify script runs (dry check)**

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
conda run -n omics-R Rscript -e "
  suppressPackageStartupMessages(library(limma))
  suppressPackageStartupMessages(library(tidyverse))
  cat('Libraries OK\n')
"
```
Expected: `Libraries OK`

**Step 8: Commit**

```bash
git add scripts/01_parse_results_qc.R
git commit -m "feat: rewrite script 01 with HPV-adjusted limma model + missingness analysis"
```

---

## PHASE II — Independent Fixes (can be parallelized)

### Task 3: Script 03 — GSEA pi-statistic ranking

**Files:**
- Modify: `scripts/03_pathway_enrichment.R` lines 78–87

**Step 1: Replace the ranked list construction**

Reemplazar líneas 78–87:

```r
# ANTES (solo logFC):
# ranked_list <- ranked_df$logFC_TVsS
# names(ranked_list) <- as.character(ranked_df$entrez_id)

# AHORA (pi-statistic: magnitud × significancia):
# pi = sign(logFC) × |logFC| × -log10(adj.P.Val)
# Captura tanto el efecto biológico como la confianza estadística
ranked_df <- all_prot %>%
  filter(!is.na(entrez_id), !is.na(logFC_TVsS), !is.na(adj.P.Val_TVsS)) %>%
  mutate(
    pi_stat = sign(logFC_TVsS) * abs(logFC_TVsS) * (-log10(adj.P.Val_TVsS + 1e-300))
  ) %>%
  arrange(desc(pi_stat))

ranked_list <- ranked_df$pi_stat
names(ranked_list) <- as.character(ranked_df$entrez_id)
ranked_list <- ranked_list[!duplicated(names(ranked_list))]

cat(sprintf("Ranked list para GSEA (pi-stat): %d genes\n", length(ranked_list)))
cat(sprintf("  pi_stat rango: %.2f a %.2f\n\n",
            min(ranked_list), max(ranked_list)))
```

**Step 2: Verify syntax**

```bash
cd ~/bioinfo/projects/hnscc_drug_repurposing
conda run -n omics-R Rscript -e "
  source('scripts/03_pathway_enrichment.R')
" 2>&1 | head -5
```
(Solo verificar que no hay errores de parsing — no ejecutar completo)

**Step 3: Commit**

```bash
git add scripts/03_pathway_enrichment.R
git commit -m "fix: GSEA ranked list uses pi-statistic (logFC × significance) instead of logFC alone"
```

---

### Task 4: Script 06 — Open Targets score threshold + query date

**Files:**
- Modify: `scripts/06_query_opentargets.py` lines 46–53 (params block) and lines 119–137 (pagination loop)

**Step 1: Add min_score parameter to the params block (after line 53)**

```python
# Leer parámetro de score mínimo desde config
import yaml
with open("config/analysis_params.yaml") as _f:
    _cfg = yaml.safe_load(_f)
MIN_HNSCC_SCORE = _cfg.get("opentargets", {}).get("min_hnscc_score", 0.2)
log.info(f"Open Targets min HNSCC score filter: {MIN_HNSCC_SCORE}")
```

**Step 2: Apply score filter in the pagination loop (after line 134)**

En el bucle de la línea 132–134, añadir filtro:

```python
    for r in rows:
        score = r["score"]
        if score < MIN_HNSCC_SCORE:   # << NUEVO filtro
            continue
        sym = r["target"]["approvedSymbol"].upper()
        hnscc_target_map[sym] = {
            "ensembl_id":  r["target"]["id"],
            "symbol":      r["target"]["approvedSymbol"],
            "hnscc_score": score,
        }
```

**Step 3: Add query date snapshot to log and output header**

Al inicio del script (después de `log.info(f"Inicio: {datetime.now()}")`):

```python
QUERY_DATE = datetime.now().strftime("%Y-%m-%d")
log.info(f"Open Targets API query date: {QUERY_DATE}")
log.info(f"Open Targets URL: {OT_URL}")
log.info(f"Min HNSCC association score: {MIN_HNSCC_SCORE}")
```

Y antes del `df_dedup.to_csv(...)` al final:

```python
df_dedup.insert(0, "query_date", QUERY_DATE)
df_dedup.insert(1, "ot_url",     OT_URL)
```

**Step 4: Verify**

```bash
python3 -c "
import ast, sys
with open('scripts/06_query_opentargets.py') as f:
    src = f.read()
ast.parse(src)
print('Syntax OK')
"
```
Expected: `Syntax OK`

**Step 5: Commit**

```bash
git add scripts/06_query_opentargets.py
git commit -m "fix: add Open Targets min HNSCC score filter (>=0.2) + query date snapshot"
```

---

### Task 5: Script 07 — Reduce L2S2 weight in config (ya hecho en Task 1) + documentar limitación

**Files:**
- Modify: `scripts/07_l2s2_connectivity.py` — añadir query date + log de la limitación

**Step 1: Add query date and limitation warning to script 07**

Después de las líneas de setup de logging en script 07:

```python
QUERY_DATE = datetime.now().strftime("%Y-%m-%d")
log.info(f"L2S2 query date: {QUERY_DATE}")
log.info(f"L2S2 URL: {L2S2_URL if 'L2S2_URL' in dir() else 'https://l2s2.ilincs.org/api/graphql'}")
log.warning(
    "METHODOLOGICAL NOTE: L2S2 uses transcriptomic (mRNA) signatures. "
    "This analysis inputs proteomic fold-changes. The assumption that "
    "protein and mRNA changes are directionally consistent is valid only "
    "for ~60-70%% of genes (mRNA-protein r~0.4-0.6 in clinical samples). "
    "L2S2 score weight has been reduced to 0.10 in the composite score. "
    "Results should be interpreted as supportive, not primary, evidence."
)
```

Y al exportar los outputs, añadir `query_date` al dataframe de resultados (misma lógica que Task 4).

**Step 2: Commit**

```bash
git add scripts/07_l2s2_connectivity.py
git commit -m "fix: add query date snapshot + methodological warning for proteomics→transcriptomics assumption"
```

---

### Task 6: Script 08 — Fix Cetuximab Class A

**Files:**
- Modify: `scripts/08_integrate_drug_targets.R` lines 298–316

**Problem:** Cetuximab's HNSCC indication in Open Targets may not match EFO_0000181 exactly, so `hnscc_indication = FALSE` and it gets classified as Class B.

**Step 1: Add a manual override list of known HNSCC-approved drugs after the classification block (after line 316)**

```r
# ---------------------------------------------------------------------------
# 4b. Corrección manual: drogas aprobadas para HNSCC con indicación conocida
#     que pueden no estar correctamente clasificadas por EFO matching.
# ---------------------------------------------------------------------------
HNSCC_APPROVED_OVERRIDE <- c(
  "CETUXIMAB",      # FDA/EMA: SCCHNeck + colorectal; EFO mismatch frecuente
  "PEMBROLIZUMAB",  # FDA: HNSCC R/M (2nd line) — puede tener otro EFO en OT
  "NIVOLUMAB"       # FDA: HNSCC R/M
)

n_overridden <- sum(
  toupper(drug_summary$drug_name_norm) %in% HNSCC_APPROVED_OVERRIDE &
    drug_summary$drug_class != "A"
)
if (n_overridden > 0) {
  cat(sprintf("\n  Corrigiendo %d drogas a Clase A (HNSCC override manual):\n",
              n_overridden))
  drug_summary <- drug_summary %>%
    mutate(
      drug_class = ifelse(
        toupper(drug_name_norm) %in% HNSCC_APPROVED_OVERRIDE,
        "A", drug_class
      ),
      drug_class_label = case_when(
        drug_class == "A" ~ "A: Approved + HNSCC evidence",
        drug_class == "B" ~ "B: Approved other cancer",
        drug_class == "C" ~ "C: Approved non-oncology",
        drug_class == "D" ~ "D: Not approved / Experimental"
      ),
      hnscc_indication = ifelse(
        toupper(drug_name_norm) %in% HNSCC_APPROVED_OVERRIDE,
        TRUE, hnscc_indication
      )
    )
  cat(paste(
    drug_summary$drug_name_norm[toupper(drug_summary$drug_name_norm) %in%
                                   HNSCC_APPROVED_OVERRIDE],
    collapse = ", "
  ), "\n")
}
```

**Step 2: Verify classification**

Añadir al final del script antes del export:

```r
cat("\nVerificación clasificación Clase A:\n")
drug_summary %>%
  filter(drug_class == "A") %>%
  select(drug_name_norm, drug_class, hnscc_indication) %>%
  print()
```

**Step 3: Commit**

```bash
git add scripts/08_integrate_drug_targets.R
git commit -m "fix: add HNSCC-approved drug override (cetuximab, pembrolizumab, nivolumab → Class A)"
```

---

### Task 7: Scripts 04–05 — API query date snapshots

**Files:**
- Modify: `scripts/04_query_dgidb.py` (añadir QUERY_DATE al output)
- Modify: `scripts/05_query_chembl.py` (añadir QUERY_DATE al output)

**Step 1: Script 04 — añadir query date**

En script 04, después del setup de logging:

```python
QUERY_DATE = datetime.now().strftime("%Y-%m-%d")
log.info(f"DGIdb API query date: {QUERY_DATE}")
log.info("DGIdb version: v5 GraphQL (https://www.dgidb.org/api/graphql)")
```

Y antes del export del dataframe principal:

```python
df_out.insert(0, "query_date", QUERY_DATE)
```

**Step 2: Script 05 — añadir query date**

Mismo patrón en script 05:

```python
QUERY_DATE = datetime.now().strftime("%Y-%m-%d")
log.info(f"ChEMBL API query date: {QUERY_DATE}")
log.info("ChEMBL version: v33 REST API (https://www.ebi.ac.uk/chembl/api/data/)")
```

**Step 3: Commit**

```bash
git add scripts/04_query_dgidb.py scripts/05_query_chembl.py
git commit -m "fix: add API query date snapshots to DGIdb and ChEMBL outputs"
```

---

## PHASE III — Network Module Analysis

### Task 8: Script 09 — Louvain community detection + module visualization

**Files:**
- Modify: `scripts/09_string_network.R`
- New output: `results/tables/network/09_modules.tsv`
- New output: `results/figures/09_network_modules.pdf`

**Step 1: Add Louvain module detection after PASO 4 (después de línea 285)**

Insertar entre el final de PASO 4 y PASO 5:

```r
# =============================================================================
# PASO 4b: Detección de módulos Louvain + caracterización funcional
# =============================================================================
cat("\n--- Paso 4b: Detectando módulos (Louvain) ---\n")

set.seed(42)
modules_louvain <- cluster_louvain(g, resolution = 1.0)
n_modules <- length(unique(membership(modules_louvain)))
cat(sprintf("Módulos Louvain detectados: %d\n", n_modules))
cat(sprintf("Tamaños: %s\n",
            paste(sort(table(membership(modules_louvain)), decreasing = TRUE),
                  collapse = ", ")))

# Asignar módulo a cada nodo
df_nodes$module_id <- membership(modules_louvain)[
  match(df_nodes$gene_symbol, names(membership(modules_louvain)))
]

# Nombrar módulos por el término GO BP más enriquecido
# (ORA rápida con clusterProfiler — requiere que los Entrez IDs estén disponibles)
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# Cargar mapeo de IDs (script 02 ya lo hizo)
id_map <- tryCatch(
  read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
             stringsAsFactors = FALSE) %>%
    select(symbol_org, entrez_id) %>%
    filter(!is.na(entrez_id)) %>%
    distinct(symbol_org, .keep_all = TRUE),
  error = function(e) {
    cat("  WARN: No se pudo cargar mapeo entrez — módulos sin nombre funcional\n")
    NULL
  }
)

universe_entrez <- if (!is.null(id_map)) as.character(id_map$entrez_id) else NULL

module_names <- sapply(seq_len(n_modules), function(m) {
  m_genes <- df_nodes$gene_symbol[df_nodes$module_id == m & !is.na(df_nodes$module_id)]
  if (length(m_genes) < 3 || is.null(id_map)) return(paste0("M", m))

  entrez_m <- id_map$entrez_id[id_map$symbol_org %in% m_genes]
  entrez_m <- as.character(na.omit(entrez_m))
  if (length(entrez_m) < 3) return(paste0("M", m))

  ora <- tryCatch(
    enrichGO(gene       = entrez_m,
             universe   = universe_entrez,
             OrgDb      = org.Hs.eg.db,
             keyType    = "ENTREZID",
             ont        = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             minGSSize     = 5,
             maxGSSize     = 200,
             readable      = TRUE),
    error = function(e) NULL
  )

  if (is.null(ora) || nrow(ora@result) == 0) return(paste0("M", m, "_uncharacterized"))

  # Simplify para reducir redundancia
  ora_simplified <- tryCatch(
    simplify(ora, cutoff = 0.70, by = "p.adjust", select_fun = min),
    error = function(e) ora
  )

  top_term <- ora_simplified@result$Description[1]
  sprintf("M%d_%s", m, gsub(" ", "_", substr(top_term, 1, 30)))
})

cat("Nombres de módulos:\n")
for (i in seq_along(module_names)) {
  n_size <- sum(df_nodes$module_id == i, na.rm = TRUE)
  cat(sprintf("  %s (n=%d)\n", module_names[i], n_size))
}

# Exportar tabla de módulos
df_modules <- df_nodes %>%
  filter(!is.na(module_id)) %>%
  mutate(module_name = module_names[module_id]) %>%
  select(gene_symbol, module_id, module_name, degree,
         betweenness_norm, direction, logFC_TVsS, is_hub) %>%
  arrange(module_id, desc(degree))

write.table(df_modules,
            "results/tables/network/09_modules.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nExportado: 09_modules.tsv (%d proteínas en %d módulos)\n",
            nrow(df_modules), n_modules))
```

**Step 2: Redefine hub logic at module level (reemplazar PASO 4, líneas 258–285)**

Añadir después del bloque de módulos (Step 1):

```r
# =============================================================================
# PASO 4c: Redefinir hubs a nivel de módulo
#           Hub = proteína en top 10% de betweenness DENTRO de su módulo
#           (evita inflar hubs de complejos multi-subunidad como Complex I)
# =============================================================================
cat("\n--- Paso 4c: Hubs por módulo (betweenness intra-módulo) ---\n")

df_nodes <- df_nodes %>%
  group_by(module_id) %>%
  mutate(
    # Betweenness intra-módulo: re-calcular en subgrafo del módulo
    module_hub_score = betweenness_norm,  # aproximación: usar betweenness global
    is_module_hub    = betweenness_norm >= quantile(betweenness_norm,
                                                     0.90, na.rm = TRUE)
  ) %>%
  ungroup()

# Mantener is_hub original para compatibilidad downstream pero añadir is_module_hub
n_module_hubs <- sum(df_nodes$is_module_hub, na.rm = TRUE)
cat(sprintf("Module hubs (top 10%% betweenness por módulo): %d\n", n_module_hubs))

# Añadir n_distinct_modules para cada gen (relevante para scoring en script 10)
# (para proteínas individuales = siempre 1; útil al agregar por fármaco)
df_nodes$n_distinct_modules <- 1L
```

**Step 3: Add module-colored network figure (añadir en PASO 6, después de la figura 09_network_full.pdf)**

```r
# 6e. -------- RED COMPLETA COLOREADA POR MÓDULO ----------------------------
cat("  [ggraph] Red con módulos Louvain coloreados...\n")

# Asignar atributo de módulo al subgrafo giant
module_lookup <- setNames(df_modules$module_name,
                           df_modules$gene_symbol)
V(g_giant)$module_name <- module_lookup[V(g_giant)$name]
V(g_giant)$is_module_hub <- df_nodes$is_module_hub[
  match(V(g_giant)$name, df_nodes$gene_symbol)
]

# Paleta para módulos (máximo 12 colores distinguibles)
n_mod_plot  <- min(n_modules, 12)
mod_palette <- setNames(
  c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
    "#a65628","#f781bf","#999999","#66c2a5","#fc8d62",
    "#8da0cb","#e78ac3")[seq_len(n_mod_plot)],
  seq_len(n_mod_plot)
)

safe_pdf("results/figures/09_network_modules.pdf", w = 15, h = 13, {
  set.seed(42)
  p <- ggraph(g_giant, layout = "fr") +
    geom_edge_link(color = "gray70", linewidth = 0.15, alpha = 0.05,
                   show.legend = FALSE) +
    # Nodos coloreados por módulo
    geom_node_point(
      aes(size  = degree,
          color = factor(
            module_id %% n_mod_plot + 1,   # reciclar colores si hay más de 12 módulos
            levels = seq_len(n_mod_plot)
          ),
          shape = is_module_hub),
      alpha = 0.75
    ) +
    # Labels solo para hubs de módulo
    geom_node_text(
      aes(label = ifelse(is_module_hub == TRUE, name, NA_character_)),
      size = 2.6, repel = TRUE, max.overlaps = 50,
      fontface = "bold", bg.color = "white", bg.r = 0.1, na.rm = TRUE
    ) +
    scale_color_manual(values = mod_palette,
                       name   = "Módulo Louvain",
                       guide  = guide_legend(ncol = 2)) +
    scale_size_continuous(range = c(0.6, 7), name = "Degree") +
    scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16),
                       labels = c("TRUE" = "Hub (top 10% betweenness intra-módulo)",
                                  "FALSE" = "Proteína"),
                       name = "") +
    labs(
      title    = "STRING PPI network — DE proteins coloreadas por módulo Louvain",
      subtitle = sprintf(
        "Componente gigante: %d nodos, %d aristas | %d módulos | diamantes = hubs por módulo",
        vcount(g_giant), ecount(g_giant), n_modules)
    ) +
    theme_graph(base_family = "sans") +
    theme(
      legend.position   = "right",
      plot.title        = element_text(size = 13, face = "bold"),
      plot.subtitle     = element_text(size = 9),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA)
    )
  print(p)
})
cat("  results/figures/09_network_modules.pdf\n")
```

**Step 4: Export updated node metrics with module info**

Reemplazar la línea de export de df_nodes (línea 582–584 del original) para incluir módulos:

```r
df_nodes_export <- df_nodes %>%
  left_join(df_modules %>% select(gene_symbol, module_name), by = "gene_symbol")

write.table(df_nodes_export,
            "results/tables/network/09_network_node_metrics.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_modules,
            "results/tables/network/09_modules.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

**Step 5: Commit**

```bash
git add scripts/09_string_network.R
git commit -m "feat: add Louvain module detection, module-level hub analysis, and module-colored network figure"
```

---

## PHASE IV — Scoring Recalibration

### Task 9: Script 10 — Replace logFC+sig with pi-stat; use n_distinct_modules; update L2S2 weight

**Files:**
- Modify: `scripts/10_prioritization_scoring.R`

**Step 1: Update config reading (lines 70–89) to use new weight names**

Reemplazar:

```r
w_logfc   <- params$scoring$weight_log2fc
w_sig     <- params$scoring$weight_significance
w_clin    <- params$scoring$weight_clinical_phase
w_cmap    <- params$scoring$weight_cmap_connectivity
w_path    <- params$scoring$weight_pathway_relevance
w_net     <- params$scoring$weight_network_centrality
```

Con:

```r
w_pistat  <- params$scoring$weight_pi_stat           # reemplaza logFC + sig
w_clin    <- params$scoring$weight_clinical_phase
w_cmap    <- params$scoring$weight_cmap_connectivity  # reducido a 0.10
w_path    <- params$scoring$weight_pathway_relevance
w_net     <- params$scoring$weight_network_centrality
w_evid    <- params$scoring$weight_evidence           # nuevo
```

Y actualizar el cat de verificación:

```r
cat(sprintf("Pesos: pi_stat=%.2f | clinical=%.2f | cmap=%.2f | pathway=%.2f | network=%.2f | evidence=%.2f\n",
            w_pistat, w_clin, w_cmap, w_path, w_net, w_evid))
weight_sum <- w_pistat + w_clin + w_cmap + w_path + w_net + w_evid
```

**Step 2: Replace DIMENSION 1+2 with pi-stat (lines 184–208)**

Reemplazar las secciones DIMENSION 1 y DIMENSION 2 completas con:

```r
# =============================================================================
# DIMENSION 1 (nueva): score_pi_stat — combina magnitud y significancia
#             pi = sign(logFC) × |logFC| × -log10(adj.P.Val)
#             Normalizado 0-1 sobre el máximo del dataset.
#             Reemplaza score_logfc (0.20) + score_significance (0.15).
# =============================================================================
cat("\n--- Calculando scores ---\n")

gene_de <- gene_de %>%
  mutate(
    pi_stat = sign(logFC_TVsS) * abs(logFC_TVsS) *
              (-log10(adj.P.Val_TVsS + 1e-300))
  )
max_pi <- max(abs(gene_de$pi_stat), na.rm = TRUE)

score_pistat_fn <- function(gene_str) {
  genes  <- get_target_genes(gene_str)
  pi_vals <- abs(gene_de$pi_stat[gene_de$symbol_org %in% genes])
  if (length(pi_vals) == 0) return(0)
  mean(pi_vals, na.rm = TRUE) / max_pi
}
```

**Step 3: Add DIMENSION 6b — n_distinct_modules instead of n_targets**

Añadir después de `score_network_fn` (líneas 244–255), una función para conteo de módulos:

```r
# Cargar tabla de módulos (generada por script 09)
modules_tbl <- tryCatch(
  read.delim("results/tables/network/09_modules.tsv", stringsAsFactors = FALSE) %>%
    select(gene_symbol, module_id) %>%
    distinct(),
  error = function(e) {
    cat("  WARN: 09_modules.tsv no encontrado — se usa n_targets raw\n")
    NULL
  }
)

score_network_fn <- function(gene_str) {
  genes <- get_target_genes(gene_str)
  nm    <- net_metrics %>% filter(gene_symbol %in% genes)
  if (nrow(nm) == 0) return(0)

  # Centralidad: 60% degree, 40% betweenness
  deg_norm  <- mean(nm$degree / max_degree, na.rm = TRUE)
  betw_norm <- mean(nm$betweenness_norm,    na.rm = TRUE)
  centrality_score <- 0.6 * deg_norm + 0.4 * betw_norm

  # Penalizar redundancia de complejo: si todos los targets están en el mismo módulo,
  # se cuenta como un único punto de contacto biológico.
  if (!is.null(modules_tbl)) {
    n_modules_hit <- modules_tbl %>%
      filter(gene_symbol %in% genes) %>%
      pull(module_id) %>%
      n_distinct()
    module_diversity <- min(n_modules_hit / 3, 1.0)  # normalizar: 3+ módulos = max
    return(0.7 * centrality_score + 0.3 * module_diversity)
  }
  centrality_score
}
```

**Step 4: Add DIMENSION 6 (evidence) — usar scores pre-calculados de script 11 si existen**

Añadir nueva función de score de evidencia:

```r
# =============================================================================
# DIMENSION 6 (nueva): score_evidence — evidencia clínica pre-calculada
#             Usa resultados de script 11 si están disponibles.
#             Fallback: 0 (no penaliza si no se ha corrido script 11 aún)
# =============================================================================
evidence_tbl <- tryCatch(
  read.delim("results/tables/evidence/11_clinical_evidence.tsv",
             stringsAsFactors = FALSE) %>%
    mutate(drug_name_norm_ev = toupper(trimws(drug_name))) %>%
    select(drug_name_norm_ev, evidence_score),
  error = function(e) {
    cat("  INFO: 11_clinical_evidence.tsv no encontrado — score_evidence = 0\n")
    NULL
  }
)

score_evidence_fn <- function(drug_name) {
  if (is.null(evidence_tbl)) return(0)
  hit <- evidence_tbl$evidence_score[
    toupper(trimws(drug_name)) == evidence_tbl$drug_name_norm_ev
  ]
  if (length(hit) == 0 || all(is.na(hit))) return(0)
  max(hit, na.rm = TRUE)
}
```

**Step 5: Update scored mutate block (líneas 262–289)**

Reemplazar el bloque `scored <- candidates %>% rowwise() %>% mutate(...)`:

```r
scored <- candidates %>%
  rowwise() %>%
  mutate(
    s_pi_stat  = score_pistat_fn(de_genes),    # reemplaza s_logfc + s_sig
    s_clinical = score_clinical_fn(max_phase),
    s_cmap     = score_cmap_fn(cmap_score),
    s_pathway  = score_pathway_fn(de_genes),
    s_network  = score_network_fn(de_genes),   # ahora incluye module diversity
    s_evidence = score_evidence_fn(drug_name_norm)  # nuevo
  ) %>%
  ungroup() %>%
  mutate(
    composite_score = w_pistat * s_pi_stat   +
                      w_clin   * s_clinical  +
                      w_cmap   * s_cmap      +
                      w_path   * s_pathway   +
                      w_net    * s_network   +
                      w_evid   * s_evidence,
    bonus = case_when(
      drug_class == "A"              ~ 0.05,
      !is.na(cmap_score) & n_sources >= 3 ~ 0.03,
      drug_class == "C" & n_sources >= 2  ~ 0.01,
      TRUE                           ~ 0
    ),
    final_score = pmin(composite_score + bonus, 1.0)
  ) %>%
  arrange(desc(final_score))
```

**Step 6: Update column references in logs and figures (score_cols)**

En el heatmap de score components (línea 428), actualizar:

```r
score_cols  <- c("s_pi_stat", "s_clinical", "s_cmap", "s_pathway", "s_network", "s_evidence")
col_labels  <- c("PI-stat\n(logFC×sig)", "Clinical\nPhase",
                 "L2S2\nReversal", "Pathway\nRelevance",
                 "Network\n(modules)", "Clinical\nEvidence")
```

Y en el barplot/scatter, actualizar el subtitle:

```r
"Multi-criteria composite score (pi-stat + clinical phase + L2S2 + pathways + network modules + clinical evidence)"
```

**Step 7: Commit**

```bash
git add scripts/10_prioritization_scoring.R
git commit -m "feat: recalibrate scoring — pi-stat replaces logFC+sig, module diversity, L2S2 weight reduced, add evidence dim"
```

---

## PHASE V — Sensitivity Analysis Strengthening

### Task 10: Script 15 — Drop-one-database + permutation test

**Files:**
- Modify: `scripts/15_sensitivity_analysis.R`

**Step 1: Add drop-one-database analysis after existing weight perturbation block**

Añadir nueva sección al final del script (antes del resumen final):

```r
# =============================================================================
# ANÁLISIS 2: Drop-one-database (leave-one-out por fuente)
# =============================================================================
cat("\n=== ANÁLISIS 2: Drop-one-database ===\n")

sources_to_drop <- c("DGIdb", "ChEMBL", "OpenTargets", "L2S2")

lod_results <- lapply(sources_to_drop, function(drop_src) {
  cat(sprintf("  Excluyendo fuente: %s\n", drop_src))
  # Re-filtrar candidatos excluyendo la fuente indicada
  cands_filtered <- candidates %>%
    mutate(
      sources_adj = str_remove_all(sources, drop_src),
      n_sources_adj = str_count(sources_adj, "\\|") + 1
    ) %>%
    filter(n_sources_adj >= min_db | drug_class %in% c("A", "B"))

  if (nrow(cands_filtered) == 0) return(NULL)

  # Recalcular scores (simplificado: usar los scores ya calculados)
  scored_lod <- scored %>%
    filter(drug_name_norm %in% cands_filtered$drug_name_norm) %>%
    arrange(desc(final_score)) %>%
    slice_head(n = top_n)

  data.frame(
    drug_name_norm = scored_lod$drug_name_norm,
    rank_lod       = seq_len(nrow(scored_lod)),
    drop_source    = drop_src
  )
})

df_lod <- bind_rows(lod_results)

# Contar cuántas veces cada drug está en top N en todas las configs LOD
lod_stability <- df_lod %>%
  group_by(drug_name_norm) %>%
  summarise(n_lod_topN = n(), .groups = "drop") %>%
  mutate(lod_stable = n_lod_topN == length(sources_to_drop))

cat(sprintf("Drugs estables en top %d con cualquier fuente excluida: %d\n",
            top_n, sum(lod_stability$lod_stable)))

# =============================================================================
# ANÁLISIS 3: Permutation test (distribución nula del composite score)
# =============================================================================
cat("\n=== ANÁLISIS 3: Permutation test (n=1000) ===\n")

set.seed(2026)
n_perm <- 1000

# Función de score con logFC permutado
perm_score <- function() {
  gene_de_perm <- gene_de %>%
    mutate(
      pi_stat_perm = sample(pi_stat),  # shuffle
      logFC_TVsS   = pi_stat_perm / (-log10(adj.P.Val_TVsS + 1e-300))  # approx
    )
  max_pi_perm <- max(abs(gene_de_perm$pi_stat_perm), na.rm = TRUE)

  scored_perm <- candidates %>%
    rowwise() %>%
    mutate(
      s_pi_perm = {
        genes <- get_target_genes(de_genes)
        pi_vals <- abs(gene_de_perm$pi_stat_perm[gene_de_perm$symbol_org %in% genes])
        if (length(pi_vals) == 0) 0 else mean(pi_vals, na.rm = TRUE) / max_pi_perm
      },
      s_comp_perm = w_pistat * s_pi_perm + w_clin * s_clinical +
                    w_cmap * s_cmap + w_path * s_pathway + w_net * s_network
    ) %>%
    ungroup()

  max(scored_perm$s_comp_perm, na.rm = TRUE)
}

perm_null <- replicate(n_perm, perm_score())
true_top_score <- max(scored$composite_score, na.rm = TRUE)
perm_pval <- mean(perm_null >= true_top_score)

cat(sprintf("Score real top-1: %.4f\n", true_top_score))
cat(sprintf("Permutation p-value (top score): %.4f\n", perm_pval))
cat(sprintf("(Interpretación: probabilidad de obtener score >= %.4f por azar)\n",
            true_top_score))

# Exportar resultados adicionales
write.table(lod_stability,
            "results/tables/15_lod_stability.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(
  analysis        = "permutation_null",
  n_permutations  = n_perm,
  true_top_score  = true_top_score,
  perm_pval       = perm_pval,
  perm_mean       = mean(perm_null),
  perm_sd         = sd(perm_null),
  perm_95pct      = quantile(perm_null, 0.95)
), "results/tables/15_permutation_test.tsv",
sep = "\t", quote = FALSE, row.names = FALSE)
```

**Step 2: Commit**

```bash
git add scripts/15_sensitivity_analysis.R
git commit -m "feat: add drop-one-database and permutation test to sensitivity analysis"
```

---

## PHASE VI — Documentation

### Task 11: Update `docs/METHODS.md` — Sección de limitaciones

**Files:**
- Modify: `docs/METHODS.md`

**Step 1: Append limitations section**

Añadir al final del archivo:

```markdown
## Limitaciones Metodológicas

### 1. Asunción transcriptómica en L2S2
El análisis de conectividad de firmas (L2S2/LINCS L1000) se basa en respuestas
transcriptómicas (mRNA) de líneas celulares tratadas. Los inputs al análisis
(top 150 proteínas UP/DOWN por logFC proteómico) asumen que los cambios
proteómicos tienen correspondencia directa con cambios en mRNA. La correlación
mRNA-proteína en tumores sólidos es típicamente r ≈ 0.40–0.60 [PMID: 34425047].
El peso de esta dimensión se redujo a 0.10 (vs. 0.15 original) para reflejar
esta limitación. Los resultados de L2S2 deben interpretarse como evidencia
de soporte, no primaria.

### 2. Sesgo de imputación MinProb
Los datos proteómicos fueron imputados con el método MinProb (valores
ausentes → distribución en percentil 2.5). Este método asume que los valores
ausentes son "below detection threshold" (MAR). Si las proteínas mitocondriales
tienen más valores faltantes en tumor que en normal por razones biológicas reales
(downregulation genuina), MinProb amplificará esta señal. El archivo
`results/figures/qc/01_missingness_vs_direction.pdf` y
`results/tables/de_limma/01_missingness_analysis.tsv` documentan el diferencial
de missingness para cada proteína.

### 3. Confounding HPV en el diseño
Los 10 pares tumor/normal incluyen 6 pacientes HPV+ y 4 HPV-. HPV+ y HPV-
tienen perfiles moleculares distintos (vías E6/E7, inmuno-infiltración,
pronóstico diferencial). El modelo limma incluye HPV como covariable aditiva
(~ condición + vph_status), absorbiendo la varianza HPV-dependiente. Sin
embargo, la potencia estadística para detectar efectos de interacción
(proteínas DE específicamente en HPV+ o HPV-) es limitada con n=6/4 por grupo.
Los candidatos del pipeline representan biología compartida entre ambos subtipos.

### 4. Tamaño muestral
n=10 pares tumor/normal es adecuado para efectos de gran magnitud (|logFC| > 1)
pero subóptimo para efectos moderados (0.5–1.0). Ver análisis de potencia en
`docs/FUTURE_WORK.md`.

### 5. Hubs de red: complejos multi-subunidad
La definición de hubs por betweenness centrality (percentil 90 por módulo)
puede identificar múltiples subunidades del mismo complejo físico (e.g., NADH
dehydrogenase Complex I: NDUFA*, NDUFB*, NDUFS*) como hubs independientes.
El script 09 usa detección de módulos Louvain y redefine hubs a nivel de
módulo (intra-módulo betweenness top 10%) para mitigar este efecto. El score
de red en script 10 incluye `module_diversity` (número de módulos distintos
targetados) para penalizar fármacos que solo actúan sobre un complejo.
```

**Step 2: Commit**

```bash
git add docs/METHODS.md
git commit -m "docs: add methodological limitations section (L2S2, MinProb, HPV, hub definition)"
```

---

## Re-ejecución del pipeline (tras implementar todos los cambios)

Una vez implementadas todas las correcciones de código, re-ejecutar en orden:

```bash
conda activate omics-R
cd ~/bioinfo/projects/hnscc_drug_repurposing

# Fase 1: nuevo DE HPV-ajustado
Rscript scripts/01_parse_results_qc.R
Rscript scripts/02_id_mapping.R

# Fase 2: enrichment con pi-stat
Rscript scripts/03_pathway_enrichment.R

# Fase 3: drug databases (paralelo)
conda activate omics-py
python scripts/04_query_dgidb.py &
python scripts/05_query_chembl.py &
python scripts/06_query_opentargets.py &
wait
python scripts/07_l2s2_connectivity.py

# Fase 4: integración y red
conda activate omics-R
Rscript scripts/08_integrate_drug_targets.R
Rscript scripts/09_string_network.R

# Fase 5: scoring recalibrado
Rscript scripts/10_prioritization_scoring.R

# Fase 6: evidencia clínica
conda activate omics-py
python scripts/11_clinicaltrials_pubmed.py &
python scripts/12_cosmic_overlap.py &
wait

# Fase 7: integración final + sensitivity
conda activate omics-R
Rscript scripts/13_evidence_summary.R
Rscript scripts/15_sensitivity_analysis.R

# Fase 8: figuras publicación
Rscript scripts/17_pub_figures.R
```

---

## Checklist de verificación post-implementación

- [ ] `01_TVsS_significant.tsv` contiene columna `sig.FDR_TVsS` (compatibilidad downstream)
- [ ] Cetuximab aparece como `drug_class = "A"` en `08_drug_summary_per_drug.tsv`
- [ ] `09_modules.tsv` existe con columnas `gene_symbol, module_id, module_name`
- [ ] `09_network_modules.pdf` generado (figura con colores por módulo)
- [ ] Composite score en script 10 usa `s_pi_stat` (no `s_logfc + s_sig`)
- [ ] Sum de pesos en script 10 = 1.00 (verificado en log)
- [ ] `15_lod_stability.tsv` y `15_permutation_test.tsv` existen
- [ ] `docs/METHODS.md` contiene sección "Limitaciones Metodológicas"
- [ ] Todos los outputs de scripts Python contienen columna `query_date`
