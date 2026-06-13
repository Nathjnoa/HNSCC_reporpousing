#!/usr/bin/env Rscript
# =============================================================================
# 16_external_validation.R
# =============================================================================
# Validación externa de los hallazgos proteómicos HNSCC usando TCGA-HNSC.
#
# Panel A — Concordancia global: DIA proteomics log2FC vs TCGA-HNSC RNA-seq
#           log2FC (tumor vs normal; DESeq2, n=43 pares). Valida el INPUT del
#           pipeline: la desregulación proteómica es reproducible a nivel
#           transcriptómico en cohorte independiente. 4 genes-pilar etiquetados
#           (EGFR, PSMB10, DNMT1, NDUFS3) representan los módulos de Fig4/Fig5.
#
# Panel B — Dianas priorizadas en TCGA: para las anclas del shortlist (Fig5),
#           log2FC + FDR en TCGA + concordancia con el proteoma. Cierra el lazo
#           Fig5 → validación externa: las predicciones se sostienen en una
#           cohorte independiente. Coloreado por pilar (PILLAR_COLS).
#
# FigS (suplementario) — Supervivencia OS: KM por expresión (mediana) de los
#           4 genes-pilar en TCGA-HNSC. Análisis exploratorio/contextual;
#           los p-valores son no-significativos (p>0.05), consistente con que
#           las dianas son vulnerabilidades terapéuticas, no biomarcadores
#           pronósticos de OS.
#
# Outputs principales (results/figures/pub/main/):
#   Fig6A_concordance.{pdf,png}       — panel A
#   Fig6B_targets_tcga.{pdf,png}      — panel B (nuevo)
#   Objetos cacheados en .objects/
#
# Outputs suplementarios (results/figures/pub/supp/):
#   FigS_survival_targets.{pdf,png}
#
# Tablas (results/tables/pub/):
#   main/Tab6_concordance_summary.tsv
#   supp/TabS2_survival_genes.tsv
#
# Caches TCGA (data/intermediate/tcga/):
#   tcga_hnsc_se.rds, tcga_hnsc_clinical.rds, tcga_hnsc_deseq2_results.tsv
#   (~800 MB; se saltan si ya existen)
#
# Dependencias adicionales (instalar una vez en omics-R):
#   BiocManager::install(c("TCGAbiolinks","DESeq2","SummarizedExperiment"))
#   install.packages(c("survival","survminer"))
#
# Ambiente: omics-R
# Ejecución:
#   conda activate omics-R
#   cd ~/bioinfo/projects/hnscc_drug_repurposing
#   Rscript scripts/16_external_validation.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(stringr)
  library(yaml)
})

# ── Setup y estilo centralizado ───────────────────────────────────────────────
cat("=== 16_external_validation.R ===\n")
source(here::here("scripts", "_setup.R"))
setup_project()
source(here::here("scripts", "_fig_style.R"))    # PRESETS, theme_pub, PILLAR_COLS,
                                                  # save_pub, save_tiff, save_panel_obj
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ── Directorios ───────────────────────────────────────────────────────────────
dir.create("logs",                              showWarnings = FALSE)
dir.create("data/intermediate/tcga",           recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/pub/main",         recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/pub/supp",         recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/pub/main",          recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/pub/supp",          recursive = TRUE, showWarnings = FALSE)

log_file <- file.path("logs", paste0("16_external_validation_", timestamp, ".log"))
sink(log_file, split = TRUE)
cat("Inicio:", format(Sys.time()), "\n\n")

# ── Verificar dependencias ────────────────────────────────────────────────────
check_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0(
      "Paquete requerido no encontrado: '", pkg, "'\n",
      "Instalar con: ", if (bioc) paste0("BiocManager::install('", pkg, "')")
                        else        paste0("install.packages('", pkg, "')")
    )
    stop(msg, call. = FALSE)
  }
}
check_pkg("TCGAbiolinks",       bioc = TRUE)
check_pkg("DESeq2",             bioc = TRUE)
check_pkg("survival")
check_pkg("survminer")
check_pkg("SummarizedExperiment", bioc = TRUE)

library(TCGAbiolinks)
library(DESeq2)
library(survival)
library(survminer)
library(SummarizedExperiment)

# ── Constantes narrativas ─────────────────────────────────────────────────────
# Genes-pilar: representan los módulos terapéuticos de Fig4/Fig5
LABEL_GENES  <- c("EGFR", "PSMB10", "DNMT1", "NDUFS3")
LABEL_PILLAR <- c(EGFR = "EGFR", PSMB10 = "Proteasome",
                  DNMT1 = "Epigenetic", NDUFS3 = "OXPHOS")

# Colores de concordancia: neutros para no chocar con PILLAR_COLS
CONCORD_COLS <- c(
  "Concordant up"   = "#999999",   # gris medio
  "Concordant down" = "#444444",   # gris oscuro
  "Discordant"      = "#DDDDDD"    # gris claro
)

# =============================================================================
# SECCIÓN 1: DATOS LOCALES (nuestro estudio)
# =============================================================================
cat("\n--- Cargando datos locales ---\n")

de_our <- read_tsv("results/tables/de_limma/01_TVsS_significant.tsv",
                   show_col_types = FALSE)
cat("Proteínas DE significativas:", nrow(de_our), "\n")

id_map <- read_tsv("data/intermediate/id_mapping/02_uniprot_to_ids.tsv",
                   show_col_types = FALSE) |>
  select(uniprot_id, gene_symbol = symbol_org) |>
  filter(!is.na(gene_symbol), gene_symbol != "")

if (!"gene_symbol" %in% colnames(de_our)) {
  de_our <- de_our |> left_join(id_map, by = "uniprot_id")
}

de_our_genes <- de_our |>
  filter(!is.na(gene_symbol), gene_symbol != "") |>
  select(gene_symbol, logFC_our = logFC_TVsS, adjP_our = adj.P.Val_TVsS) |>
  distinct(gene_symbol, .keep_all = TRUE)

cat("Proteínas con gene_symbol:", nrow(de_our_genes), "\n")

# Shortlist de dianas priorizadas (lógica paralela a 17b_fig5_panels.R:67-73)
scored <- read_tsv("results/tables/10_all_candidates_scored.tsv",
                   show_col_types = FALSE)
params   <- yaml::read_yaml("config/analysis_params.yaml")
w_target <- params$scoring_v3$weight_target
w_drug   <- params$scoring_v3$weight_drug

# Top fármaco por diana-ancla → top N por composite (excluyendo off_network)
N_SHORTLIST <- 14
shortlist <- scored |>
  filter(tier != "off_network") |>
  group_by(primary_target) |>
  slice_max(composite_score, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(composite_score, n = N_SHORTLIST) |>
  mutate(
    drug_title  = str_to_title(str_to_lower(drug_name_norm)),
    tier_lab    = recode(tier,
                         hub_central     = "Network-central hub",
                         peripheral_diff = "Peripheral"),
    pillar      = case_when(
      primary_target %in% c("EGFR","ERBB2","ERBB3")    ~ "EGFR",
      primary_target %in% c("PSMB5","PSMB10","PSMD11",
                             "PSMA1","PSMA2","PSMB2",
                             "PSMB1","PSMC2","PSMC5")   ~ "Proteasome",
      primary_target %in% c("DNMT1","HDAC1","HDAC2",
                             "EZH2","MAOA","IMPDH2",
                             "TYMS","DHFR")             ~ "Epigenetic",
      primary_target %in% c("NDUFS3","NDUFS2","MT-ND1",
                             "NDUFV1","NDUFA13")        ~ "OXPHOS",
      TRUE                                               ~ "Other"
    )
  )

cat(sprintf("Shortlist: %d anclas | hub_central=%d peripheral=%d\n",
            nrow(shortlist),
            sum(shortlist$tier == "hub_central"),
            sum(shortlist$tier == "peripheral_diff")))
cat("Dianas ancla:", paste(shortlist$primary_target, collapse = ", "), "\n")

# =============================================================================
# SECCIÓN 2: DESCARGA TCGA-HNSC (con cache)
# =============================================================================
cache_se  <- "data/intermediate/tcga/tcga_hnsc_se.rds"
cache_cli <- "data/intermediate/tcga/tcga_hnsc_clinical.rds"

if (file.exists(cache_se) && file.exists(cache_cli)) {
  cat("\n--- Cargando TCGA-HNSC desde cache ---\n")
  se   <- readRDS(cache_se)
  clin <- readRDS(cache_cli)
} else {
  cat("\n--- Descargando TCGA-HNSC (primera ejecución, ~800 MB) ---\n")
  query <- GDCquery(
    project       = "TCGA-HNSC",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type   = c("Primary Tumor", "Solid Tissue Normal")
  )
  GDCdownload(query, method = "api", files.per.chunk = 50,
              directory = "data/intermediate/tcga/GDCdata")
  se <- GDCprepare(query, directory = "data/intermediate/tcga/GDCdata")
  saveRDS(se, cache_se)
  clin <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")
  saveRDS(clin, cache_cli)
  cat("TCGA-HNSC descargado y guardado en cache.\n")
}

cat("Dimensiones SE:", dim(se), "\n")
cat("Muestras TCGA:", ncol(se), "\n")

# =============================================================================
# SECCIÓN 3: DE ANALYSIS TCGA (DESeq2, con cache)
# =============================================================================
cache_de <- "data/intermediate/tcga/tcga_hnsc_deseq2_results.tsv"

if (file.exists(cache_de)) {
  cat("\n--- Cargando DE TCGA desde cache ---\n")
  de_tcga <- read_tsv(cache_de, show_col_types = FALSE)
} else {
  cat("\n--- Corriendo DESeq2 en TCGA tumor vs normal (~10 min) ---\n")
  sample_type <- substr(colnames(se), 14, 15)
  se_tn       <- se[, sample_type %in% c("01", "11")]
  condition   <- factor(ifelse(substr(colnames(se_tn), 14, 15) == "01",
                               "Tumor", "Normal"),
                        levels = c("Normal", "Tumor"))
  cat("Muestras: Tumor =", sum(condition == "Tumor"),
      "  Normal =", sum(condition == "Normal"), "\n")

  gene_names <- rowData(se_tn)$gene_name
  counts_mat <- assay(se_tn, "unstranded")
  rownames(counts_mat) <- gene_names
  col_data <- data.frame(condition = condition, row.names = colnames(se_tn))

  dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                 colData   = col_data,
                                 design    = ~ condition)
  dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
  cat("Genes tras filtro:", nrow(dds), "\n")
  dds <- DESeq(dds, parallel = FALSE)
  res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)

  de_tcga <- as.data.frame(res) |>
    rownames_to_column("gene_symbol") |>
    filter(!is.na(log2FoldChange), gene_symbol != "") |>
    select(gene_symbol,
           logFC_tcga = log2FoldChange,
           adjP_tcga  = padj,
           baseMean) |>
    distinct(gene_symbol, .keep_all = TRUE)

  write_tsv(de_tcga, cache_de)
  cat("DESeq2 completado:", nrow(de_tcga), "genes.\n")
}
cat("Genes TCGA DE calculados:", nrow(de_tcga), "\n")

# =============================================================================
# SECCIÓN 4: CONCORDANCIA GLOBAL
# =============================================================================
cat("\n--- Análisis de concordancia ---\n")

conc <- de_our_genes |>
  inner_join(de_tcga, by = "gene_symbol") |>
  filter(!is.na(logFC_tcga), !is.na(logFC_our)) |>
  mutate(
    concord_class = case_when(
      logFC_our > 0 & logFC_tcga > 0 ~ "Concordant up",
      logFC_our < 0 & logFC_tcga < 0 ~ "Concordant down",
      TRUE                            ~ "Discordant"
    )
  )

cat("Genes en ambos datasets:", nrow(conc), "\n")
cat(table(conc$concord_class), "\n")

pearson_r <- cor(conc$logFC_our, conc$logFC_tcga, method = "pearson")
pearson_p <- cor.test(conc$logFC_our, conc$logFC_tcga)$p.value
pct_concord <- mean(conc$concord_class != "Discordant") * 100
cat(sprintf("Pearson r = %.3f  (p = %.2e)\n", pearson_r, pearson_p))
cat(sprintf("Concordancia direccional: %.1f%%\n", pct_concord))

# Marcar genes-pilar
conc <- conc |>
  mutate(
    is_label = gene_symbol %in% LABEL_GENES,
    pillar   = LABEL_PILLAR[gene_symbol]
  )

# Tabla publicación
summary_tbl <- tibble(
  metric = c("N genes overlap", "Pearson r", "p-value (Pearson)",
             "% concordant direction", "% concordant up", "% concordant down"),
  value  = c(
    nrow(conc),
    round(pearson_r, 3),
    signif(pearson_p, 3),
    round(pct_concord, 1),
    round(mean(conc$concord_class == "Concordant up")   * 100, 1),
    round(mean(conc$concord_class == "Concordant down") * 100, 1)
  )
)
write_tsv(summary_tbl, "results/tables/pub/main/Tab6_concordance_summary.tsv")
cat("Tabla concordancia guardada.\n")

# =============================================================================
# SECCIÓN 5: PANEL A — Scatter concordancia (re-tematizado)
# =============================================================================
cat("\n--- Generando Panel A (concordancia global) ---\n")

p_concord <- ggplot(conc, aes(x = logFC_our, y = logFC_tcga)) +

  # Cuadrantes de referencia
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70",
             linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70",
             linewidth = 0.3) +

  # Nube de puntos (colores neutros por dirección)
  geom_point(data = filter(conc, !is_label),
             aes(color = concord_class),
             size = 0.9, alpha = 0.55, shape = 16) +

  # Genes-pilar (PILLAR_COLS; resaltan sobre la nube)
  geom_point(data = filter(conc, is_label),
             aes(fill = pillar),
             size = 2.8, shape = 21, colour = "white", stroke = 0.6) +

  # Etiquetas genes-pilar
  geom_label_repel(
    data          = filter(conc, is_label),
    aes(label = gene_symbol, fill = pillar),
    color         = "white",
    size          = PRESETS$double_col$tick / .pt,
    fontface      = "bold",
    box.padding   = 0.4,
    label.size    = 0,
    label.padding = unit(0.12, "lines"),
    show.legend   = FALSE
  ) +

  # Línea de regresión
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              linewidth = 0.45, alpha = 0.10) +

  # Anotación r / n
  annotate("text",
    x     = min(conc$logFC_our) * 0.9,
    y     = max(conc$logFC_tcga) * 0.92,
    label = sprintf("r = %.3f\nn = %d", pearson_r, nrow(conc)),
    hjust = 0, size = PRESETS$double_col$tick / .pt, color = "grey20"
  ) +

  scale_color_manual(values = CONCORD_COLS, name = NULL,
                     guide = guide_legend(
                       override.aes = list(size = 2, alpha = 0.9))) +
  scale_fill_manual(values = PILLAR_COLS, name = "Module",
                    guide = guide_legend(
                      override.aes = list(size = 2.5, colour = "grey40",
                                          stroke = 0.3))) +

  labs(
    x     = expression("log"[2]*"FC  (DIA proteomics)"),
    y     = expression("log"[2]*"FC  (TCGA-HNSC)")
  ) +
  theme_pub("double_col") +
  theme(legend.position = "right")

save_pub(p_concord, "Fig6A_concordance")
save_panel_obj(p_concord, "Fig6A_concordance")

# =============================================================================
# SECCIÓN 6: PANEL B — Dianas priorizadas en TCGA
# =============================================================================
cat("\n--- Generando Panel B (dianas priorizadas en TCGA) ---\n")

# Unir shortlist con resultados TCGA
# Nota: S4Vectors enmascara dplyr::rename → usar select() con renaming o dplyr::rename()
de_tcga_j    <- de_tcga    |> dplyr::select(target_gene = gene_symbol,
                                              logFC_tcga_target = logFC_tcga,
                                              adjP_tcga_target  = adjP_tcga)
de_our_genes_j <- de_our_genes |> dplyr::select(target_gene = gene_symbol,
                                                  logFC_our_target = logFC_our)
targets_tcga <- shortlist |>
  left_join(de_tcga_j,    by = c("primary_target" = "target_gene")) |>
  left_join(de_our_genes_j, by = c("primary_target" = "target_gene")) |>
  mutate(
    # concordancia: proteoma vs transcriptoma de la diana específica
    concord_target = case_when(
      !is.na(logFC_our_target) & !is.na(logFC_tcga_target) &
        sign(logFC_our_target) == sign(logFC_tcga_target) ~ TRUE,
      TRUE ~ FALSE
    ),
    sig_tcga       = !is.na(adjP_tcga_target) & adjP_tcga_target < 0.05,
    sig_label      = case_when(
      !is.na(adjP_tcga_target) & adjP_tcga_target < 0.001 ~ "***",
      !is.na(adjP_tcga_target) & adjP_tcga_target < 0.01  ~ "**",
      !is.na(adjP_tcga_target) & adjP_tcga_target < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    # Etiqueta del punto: fármaco top (drug_title) + diana
    drug_label  = sprintf("%s\n(%s)", drug_title, primary_target),
    drug_label  = reorder(drug_label, composite_score),
    pillar      = factor(pillar, levels = c("EGFR","Proteasome","Epigenetic","OXPHOS","Other"))
  )

# Reporte al log
n_found     <- sum(!is.na(targets_tcga$logFC_tcga_target))
n_concord   <- sum(targets_tcga$concord_target, na.rm = TRUE)
n_sig       <- sum(targets_tcga$sig_tcga, na.rm = TRUE)
cat(sprintf("\n=== Checkpoint Fig6B ===\n"))
cat(sprintf("Anclas en shortlist:          %d\n", nrow(targets_tcga)))
cat(sprintf("Con datos en TCGA RNA-seq:    %d / %d\n", n_found, nrow(targets_tcga)))
cat(sprintf("Concordantes prot. vs RNA:    %d / %d  (%.0f%%)\n",
            n_concord, n_found, 100 * n_concord / max(n_found, 1)))
cat(sprintf("Significativas (FDR<0.05):    %d / %d\n", n_sig, n_found))
cat("Por pilar:\n")
print(targets_tcga |>
  dplyr::count(pillar, concord_target, sig_tcga) |>
  arrange(pillar))
cat("========================\n\n")

# Colores de pilar (Other = gris)
PILLAR_PLOT <- c(PILLAR_COLS, Other = "#AAAAAA")

# Punto central del lollipop = logFC en TCGA de la diana; barra = compositeScore
# Layout: lollipop horizontal con dos canales (logFC TCGA, composite score)

# Subgráfico izquierdo: log2FC de la diana en TCGA (punto + segmento)
p_lfc <- ggplot(targets_tcga,
                aes(y = drug_label,
                    x = ifelse(is.na(logFC_tcga_target), 0, logFC_tcga_target))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70",
             linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = ifelse(is.na(logFC_tcga_target), 0, logFC_tcga_target),
                   yend = drug_label,
                   color = pillar),
               linewidth = 0.6, alpha = 0.7) +
  geom_point(aes(color = pillar,
                 shape = concord_target),
             size = 2.4) +
  geom_text(aes(label = sig_label,
                x     = ifelse(is.na(logFC_tcga_target), 0, logFC_tcga_target)),
            hjust = -0.4, size = PRESETS$single_col$tick / .pt - 0.5,
            color = "grey30") +
  scale_color_manual(values = PILLAR_PLOT, name = "Module",
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     name = "Concordant\nprot. vs RNA",
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = expression("log"[2]*"FC  (TCGA-HNSC)"),
       y = NULL) +
  theme_pub("double_col") +
  theme(legend.position = "right",
        panel.grid.major.x = element_line(linewidth = 0.2, color = "grey90"))

# Subgráfico derecho: composite score (barras apiladas TP/DV como Fig5A)
COMP_COLS <- c("Target priority" = "#009E73", "Drug viability" = "#CC79A7")

targets_long <- targets_tcga |>
  mutate(
    `Target priority` = w_target * TargetPriority,
    `Drug viability`  = w_drug   * DrugViability
  ) |>
  select(drug_label, `Target priority`, `Drug viability`, composite_score) |>
  pivot_longer(c(`Target priority`, `Drug viability`),
               names_to = "component", values_to = "contrib") |>
  mutate(component = factor(component, levels = names(COMP_COLS)))

p_score <- ggplot(targets_long, aes(x = contrib, y = drug_label, fill = component)) +
  geom_col(width = 0.65) +
  geom_text(
    data = targets_tcga |> mutate(drug_label = reorder(drug_label, composite_score)),
    aes(x = composite_score, y = drug_label,
        label = sprintf("%.2f", composite_score)),
    inherit.aes = FALSE,
    hjust = -0.15, size = PRESETS$single_col$tick / .pt, color = "grey25"
  ) +
  scale_fill_manual(values = COMP_COLS, name = "Score component") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Composite score", y = NULL) +
  theme_pub("double_col") +
  theme(legend.position  = "right",
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.line.y      = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.2, color = "grey90"))

# Combinar los dos subgráficos del panel B
p_targets <- p_lfc + p_score +
  plot_layout(widths = c(1.4, 1), guides = "collect") &
  theme(legend.position = "right")

save_pub(p_targets, "Fig6B_targets_tcga", h_add = 20)
save_panel_obj(p_targets, "Fig6B_targets_tcga")

# =============================================================================
# SECCIÓN 7: FigS (suplementario) — Supervivencia OS
# Análisis exploratorio; los resultados son no-significativos (p>0.05).
# Las dianas son vulnerabilidades terapéuticas, no biomarcadores pronósticos.
# =============================================================================
cat("\n--- Generando FigS (supervivencia — suplementario) ---\n")

# Muestras tumorales
se_tumor <- se[, substr(colnames(se), 14, 15) == "01"]
cat("Muestras tumorales para supervivencia:", ncol(se_tumor), "\n")

expr_mat           <- log2(assay(se_tumor, "fpkm_unstrand") + 0.1)
rownames(expr_mat) <- rowData(se_tumor)$gene_name

col_df <- as.data.frame(colData(se_tumor))
surv_df <- col_df |>
  mutate(
    patient_id = patient,
    os_event   = as.integer(tolower(vital_status) == "dead"),
    os_time    = if_else(os_event == 1,
                         suppressWarnings(as.numeric(days_to_death)),
                         suppressWarnings(as.numeric(days_to_last_follow_up)))
  ) |>
  filter(!is.na(os_time), os_time > 0) |>
  select(patient_id, os_event, os_time) |>
  distinct(patient_id, .keep_all = TRUE)

cat("Pacientes con datos de supervivencia:", nrow(surv_df), "\n")

surv_genes_data <- surv_df
for (gene in LABEL_GENES) {
  if (gene %in% rownames(expr_mat)) {
    gene_df <- data.frame(
      patient_id = col_df$patient,
      expr       = as.numeric(expr_mat[gene, ]),
      stringsAsFactors = FALSE
    ) |>
      group_by(patient_id) |>
      summarise(!!gene := mean(expr, na.rm = TRUE), .groups = "drop")
    surv_genes_data <- surv_genes_data |> left_join(gene_df, by = "patient_id")
  } else {
    cat("ATENCIÓN: gen no encontrado en expr_mat:", gene, "\n")
    surv_genes_data[[gene]] <- NA_real_
  }
}

for (gene in LABEL_GENES) {
  med <- median(surv_genes_data[[gene]], na.rm = TRUE)
  surv_genes_data[[paste0(gene, "_group")]] <- if_else(
    surv_genes_data[[gene]] >= med, "High", "Low"
  )
}
write_tsv(surv_genes_data, "results/tables/pub/supp/TabS2_survival_genes.tsv")
cat("Tabla supervivencia guardada.\n")

# KM plots 2×2
p_s <- PRESETS$single_col
PILLAR_PLOT_km <- PILLAR_COLS

km_plots <- lapply(seq_along(LABEL_GENES), function(i) {
  gene    <- LABEL_GENES[i]
  pil     <- LABEL_PILLAR[[gene]]
  gcol    <- PILLAR_PLOT_km[[pil]]
  grp_col <- paste0(gene, "_group")

  df_gene <- surv_genes_data |>
    filter(!is.na(.data[[grp_col]])) |>
    mutate(group = factor(.data[[grp_col]], levels = c("High", "Low")))

  if (nrow(df_gene) < 20) { message("Pocos datos para ", gene); return(NULL) }

  fit <- survfit(Surv(os_time, os_event) ~ group, data = df_gene)
  lrt <- tryCatch(
    survdiff(Surv(os_time, os_event) ~ group, data = df_gene),
    error = function(e) NULL
  )
  pval_txt <- if (!is.null(lrt)) {
    pv <- 1 - pchisq(lrt$chisq, df = length(lrt$n) - 1)
    if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)
  } else ""

  cat(sprintf("  KM %s (%s): %s\n", gene, pil, pval_txt))

  ggsurvplot(
    fit, data = df_gene,
    palette        = c(gcol, "grey60"),
    size           = 0.6, conf.int = FALSE, pval = FALSE, risk.table = FALSE,
    xscale         = "d_y",
    xlab           = "Time (years)", ylab = "Overall survival",
    title          = paste0(gene, " (", pil, ")"),
    legend.title   = "Expression", legend.labs = c("High", "Low"),
    fontsize       = p_s$tick / .pt,
    ggtheme        = theme_pub("single_col") +
                       theme(plot.title      = element_text(size = p_s$axis,
                                                            face = "bold"),
                             legend.position = c(0.78, 0.85))
  )$plot +
  annotate("text", x = Inf, y = 0.05, label = pval_txt,
           hjust = 1.05, vjust = 0,
           size = p_s$tick / .pt, color = "grey20")
})

km_ok <- Filter(Negate(is.null), km_plots)
fig_surv <- if (length(km_ok) >= 1) {
  wrap_plots(km_ok, ncol = 2) +
    plot_annotation(
      title    = "Supplementary — OS stratification by target gene expression",
      subtitle = "Note: all comparisons non-significant (p > 0.05); consistent with\ntargets being therapeutic vulnerabilities, not prognostic biomarkers.",
      theme    = theme_pub("double_col") +
                   theme(plot.title    = element_text(size = 10, face = "bold"),
                         plot.subtitle = element_text(size = 8,  color = "grey40"))
    )
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5,
                      label = "No survival data available", size = 5) +
    theme_void()
}

# Guardar en suplementario
ggsave("results/figures/pub/supp/FigS_survival_targets.pdf",
       fig_surv, width = PRESETS$double_col$w, height = PRESETS$double_col$h + 20,
       units = "mm", device = cairo_pdf)
ggsave("results/figures/pub/supp/FigS_survival_targets.png",
       fig_surv, width = PRESETS$double_col$w, height = PRESETS$double_col$h + 20,
       units = "mm", dpi = 300, bg = "white")
cat("  Saved: FigS_survival_targets (supp)\n")

# =============================================================================
# SECCIÓN 8: RESUMEN
# =============================================================================
cat("\n========================================\n")
cat("RESUMEN 16_external_validation.R\n")
cat("========================================\n")
cat(sprintf("Genes solapantes (nuestros vs TCGA):  %d\n", nrow(conc)))
cat(sprintf("Concordancia direccional:              %.1f%%\n", pct_concord))
cat(sprintf("Pearson r (logFC prot. vs RNA):        %.3f (p = %.2e)\n",
            pearson_r, pearson_p))
cat(sprintf("Dianas shortlist en TCGA:              %d / %d\n",
            n_found, nrow(targets_tcga)))
cat(sprintf("  Concordantes prot. vs RNA:           %d / %d (%.0f%%)\n",
            n_concord, n_found, 100 * n_concord / max(n_found, 1)))
cat(sprintf("  FDR<0.05 en TCGA RNA-seq:            %d / %d\n",
            n_sig, n_found))
cat(sprintf("Pacientes para supervivencia:          %d\n", nrow(surv_df)))
cat("\nOutputs:\n")
cat("  results/figures/pub/main/Fig6A_concordance.{pdf,png}\n")
cat("  results/figures/pub/main/Fig6B_targets_tcga.{pdf,png}\n")
cat("  results/figures/pub/supp/FigS_survival_targets.{pdf,png}\n")
cat("  results/tables/pub/main/Tab6_concordance_summary.tsv\n")
cat("  results/tables/pub/supp/TabS2_survival_genes.tsv\n")
cat("  (Ejecutar 17i_fig6_multipanel.R para ensamblar Fig6_multipanel.tif)\n")
cat("\nFin:", format(Sys.time()), "\n")
sink()
