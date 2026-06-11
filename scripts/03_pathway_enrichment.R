#!/usr/bin/env Rscript
# =============================================================================
# Script 03: Pathway Enrichment Analysis
# HNSCC Drug Repurposing — Fase 2
# =============================================================================
# Input:  results/tables/de_limma/02_TVsS_significant_with_ids.tsv
#         results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv
#         config/analysis_params.yaml
# Output: results/tables/pathway_enrichment/03_*.tsv
#         results/figures/pathway_enrichment/03_*.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(msigdbr)
  library(ggplot2)
  library(yaml)
  library(dplyr)
})

# --- Log setup ----------------------------------------------------------------
log_file <- file.path("logs", paste0("03_pathway_enrichment_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
dir.create("logs", showWarnings = FALSE)
con_log <- file(log_file, open = "wt")
sink(con_log, type = "output")
sink(con_log, type = "message", append = TRUE)

cat("=== 03_pathway_enrichment.R ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# --- Paths -------------------------------------------------------------------
dir_tables  <- "results/tables/pathway_enrichment"
dir_figures <- "results/figures/pathway_enrichment"
dir.create(dir_tables,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

# --- Parametros desde config -------------------------------------------------
params <- yaml::read_yaml("config/analysis_params.yaml")
pval_cutoff    <- params$pathway$pval_cutoff      # 0.05
min_gs_size    <- params$pathway$min_gs_size      # 10
max_gs_size    <- params$pathway$max_gs_size      # 500

cat(sprintf("Parametros: pval=%.2f | min_gs=%d | max_gs=%d\n\n",
            pval_cutoff, min_gs_size, max_gs_size))

# --- Cargar datos ------------------------------------------------------------
cat("--- Cargando datos ---\n")
sig  <- read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
                   stringsAsFactors = FALSE)
all_prot <- read.delim("results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv",
                       stringsAsFactors = FALSE)

cat(sprintf("Proteinas significativas: %d\n", nrow(sig)))
cat(sprintf("Proteinas totales:        %d\n", nrow(all_prot)))

# GSEA ranked list: pi-statistic (magnitud × significancia)
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

# Semilla fija para reproducibilidad de GSEA (fgsea usa permutaciones aleatorias).
# Se reaplica con set.seed() antes de cada llamada GSEA + seed=TRUE.
GSEA_SEED <- 42

# =============================================================================
# BLOQUE 4: GSEA — Hallmarks MSigDB
# =============================================================================
cat("\n=== GSEA: Hallmarks MSigDB ===\n")

hallmarks_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

cat(sprintf("  Gene sets Hallmarks: %d\n", length(unique(hallmarks_sets$gs_name))))

gsea_hallmarks <- tryCatch({
  set.seed(GSEA_SEED)
  res <- GSEA(
    geneList     = ranked_list,
    TERM2GENE    = hallmarks_sets,
    pvalueCutoff = pval_cutoff,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    verbose      = FALSE,
    eps          = 0,
    seed         = TRUE
  )
  cat(sprintf("  Hallmarks GSEA: %d gene sets significativos\n", nrow(as.data.frame(res))))
  res
}, error = function(e) {
  cat(sprintf("  Hallmarks GSEA ERROR: %s\n", e$message))
  NULL
})

# =============================================================================
# BLOQUE 5: GSEA — GO Biological Process
# =============================================================================
cat("\n=== GSEA: GO Biological Process ===\n")
gsea_go_bp <- tryCatch({
  set.seed(GSEA_SEED)
  res <- gseGO(
    geneList     = ranked_list,
    OrgDb        = org.Hs.eg.db,
    keyType      = "ENTREZID",
    ont          = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    verbose      = FALSE,
    eps          = 0,
    seed         = TRUE
  )
  cat(sprintf("  GSEA GO BP: %d terminos\n", nrow(as.data.frame(res))))
  res
}, error = function(e) {
  cat(sprintf("  GSEA GO BP ERROR: %s\n", e$message))
  NULL
})

# =============================================================================
# BLOQUE 5b: GSEA — KEGG
# =============================================================================
cat("\n=== GSEA: KEGG ===\n")
gsea_kegg <- tryCatch({
  set.seed(GSEA_SEED)
  res <- gseKEGG(
    geneList     = ranked_list,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    verbose      = FALSE,
    eps          = 0,
    seed         = TRUE
  )
  res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  cat(sprintf("  GSEA KEGG: %d rutas\n", nrow(as.data.frame(res))))
  res
}, error = function(e) {
  cat(sprintf("  GSEA KEGG ERROR: %s\n", e$message))
  NULL
})

# =============================================================================
# BLOQUE 5c: GSEA — Reactome
# =============================================================================
cat("\n=== GSEA: Reactome ===\n")
gsea_reactome <- tryCatch({
  set.seed(GSEA_SEED)
  res <- gsePathway(
    geneList     = ranked_list,
    organism     = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    verbose      = FALSE,
    eps          = 0,
    seed         = TRUE
  )
  res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  cat(sprintf("  GSEA Reactome: %d rutas\n", nrow(as.data.frame(res))))
  res
}, error = function(e) {
  cat(sprintf("  GSEA Reactome ERROR: %s\n", e$message))
  NULL
})

# =============================================================================
# BLOQUE 6: FIGURAS
# =============================================================================
cat("\n=== Generando figuras ===\n")

# Helper: dotplot GSEA balanceado (top n_each por p.adjust por dirección NES)
# Garantiza representación simétrica de términos activados y suprimidos.
gsea_balanced_dotplot <- function(gsea_res, n_each = 10, title = "",
                                   label_wrap = NULL, base_size = 18) {
  df <- as.data.frame(gsea_res)
  top_act  <- df %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice_head(n = n_each)
  top_supp <- df %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice_head(n = n_each)
  plot_df  <- bind_rows(top_act, top_supp) %>% arrange(NES)
  if (!is.null(label_wrap))
    plot_df$Description <- stringr::str_wrap(plot_df$Description, width = label_wrap)
  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)
  ggplot(plot_df, aes(x = NES, y = Description,
                      size = setSize, color = p.adjust)) +
    geom_point() +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    scale_color_gradient(low = "#D55E00", high = "#56B4E9", name = "p.adjust",
                         trans = "log10",
                         guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(name = "Set size", range = c(2, 8)) +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw(base_size = base_size)
}

safe_pdf_plot <- function(file, width, height, expr) {
  tryCatch({
    pdf(file, width = width, height = height)
    force(expr)
    dev.off()
    cat(sprintf("  Guardado: %s\n", basename(file)))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("  FIGURA OMITIDA (%s): %s\n", basename(file), e$message))
  })
}

# --- Hallmarks GSEA dotplot ---
if (!is.null(gsea_hallmarks) && nrow(as.data.frame(gsea_hallmarks)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_Hallmarks_GSEA_dotplot.pdf"), 12, 10, {
    p <- gsea_balanced_dotplot(gsea_hallmarks, n_each = 10,
                               title = "MSigDB Hallmarks GSEA — TVsS", base_size = 18)
    print(p)
  })
  # Ridge plot de distribuciones
  safe_pdf_plot(file.path(dir_figures, "03_Hallmarks_GSEA_ridgeplot.pdf"), 11, 12, {
    p <- ridgeplot(gsea_hallmarks, showCategory = 20, fill = "p.adjust") +
      labs(title = "Hallmarks GSEA — distribución de ranks") +
      theme_bw(base_size = 18)
    print(p)
  })
}

# --- GSEA GO BP dotplot ---
if (!is.null(gsea_go_bp) && nrow(as.data.frame(gsea_go_bp)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_GO_BP_GSEA_dotplot.pdf"), 12, 14, {
    p <- gsea_balanced_dotplot(gsea_go_bp, n_each = 10, label_wrap = 35,
                               title = "GO BP GSEA — TVsS", base_size = 18)
    print(p)
  })
}

# --- GSEA KEGG dotplot ---
if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_KEGG_GSEA_dotplot.pdf"), 12, 10, {
    p <- gsea_balanced_dotplot(gsea_kegg, n_each = 10,
                               title = "KEGG GSEA — TVsS", base_size = 18)
    print(p)
  })
}

# --- GSEA Reactome dotplot ---
if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_Reactome_GSEA_dotplot.pdf"), 12, 10, {
    p <- gsea_balanced_dotplot(gsea_reactome, n_each = 10,
                               title = "Reactome GSEA — TVsS", base_size = 18)
    print(p)
  })
}

# =============================================================================
# BLOQUE 7: EXPORTAR TABLAS
# =============================================================================
cat("\n=== Exportando tablas ===\n")

save_enrichment_tsv <- function(obj, filename) {
  if (!is.null(obj)) {
    df <- as.data.frame(obj)
    if (nrow(df) > 0) {
      write.table(df, file.path(dir_tables, filename),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat(sprintf("  %s (%d filas)\n", filename, nrow(df)))
    } else {
      cat(sprintf("  %s — sin resultados significativos\n", filename))
    }
  } else {
    cat(sprintf("  %s — objeto NULL (error en analisis)\n", filename))
  }
}

save_enrichment_tsv(gsea_hallmarks, "03_Hallmarks_GSEA.tsv")
save_enrichment_tsv(gsea_go_bp,     "03_GO_BP_GSEA.tsv")
save_enrichment_tsv(gsea_kegg,      "03_KEGG_GSEA.tsv")
save_enrichment_tsv(gsea_reactome,  "03_Reactome_GSEA.tsv")

# Tabla resumen: top 10 de cada analisis para vista rapida
make_top10 <- function(obj, source_name) {
  if (!is.null(obj)) {
    df <- as.data.frame(obj)
    if (nrow(df) > 0) {
      df %>%
        slice_head(n = 10) %>%
        mutate(source = source_name) %>%
        select(source, ID, Description, pvalue, p.adjust, everything())
    }
  }
}

top_list <- list(
  make_top10(gsea_hallmarks, "Hallmarks_GSEA"),
  make_top10(gsea_go_bp,     "GO_BP_GSEA"),
  make_top10(gsea_kegg,      "KEGG_GSEA"),
  make_top10(gsea_reactome,  "Reactome_GSEA")
)
top_list <- Filter(Negate(is.null), top_list)

if (length(top_list) > 0) {
  # Columnas comunes
  common_cols <- Reduce(intersect, lapply(top_list, colnames))
  top_summary <- bind_rows(lapply(top_list, function(x) x[, common_cols, drop = FALSE]))
  write.table(top_summary, file.path(dir_tables, "03_summary_top10_all.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  03_summary_top10_all.tsv (%d filas)\n", nrow(top_summary)))
}

# =============================================================================
# RESUMEN FINAL
# =============================================================================
cat("\n=== RESUMEN ===\n")

n_results <- function(obj) if (!is.null(obj)) nrow(as.data.frame(obj)) else 0

cat(sprintf("  Hallmarks GSEA:         %d gene sets\n",  n_results(gsea_hallmarks)))
cat(sprintf("  GO BP GSEA:             %d terminos\n",   n_results(gsea_go_bp)))
cat(sprintf("  KEGG GSEA:              %d rutas\n",      n_results(gsea_kegg)))
cat(sprintf("  Reactome GSEA:          %d rutas\n",      n_results(gsea_reactome)))

cat(sprintf("\nFiguras en:  %s\n", dir_figures))
cat(sprintf("Tablas en:   %s\n", dir_tables))
cat(sprintf("\nSiguiente:   scripts/04_query_dgidb.py\n"))
cat("\nFin:", format(Sys.time()), "\n")

sink(type = "message")
sink()
close(con_log)
cat("Script completado. Log en:", log_file, "\n")
