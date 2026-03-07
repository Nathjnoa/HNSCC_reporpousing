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
qval_cutoff    <- params$pathway$qval_cutoff      # 0.20
min_gs_size    <- params$pathway$min_gs_size      # 10
max_gs_size    <- params$pathway$max_gs_size      # 500
simplify_co    <- params$pathway$simplify_cutoff  # 0.70

cat(sprintf("Parametros: pval=%.2f | qval=%.2f | min_gs=%d | max_gs=%d | simplify=%.2f\n\n",
            pval_cutoff, qval_cutoff, min_gs_size, max_gs_size, simplify_co))

# --- Cargar datos ------------------------------------------------------------
cat("--- Cargando datos ---\n")
sig  <- read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
                   stringsAsFactors = FALSE)
all_prot <- read.delim("results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv",
                       stringsAsFactors = FALSE)

cat(sprintf("Proteinas significativas: %d\n", nrow(sig)))
cat(sprintf("Proteinas totales:        %d\n", nrow(all_prot)))

# Entrez IDs para ORA (solo con ID valido, valor exactamente entero)
sig_entrez <- sig %>%
  filter(!is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  as.character() %>%
  unique()

# Universo de background = todas las proteinas cuantificadas con Entrez
universe_entrez <- all_prot %>%
  filter(!is.na(entrez_id)) %>%
  pull(entrez_id) %>%
  as.character() %>%
  unique()

cat(sprintf("Entrez IDs para ORA:      %d / %d con universo=%d\n\n",
            length(sig_entrez), nrow(sig), length(universe_entrez)))

# Lista rankeada para GSEA (logFC, ordenada descendente, sin NA)
ranked_df <- all_prot %>%
  filter(!is.na(entrez_id), !is.na(logFC_TVsS)) %>%
  arrange(desc(logFC_TVsS))
ranked_list <- ranked_df$logFC_TVsS
names(ranked_list) <- as.character(ranked_df$entrez_id)
# En caso de duplicados, conservar el primero (mayor |logFC|)
ranked_list <- ranked_list[!duplicated(names(ranked_list))]
cat(sprintf("Ranked list para GSEA: %d genes (rango logFC: %.2f a %.2f)\n\n",
            length(ranked_list), min(ranked_list), max(ranked_list)))

# =============================================================================
# BLOQUE 1: ORA — Gene Ontology
# =============================================================================
cat("=== ORA: Gene Ontology ===\n")

run_go_ora <- function(ont) {
  cat(sprintf("  GO %s ... ", ont))
  tryCatch({
    res <- enrichGO(
      gene          = sig_entrez,
      universe      = universe_entrez,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = pval_cutoff,
      qvalueCutoff  = qval_cutoff,
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      readable      = TRUE
    )
    n <- nrow(as.data.frame(res))
    cat(sprintf("%d terminos\n", n))
    res
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    NULL
  })
}

go_bp <- run_go_ora("BP")
go_mf <- run_go_ora("MF")
go_cc <- run_go_ora("CC")

# Simplificar GO BP (elimina redundancia semantica)
go_bp_simple <- NULL
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  go_bp_simple <- simplify(go_bp, cutoff = simplify_co, by = "p.adjust", select_fun = min)
  cat(sprintf("  GO BP simplificado: %d terminos\n", nrow(as.data.frame(go_bp_simple))))
}

# =============================================================================
# BLOQUE 2: ORA — KEGG
# =============================================================================
cat("\n=== ORA: KEGG ===\n")
kegg_ora <- tryCatch({
  res <- enrichKEGG(
    gene          = sig_entrez,
    universe      = universe_entrez,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = pval_cutoff,
    qvalueCutoff  = qval_cutoff,
    minGSSize     = min_gs_size,
    maxGSSize     = max_gs_size
  )
  res <- setReadable(res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  cat(sprintf("  KEGG: %d rutas\n", nrow(as.data.frame(res))))
  res
}, error = function(e) {
  cat(sprintf("  KEGG ERROR: %s\n", e$message))
  NULL
})

# =============================================================================
# BLOQUE 3: ORA — Reactome
# =============================================================================
cat("\n=== ORA: Reactome ===\n")
reactome_ora <- tryCatch({
  res <- enrichPathway(
    gene          = sig_entrez,
    universe      = universe_entrez,
    organism      = "human",
    pAdjustMethod = "BH",
    pvalueCutoff  = pval_cutoff,
    qvalueCutoff  = qval_cutoff,
    minGSSize     = min_gs_size,
    maxGSSize     = max_gs_size,
    readable      = TRUE
  )
  cat(sprintf("  Reactome: %d rutas\n", nrow(as.data.frame(res))))
  res
}, error = function(e) {
  cat(sprintf("  Reactome ERROR: %s\n", e$message))
  NULL
})

# =============================================================================
# BLOQUE 4: GSEA — Hallmarks MSigDB
# =============================================================================
cat("\n=== GSEA: Hallmarks MSigDB ===\n")

hallmarks_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene) %>%
  mutate(entrez_gene = as.character(entrez_gene))

cat(sprintf("  Gene sets Hallmarks: %d\n", length(unique(hallmarks_sets$gs_name))))

gsea_hallmarks <- tryCatch({
  res <- GSEA(
    geneList     = ranked_list,
    TERM2GENE    = hallmarks_sets,
    pvalueCutoff = pval_cutoff,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    verbose      = FALSE,
    eps          = 0
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
    eps          = 0
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
  res <- gseKEGG(
    geneList     = ranked_list,
    organism     = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    verbose      = FALSE,
    eps          = 0
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
  res <- gsePathway(
    geneList     = ranked_list,
    organism     = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    minGSSize    = min_gs_size,
    maxGSSize    = max_gs_size,
    verbose      = FALSE,
    eps          = 0
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

# --- GO BP dotplot ---
if (!is.null(go_bp_simple) && nrow(as.data.frame(go_bp_simple)) >= 5) {
  safe_pdf_plot(file.path(dir_figures, "03_GO_BP_dotplot.pdf"), 10, 8, {
    p <- dotplot(go_bp_simple, showCategory = 20, font.size = 9) +
      labs(title = "GO Biological Process ORA — TVsS",
           subtitle = sprintf("Simplificado (cutoff=%.2f) | %d términos", simplify_co,
                              nrow(as.data.frame(go_bp_simple)))) +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- GO MF dotplot ---
if (!is.null(go_mf) && nrow(as.data.frame(go_mf)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_GO_MF_dotplot.pdf"), 10, 7, {
    p <- dotplot(go_mf, showCategory = 15, font.size = 9) +
      labs(title = "GO Molecular Function ORA — TVsS") +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- GO CC dotplot ---
if (!is.null(go_cc) && nrow(as.data.frame(go_cc)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_GO_CC_dotplot.pdf"), 10, 7, {
    p <- dotplot(go_cc, showCategory = 15, font.size = 9) +
      labs(title = "GO Cellular Component ORA — TVsS") +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- KEGG barplot ---
if (!is.null(kegg_ora) && nrow(as.data.frame(kegg_ora)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_KEGG_barplot.pdf"), 10, 7, {
    p <- barplot(kegg_ora, showCategory = 15, font.size = 9) +
      labs(title = "KEGG Pathway ORA — TVsS") +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- Reactome dotplot ---
if (!is.null(reactome_ora) && nrow(as.data.frame(reactome_ora)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_Reactome_dotplot.pdf"), 10, 7, {
    p <- dotplot(reactome_ora, showCategory = 15, font.size = 9) +
      labs(title = "Reactome Pathway ORA — TVsS") +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- GO BP emapplot (red de terminos) ---
# Requiere >= 10 terminos para que la red sea interpretable
n_bp_simple <- if (!is.null(go_bp_simple)) nrow(as.data.frame(go_bp_simple)) else 0
if (n_bp_simple >= 10) {
  safe_pdf_plot(file.path(dir_figures, "03_GO_BP_emapplot.pdf"), 11, 10, {
    go_bp_sim2 <- pairwise_termsim(go_bp_simple)
    p <- emapplot(go_bp_sim2, showCategory = min(30, n_bp_simple)) +
      labs(title = "GO BP — Red de similitud semántica",
           subtitle = sprintf("%d términos (simplificado)", n_bp_simple)) +
      theme(plot.title = element_text(size = 10))
    print(p)
  })
} else {
  cat(sprintf("  emapplot omitido: solo %d terminos (requiere >= 10)\n", n_bp_simple))
}

# --- GO BP cnetplot (genes por termino) ---
if (!is.null(go_bp_simple) && nrow(as.data.frame(go_bp_simple)) >= 3) {
  # Pasar fold change para colorear los nodos gen
  fc_vec <- sig$logFC_TVsS
  names(fc_vec) <- sig$symbol_org
  safe_pdf_plot(file.path(dir_figures, "03_GO_BP_cnetplot.pdf"), 13, 11, {
    p <- cnetplot(go_bp_simple, showCategory = 8,
                  foldChange = fc_vec,
                  cex.params = list(category_label = 0.8, gene_label = 0.6)) +
      scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                            midpoint = 0, name = "log2FC") +
      labs(title = "GO BP — Red gen-concepto (top 8 términos)") +
      theme(plot.title = element_text(size = 10))
    print(p)
  })
}

# --- Hallmarks GSEA dotplot ---
if (!is.null(gsea_hallmarks) && nrow(as.data.frame(gsea_hallmarks)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_Hallmarks_GSEA_dotplot.pdf"), 10, 8, {
    p <- dotplot(gsea_hallmarks, showCategory = 20, split = ".sign",
                 font.size = 9) +
      facet_grid(. ~ .sign) +
      labs(title = "MSigDB Hallmarks GSEA — TVsS") +
      theme_bw(base_size = 9)
    print(p)
  })
  # Ridge plot de distribuciones
  safe_pdf_plot(file.path(dir_figures, "03_Hallmarks_GSEA_ridgeplot.pdf"), 9, 10, {
    p <- ridgeplot(gsea_hallmarks, showCategory = 20, fill = "p.adjust") +
      labs(title = "Hallmarks GSEA — distribución de ranks") +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- GSEA GO BP dotplot ---
if (!is.null(gsea_go_bp) && nrow(as.data.frame(gsea_go_bp)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_GO_BP_GSEA_dotplot.pdf"), 10, 8, {
    p <- dotplot(gsea_go_bp, showCategory = 20, split = ".sign",
                 font.size = 9) +
      facet_grid(. ~ .sign) +
      labs(title = "GO BP GSEA — TVsS") +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- GSEA KEGG dotplot ---
if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_KEGG_GSEA_dotplot.pdf"), 10, 8, {
    p <- dotplot(gsea_kegg, showCategory = 20, split = ".sign",
                 font.size = 9) +
      facet_grid(. ~ .sign) +
      labs(title = "KEGG GSEA — TVsS") +
      theme_bw(base_size = 9)
    print(p)
  })
}

# --- GSEA Reactome dotplot ---
if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) >= 3) {
  safe_pdf_plot(file.path(dir_figures, "03_Reactome_GSEA_dotplot.pdf"), 10, 8, {
    p <- dotplot(gsea_reactome, showCategory = 20, split = ".sign",
                 font.size = 9) +
      facet_grid(. ~ .sign) +
      labs(title = "Reactome GSEA — TVsS") +
      theme_bw(base_size = 9)
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

save_enrichment_tsv(go_bp,          "03_GO_BP_ORA_full.tsv")
save_enrichment_tsv(go_bp_simple,   "03_GO_BP_ORA_simplified.tsv")
save_enrichment_tsv(go_mf,          "03_GO_MF_ORA.tsv")
save_enrichment_tsv(go_cc,          "03_GO_CC_ORA.tsv")
save_enrichment_tsv(kegg_ora,       "03_KEGG_ORA.tsv")
save_enrichment_tsv(reactome_ora,   "03_Reactome_ORA.tsv")
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
  make_top10(go_bp_simple,   "GO_BP_ORA"),
  make_top10(go_mf,          "GO_MF_ORA"),
  make_top10(go_cc,          "GO_CC_ORA"),
  make_top10(kegg_ora,       "KEGG_ORA"),
  make_top10(reactome_ora,   "Reactome_ORA"),
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

cat(sprintf("  GO BP ORA (full):       %d terminos\n",   n_results(go_bp)))
cat(sprintf("  GO BP ORA (simplif.):   %d terminos\n",   n_results(go_bp_simple)))
cat(sprintf("  GO MF ORA:              %d terminos\n",   n_results(go_mf)))
cat(sprintf("  GO CC ORA:              %d terminos\n",   n_results(go_cc)))
cat(sprintf("  KEGG ORA:               %d rutas\n",      n_results(kegg_ora)))
cat(sprintf("  Reactome ORA:           %d rutas\n",      n_results(reactome_ora)))
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
