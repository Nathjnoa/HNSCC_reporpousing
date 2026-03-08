#!/usr/bin/env Rscript
# =============================================================================
# Script 09: Red PPI con STRING — via REST API + igraph
# HNSCC Drug Repurposing — Fase 3
# =============================================================================
# Objetivo: Construir red PPI de genes DE usando STRING API v12. Calcular
#           metricas de centralidad e identificar hubs druggables.
#
# Input:  results/tables/de_limma/02_TVsS_significant_with_ids.tsv
#         results/tables/drug_targets/08_multi_source_candidates.tsv
#         config/analysis_params.yaml
#
# Output: results/tables/network/09_network_node_metrics.tsv
#         results/tables/network/09_network_edges.tsv
#         results/tables/network/09_druggable_hubs.tsv
#         results/figures/09_network_degree_dist.pdf
#         results/figures/09_top25_hubs_degree.pdf
#         results/figures/09_hub_druggable_scatter.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(igraph)
  library(ggraph)
  library(dplyr)
  library(ggplot2)
  library(yaml)
})

# --- Establecer working directory (desde ruta del script) -------------------
args <- commandArgs(trailingOnly = FALSE)
script_flag <- args[grep("^--file=", args)]
if (length(script_flag) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_flag))
  proj_dir    <- dirname(dirname(script_path))   # scripts/ -> proyecto/
  if (file.exists(file.path(proj_dir, "config/analysis_params.yaml"))) {
    setwd(proj_dir)
  }
}
cat("Working dir:", getwd(), "\n")

# --- Log setup ---------------------------------------------------------------
dir.create("logs",                   showWarnings = FALSE)
dir.create("results/tables/network", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures",        showWarnings = FALSE, recursive = TRUE)

log_file <- paste0("logs/09_string_network_",
                   format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
con_log <- file(log_file, open = "wt")
sink(con_log, type = "output")
sink(con_log, type = "message", append = TRUE)

cat("=== 09_string_network.R ===\n")
cat("Inicio:", format(Sys.time()), "\n\n")

# --- Parametros --------------------------------------------------------------
params         <- yaml::read_yaml("config/analysis_params.yaml")
score_thr      <- params$network$string_score_threshold   # 700
hub_percentile <- params$network$hub_percentile           # 90
STRING_BASE    <- "https://string-db.org/api/json"
SPECIES        <- 9606

cat(sprintf("STRING score threshold:  %d\n", score_thr))
cat(sprintf("Hub percentile: %d%% (top %d%%)\n", hub_percentile, 100 - hub_percentile))

# =============================================================================
# HELPER: STRING REST API
# STRING devuelve Content-Type: text/json — parsear con jsonlite
# =============================================================================
string_post <- function(endpoint, body_list, retries = 2, wait = 5) {
  url <- paste0(STRING_BASE, "/", endpoint)
  for (attempt in seq_len(retries + 1)) {
    tryCatch({
      resp <- request(url) %>%
        req_body_form(!!!body_list) %>%
        req_timeout(120) %>%
        req_perform()
      raw <- resp_body_string(resp)
      return(fromJSON(raw))
    }, error = function(e) {
      if (attempt <= retries) {
        cat(sprintf("  Reintento %d/%d: %s\n", attempt, retries, e$message))
        Sys.sleep(wait)
      } else {
        cat(sprintf("  ERROR en %s: %s\n", endpoint, e$message))
      }
    })
  }
  NULL
}

# --- Cargar genes DE ---------------------------------------------------------
cat("\n--- Cargando datos ---\n")
sig <- read.delim("results/tables/de_limma/02_TVsS_significant_with_ids.tsv",
                  stringsAsFactors = FALSE)
gene_meta <- sig %>%
  select(symbol_org, uniprot_id, logFC_TVsS, adj.P.Val_TVsS) %>%
  mutate(direction = ifelse(logFC_TVsS > 0, "up", "down")) %>%
  distinct(symbol_org, .keep_all = TRUE)

gene_symbols <- unique(na.omit(sig$symbol_org))
cat(sprintf("Genes DE cargados: %d\n", length(gene_symbols)))

# Candidatos de farmaco (script 08)
drug_cands <- read.delim("results/tables/drug_targets/08_multi_source_candidates.tsv",
                         stringsAsFactors = FALSE)
cat(sprintf("Candidatos de farmaco: %d\n", nrow(drug_cands)))

# =============================================================================
# PASO 1: Mapear genes → STRING IDs
# =============================================================================
cat("\n--- Paso 1: Mapeando genes a STRING IDs ---\n")

# STRING acepta genes separados por \r (carriage return en el body)
map_result <- string_post("get_string_ids", list(
  identifiers     = paste(gene_symbols, collapse = "\r"),
  species         = as.character(SPECIES),
  limit           = "1",
  echo_query      = "1",
  caller_identity = "hnscc_drug_repurposing"
))

if (is.null(map_result) || nrow(as.data.frame(map_result)) == 0) {
  cat("ERROR: Mapeo de genes fallido\n")
  sink(type = "message"); sink(); close(con_log); quit(status = 1)
}

df_map <- as.data.frame(map_result)
cat(sprintf("Genes mapeados: %d / %d\n", nrow(df_map), length(gene_symbols)))
cat(sprintf("Columnas: %s\n", paste(colnames(df_map), collapse = ", ")))

string_ids  <- df_map$stringId
mapped_syms <- df_map$preferredName
cat(sprintf("Ejemplo: %s -> %s\n", mapped_syms[1], string_ids[1]))

# =============================================================================
# PASO 2: Obtener interacciones STRING
# =============================================================================
cat(sprintf("\n--- Paso 2: Obteniendo interacciones (score >= %d) ---\n", score_thr))
cat("(Puede tardar 1-2 minutos para 500+ proteinas...)\n")

net_result <- string_post("network", list(
  identifiers     = paste(string_ids, collapse = "\r"),
  species         = as.character(SPECIES),
  required_score  = as.character(score_thr),
  caller_identity = "hnscc_drug_repurposing"
))

if (is.null(net_result) || nrow(as.data.frame(net_result)) == 0) {
  cat("Sin interacciones con score >= ", score_thr, ". Reintentando con 400...\n")
  net_result <- string_post("network", list(
    identifiers    = paste(string_ids, collapse = "\r"),
    species        = as.character(SPECIES),
    required_score = "400",
    caller_identity = "hnscc_drug_repurposing"
  ))
  if (!is.null(net_result) && nrow(as.data.frame(net_result)) > 0) {
    score_thr <- 400L  # actualizar para que filtros posteriores sean consistentes
    cat("Fallback exitoso. score_thr actualizado a 400.\n")
  }
}

if (is.null(net_result)) {
  cat("ERROR: No se obtuvieron interacciones\n")
  sink(type = "message"); sink(); close(con_log); quit(status = 1)
}

df_edges_raw <- as.data.frame(net_result)
cat(sprintf("Interacciones obtenidas: %d\n", nrow(df_edges_raw)))
cat(sprintf("Columnas: %s\n", paste(colnames(df_edges_raw), collapse = ", ")))

# Detectar columnas (STRING puede variar entre versiones)
score_col <- intersect(c("score", "combined_score"), colnames(df_edges_raw))[1]
nameA_col <- intersect(c("preferredName_A", "preferredName.A"), colnames(df_edges_raw))[1]
nameB_col <- intersect(c("preferredName_B", "preferredName.B"), colnames(df_edges_raw))[1]

cat(sprintf("Score col: %s | Names: %s, %s\n", score_col, nameA_col, nameB_col))

# Score puede estar en [0,1] o [0,1000] segun version de API
score_vals <- as.numeric(df_edges_raw[[score_col]])
score_max  <- max(score_vals, na.rm = TRUE)
cat(sprintf("Score range: %.4f - %.4f\n", min(score_vals, na.rm=TRUE), score_max))

# Normalizar y filtrar
if (score_max > 1) {
  # Score en escala 0-1000 (API mas antigua)
  df_edges <- df_edges_raw %>%
    mutate(
      gene_A = .data[[nameA_col]],
      gene_B = .data[[nameB_col]],
      score  = as.numeric(.data[[score_col]])
    ) %>%
    filter(score >= score_thr) %>%
    select(gene_A, gene_B, score) %>%
    distinct()
} else {
  # Score en escala 0-1 (API v12)
  df_edges <- df_edges_raw %>%
    mutate(
      gene_A = .data[[nameA_col]],
      gene_B = .data[[nameB_col]],
      score  = as.numeric(.data[[score_col]])
    ) %>%
    filter(score >= score_thr / 1000) %>%
    select(gene_A, gene_B, score) %>%
    distinct()
}

cat(sprintf("Interacciones tras filtro: %d\n", nrow(df_edges)))

write.table(df_edges,
            "results/tables/network/09_network_edges.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# PASO 3: Construir grafo igraph
# =============================================================================
cat("\n--- Paso 3: Construyendo grafo ---\n")

g <- graph_from_data_frame(df_edges, directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

n_nodes <- vcount(g)
n_edges <- ecount(g)
cat(sprintf("Nodos: %d | Aristas: %d | Densidad: %.5f\n",
            n_nodes, n_edges, edge_density(g)))

comps <- components(g)
cat(sprintf("Componentes: %d | Mayor componente: %d nodos\n",
            comps$no, max(comps$csize)))

# =============================================================================
# PASO 4: Metricas de centralidad
# =============================================================================
cat("\n--- Paso 4: Metricas de centralidad ---\n")

deg   <- degree(g)
betw  <- betweenness(g, normalized = TRUE)
eig   <- tryCatch(
  eigen_centrality(g)$vector,
  error = function(e) { cat("  eigenvector centrality: NA (grafo desconectado)\n")
                        rep(NA_real_, n_nodes) }
)
clust <- transitivity(g, type = "local", isolates = "zero")

df_nodes <- data.frame(
  gene_symbol      = V(g)$name,
  degree           = deg,
  betweenness_norm = betw,
  eigenvector      = eig,
  clustering_coeff = clust,
  stringsAsFactors = FALSE
) %>%
  left_join(gene_meta, by = c("gene_symbol" = "symbol_org")) %>%
  arrange(desc(degree))

# Hub criteria: union de degree >= p90 OR betweenness >= p90
# Captura "party hubs" (grado alto) y "bottleneck hubs" (intermediacion alta)
hub_thr_deg  <- quantile(deg,  probs = hub_percentile / 100, na.rm = TRUE)
hub_thr_betw <- quantile(betw, probs = hub_percentile / 100, na.rm = TRUE)
hub_thr      <- hub_thr_deg  # referencia para plots de distribucion de grado

df_nodes <- df_nodes %>%
  mutate(
    # Score compuesto: media de rangos normalizados (robusto ante escalas distintas)
    hub_score = (rank(degree) + rank(betweenness_norm) +
                 rank(coalesce(eigenvector, 0))) / (3 * n()),
    is_hub    = degree >= hub_thr_deg | betweenness_norm >= hub_thr_betw
  )

n_hubs <- sum(df_nodes$is_hub)
cat(sprintf("Hub threshold: degree >= %.0f (p%d) OR betweenness >= %.4f (p%d)\n",
            hub_thr_deg, hub_percentile, hub_thr_betw, hub_percentile))
cat(sprintf("Hubs identificados: %d (union criterio dual)\n", n_hubs))

cat("\nTop 20 hubs:\n")
top20 <- head(df_nodes, 20)
for (i in seq_len(nrow(top20))) {
  r <- top20[i, ]
  cat(sprintf("  %2d. %-15s deg=%3d  betw=%.4f  %s  logFC=%s\n",
              i, r$gene_symbol, r$degree, r$betweenness_norm,
              coalesce(r$direction, "?"),
              ifelse(is.na(r$logFC_TVsS), "NA", sprintf("%.2f", r$logFC_TVsS))))
}

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
dir.create("results/tables/network", showWarnings = FALSE, recursive = TRUE)
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

# =============================================================================
# PASO 4c: Redefinir hubs a nivel de módulo
#           Hub = proteína en top 10% de betweenness DENTRO de su módulo
# =============================================================================
cat("\n--- Paso 4c: Hubs por módulo (betweenness intra-módulo) ---\n")

df_nodes <- df_nodes %>%
  group_by(module_id) %>%
  mutate(
    module_hub_score = betweenness_norm,
    is_module_hub    = betweenness_norm >= quantile(betweenness_norm,
                                                     0.90, na.rm = TRUE)
  ) %>%
  ungroup()

n_module_hubs <- sum(df_nodes$is_module_hub, na.rm = TRUE)
cat(sprintf("Module hubs (top 10%% betweenness por módulo): %d\n", n_module_hubs))

df_nodes$n_distinct_modules <- 1L

# =============================================================================
# PASO 5: Hubs druggables
# =============================================================================
cat("\n--- Paso 5: Cruzando hubs con candidatos de farmaco ---\n")

hubs <- df_nodes %>% filter(is_hub)

hub_drug_list <- lapply(hubs$gene_symbol, function(gene) {
  hits <- drug_cands %>%
    filter(grepl(paste0("(^|\\|)", gene, "(\\||$)"), de_genes))
  if (nrow(hits) > 0) {
    data.frame(
      gene_symbol  = gene,
      n_drugs      = nrow(hits),
      drug_classes = paste(sort(unique(hits$drug_class)), collapse = "|"),
      top_drugs    = paste(head(hits$drug_name_norm, 5), collapse = "|"),
      min_class    = min(hits$drug_class),
      stringsAsFactors = FALSE
    )
  } else NULL
})

df_hub_drugs <- bind_rows(hub_drug_list)

df_druggable <- hubs %>%
  left_join(df_hub_drugs, by = "gene_symbol") %>%
  mutate(is_druggable = !is.na(n_drugs)) %>%
  arrange(desc(is_druggable), desc(degree))

cat(sprintf("Hubs druggables: %d / %d\n",
            sum(df_druggable$is_druggable), n_hubs))

cat("\nHubs druggables (clase A o B preferentemente):\n")
df_druggable %>%
  filter(is_druggable) %>%
  arrange(min_class, desc(degree)) %>%
  { for (i in seq_len(nrow(.))) {
      r <- .[i, ]
      cat(sprintf("  %-15s deg=%3d  %d drugs  clases=%s  | %s\n",
                  r$gene_symbol, r$degree, r$n_drugs,
                  r$drug_classes,
                  substr(r$top_drugs, 1, 60)))
  }; invisible(.) }

# =============================================================================
# PASO 6: Figuras
# =============================================================================
cat("\n--- Paso 6: Figuras ---\n")

safe_pdf <- function(file, expr, w = 8, h = 6) {
  tryCatch({
    pdf(file, width = w, height = h)
    force(expr)
    dev.off()
    cat(sprintf("  %s\n", file))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat(sprintf("  WARN %s: %s\n", file, e$message))
  })
}

# 6a. Distribucion de grados
safe_pdf("results/figures/09_network_degree_dist.pdf", {
  p <- data.frame(degree = deg) %>%
    ggplot(aes(x = degree)) +
    geom_histogram(binwidth = 2, fill = "#2c7bb6", color = "white", alpha = 0.85) +
    geom_vline(xintercept = hub_thr, linetype = "dashed", color = "red", linewidth = 0.9) +
    annotate("text", x = hub_thr + 0.5, y = Inf,
             label = paste0("Hub\n(>=", hub_thr, ")"),
             hjust = 0, vjust = 1.5, color = "red", size = 3.5) +
    labs(title = "STRING PPI network — degree distribution",
         subtitle = sprintf("%d nodes, %d edges | score >= %d",
                            n_nodes, n_edges, score_thr),
         x = "Degree", y = "Count") +
    theme_bw(base_size = 13)
  print(p)
})

# 6b. Top 25 hubs barplot
top25_hubs <- df_nodes %>%
  filter(is_hub) %>%
  slice_head(n = 25) %>%
  mutate(direction   = coalesce(direction, "unknown"),
         gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)))

safe_pdf("results/figures/09_top25_hubs_degree.pdf", w = 9, h = 7, {
  p <- ggplot(top25_hubs, aes(x = gene_symbol, y = degree, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(up = "#e31a1c", down = "#1f78b4",
                                 unknown = "#aaaaaa")) +
    coord_flip() +
    labs(title = "Top 25 hub proteins — STRING PPI network",
         subtitle = sprintf("DE genes | score >= %d | color = expression vs Normal",
                            score_thr),
         x = NULL, y = "Degree", fill = "Direction") +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank())
  print(p)
})

# 6c. Hubs druggables scatter
if (nrow(df_druggable) >= 3) {
  safe_pdf("results/figures/09_hub_druggable_scatter.pdf", w = 9, h = 6, {
    p <- df_druggable %>%
      mutate(label = ifelse(is_druggable, gene_symbol, NA)) %>%
      ggplot(aes(x = degree, y = betweenness_norm,
                 color = is_druggable, size = is_druggable)) +
      geom_point(alpha = 0.8) +
      ggrepel::geom_text_repel(aes(label = label),
                               size = 3, max.overlaps = 25, na.rm = TRUE) +
      scale_color_manual(values = c("TRUE" = "#d62728", "FALSE" = "#aaaaaa"),
                         labels = c("TRUE" = "Has drug candidate", "FALSE" = "No drug")) +
      scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.5)) +
      labs(title = "Hub proteins — druggability",
           subtitle = "Red labels = hubs with drug candidate in integrated DB",
           x = "Degree", y = "Betweenness centrality (normalized)",
           color = NULL) +
      guides(size = "none") +
      theme_bw(base_size = 12)
    print(p)
  })
}

# 6d. -------- NETWORK PLOTS con ggraph ----------------------------------------
# Plot 1: Red completa del componente mayor — nodos coloreados por direction
#         tamaño = degree, hubs etiquetados
cat("  [ggraph] Red completa — componente mayor...\n")

# Trabajar solo sobre el componente mas grande
giant_nodes <- names(which(comps$membership == which.max(comps$csize)))
g_giant     <- induced_subgraph(g, vids = giant_nodes)

# Atributos de nodos para ggraph
node_attrs <- df_nodes %>%
  filter(gene_symbol %in% V(g_giant)$name) %>%
  select(gene_symbol, degree, betweenness_norm, direction,
         logFC_TVsS, is_hub, hub_score) %>%
  mutate(direction = coalesce(direction, "unknown"))

V(g_giant)$degree        <- node_attrs$degree[match(V(g_giant)$name, node_attrs$gene_symbol)]
V(g_giant)$betweenness   <- node_attrs$betweenness_norm[match(V(g_giant)$name, node_attrs$gene_symbol)]
V(g_giant)$direction     <- node_attrs$direction[match(V(g_giant)$name, node_attrs$gene_symbol)]
V(g_giant)$logFC         <- node_attrs$logFC_TVsS[match(V(g_giant)$name, node_attrs$gene_symbol)]
V(g_giant)$is_hub        <- node_attrs$is_hub[match(V(g_giant)$name, node_attrs$gene_symbol)]
V(g_giant)$hub_score     <- node_attrs$hub_score[match(V(g_giant)$name, node_attrs$gene_symbol)]
V(g_giant)$label         <- ifelse(V(g_giant)$is_hub == TRUE, V(g_giant)$name, NA_character_)

# Asignar score de aristas (match insensible al orden A-B / B-A)
edge_ends    <- ends(g_giant, E(g_giant), names = TRUE)
score_lookup <- setNames(df_edges$score,
                         paste(df_edges$gene_A, df_edges$gene_B))
score_lookup_rev <- setNames(df_edges$score,
                             paste(df_edges$gene_B, df_edges$gene_A))
e_keys       <- paste(edge_ends[,1], edge_ends[,2])
e_scores     <- score_lookup[e_keys]
na_mask      <- is.na(e_scores)
e_scores[na_mask] <- score_lookup_rev[e_keys[na_mask]]
E(g_giant)$score <- as.numeric(e_scores)

# Propagar score al subgrafo de hubs (se creara despues, pero lo necesitamos)
# -> se hara tras crear g_hubs (linea siguiente en el bloque de Plot 2)

# Red minimalista: nodos pequenos para proteinas comunes, grandes para hubs.
# Solo se renderizan edges; el alpha muy bajo evita el efecto "hairball"
# manteniendo la estructura topologica visible.
safe_pdf("results/figures/09_network_full.pdf", w = 14, h = 12, {
  set.seed(42)
  p <- ggraph(g_giant, layout = "fr") +
    # Aristas: alpha bajo — muestra estructura sin crear hairball
    geom_edge_link(color = "gray70", linewidth = 0.2, alpha = 0.06,
                   show.legend = FALSE) +
    # No-hubs: puntos pequenos y semitransparentes (contexto)
    geom_node_point(aes(size = degree, color = direction),
                    data = function(d) d[d$is_hub == FALSE, ],
                    alpha = 0.45, shape = 16) +
    # Hubs: diamantes grandes, completamente visibles (protagonistas)
    geom_node_point(aes(size = degree, color = direction),
                    data = function(d) d[d$is_hub == TRUE, ],
                    alpha = 1.0, shape = 18) +
    # Labels solo en hubs
    geom_node_text(aes(label = label),
                   size = 2.8, repel = TRUE,
                   max.overlaps = 60, na.rm = TRUE,
                   fontface = "bold", bg.color = "white", bg.r = 0.1) +
    scale_color_manual(values = c(up = "#d62728", down = "#1f78b4",
                                  unknown = "#999999"),
                       name = "Expression\nvs Normal") +
    scale_size_continuous(range = c(0.8, 8), name = "Degree") +
    labs(title = "STRING PPI network — DE proteins in HNSCC",
         subtitle = sprintf(
           "Giant component: %d nodes, %d edges | score \u2265 %d | diamonds = hubs (degree OR betweenness \u2265 p%d)",
           vcount(g_giant), ecount(g_giant), score_thr, hub_percentile)) +
    theme_graph(base_family = "sans") +
    theme(legend.position    = "right",
          plot.title         = element_text(size = 14, face = "bold"),
          plot.subtitle      = element_text(size = 9),
          plot.background    = element_rect(fill = "white", color = NA),
          panel.background   = element_rect(fill = "white", color = NA))
  print(p)
})

# 6e. -------- RED COMPLETA COLOREADA POR MÓDULO ----------------------------
cat("  [ggraph] Red con módulos Louvain coloreados...\n")

module_lookup <- setNames(df_modules$module_name, df_modules$gene_symbol)
V(g_giant)$module_id_v <- df_nodes$module_id[
  match(V(g_giant)$name, df_nodes$gene_symbol)
]
V(g_giant)$is_module_hub <- df_nodes$is_module_hub[
  match(V(g_giant)$name, df_nodes$gene_symbol)
]

n_mod_plot  <- min(n_modules, 12)
mod_palette <- setNames(
  c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
    "#a65628","#f781bf","#999999","#66c2a5","#fc8d62",
    "#8da0cb","#e78ac3")[seq_len(n_mod_plot)],
  as.character(seq_len(n_mod_plot))
)

safe_pdf("results/figures/09_network_modules.pdf", w = 15, h = 13, {
  set.seed(42)
  p_mod <- ggraph(g_giant, layout = "fr") +
    geom_edge_link(color = "gray70", linewidth = 0.15, alpha = 0.05,
                   show.legend = FALSE) +
    geom_node_point(
      aes(size  = degree,
          color = factor(V(g_giant)$module_id_v %% n_mod_plot + 1),
          shape = V(g_giant)$is_module_hub),
      alpha = 0.75
    ) +
    geom_node_text(
      aes(label = ifelse(V(g_giant)$is_module_hub == TRUE, name, NA_character_)),
      size = 2.6, repel = TRUE, max.overlaps = 50,
      fontface = "bold", na.rm = TRUE
    ) +
    scale_color_manual(values = mod_palette, name = "Módulo Louvain",
                       guide  = guide_legend(ncol = 2)) +
    scale_size_continuous(range = c(0.6, 7), name = "Degree") +
    scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16),
                       labels = c("TRUE" = "Hub (top 10% betweenness)",
                                  "FALSE" = "Proteína"),
                       name = "") +
    labs(
      title    = "STRING PPI network — DE proteins coloreadas por módulo Louvain",
      subtitle = sprintf(
        "Componente gigante: %d nodos, %d aristas | %d módulos | diamantes = hubs por módulo",
        vcount(g_giant), ecount(g_giant), n_modules)
    ) +
    theme_graph(base_family = "sans") +
    theme(legend.position   = "right",
          plot.title        = element_text(size = 13, face = "bold"),
          plot.subtitle     = element_text(size = 9),
          plot.background   = element_rect(fill = "white", color = NA),
          panel.background  = element_rect(fill = "white", color = NA))
  print(p_mod)
})

# Plot 2: Subred de hubs solamente — layout KK para grafos medianos (~50-150 nodos)
cat("  [ggraph] Subred de hubs...\n")
hub_names <- df_nodes$gene_symbol[df_nodes$is_hub]
g_hubs    <- induced_subgraph(g, vids = hub_names)

# Propagar score de aristas desde g al subgrafo
hub_edge_ends <- ends(g_hubs, E(g_hubs), names = TRUE)
hub_e_keys    <- paste(hub_edge_ends[,1], hub_edge_ends[,2])
hub_e_scores  <- score_lookup[hub_e_keys]
hub_na        <- is.na(hub_e_scores)
hub_e_scores[hub_na] <- score_lookup_rev[hub_e_keys[hub_na]]
E(g_hubs)$score <- as.numeric(hub_e_scores)

# Grado en el subgrafo de hubs (no el grado en la red completa)
deg_hubs <- degree(g_hubs)
V(g_hubs)$degree       <- deg_hubs
V(g_hubs)$degree_full  <- deg[match(V(g_hubs)$name, names(deg))]  # grado en la red completa
V(g_hubs)$direction    <- node_attrs$direction[match(V(g_hubs)$name, node_attrs$gene_symbol)]
V(g_hubs)$hub_score    <- node_attrs$hub_score[match(V(g_hubs)$name, node_attrs$gene_symbol)]
V(g_hubs)$is_druggable <- V(g_hubs)$name %in%
  df_druggable$gene_symbol[df_druggable$is_druggable]

# Si el subgrafo es muy grande (>150 nodos), usar FR; si no, KK (mas estetico)
hub_layout <- if (vcount(g_hubs) <= 150) "kk" else "fr"
cat(sprintf("    Hub subgraph: %d nodos, %d aristas | layout = %s\n",
            vcount(g_hubs), ecount(g_hubs), hub_layout))

# Normalizar scores de aristas a [0,1] para alpha; fallback a 0.5 si todos NA
hub_scores_raw <- E(g_hubs)$score
if (all(is.na(hub_scores_raw))) {
  E(g_hubs)$score_norm <- rep(0.5, ecount(g_hubs))
} else {
  s_min <- min(hub_scores_raw, na.rm = TRUE)
  s_max <- max(hub_scores_raw, na.rm = TRUE)
  E(g_hubs)$score_norm <- ifelse(
    s_max > s_min,
    (hub_scores_raw - s_min) / (s_max - s_min),
    0.5
  )
}

safe_pdf("results/figures/09_network_hubs_only.pdf", w = 11, h = 10, {
  set.seed(42)
  p <- ggraph(g_hubs, layout = hub_layout) +
    geom_edge_link(aes(edge_alpha = score_norm),
                   color = "gray40", linewidth = 0.6,
                   show.legend = FALSE) +
    scale_edge_alpha_continuous(range = c(0.2, 0.85)) +
    # Nodos — tamaño por hub_score compuesto (mejor ranking que grado simple)
    geom_node_point(aes(size  = hub_score,
                        color = direction,
                        shape = is_druggable),
                    alpha = 0.92) +
    # Labels en todos los hubs (grafo pequeno, todos son relevantes)
    geom_node_text(aes(label = name),
                   size = 3.2, repel = TRUE,
                   max.overlaps = 80, fontface = "bold",
                   bg.color = "white", bg.r = 0.12,
                   force = 1.5) +
    scale_color_manual(values = c(up = "#d62728", down = "#1f78b4",
                                  unknown = "#999999"),
                       name = "Expression\nvs Normal") +
    scale_size_continuous(range = c(3, 11), name = "Hub score\n(composite)") +
    scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16),
                       labels = c("TRUE" = "Has drug candidate",
                                  "FALSE" = "No drug candidate"),
                       name = "Druggable") +
    labs(title = "STRING PPI — Hub protein subnetwork",
         subtitle = sprintf(
           "%d hub proteins (degree OR betweenness \u2265 p%d) | diamonds = druggable | edge alpha \u221d STRING score",
           vcount(g_hubs), hub_percentile)) +
    theme_graph(base_family = "sans") +
    theme(legend.position   = "right",
          plot.title        = element_text(size = 14, face = "bold"),
          plot.subtitle     = element_text(size = 9),
          plot.background   = element_rect(fill = "white", color = NA),
          panel.background  = element_rect(fill = "white", color = NA))
  print(p)
})

# Exportar tambien en formato GraphML (para Cytoscape)
cat("  Exportando GraphML para Cytoscape...\n")
tryCatch({
  write_graph(g_giant, "results/tables/network/09_network_giant.graphml",
              format = "graphml")
  cat("  results/tables/network/09_network_giant.graphml\n")
}, error = function(e) {
  cat("  WARN GraphML:", e$message, "\n")
})

# =============================================================================
# EXPORTAR
# =============================================================================
cat("\n--- Exportando ---\n")
df_nodes_export <- df_nodes %>%
  left_join(df_modules %>% select(gene_symbol, module_name), by = "gene_symbol")

write.table(df_nodes_export,
            "results/tables/network/09_network_node_metrics.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_modules,
            "results/tables/network/09_modules.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_druggable,
            "results/tables/network/09_druggable_hubs.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== RESUMEN ===\n")
cat(sprintf("  Genes consultados en STRING: %d\n", length(gene_symbols)))
cat(sprintf("  Nodos en red:               %d\n", n_nodes))
cat(sprintf("  Interacciones:              %d\n", n_edges))
cat(sprintf("  Hubs (top %d%%):            %d\n", 100 - hub_percentile, n_hubs))
cat(sprintf("  Hubs druggables:            %d\n", sum(df_druggable$is_druggable)))
cat(sprintf("\nSiguiente: scripts/10_prioritization_scoring.R\n"))
cat("Fin:", format(Sys.time()), "\n")

sink(type = "message"); sink(); close(con_log)
cat("Script completado. Log en:", log_file, "\n")
