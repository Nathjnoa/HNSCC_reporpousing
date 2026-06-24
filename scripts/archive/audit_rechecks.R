# audit_rechecks.R — Re-verificaciones empíricas para auditoría pre-manuscrito
# Ambiente: omics-R
# Replica fielmente: pi_stat ranking (script 03), Louvain seed=42 + centralidades (script 09)
# Salidas: results/audit/
# NO modifica ningún resultado canónico.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})
set.seed(2026)
root <- "~/bioinfo/projects/hnscc_drug_repurposing"
setwd(path.expand(root))
dir.create("results/audit", showWarnings = FALSE, recursive = TRUE)
log_lines <- c()
LOG <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines[[length(log_lines)+1]] <<- s }

# ---------------------------------------------------------------------------
# Marcadores de composición tisular (lista curada; ver AUDIT para citas)
# Músculo esquelético/cardíaco contráctil + sarcómero + glucólisis muscular
muscle_exact <- c("SLC25A4","ATP2A1","ATP2A2","PYGM","TCAP","LDB3","FHL1","CASQ1",
                  "CASQ2","CKM","MB","DES","NEB","TTN","PGM1","ALDOA","ENO3",
                  "PVALB","SMPX","MYOT","XIRP2","SYPL2","PHKA1")
muscle_prefix <- c("MYH","MYL","MYBP","TNN","TNNT","TNNI","TNNC","ACTN","ACTA",
                   "TPM","TMOD","MYOM","MYOZ","CKMT","KLHL")
# Mucosa respiratoria / glándula salival / secretor de vía aérea
airway_exact <- c("BPIFA1","BPIFA2","BPIFB1","BPIFB2","LPO","LTF","STATH","HTN1",
                  "HTN3","PRB1","PRB2","PRH1","PRH2","ZG16B","DMBT1","SMR3B")
airway_prefix <- c("MUC","SCGB","SCGB1","SCGB2","SCGB3","PIP","CST","SMR")
is_tissue_marker <- function(g) {
  g <- toupper(g)
  hit_m <- g %in% muscle_exact | sapply(g, function(x) any(startsWith(x, muscle_prefix)))
  hit_a <- g %in% airway_exact | sapply(g, function(x) any(startsWith(x, airway_prefix)))
  data.frame(gene = g, muscle = hit_m, airway = hit_a,
             tissue_marker = hit_m | hit_a, stringsAsFactors = FALSE)
}

# ===========================================================================
# RE-CHECK 1 — Composición tisular en el eje down
# ===========================================================================
LOG("===== RE-CHECK 1: COMPOSICIÓN TISULAR (eje down) =====")
de <- read_tsv("results/tables/de_limma/01_TVsS_all_proteins.tsv", show_col_types = FALSE)
sig <- de %>% filter(adj.P.Val_TVsS < 0.05, abs(logFC_TVsS) > 1) %>%
  mutate(dir = ifelse(logFC_TVsS > 0, "up", "down"))
mk <- is_tissue_marker(sig$gene_symbol)
sig <- bind_cols(sig, mk[, c("muscle","airway","tissue_marker")])
n_down <- sum(sig$dir == "down"); n_up <- sum(sig$dir == "up")
down_mk <- sig %>% filter(dir == "down")
up_mk   <- sig %>% filter(dir == "up")
LOG("DE significativas: %d (up=%d, down=%d)", nrow(sig), n_up, n_down)
LOG("Marcadores tisulares en DOWN: %d/%d (%.1f%%)  [músculo=%d, airway=%d]",
    sum(down_mk$tissue_marker), n_down, 100*mean(down_mk$tissue_marker),
    sum(down_mk$muscle), sum(down_mk$airway))
LOG("Marcadores tisulares en UP:   %d/%d (%.1f%%)",
    sum(up_mk$tissue_marker), n_up, 100*mean(up_mk$tissue_marker))
# Distribución en el ranking down (más negativo = top)
down_sorted <- down_mk %>% arrange(logFC_TVsS)
LOG("De los 20 down más extremos: %d son marcadores tisulares (%.0f%%)",
    sum(down_sorted$tissue_marker[1:20]), 100*mean(down_sorted$tissue_marker[1:20]))
LOG("De los 50 down más extremos: %d son marcadores tisulares (%.0f%%)",
    sum(down_sorted$tissue_marker[1:min(50,n_down)]),
    100*mean(down_sorted$tissue_marker[1:min(50,n_down)]))
write_tsv(sig %>% select(gene_symbol, logFC_TVsS, adj.P.Val_TVsS, dir,
                         muscle, airway, tissue_marker) %>% arrange(logFC_TVsS),
          "results/audit/recheck1_tissue_markers.tsv")

# ===========================================================================
# RE-CHECK 2 — Sensibilidad GSEA al excluir contaminantes
# ===========================================================================
LOG("\n===== RE-CHECK 2: SENSIBILIDAD GSEA HALLMARKS =====")
gsea_ok <- tryCatch({
  suppressPackageStartupMessages({ library(clusterProfiler); library(msigdbr) })
  TRUE }, error = function(e) FALSE)
if (gsea_ok) {
  # Replica pi_stat ranking (script 03 línea 65)
  rk <- de %>% filter(!is.na(logFC_TVsS), !is.na(adj.P.Val_TVsS)) %>%
    mutate(pi = sign(logFC_TVsS)*abs(logFC_TVsS)*(-log10(adj.P.Val_TVsS + 1e-300))) %>%
    arrange(desc(pi))
  msig <- tryCatch(msigdbr(species = "Homo sapiens", category = "H"),
                   error = function(e) msigdbr(species = "Homo sapiens", collection = "H"))
  t2g <- msig %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
  run_gsea <- function(genes_vec, pi_vec) {
    gl <- pi_vec; names(gl) <- genes_vec; gl <- sort(gl, decreasing = TRUE)
    set.seed(42)
    suppressWarnings(GSEA(gl, TERM2GENE = t2g, pvalueCutoff = 1.0,
         minGSSize = 10, maxGSSize = 500, seed = TRUE, eps = 0, verbose = FALSE))
  }
  res_full <- run_gsea(rk$gene_symbol, rk$pi)
  # Excluir contaminantes del ranking
  mk_all <- is_tissue_marker(rk$gene_symbol)
  rk_clean <- rk[!mk_all$tissue_marker, ]
  res_clean <- run_gsea(rk_clean$gene_symbol, rk_clean$pi)
  foc <- c("HALLMARK_MYOGENESIS","HALLMARK_ADIPOGENESIS","HALLMARK_HEME_METABOLISM",
           "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_E2F_TARGETS",
           "HALLMARK_G2M_CHECKPOINT","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_FATTY_ACID_METABOLISM","HALLMARK_INTERFERON_GAMMA_RESPONSE")
  gx <- function(r) as.data.frame(r) %>% filter(ID %in% foc) %>%
    select(ID, NES, p.adjust)
  cmp <- full_join(gx(res_full) %>% rename(NES_full=NES, padj_full=p.adjust),
                   gx(res_clean) %>% rename(NES_clean=NES, padj_clean=p.adjust),
                   by = "ID") %>% arrange(ID)
  print(cmp)
  write_tsv(cmp, "results/audit/recheck2_gsea_sensitivity.tsv")
  for (i in seq_len(nrow(cmp))) {
    r <- cmp[i,]
    sig_full <- !is.na(r$padj_full) && r$padj_full < 0.05
    sig_clean <- !is.na(r$padj_clean) && r$padj_clean < 0.05
    LOG("%s: full %s (NES=%.2f) -> clean %s (NES=%s)",
        gsub("HALLMARK_","",r$ID),
        ifelse(sig_full,"SIG","ns"), ifelse(is.na(r$NES_full),NA,r$NES_full),
        ifelse(sig_clean,"SIG","ns"),
        ifelse(is.na(r$NES_clean),"NA",sprintf("%.2f",r$NES_clean)))
  }
} else { LOG("GSEA omitido (paquetes no disponibles)") }

# ===========================================================================
# RE-CHECK 3 — Sensibilidad de la red al excluir contaminantes
# ===========================================================================
LOG("\n===== RE-CHECK 3: SENSIBILIDAD DE RED (Louvain + centralidad) =====")
library(igraph)
edges <- read_tsv("results/tables/network/09_network_edges.tsv", show_col_types = FALSE)
nm_orig <- read_tsv("results/tables/network/09_network_node_metrics.tsv", show_col_types = FALSE)
build_metrics <- function(edge_df) {
  g <- graph_from_data_frame(edge_df[,c("gene_A","gene_B")], directed = FALSE)
  comp <- components(g); g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  deg <- degree(g); betw <- betweenness(g, normalized = TRUE)
  eig <- tryCatch(eigen_centrality(g)$vector, error=function(e) rep(NA, vcount(g)))
  set.seed(42); lou <- cluster_louvain(g, resolution = 1.0)
  data.frame(gene_symbol = names(deg), degree = deg, betweenness_norm = betw,
             eigenvector = eig, module_id = membership(lou), row.names = NULL)
}
m_full <- build_metrics(edges)
contam <- is_tissue_marker(unique(c(edges$gene_A, edges$gene_B)))
contam_genes <- contam$gene[contam$tissue_marker]
edges_clean <- edges %>% filter(!(toupper(gene_A) %in% contam_genes |
                                  toupper(gene_B) %in% contam_genes))
m_clean <- build_metrics(edges_clean)
LOG("Red original: %d nodos, %d aristas; contaminantes en red: %d",
    nrow(m_full), nrow(edges), length(intersect(toupper(m_full$gene_symbol), contam_genes)))
LOG("Red limpia:   %d nodos, %d aristas", nrow(m_clean), nrow(edges_clean))
hub_status <- function(m, genes) {
  thr_d <- quantile(m$degree, 0.90); thr_b <- quantile(m$betweenness_norm, 0.90)
  m %>% filter(gene_symbol %in% genes) %>%
    mutate(is_hub = degree >= thr_d | betweenness_norm >= thr_b,
           deg_pctl = round(100*rank(degree)/nrow(m))) %>%
    select(gene_symbol, degree, betweenness_norm, is_hub)
}
focus_genes <- c("EGFR","NDUFS2","NDUFS3","NDUFA9","PSMA2","PSMB10","PSMB5","DNMT1","RPS11")
hs_full <- hub_status(m_full, focus_genes) %>% rename(deg_full=degree, betw_full=betweenness_norm, hub_full=is_hub)
hs_clean<- hub_status(m_clean,focus_genes) %>% rename(deg_clean=degree, betw_clean=betweenness_norm, hub_clean=is_hub)
hs <- full_join(hs_full, hs_clean, by="gene_symbol")
print(hs)
write_tsv(hs, "results/audit/recheck3_network_hub_sensitivity.tsv")
for (i in seq_len(nrow(hs))) {
  r <- hs[i,]
  LOG("%s: hub_full=%s (deg=%s) -> hub_clean=%s (deg=%s)", r$gene_symbol,
      r$hub_full, ifelse(is.na(r$deg_full),"NA",r$deg_full),
      ifelse(is.na(r$hub_clean),"NA(out)",as.character(r$hub_clean)),
      ifelse(is.na(r$deg_clean),"NA",r$deg_clean))
}

# ===========================================================================
# RE-CHECK 5 — Fig6C por gen (CPTAC + TCGA por ancla)
# ===========================================================================
LOG("\n===== RE-CHECK 5: VALIDACIÓN POR GEN (Fig6C) =====")
scored <- read_tsv("results/tables/10_all_candidates_scored.tsv", show_col_types = FALSE)
shortlist <- scored %>% filter(tier != "off_network", !is.na(primary_target)) %>%
  group_by(primary_target) %>% slice_max(composite_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>% slice_max(composite_score, n = 14) %>%
  select(primary_target, drug_name_norm, composite_score, tier, primary_pi)
cptac <- read_tsv("data/intermediate/cptac/16c_cptac_hnscc_de.tsv", show_col_types = FALSE)
tcga  <- read_tsv("data/intermediate/tcga/tcga_hnsc_deseq2_results.tsv", show_col_types = FALSE)
disc  <- de %>% select(gene_symbol, logFC_disc = logFC_TVsS, adjP_disc = adj.P.Val_TVsS)
fig6c <- shortlist %>%
  left_join(disc, by = c("primary_target"="gene_symbol")) %>%
  left_join(cptac, by = c("primary_target"="gene_symbol")) %>%
  left_join(tcga,  by = c("primary_target"="gene_symbol")) %>%
  mutate(conc_cptac = sign(logFC_disc) == sign(logFC_cptac),
         conc_tcga  = sign(logFC_disc) == sign(logFC_tcga),
         sig_cptac  = !is.na(adjP_cptac) & adjP_cptac < 0.05,
         sig_tcga   = !is.na(adjP_tcga)  & adjP_tcga  < 0.05)
print(as.data.frame(fig6c %>% select(primary_target, drug_name_norm, logFC_disc,
      logFC_cptac, adjP_cptac, conc_cptac, logFC_tcga, adjP_tcga, conc_tcga)))
write_tsv(fig6c, "results/audit/recheck5_fig6c_per_gene.tsv")
LOG("CPTAC: con datos %d/14, concordantes %d, FDR<0.05 %d, discordantes: %s",
    sum(!is.na(fig6c$logFC_cptac)), sum(fig6c$conc_cptac, na.rm=TRUE),
    sum(fig6c$sig_cptac, na.rm=TRUE),
    paste(fig6c$primary_target[which(!fig6c$conc_cptac)], collapse=","))
LOG("TCGA:  con datos %d/14, concordantes %d, FDR<0.05 %d, discordantes: %s",
    sum(!is.na(fig6c$logFC_tcga)), sum(fig6c$conc_tcga, na.rm=TRUE),
    sum(fig6c$sig_tcga, na.rm=TRUE),
    paste(fig6c$primary_target[which(!fig6c$conc_tcga)], collapse=","))

# ===========================================================================
# RE-CHECK 6 — Ablación de baseline (¿el composite aporta sobre centralidad/π?)
# ===========================================================================
LOG("\n===== RE-CHECK 6: ABLACIÓN BASELINE =====")
cand <- scored %>% filter(tier != "off_network")
# Rankings alternativos
cand <- cand %>% mutate(
  rank_composite = rank(-composite_score, ties.method="min"),
  rank_cent = rank(-primary_target_centrality, ties.method="min"),
  rank_pi   = rank(-primary_pi, ties.method="min"),
  rank_tp   = rank(-TargetPriority, ties.method="min"))
sp <- function(a,b) suppressWarnings(cor(a,b,method="spearman",use="complete.obs"))
LOG("Spearman composite vs centrality-only: %.3f", sp(cand$composite_score, cand$primary_target_centrality))
LOG("Spearman composite vs pi-only:         %.3f", sp(cand$composite_score, cand$primary_pi))
LOG("Spearman composite vs TargetPriority:  %.3f", sp(cand$composite_score, cand$TargetPriority))
top14_comp <- cand %>% slice_min(rank_composite, n=14) %>% pull(drug_name_norm)
top14_cent <- cand %>% slice_min(rank_cent, n=14) %>% pull(drug_name_norm)
top14_pi   <- cand %>% slice_min(rank_pi, n=14) %>% pull(drug_name_norm)
LOG("Solapamiento top14 composite vs centrality-only: %d/14", length(intersect(top14_comp, top14_cent)))
LOG("Solapamiento top14 composite vs pi-only:         %d/14", length(intersect(top14_comp, top14_pi)))
write_tsv(cand %>% select(drug_name_norm, primary_target, composite_score, TargetPriority,
          DrugViability, primary_target_centrality, primary_pi,
          rank_composite, rank_cent, rank_pi) %>% arrange(rank_composite),
          "results/audit/recheck6_baseline_ablation.tsv")

# ===========================================================================
# RE-CHECK 7 — Manejo de NA en L2S2 (sesgo contra anticuerpos)
# ===========================================================================
LOG("\n===== RE-CHECK 7: NA-HANDLING L2S2 (anticuerpos) =====")
# Pesos v3: composite = 0.60*TP + 0.40*DV ; DV = 0.40*rev + 0.40*reg + 0.20*breadth
w_target<-0.60; w_drug<-0.40; w_rev<-0.40; w_reg<-0.40; w_brd<-0.20
ab <- scored %>% mutate(
  cmap_na = is.na(cmap_score),
  # DV original (NA rev -> 0, ya está en s_reversal)
  DV_orig = DrugViability,
  # DV alternativo: si no hay reversal (cmap NA), renormalizar reg+breadth
  DV_alt = ifelse(cmap_na,
                  (w_reg*s_regulatory + w_brd*s_breadth)/(w_reg+w_brd),
                  DrugViability),
  composite_alt = w_target*TargetPriority + w_drug*DV_alt,
  rank_orig = rank(-composite_score, ties.method="min"),
  rank_alt  = rank(-composite_alt,  ties.method="min"))
mabs <- ab %>% filter(str_detect(toupper(drug_name_norm),
        "CETUXIMAB|PANITUMUMAB|NIMOTUZUMAB|PEMBROLIZUMAB|NIVOLUMAB|TRASTUZUMAB|BEVACIZUMAB|DURVALUMAB"))
print(as.data.frame(mabs %>% select(drug_name_norm, primary_target, cmap_na,
      composite_score, composite_alt, rank_orig, rank_alt)))
write_tsv(ab %>% filter(cmap_na) %>% select(drug_name_norm, primary_target, drug_class,
          composite_score, composite_alt, rank_orig, rank_alt) %>% arrange(rank_alt),
          "results/audit/recheck7_l2s2_na_handling.tsv")
LOG("Drogas con cmap NA (sin firma L2S2): %d; de ellas anticuerpos mostrados arriba", sum(ab$cmap_na))
for (i in seq_len(nrow(mabs))) { r <- mabs[i,]
  LOG("%s: composite %.3f (rank %d) -> NA-excluido %.3f (rank %d)", r$drug_name_norm,
      r$composite_score, r$rank_orig, r$composite_alt, r$rank_alt) }

writeLines(unlist(log_lines), "results/audit/rechecks_log.txt")
LOG("\n===== FIN. Salidas en results/audit/ =====")
