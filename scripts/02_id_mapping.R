## ============================================================
## Script 02: ID Mapping — UniProt → Gene Symbol + Entrez ID
## Proyecto: HNSCC Drug Repurposing
## Input:  results/tables/de_limma/01_TVsS_significant.tsv
##         results/tables/de_limma/01_TVsS_all_proteins.tsv
## Output: results/tables/de_limma/02_TVsS_significant_with_ids.tsv
##         results/tables/de_limma/02_TVsS_all_proteins_with_ids.tsv
##         data/intermediate/id_mapping/02_uniprot_to_ids.tsv
## ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(here)
})

## ── Paths ────────────────────────────────────────────────────
proj_dir   <- here::here()
de_dir     <- file.path(proj_dir, "results/tables/de_limma")
map_dir    <- file.path(proj_dir, "data/intermediate/id_mapping")
log_dir    <- file.path(proj_dir, "logs")

dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("02_id_mapping_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
con <- file(log_file, open = "wt")
sink(con, type = "message")
sink(con, type = "output", append = TRUE)

cat("=== Script 02: ID Mapping ===\n")
cat("Started:", format(Sys.time()), "\n\n")

## ── Load tables ──────────────────────────────────────────────
sig_tbl <- read.delim(file.path(de_dir, "01_TVsS_significant.tsv"),
                      stringsAsFactors = FALSE)
all_tbl <- read.delim(file.path(de_dir, "01_TVsS_all_proteins.tsv"),
                      stringsAsFactors = FALSE)

cat("Proteins loaded — significant:", nrow(sig_tbl),
    " | all:", nrow(all_tbl), "\n")

## ── Map UniProt → Entrez ID + canonical Symbol ───────────────
uniprot_ids <- unique(all_tbl$uniprot_id)

# Some entries may be isoform IDs (e.g. Q9Y4L1-2); strip suffix for mapping
uniprot_clean <- sub("-\\d+$", "", uniprot_ids)

entrez_map <- mapIds(org.Hs.eg.db,
                     keys    = uniprot_clean,
                     keytype = "UNIPROT",
                     column  = "ENTREZID",
                     multiVals = "first")

symbol_map <- mapIds(org.Hs.eg.db,
                     keys    = uniprot_clean,
                     keytype = "UNIPROT",
                     column  = "SYMBOL",
                     multiVals = "first")

id_map_df <- data.frame(
  uniprot_id     = uniprot_ids,
  uniprot_clean  = uniprot_clean,
  symbol_org     = unname(symbol_map[uniprot_clean]),
  entrez_id      = unname(entrez_map[uniprot_clean]),
  stringsAsFactors = FALSE
)

## ── Mapping statistics ───────────────────────────────────────
n_entrez_mapped   <- sum(!is.na(id_map_df$entrez_id))
n_entrez_unmapped <- sum(is.na(id_map_df$entrez_id))

cat("\n--- Mapping statistics ---\n")
cat("Total UniProt IDs:    ", nrow(id_map_df), "\n")
cat("Entrez mapped:        ", n_entrez_mapped,
    sprintf("(%.1f%%)\n", 100 * n_entrez_mapped / nrow(id_map_df)))
cat("Entrez unmapped:      ", n_entrez_unmapped,
    sprintf("(%.1f%%)\n", 100 * n_entrez_unmapped / nrow(id_map_df)))

# Symbol consistency check
id_map_df$symbol_match <- id_map_df$symbol_org == id_map_df$uniprot_clean |
                          is.na(id_map_df$symbol_org) # unused — candidate for removal
# Compare with gene_symbol from proteoDA output
all_merged_check <- left_join(all_tbl, id_map_df, by = "uniprot_id")
symbol_discrepancies <- sum(
  !is.na(all_merged_check$symbol_org) &
  all_merged_check$gene_symbol != all_merged_check$symbol_org,
  na.rm = TRUE
)
cat("Symbol discrepancies vs proteoDA: ", symbol_discrepancies,
    "(using proteoDA gene_symbol as canonical)\n\n")

## ── Merge IDs into tables ────────────────────────────────────
add_ids <- function(tbl, map_df) {
  tbl %>%
    left_join(map_df %>% select(uniprot_id, entrez_id, symbol_org),
              by = "uniprot_id") %>%
    relocate(gene_symbol, uniprot_id, entry_name, symbol_org, entrez_id)
}

sig_with_ids <- add_ids(sig_tbl, id_map_df)
all_with_ids <- add_ids(all_tbl, id_map_df)

## Check significant proteins mapping
n_sig_entrez <- sum(!is.na(sig_with_ids$entrez_id))
cat("Significant proteins with Entrez ID:", n_sig_entrez, "/", nrow(sig_with_ids), "\n")
cat("Significant proteins WITHOUT Entrez ID:",
    nrow(sig_with_ids) - n_sig_entrez, "\n")

# Report unmapped significant proteins
unmapped_sig <- sig_with_ids %>%
  filter(is.na(entrez_id)) %>%
  select(gene_symbol, uniprot_id, logFC_TVsS, adj.P.Val_TVsS)

if (nrow(unmapped_sig) > 0) {
  cat("\nUnmapped significant proteins:\n")
  print(unmapped_sig)
}

## ── Export ───────────────────────────────────────────────────

# 1. Significant proteins with IDs
write.table(sig_with_ids,
            file.path(de_dir, "02_TVsS_significant_with_ids.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 2. All proteins with IDs (for GSEA ranked list in script 03)
write.table(all_with_ids,
            file.path(de_dir, "02_TVsS_all_proteins_with_ids.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 3. Mapping reference table
write.table(id_map_df,
            file.path(map_dir, "02_uniprot_to_ids.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n--- Outputs exported ---\n")
cat("02_TVsS_significant_with_ids.tsv  — ", nrow(sig_with_ids), "rows\n")
cat("02_TVsS_all_proteins_with_ids.tsv — ", nrow(all_with_ids), "rows\n")
cat("02_uniprot_to_ids.tsv             — ", nrow(id_map_df), "rows\n")

cat("\nFinished:", format(Sys.time()), "\n")

sink(type = "output")
sink(type = "message")
close(con)
cat("Log saved:", log_file, "\n")
