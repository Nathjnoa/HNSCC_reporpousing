# tests/test_pathway_scoring.R
# Reproduce bug: s_pathway colapsa a 0 porque el script 10 leía archivos ORA
# que fueron eliminados en la limpieza (9f0653c). Se migró a GSEA leading-edge
# sobre GO BP + KEGG + Reactome (bases con vías de señalización donde EGFR sí
# aparece; Hallmarks no lo contiene por diseño y se reserva para la figura).
#
# Ambiente: omics-R
# Ejecutar desde la raíz del proyecto:
#   Rscript tests/test_pathway_scoring.R
#
# Pre-fix  -> FALLA (s_pathway == 0 para fármacos EGFR; pathway pool vacío)
# Post-fix -> PASA  (EGFR en leading-edge GSEA GO/KEGG/Reactome; s_pathway == 1)

suppressPackageStartupMessages(library(here))
setwd(here::here())

top <- read.delim("results/tables/10_top20_candidates.tsv",
                  stringsAsFactors = FALSE)

# Gefitinib (inhibidor EGFR, target único EGFR) debe tener s_pathway > 0:
# EGFR aparece en el core_enrichment (leading-edge) de Hallmarks GSEA.
gef <- top[top$drug_name_norm == "GEFITINIB", ]

stopifnot(
  "GEFITINIB ausente del top20"                      = nrow(gef) == 1,
  "bug: s_pathway de GEFITINIB es 0 (pool vacío)"    = gef$s_pathway[1] > 0,
  "bug: s_pathway debería ser 1 (única diana EGFR en vía)" =
    abs(gef$s_pathway[1] - 1) < 1e-6
)

cat("TEST PASSED — s_pathway(GEFITINIB) =", gef$s_pathway[1], "\n")
