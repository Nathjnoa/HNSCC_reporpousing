# =============================================================================
# scripts/_setup.R — Bootstrap compartido para los scripts R del pipeline
# hnscc_drug_repurposing.
#
# Reemplaza el bloque duplicado de resolución de working directory que estaba
# copiado en ~9 scripts. Usa el paquete {here}, que ancla la raíz del proyecto
# en el directorio que contiene .git (o .here / .Rproj), de forma independiente
# a la profundidad del script. Esto funciona igual para scripts en `scripts/`
# y en `scripts/supp/`.
#
# Uso (al inicio de cada script, tras cargar las librerías):
#     source(here::here("scripts", "_setup.R"))
#     setup_project()
#
# Requisito: ejecutar con el working directory dentro del repo
# (p. ej. `Rscript scripts/10_prioritization_scoring.R` desde la raíz).
# =============================================================================

suppressPackageStartupMessages(library(here))

# Fija el working directory en la raíz del proyecto y lo reporta.
# Devuelve (invisible) la ruta de la raíz.
setup_project <- function() {
  root <- here::here()
  setwd(root)
  cat("Working directory:", getwd(), "\n")
  invisible(root)
}
