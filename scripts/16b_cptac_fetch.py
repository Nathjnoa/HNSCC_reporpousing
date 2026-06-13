#!/usr/bin/env python3
# =============================================================================
# 16b_cptac_fetch.py — Descarga datos proteómicos CPTAC-HNSCC
# =============================================================================
# Obtiene la proteómica cuantitativa (log2-ratios TMT) del dataset CPTAC-HNSCC
# (Huang et al. 2021, Cancer Cell) vía el paquete `cptac` (Payne lab, BYU).
#
# SOLO descarga y exporta la matriz y metadata de muestras. No hace estadística:
# el DE se realiza en R con limma (script 16c_cptac_de.R) para mantener el mismo
# método que nuestra proteómica de descubrimiento (proteoDA/limma).
#
# Outputs (data/intermediate/cptac/):
#   cptac_hnscc_proteomics.tsv — genes × muestras (log2-ratios TMT, umich)
#   cptac_hnscc_samples.tsv   — sample_id, condition, patient_id
#
# Ambiente: omics-py
# Dependencia: pip install cptac==1.5.14
# Ejecución (desde raíz del proyecto):
#   conda run -n omics-py python scripts/16b_cptac_fetch.py
# Ejecutar ANTES de 16c_cptac_de.R
# =============================================================================

import os
import sys
import warnings
import logging
from datetime import datetime

import pandas as pd

# Silenciar warnings menores del paquete cptac
warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURACIÓN
# =============================================================================
CPTAC_SOURCE    = "umich"       # Preferido; bcm como fallback
CPTAC_SOURCE_FB = "bcm"
OUT_DIR         = "data/intermediate/cptac"
LOG_DIR         = "logs"
PROJECT_ROOT    = os.getcwd()

# =============================================================================
# LOGGING
# =============================================================================
os.makedirs(LOG_DIR, exist_ok=True)
log_path = os.path.join(LOG_DIR, f"16b_cptac_fetch_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(log_path),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

log.info("=== 16b_cptac_fetch.py ===")
log.info(f"Directorio de trabajo: {PROJECT_ROOT}")
log.info(f"Fuente CPTAC:          {CPTAC_SOURCE} (fallback: {CPTAC_SOURCE_FB})")

# =============================================================================
# SECCIÓN 1: CARGA DEL DATASET CPTAC-HNSCC
# =============================================================================
log.info("\n--- Cargando dataset CPTAC-HNSCC ---")

try:
    import cptac
    log.info(f"cptac versión: {getattr(cptac, '__version__', '1.5.14')}")
except ImportError:
    log.error("Paquete cptac no instalado. Ejecutar: pip install cptac")
    sys.exit(1)

# La primera ejecución descarga los datos (~varios GB) y los cachea localmente.
# Las ejecuciones posteriores leen del caché.
hnscc = cptac.Hnscc()
log.info("Dataset CPTAC-HNSCC cargado (caché local)")

# =============================================================================
# SECCIÓN 2: DATOS PROTEÓMICOS
# =============================================================================
log.info("\n--- Obteniendo datos proteómicos ---")

prot = None
used_source = None
for source in [CPTAC_SOURCE, CPTAC_SOURCE_FB]:
    try:
        prot = hnscc.get_proteomics(source=source)
        used_source = source
        log.info(f"  Fuente usada: {source}")
        break
    except Exception as e:
        log.warning(f"  Fuente {source} no disponible: {e}")

if prot is None:
    log.error("Ninguna fuente proteómica disponible. Abortando.")
    sys.exit(1)

log.info(f"  Dimensiones brutas: {prot.shape[0]} muestras × {prot.shape[1]} proteínas")

# Colapsar MultiIndex de columnas (Name, Database_ID) → gene_symbol (primer nivel)
if isinstance(prot.columns, pd.MultiIndex):
    prot.columns = prot.columns.get_level_values(0)
    log.info("  MultiIndex de columnas colapsado a gene_symbol")

# Eliminar columnas duplicadas de gene_symbol (mantener primera ocurrencia)
n_before = prot.shape[1]
prot = prot.loc[:, ~prot.columns.duplicated()]
n_after = prot.shape[1]
if n_before != n_after:
    log.info(f"  Columnas duplicadas eliminadas: {n_before} → {n_after}")

# =============================================================================
# SECCIÓN 3: METADATA DE MUESTRAS — tumor/normal y emparejamiento
# =============================================================================
log.info("\n--- Construyendo metadata de muestras ---")

# En CPTAC-HNSCC los normales tienen sufijo '.N' en el sample_id (índice)
sample_ids  = list(prot.index)
conditions  = ["Normal" if s.endswith(".N") else "Tumor" for s in sample_ids]
patient_ids = [s[:-2] if s.endswith(".N") else s for s in sample_ids]

samples = pd.DataFrame({
    "sample_id":  sample_ids,
    "condition":  conditions,
    "patient_id": patient_ids,
})

n_tumor  = (samples["condition"] == "Tumor").sum()
n_normal = (samples["condition"] == "Normal").sum()
paired   = set(samples[samples["condition"] == "Normal"]["patient_id"]) & \
           set(samples[samples["condition"] == "Tumor"]["patient_id"])
n_paired = len(paired)

log.info(f"  Tumores:    {n_tumor}")
log.info(f"  Normales:   {n_normal}")
log.info(f"  Pares T/N:  {n_paired}  (mismo patient_id en ambas condiciones)")
log.info(f"  NOTA: diseño PAREADO — usar duplicateCorrelation(block=patient_id) en limma")

# =============================================================================
# SECCIÓN 4: EXPORTAR OUTPUTS
# =============================================================================
log.info("\n--- Exportando TSV ---")
os.makedirs(OUT_DIR, exist_ok=True)

# 4a. Matriz proteómica (samples = filas, genes = columnas)
# Transponer a formato genes × muestras para que R lo lea cómodamente
prot_out = prot.T
prot_out.index.name = "gene_symbol"
prot_path = os.path.join(OUT_DIR, "cptac_hnscc_proteomics.tsv")
prot_out.to_csv(prot_path, sep="\t")
log.info(f"  {prot_path}  ({prot_out.shape[0]} genes × {prot_out.shape[1]} muestras)")

# 4b. Metadata de muestras
samples_path = os.path.join(OUT_DIR, "cptac_hnscc_samples.tsv")
samples.to_csv(samples_path, sep="\t", index=False)
log.info(f"  {samples_path}  ({len(samples)} muestras)")

# =============================================================================
# RESUMEN
# =============================================================================
log.info("\n========================================")
log.info("RESUMEN 16b_cptac_fetch.py")
log.info("========================================")
log.info(f"Fuente:       CPTAC-HNSCC ({used_source})")
log.info(f"Tumores:      {n_tumor}")
log.info(f"Normales:     {n_normal}")
log.info(f"Pares T/N:    {n_paired}")
log.info(f"Genes:        {prot_out.shape[0]}")
log.info("Diseño:       Pareado — usar duplicateCorrelation en limma (script 16c)")
log.info("\nOutputs:")
log.info(f"  {prot_path}")
log.info(f"  {samples_path}")
log.info(f"\nSiguiente paso: Rscript scripts/16c_cptac_de.R")
log.info(f"\nFin: {datetime.now()}")
