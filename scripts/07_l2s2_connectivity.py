#!/usr/bin/env python3
"""
Script 07: L2S2 Connectivity Analysis (reemplaza CMap2)
HNSCC Drug Repurposing — Fase 3

Consulta LINCS L1000 Signature Search (L2S2) via GraphQL API para encontrar
fármacos aprobados cuya firma transcriptómica revierte el perfil proteómico tumoral.

Método:
  - Query A: genes UP en tumor → busca firmas DOWN del fármaco (reversor)
  - Query B: genes DOWN en tumor → busca firmas UP del fármaco (reversor)
  - Agrega por fármaco: n_sigs, OR, pvalue → reversal_score normalizado a [-1, 0]

Ref: https://l2s2.maayanlab.cloud | Evangelista et al., NAR 2025
     https://doi.org/10.1093/nar/gkae1059

Input:  results/tables/de_limma/02_TVsS_significant_with_ids.tsv
Output: results/tables/drug_targets/07_l2s2_results.tsv
        results/tables/drug_targets/07_l2s2_top_reversors.tsv
"""

import sys
import re
import time
import math
import logging
import os
from datetime import datetime
from collections import defaultdict

import requests
import pandas as pd

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
os.makedirs("logs", exist_ok=True)
os.makedirs("results/tables/drug_targets", exist_ok=True)

log_file = f"logs/07_l2s2_connectivity_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

QUERY_DATE = datetime.now().strftime("%Y-%m-%d")
log.info(f"L2S2 query date: {QUERY_DATE}")
log.warning(
    "METHODOLOGICAL NOTE: L2S2 uses transcriptomic (mRNA) signatures. "
    "This analysis inputs proteomic fold-changes. The assumption that "
    "protein and mRNA changes are directionally consistent is valid only "
    "for ~60-70%% of genes (mRNA-protein r~0.4-0.6 in clinical samples). "
    "L2S2 score weight has been reduced to 0.10 in the composite score. "
    "Results should be interpreted as supportive, not primary, evidence."
)

log.info("=== 07_l2s2_connectivity.py ===")
log.info(f"Inicio: {datetime.now()}")

# ---------------------------------------------------------------------------
# Parámetros
# ---------------------------------------------------------------------------
L2S2_URL      = "https://l2s2.maayanlab.cloud/graphql"
TOP_N_GENES   = 150       # genes per direction (up/down) sent to L2S2 query
PVALUE_LE     = 0.001     # leading-edge p-value filter for L2S2 results
BATCH_SIZE    = 500       # nodos por request de paginación
MAX_NODES     = 50_000    # safety cap to avoid memory issues on large graphs
SLEEP_SEC     = 0.3       # delay between L2S2 API requests (rate limit)
TIMEOUT       = 90        # segundos timeout por request

# ---------------------------------------------------------------------------
# Query GraphQL
# ---------------------------------------------------------------------------
ENRICH_QUERY = """
query($genes: [String]!, $topN: Int, $pvalueLe: Float,
      $offset: Int, $first: Int, $filterFda: Boolean) {
  currentBackground {
    enrich(genes: $genes, topN: $topN, pvalueLe: $pvalueLe,
           offset: $offset, first: $first, filterFda: $filterFda) {
      totalCount
      nodes {
        pvalue
        adjPvalue
        oddsRatio
        nOverlap
        geneSets {
          nodes {
            term
            nGeneIds
          }
        }
      }
    }
  }
}
"""

# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------
CONC_RE = re.compile(r'^\d+\.?\d*[nuMmpkg/]', re.IGNORECASE)

def parse_drug_term(term):
    """
    Parsea un término LINCS L1000 para extraer nombre del fármaco, línea celular
    y dirección.
    Formato típico: PLATE_CELL_TIME_WELL_DRUG[_CONC] up/down
    Ejemplos:
      LJP009_SHSY5Y_24H_I22_alectinib_0.37uM down  → ('alectinib', 'SHSY5Y', 'down')
      CPC006_A375_24H_X1_imatinib_1uM down          → ('imatinib', 'A375', 'down')
      XPR023_A375.311_96H_E12_RTP5 up               → CRISPR KO, se omite
    """
    if term.startswith('XPR') or term.startswith('ENR'):
        return None, None, None  # CRISPR KO o enriquecimiento — omitir

    if term.endswith(' up'):
        direction = 'up'
        base = term[:-3].strip()
    elif term.endswith(' down'):
        direction = 'down'
        base = term[:-5].strip()
    else:
        return None, None, None

    parts = base.split('_')
    if len(parts) < 4:
        return None, None, None

    cell_line = parts[1] if len(parts) > 1 else 'unknown'

    # Buscar índice de concentración (último elemento que parece concentración)
    conc_idx = None
    for i in range(len(parts) - 1, 3, -1):  # minimum parts to consider a valid drug name token
        if CONC_RE.match(parts[i]):
            conc_idx = i
            break

    if conc_idx is not None and conc_idx > 4:
        drug = '_'.join(parts[4:conc_idx])
    elif len(parts) > 4:
        drug = '_'.join(parts[4:])
    else:
        drug = parts[-1]

    drug = drug.strip().lower()
    if not drug or drug.startswith('brd-') and len(drug) < 8:
        return None, None, None

    return drug, cell_line, direction


def query_l2s2_paginated(genes, target_direction, label=""):
    """
    Consulta L2S2 con una lista de genes y retorna todos los nodos con
    firmas en target_direction (up o down) para fármacos FDA-aprobados.
    """
    all_drug_hits = defaultdict(list)
    offset = 0
    total_retrieved = 0
    total_count = None

    log.info(f"  [{label}] Iniciando query → firmas '{target_direction}' | genes={len(genes)}")

    while True:
        payload = {
            "query": ENRICH_QUERY,
            "variables": {
                "genes":      genes,
                "topN":       MAX_NODES,
                "pvalueLe":   PVALUE_LE,
                "offset":     offset,
                "first":      BATCH_SIZE,
                "filterFda":  True,
            }
        }

        try:
            r = requests.post(L2S2_URL, json=payload, timeout=TIMEOUT)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            log.error(f"  [{label}] Error en request (offset={offset}): {e}")
            break

        if "errors" in data:
            log.error(f"  [{label}] GraphQL errors: {data['errors']}")
            break

        result    = data["data"]["currentBackground"]["enrich"]
        nodes     = result["nodes"]
        if total_count is None:
            total_count = result["totalCount"]
            log.info(f"  [{label}] Total nodos disponibles (pvalue<{PVALUE_LE}, FDA): {total_count}")

        if not nodes:
            break

        # Procesar nodos
        for node in nodes:
            p    = node["pvalue"]
            adjp = node["adjPvalue"]
            OR   = node["oddsRatio"]
            noverlap = node["nOverlap"]

            for gs in node["geneSets"]["nodes"]:
                term = gs["term"]
                drug, cell_line, direction = parse_drug_term(term)

                if drug is None or direction != target_direction:
                    continue

                all_drug_hits[drug].append({
                    "term":       term,
                    "cell_line":  cell_line,
                    "direction":  direction,
                    "pvalue":     p,
                    "adjPvalue":  adjp,
                    "oddsRatio":  OR,
                    "nOverlap":   noverlap,
                })

        total_retrieved += len(nodes)
        log.info(f"  [{label}] Procesados {total_retrieved}/{total_count} nodos | "
                 f"drugs únicos hasta ahora: {len(all_drug_hits)}")

        if total_retrieved >= total_count or total_retrieved >= MAX_NODES:
            break

        offset += BATCH_SIZE
        time.sleep(SLEEP_SEC)

    log.info(f"  [{label}] Completado: {len(all_drug_hits)} drugs con firmas '{target_direction}'")
    return all_drug_hits


# ---------------------------------------------------------------------------
# 1. Cargar genes DE
# ---------------------------------------------------------------------------
DE_FILE = "results/tables/de_limma/02_TVsS_significant_with_ids.tsv"
log.info(f"Cargando {DE_FILE}")

de = pd.read_csv(DE_FILE, sep="\t")
log.info(f"  {de.shape[0]} proteínas significativas")

up_genes   = (de[de["logFC_TVsS"] > 0]
              .nlargest(TOP_N_GENES, "logFC_TVsS")["gene_symbol"].tolist())
down_genes = (de[de["logFC_TVsS"] < 0]
              .nsmallest(TOP_N_GENES, "logFC_TVsS")["gene_symbol"].tolist())

log.info(f"  Top {len(up_genes)} UP genes:   {up_genes[:5]} ...")
log.info(f"  Top {len(down_genes)} DOWN genes: {down_genes[:5]} ...")

# ---------------------------------------------------------------------------
# 2. Queries a L2S2
# ---------------------------------------------------------------------------
# Query A: genes UP → busca firmas DOWN del fármaco (reversor: drug baja lo que tumor sube)
log.info("\n--- Query A: UP genes → drug DOWN signatures ---")
hits_down = query_l2s2_paginated(up_genes,   target_direction="down", label="UP→DOWN")

# Query B: genes DOWN → busca firmas UP del fármaco (reversor: drug sube lo que tumor baja)
log.info("\n--- Query B: DOWN genes → drug UP signatures ---")
hits_up   = query_l2s2_paginated(down_genes, target_direction="up",   label="DOWN→UP")

# ---------------------------------------------------------------------------
# 3. Agregar por fármaco
# ---------------------------------------------------------------------------
log.info("\nAgregando resultados por fármaco...")

all_drugs = set(hits_down.keys()) | set(hits_up.keys())
log.info(f"  Drugs totales con ≥1 hit reversal: {len(all_drugs)}")

records = []
for drug in all_drugs:
    h_down = hits_down.get(drug, [])   # firmas DOWN (reversor de UP)
    h_up   = hits_up.get(drug, [])     # firmas UP   (reversor de DOWN)
    h_all  = h_down + h_up

    n_down      = len(h_down)
    n_up        = len(h_up)
    n_total     = n_down + n_up
    consistent  = (n_down > 0) and (n_up > 0)  # revierte AMBAS direcciones

    best_pvalue = min(h["pvalue"] for h in h_all)
    or_vals     = [h["oddsRatio"] for h in h_all if h["oddsRatio"] > 0]
    mean_log2OR = (sum(math.log2(v) for v in or_vals) / len(or_vals)) if or_vals else 0.0

    cell_lines  = list({h["cell_line"] for h in h_all})
    n_cell_lines = len(cell_lines)

    # Score de reversión combinado
    # Componentes: -log10(best_p) * mean_log2(OR) * log2(n_total+1)
    # Bonus ×1.5 si hay reversión consistente en ambas direcciones
    log10p = -math.log10(best_pvalue) if best_pvalue > 0 else 30.0
    reversal_raw = log10p * max(mean_log2OR, 0) * math.log2(n_total + 1)
    if consistent:
        reversal_raw *= 1.5  # bonus multiplier for consistent directional reversal across datasets

    records.append({
        "drug_name":       drug,
        "pert":            drug.upper(),
        "type":            "trt_cp",          # campo requerido por script 08
        "n_down_sigs":     n_down,
        "n_up_sigs":       n_up,
        "n_total_sigs":    n_total,
        "consistent_rev":  consistent,
        "best_pvalue":     best_pvalue,
        "mean_log2_OR":    round(mean_log2OR, 4),
        "n_cell_lines":    n_cell_lines,
        "cell_lines":      ";".join(sorted(cell_lines)[:10]),
        "reversal_raw":    reversal_raw,
    })

df_results = pd.DataFrame(records)

# ---------------------------------------------------------------------------
# 4. Normalizar score → scaled_score ∈ [-1, 0] (negativo = reversor, como CMap2)
# ---------------------------------------------------------------------------
max_raw = df_results["reversal_raw"].max()
if max_raw > 0:
    df_results["scaled_score"] = -(df_results["reversal_raw"] / max_raw)
else:
    df_results["scaled_score"] = 0.0

df_results = df_results.sort_values("scaled_score")  # más negativo = mejor reversor

log.info(f"  Score máximo (raw):  {max_raw:.2f}")
log.info(f"  scaled_score range:  [{df_results['scaled_score'].min():.3f}, "
         f"{df_results['scaled_score'].max():.3f}]")
log.info(f"  Consistent reversors (ambas dir): "
         f"{df_results['consistent_rev'].sum()} / {len(df_results)}")

# ---------------------------------------------------------------------------
# 5. Exportar
# ---------------------------------------------------------------------------
df_results.insert(0, "query_date", QUERY_DATE)

# Archivo completo
out_all = "results/tables/drug_targets/07_l2s2_results.tsv"
df_results.to_csv(out_all, sep="\t", index=False)
log.info(f"\nGuardado: {out_all}  ({len(df_results)} drugs)")

# Top reversores (para integración con script 08)
df_top = df_results[df_results["scaled_score"] < 0].copy()
out_top = "results/tables/drug_targets/07_l2s2_top_reversors.tsv"
df_top.to_csv(out_top, sep="\t", index=False)
log.info(f"Guardado: {out_top}  ({len(df_top)} reversores)")

# Summary en log
log.info("\n=== TOP 20 REVERSORES L2S2 ===")
log.info(f"{'Drug':40s}  {'score':>8}  {'n_sigs':>6}  {'n_cell':>6}  {'consistent':>10}  {'best_p':>12}")
for _, row in df_top.head(20).iterrows():
    log.info(f"{row['drug_name']:40s}  {row['scaled_score']:8.4f}  "
             f"{row['n_total_sigs']:6d}  {row['n_cell_lines']:6d}  "
             f"{'YES' if row['consistent_rev'] else 'no':>10}  "
             f"{row['best_pvalue']:12.3e}")

log.info(f"\nFin: {datetime.now()}")
log.info("=== 07_l2s2_connectivity.py COMPLETADO ===")
