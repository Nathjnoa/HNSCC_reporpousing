#!/usr/bin/env python3
"""
Script 06: Query Open Targets Platform GraphQL API
HNSCC Drug Repurposing — Fase 3

Estrategia:
  1. Paginar disease.associatedTargets(EFO_0000181) -> mapa symbol->ensembl_id+score_HNSCC
  2. Filtrar targets que son nuestros genes DE (overlap por symbol)
  3. Para cada gen matched: query target.knownDrugs (evidencia HNSCC + reposicionamiento)

Input:  results/tables/de_limma/02_TVsS_significant_with_ids.tsv
Output: results/tables/drug_targets/06_opentargets_gene_drugs.tsv
        results/tables/drug_targets/06_opentargets_hnscc_scores.tsv
"""

import sys
import time
import logging
import os
from datetime import datetime

import requests
import pandas as pd

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
os.makedirs("logs", exist_ok=True)
os.makedirs("results/tables/drug_targets", exist_ok=True)

log_file = f"logs/06_query_opentargets_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

log.info("=== 06_query_opentargets.py ===")
log.info(f"Inicio: {datetime.now()}")

# ---------------------------------------------------------------------------
# Parametros
# ---------------------------------------------------------------------------
OT_URL      = "https://api.platform.opentargets.org/api/v4/graphql"
HNSCC_EFO   = "EFO_0000181"
PAGE_SIZE   = 500   # targets por pagina (max ~500 para evitar timeout)
SLEEP_SEC   = 0.5
TIMEOUT     = 60

# ---------------------------------------------------------------------------
# Helper de consulta
# ---------------------------------------------------------------------------
def gql(query: str, variables: dict = None) -> dict | None:
    payload = {"query": query}
    if variables:
        payload["variables"] = variables
    try:
        resp = requests.post(
            OT_URL,
            json=payload,
            timeout=TIMEOUT,
            headers={"Content-Type": "application/json"},
        )
        resp.raise_for_status()
        data = resp.json()
        if "errors" in data:
            log.warning(f"  GQL errors: {data['errors'][:1]}")
        return data.get("data")
    except requests.exceptions.RequestException as e:
        log.warning(f"  Request error: {e}")
        return None

# ---------------------------------------------------------------------------
# Cargar genes DE
# ---------------------------------------------------------------------------
sig = pd.read_csv(
    "results/tables/de_limma/02_TVsS_significant_with_ids.tsv", sep="\t"
)
gene_meta = sig[["symbol_org", "gene_symbol", "uniprot_id", "logFC_TVsS", "adj.P.Val_TVsS"]].copy()
gene_meta["direction"] = gene_meta["logFC_TVsS"].apply(lambda x: "up" if x > 0 else "down")

# Conjunto de simbolos (upper) para matching
de_symbols_upper = {str(s).upper() for s in sig["symbol_org"].dropna()}
de_symbols_upper |= {str(s).upper() for s in sig["gene_symbol"].dropna()}
log.info(f"Genes DE a buscar: {len(sig)} proteinas | {len(de_symbols_upper)} simbolos unicos")

# ---------------------------------------------------------------------------
# Paso 1: Paginar HNSCC associated targets -> mapa symbol -> {ensembl, score}
# ---------------------------------------------------------------------------
log.info(f"\n--- Paso 1: Paginando targets asociados a {HNSCC_EFO} (HNSCC) ---")

QUERY_ASSOC = """
query HNSCCTargets($efo: String!, $page: Int!, $size: Int!) {
  disease(efoId: $efo) {
    associatedTargets(page: {index: $page, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
        }
        score
      }
    }
  }
}
"""

# Obtener total primero
data0 = gql(QUERY_ASSOC, {"efo": HNSCC_EFO, "page": 0, "size": 1})
total = data0["disease"]["associatedTargets"]["count"]
n_pages = (total // PAGE_SIZE) + 1
log.info(f"Total targets HNSCC: {total} | paginas: {n_pages}")

hnscc_target_map: dict[str, dict] = {}  # symbol_upper -> {ensembl_id, score}

for page in range(n_pages):
    data = gql(QUERY_ASSOC, {"efo": HNSCC_EFO, "page": page, "size": PAGE_SIZE})
    if not data:
        log.warning(f"  Pagina {page} fallida")
        time.sleep(SLEEP_SEC)
        continue
    rows = data["disease"]["associatedTargets"]["rows"]
    for r in rows:
        sym = r["target"]["approvedSymbol"].upper()
        hnscc_target_map[sym] = {
            "ensembl_id":   r["target"]["id"],
            "symbol":       r["target"]["approvedSymbol"],
            "hnscc_score":  r["score"],
        }
    if (page + 1) % 5 == 0:
        log.info(f"  Pagina {page+1}/{n_pages} | targets en mapa: {len(hnscc_target_map)}")
    time.sleep(SLEEP_SEC)

log.info(f"Targets HNSCC en mapa: {len(hnscc_target_map)}")

# ---------------------------------------------------------------------------
# Paso 2: Intersectar con nuestros genes DE
# ---------------------------------------------------------------------------
log.info("\n--- Paso 2: Interseccion con genes DE ---")

matched: list[dict] = []
for _, row in sig.iterrows():
    sym_up = str(row["symbol_org"]).upper() if pd.notna(row["symbol_org"]) else ""
    sym_up2 = str(row["gene_symbol"]).upper() if pd.notna(row["gene_symbol"]) else ""
    hit = hnscc_target_map.get(sym_up) or hnscc_target_map.get(sym_up2)
    if hit:
        matched.append({
            "uniprot_id":   row["uniprot_id"],
            "symbol_org":   row["symbol_org"],
            "gene_symbol":  row["gene_symbol"],
            "logFC_TVsS":   row["logFC_TVsS"],
            "adj_pval":     row["adj.P.Val_TVsS"],
            "direction":    "up" if row["logFC_TVsS"] > 0 else "down",
            "ensembl_id":   hit["ensembl_id"],
            "hnscc_score":  hit["hnscc_score"],
        })

df_matched = pd.DataFrame(matched)
log.info(f"Genes DE con evidencia HNSCC en Open Targets: {len(df_matched)} / {len(sig)}")

# Exportar scores HNSCC
df_matched.to_csv(
    "results/tables/drug_targets/06_opentargets_hnscc_scores.tsv",
    sep="\t", index=False,
)
log.info(f"Exportado: 06_opentargets_hnscc_scores.tsv")

# ---------------------------------------------------------------------------
# Paso 3: knownDrugs para TODOS los genes DE (con Ensembl ID)
# — Incluye genes con y sin score HNSCC (para capturar reposicionamiento novel)
# ---------------------------------------------------------------------------
log.info("\n--- Paso 3: Consultando knownDrugs para genes DE ---")

# Para genes sin match en HNSCC, buscar Ensembl ID via search
QUERY_SEARCH = """
query SearchGene($q: String!) {
  search(queryString: $q, entityNames: ["target"]) {
    hits {
      id
      object {
        ... on Target {
          approvedSymbol
        }
      }
    }
  }
}
"""

# Construir mapa completo symbol -> ensembl_id (matched + search para no-matched)
full_ensembl_map: dict[str, str] = {
    r["symbol_org"]: r["ensembl_id"] for _, r in df_matched.iterrows()
    if pd.notna(r.get("symbol_org"))
}

# Para genes no matched, intentar busqueda directa
unmatched_genes = [
    (row["symbol_org"] or row["gene_symbol"])
    for _, row in sig.iterrows()
    if str(row.get("symbol_org", "")).upper() not in {k.upper() for k in full_ensembl_map}
    and str(row.get("gene_symbol", "")).upper() not in {k.upper() for k in full_ensembl_map}
    and pd.notna(row.get("symbol_org") or row.get("gene_symbol"))
]
unmatched_genes = list(dict.fromkeys(unmatched_genes))  # unique, ordered

log.info(f"Genes ya con Ensembl ID: {len(full_ensembl_map)}")
log.info(f"Genes a buscar Ensembl via search: {len(unmatched_genes)}")

for i, sym in enumerate(unmatched_genes):
    data = gql(QUERY_SEARCH, {"q": sym})
    if data and data.get("search", {}).get("hits"):
        for hit in data["search"]["hits"]:
            obj = hit.get("object", {})
            if obj.get("approvedSymbol", "").upper() == sym.upper():
                full_ensembl_map[sym] = hit["id"]
                break
    time.sleep(0.2)
    if (i + 1) % 50 == 0:
        log.info(f"  Busqueda Ensembl: {i+1}/{len(unmatched_genes)}")

log.info(f"Total genes con Ensembl ID: {len(full_ensembl_map)}")

# --- Consultar knownDrugs ---
QUERY_DRUGS = """
query TargetDrugs($ensembl: String!) {
  target(ensemblId: $ensembl) {
    id
    approvedSymbol
    knownDrugs {
      rows {
        drug {
          id
          name
          maximumClinicalTrialPhase
          isApproved
          drugType
          mechanismsOfAction {
            rows {
              actionType
              mechanismOfAction
            }
          }
        }
        disease {
          id
          name
        }
        phase
        mechanismOfAction
        urls {
          url
        }
      }
    }
  }
}
"""

all_drug_rows: list[dict] = []
genes_with_drugs = 0
genes_no_drugs = 0

symbol_to_meta = {
    row["symbol_org"]: row for _, row in sig.iterrows()
}

items = list(full_ensembl_map.items())
log.info(f"Consultando knownDrugs para {len(items)} genes ...")

for i, (sym, ensembl) in enumerate(items):
    data = gql(QUERY_DRUGS, {"ensembl": ensembl})
    if not data or not data.get("target"):
        time.sleep(SLEEP_SEC)
        continue

    tgt = data["target"]
    known = tgt.get("knownDrugs") or {}
    rows = known.get("rows") or []

    if not rows:
        genes_no_drugs += 1
    else:
        genes_with_drugs += 1

    # Metadatos expresion
    meta = symbol_to_meta.get(sym, {})
    logfc    = meta.get("logFC_TVsS", "")
    adj_pval = meta.get("adj.P.Val_TVsS", "")
    direction = "up" if (isinstance(logfc, float) and logfc > 0) else "down"
    hnscc_score = ""
    if sym.upper() in {k.upper() for k in df_matched["symbol_org"].dropna()}:
        match_row = df_matched[df_matched["symbol_org"].str.upper() == sym.upper()]
        if not match_row.empty:
            hnscc_score = match_row.iloc[0]["hnscc_score"]

    for r in rows:
        drug = r.get("drug") or {}
        disease = r.get("disease") or {}
        is_hnscc = disease.get("id", "") == HNSCC_EFO
        all_drug_rows.append({
            "symbol_org":           sym,
            "ensembl_id":           ensembl,
            "logFC_TVsS":           logfc,
            "adj_pval":             adj_pval,
            "direction":            direction,
            "hnscc_ot_score":       hnscc_score,
            "drug_id":              drug.get("id", ""),
            "drug_name":            (drug.get("name") or "").strip().upper(),
            "max_phase_global":     drug.get("maximumClinicalTrialPhase", ""),
            "is_approved":          drug.get("isApproved", ""),
            "drug_type":            drug.get("drugType", ""),
            "trial_phase":          r.get("phase", ""),
            "indication_efo":       disease.get("id", ""),
            "indication_name":      disease.get("name", ""),
            "is_hnscc_indication":  is_hnscc,
            "mechanism_of_action":  r.get("mechanismOfAction", ""),
        })

    time.sleep(SLEEP_SEC)
    if (i + 1) % 50 == 0:
        log.info(
            f"  Procesados {i+1}/{len(items)} | "
            f"con drugs: {genes_with_drugs} | sin drugs: {genes_no_drugs}"
        )

log.info(f"\nGenes con drogas: {genes_with_drugs}")
log.info(f"Genes sin drogas: {genes_no_drugs}")
log.info(f"Pares gen-farmaco totales: {len(all_drug_rows)}")

# ---------------------------------------------------------------------------
# Construir DataFrame y exportar
# ---------------------------------------------------------------------------
if not all_drug_rows:
    log.error("Sin resultados de knownDrugs. Revisar conectividad.")
    sys.exit(1)

df_drugs = pd.DataFrame(all_drug_rows)
df_drugs["max_phase_global"] = pd.to_numeric(df_drugs["max_phase_global"], errors="coerce")
df_drugs["trial_phase"]      = pd.to_numeric(df_drugs["trial_phase"],      errors="coerce")

# Deduplicar por gen+farmaco+indicacion
df_dedup = (
    df_drugs
    .sort_values("max_phase_global", ascending=False)
    .drop_duplicates(subset=["symbol_org", "drug_id", "indication_efo"])
    .reset_index(drop=True)
)

log.info(f"Filas dedup (gen+drug+indicacion): {len(df_dedup)}")
log.info(f"Farmacos unicos: {df_dedup['drug_id'].nunique()}")
log.info(f"Genes unicos:    {df_dedup['symbol_org'].nunique()}")

# --- Estadisticas clave ---
hnscc_direct = df_dedup[df_dedup["is_hnscc_indication"] == True]
log.info(f"\nPares gen-farmaco con indicacion HNSCC directa: {len(hnscc_direct)}")
log.info(f"Farmacos con indicacion HNSCC directa: {hnscc_direct['drug_id'].nunique()}")

approved = df_dedup[df_dedup["is_approved"] == True]
log.info(f"Farmacos aprobados (cualquier indicacion): {approved['drug_id'].nunique()}")

novel = df_dedup[
    (df_dedup["is_approved"] == True) &
    (~df_dedup["drug_id"].isin(hnscc_direct["drug_id"]))
]
log.info(f"Candidatos reposicionamiento novel (aprobados, no-HNSCC): {novel['drug_id'].nunique()}")

log.info(f"\nTop 10 farmacos aprobados por n_genes_DE:")
top_drugs = (
    df_dedup[df_dedup["is_approved"] == True]
    .groupby(["drug_name", "drug_id"])["symbol_org"]
    .nunique()
    .reset_index()
    .rename(columns={"symbol_org": "n_de_genes"})
    .sort_values("n_de_genes", ascending=False)
    .head(10)
)
for _, r in top_drugs.iterrows():
    hnscc_flag = " [HNSCC]" if r["drug_id"] in hnscc_direct["drug_id"].values else ""
    log.info(f"  {r['drug_name']:<40} n_genes={r['n_de_genes']}{hnscc_flag}")

# Exportar
out_drugs = "results/tables/drug_targets/06_opentargets_gene_drugs.tsv"
df_dedup.to_csv(out_drugs, sep="\t", index=False)
log.info(f"\nExportado: {out_drugs}  ({len(df_dedup)} filas)")
log.info(f"\nSiguiente: scripts/07_cmap_connectivity.R")
log.info(f"Fin: {datetime.now()}")
