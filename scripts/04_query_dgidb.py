#!/usr/bin/env python3
"""
Script 04: Query DGIdb GraphQL API (v5)
HNSCC Drug Repurposing — Fase 3

Input:  results/tables/de_limma/02_TVsS_significant_with_ids.tsv
Output: results/tables/drug_targets/04_dgidb_raw.tsv
        results/tables/drug_targets/04_dgidb_summary.tsv
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

log_file = f"logs/04_query_dgidb_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

log.info("=== 04_query_dgidb.py ===")
log.info(f"Inicio: {datetime.now()}")

# ---------------------------------------------------------------------------
# Parametros
# ---------------------------------------------------------------------------
DGIDB_GRAPHQL = "https://dgidb.org/api/graphql"
BATCH_SIZE    = 50   # genes por query (GraphQL)
SLEEP_SEC     = 1.0
TIMEOUT       = 60

# ---------------------------------------------------------------------------
# Cargar genes
# ---------------------------------------------------------------------------
sig = pd.read_csv(
    "results/tables/de_limma/02_TVsS_significant_with_ids.tsv", sep="\t"
)
genes = (
    sig["symbol_org"]
    .fillna(sig["gene_symbol"])
    .dropna()
    .unique()
    .tolist()
)
log.info(f"Genes a consultar: {len(genes)}")

gene_meta = (
    sig[["symbol_org", "logFC_TVsS", "adj.P.Val_TVsS"]]
    .rename(columns={
        "symbol_org":     "query_gene",
        "logFC_TVsS":     "logFC",
        "adj.P.Val_TVsS": "adj_pval",
    })
    .assign(direction=lambda df: df["logFC"].apply(lambda x: "up" if x > 0 else "down"))
)

# ---------------------------------------------------------------------------
# GraphQL query
# ---------------------------------------------------------------------------
GQL_QUERY = """
query GeneInteractions($names: [String!]!) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        drug {
          name
          conceptId
        }
        interactionScore
        interactionTypes {
          type
          directionality
        }
        sources {
          sourceDbName
        }
        publications {
          pmid
        }
      }
    }
  }
}
"""


def query_dgidb_batch(gene_list: list[str]) -> list[dict] | None:
    """GraphQL POST para una lista de genes. Retorna lista de nodos o None."""
    payload = {
        "query": GQL_QUERY,
        "variables": {"names": gene_list},
    }
    try:
        resp = requests.post(
            DGIDB_GRAPHQL,
            json=payload,
            timeout=TIMEOUT,
            headers={"Content-Type": "application/json"},
        )
        resp.raise_for_status()
        data = resp.json()
        if "errors" in data:
            log.warning(f"  GraphQL errors: {data['errors']}")
        return data.get("data", {}).get("genes", {}).get("nodes", [])
    except requests.exceptions.RequestException as e:
        log.warning(f"  Request error: {e}")
        return None


# ---------------------------------------------------------------------------
# Consulta en batches
# ---------------------------------------------------------------------------
batches = [genes[i:i + BATCH_SIZE] for i in range(0, len(genes), BATCH_SIZE)]
log.info(f"Batches: {len(batches)} x ~{BATCH_SIZE} genes")

all_interactions: list[dict] = []
matched_genes: set[str] = set()
empty_genes: list[str] = []

for i, batch in enumerate(batches, 1):
    log.info(f"  Batch {i}/{len(batches)} ({len(batch)} genes) ...")
    nodes = query_dgidb_batch(batch)

    if nodes is None:
        log.warning(f"  Batch {i} fallido — omitido")
        time.sleep(SLEEP_SEC)
        continue

    # genes retornados sin interacciones
    returned_names = {n["name"] for n in nodes}
    for g in batch:
        if g.upper() not in {r.upper() for r in returned_names}:
            empty_genes.append(g)

    for node in nodes:
        gene_name = node["name"]
        ixs = node.get("interactions") or []
        if not ixs:
            empty_genes.append(gene_name)
            continue
        matched_genes.add(gene_name)
        for ix in ixs:
            drug = ix.get("drug") or {}
            drug_name = (drug.get("name") or "").strip().upper()
            if not drug_name:
                continue

            concept_id = drug.get("conceptId", "")
            # extraer ChEMBL ID si está disponible en conceptId
            chembl_id = concept_id if concept_id.startswith("chembl:") else ""

            itypes = "|".join(
                t["type"] for t in (ix.get("interactionTypes") or []) if t.get("type")
            )
            idirs = "|".join(
                t["directionality"]
                for t in (ix.get("interactionTypes") or [])
                if t.get("directionality")
            )
            sources = "|".join(
                s["sourceDbName"] for s in (ix.get("sources") or []) if s.get("sourceDbName")
            )
            pmids = "|".join(
                str(p["pmid"]) for p in (ix.get("publications") or []) if p.get("pmid")
            )

            all_interactions.append({
                "query_gene":        gene_name,
                "drug_name":         drug_name,
                "drug_concept_id":   concept_id,
                "drug_chembl_id":    chembl_id,
                "interaction_score": ix.get("interactionScore", ""),
                "interaction_types": itypes,
                "interaction_dirs":  idirs,
                "sources":           sources,
                "n_sources":         len([s for s in sources.split("|") if s]),
                "pmids":             pmids,
            })

    time.sleep(SLEEP_SEC)

log.info(f"\nGenes con interacciones:       {len(matched_genes)}")
log.info(f"Genes sin interacciones/mapeo: {len(set(empty_genes))}")
log.info(f"Interacciones totales crudas:  {len(all_interactions)}")

# ---------------------------------------------------------------------------
# DataFrame principal
# ---------------------------------------------------------------------------
if not all_interactions:
    log.error("No se recuperaron interacciones. Revisar conectividad.")
    sys.exit(1)

df = pd.DataFrame(all_interactions)

# Join metadata de expresion
df = df.merge(gene_meta, left_on="query_gene", right_on="query_gene", how="left")

# Deduplicar por gen+farmaco (conservar fila con mayor interaction_score)
df["interaction_score"] = pd.to_numeric(df["interaction_score"], errors="coerce").fillna(0)
df_dedup = (
    df.sort_values("interaction_score", ascending=False)
    .drop_duplicates(subset=["query_gene", "drug_name"])
    .reset_index(drop=True)
)

log.info(f"Interacciones dedup (gene+drug): {len(df_dedup)}")
log.info(f"Farmacos unicos:                 {df_dedup['drug_name'].nunique()}")
log.info(f"Genes con interacciones:         {df_dedup['query_gene'].nunique()}")

# ---------------------------------------------------------------------------
# Tabla resumen por farmaco
# ---------------------------------------------------------------------------
summary = (
    df_dedup.groupby("drug_name")
    .agg(
        n_targets=(      "query_gene",       "nunique"),
        target_genes=(   "query_gene",       lambda x: "|".join(sorted(x.unique()))),
        n_up_targets=(   "direction",        lambda x: (x == "up").sum()),
        n_down_targets=( "direction",        lambda x: (x == "down").sum()),
        mean_logFC=(     "logFC",            "mean"),
        interaction_types=("interaction_types",
                           lambda x: "|".join(sorted({t for v in x for t in v.split("|") if t}))),
        max_interaction_score=("interaction_score", "max"),
        drug_concept_id=("drug_concept_id",  lambda x: next((v for v in x if v), "")),
        drug_chembl_id=( "drug_chembl_id",   lambda x: next((v for v in x if v), "")),
        max_n_sources=(  "n_sources",        "max"),
    )
    .reset_index()
    .sort_values("n_targets", ascending=False)
)

log.info(f"\nTop 10 farmacos por numero de targets:")
for _, r in summary.head(10).iterrows():
    log.info(
        f"  {r['drug_name']:<35} targets={r['n_targets']:3d}  "
        f"up={r['n_up_targets']:3d}  down={r['n_down_targets']:3d}  "
        f"score={r['max_interaction_score']:.3f}"
    )

# ---------------------------------------------------------------------------
# Exportar
# ---------------------------------------------------------------------------
out_raw     = "results/tables/drug_targets/04_dgidb_raw.tsv"
out_summary = "results/tables/drug_targets/04_dgidb_summary.tsv"
out_empty   = "results/tables/drug_targets/04_dgidb_no_interactions.txt"

df_dedup.to_csv(out_raw,     sep="\t", index=False)
summary.to_csv(out_summary,  sep="\t", index=False)
with open(out_empty, "w") as f:
    f.write("\n".join(sorted(set(empty_genes))))

log.info(f"\nExportado: {out_raw}  ({len(df_dedup)} filas)")
log.info(f"Exportado: {out_summary}  ({len(summary)} farmacos unicos)")
log.info(f"Exportado: {out_empty}")
log.info(f"\nSiguiente: scripts/05_query_chembl.py")
log.info(f"Fin: {datetime.now()}")
