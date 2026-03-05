#!/usr/bin/env python3
"""
Script 05: Query ChEMBL REST API — farmacos fase >= 3 por target proteomico
HNSCC Drug Repurposing — Fase 3

Flujo (REST puro, sin libreria chembl):
  UniProt IDs
    -> ChEMBL target ID (target_components__accession, SINGLE PROTEIN, Homo sapiens)
    -> Drug mechanisms (target -> mechanism -> molecule)
    -> Molecule details (max_phase, nombre, tipo)
    -> Filtro: max_phase >= 3

Input:  results/tables/de_limma/02_TVsS_significant_with_ids.tsv
Output: results/tables/drug_targets/05_chembl_drugs.tsv
        results/tables/drug_targets/05_chembl_summary.tsv
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

log_file = f"logs/05_query_chembl_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

log.info("=== 05_query_chembl.py (REST, sin libreria) ===")
log.info(f"Inicio: {datetime.now()}")

# ---------------------------------------------------------------------------
# Parametros
# ---------------------------------------------------------------------------
BASE   = "https://www.ebi.ac.uk/chembl/api/data"
MIN_PHASE = 3
TIMEOUT   = 20    # segundos
SLEEP_SEC = 0.3
BATCH     = 50    # moleculas por request

session = requests.Session()
session.headers.update({"Accept": "application/json"})

def get_json(url: str, params: dict = None, retries: int = 2) -> dict | None:
    for attempt in range(retries + 1):
        try:
            r = session.get(url, params=params, timeout=TIMEOUT)
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if attempt < retries:
                time.sleep(1)
            else:
                log.debug(f"    GET fallido {url}: {e}")
    return None

# ---------------------------------------------------------------------------
# Cargar datos
# ---------------------------------------------------------------------------
sig = pd.read_csv(
    "results/tables/de_limma/02_TVsS_significant_with_ids.tsv", sep="\t"
)
uniprot_ids = sig["uniprot_id"].dropna().unique().tolist()
log.info(f"Proteinas significativas: {len(sig)} | UniProt IDs: {len(uniprot_ids)}")

gene_meta = (
    sig[["uniprot_id", "symbol_org", "gene_symbol", "logFC_TVsS", "adj.P.Val_TVsS"]]
    .assign(direction=lambda df: df["logFC_TVsS"].apply(lambda x: "up" if x > 0 else "down"))
)

# ---------------------------------------------------------------------------
# Paso 1: UniProt → ChEMBL target (REST)
# ---------------------------------------------------------------------------
log.info("\n--- Paso 1: UniProt -> ChEMBL target ---")

uniprot_to_target: list[dict] = []
no_target: list[str] = []

for i, uniprot in enumerate(uniprot_ids):
    params = {
        "target_components__accession": uniprot,
        "target_type": "SINGLE PROTEIN",
        "organism": "Homo sapiens",
        "format": "json",
        "limit": 5,
    }
    data = get_json(f"{BASE}/target.json", params)
    if data:
        for t in data.get("targets", []):
            uniprot_to_target.append({
                "uniprot_id":       uniprot,
                "target_chembl_id": t["target_chembl_id"],
                "target_pref_name": t.get("pref_name", ""),
            })
    else:
        no_target.append(uniprot)

    time.sleep(SLEEP_SEC)
    if (i + 1) % 50 == 0:
        log.info(f"  {i+1}/{len(uniprot_ids)} | targets encontrados: {len(uniprot_to_target)}")

df_targets = pd.DataFrame(uniprot_to_target) if uniprot_to_target else pd.DataFrame(
    columns=["uniprot_id", "target_chembl_id", "target_pref_name"]
)
target_ids = df_targets["target_chembl_id"].unique().tolist()

log.info(f"UniProt con target ChEMBL: {df_targets['uniprot_id'].nunique()} / {len(uniprot_ids)}")
log.info(f"Targets ChEMBL unicos:     {len(target_ids)}")

# ---------------------------------------------------------------------------
# Paso 2: Targets → Drug mechanisms
# ---------------------------------------------------------------------------
log.info("\n--- Paso 2: Targets -> mecanismos ---")

all_mechs: list[dict] = []

for i, target_id in enumerate(target_ids):
    params = {
        "target_chembl_id": target_id,
        "format": "json",
        "limit": 100,
    }
    data = get_json(f"{BASE}/mechanism.json", params)
    if data:
        for m in data.get("mechanisms", []):
            mol_id = m.get("molecule_chembl_id", "")
            if mol_id:
                all_mechs.append({
                    "target_chembl_id":   target_id,
                    "molecule_chembl_id": mol_id,
                    "mechanism_of_action": m.get("mechanism_of_action", ""),
                    "action_type":         m.get("action_type", ""),
                })
    time.sleep(SLEEP_SEC)
    if (i + 1) % 50 == 0:
        log.info(f"  {i+1}/{len(target_ids)} targets | mechs: {len(all_mechs)}")

df_mechs = pd.DataFrame(all_mechs) if all_mechs else pd.DataFrame(
    columns=["target_chembl_id", "molecule_chembl_id", "mechanism_of_action", "action_type"]
)
log.info(f"Mecanismos totales: {len(df_mechs)} | moleculas unicas: {df_mechs['molecule_chembl_id'].nunique()}")

# ---------------------------------------------------------------------------
# Paso 3: Molecule details (batch por __in)
# ---------------------------------------------------------------------------
log.info("\n--- Paso 3: Detalles de moleculas ---")

mol_ids = df_mechs["molecule_chembl_id"].dropna().unique().tolist()
log.info(f"Moleculas a consultar: {len(mol_ids)}")

mol_records: list[dict] = []

for i in range(0, len(mol_ids), BATCH):
    batch = mol_ids[i:i + BATCH]
    ids_str = ",".join(batch)
    params = {
        "molecule_chembl_id__in": ids_str,
        "format": "json",
        "limit": BATCH,
    }
    data = get_json(f"{BASE}/molecule.json", params)
    if data:
        for m in data.get("molecules", []):
            props = m.get("molecule_properties") or {}
            mol_records.append({
                "molecule_chembl_id": m.get("molecule_chembl_id", ""),
                "pref_name":          m.get("pref_name", ""),
                "max_phase":          m.get("max_phase"),
                "molecule_type":      m.get("molecule_type", ""),
                "first_approval":     m.get("first_approval", ""),
                "oral":               m.get("oral", ""),
                "indication_class":   m.get("indication_class", ""),
                "mw_freebase":        props.get("mw_freebase", ""),
                "alogp":              props.get("alogp", ""),
            })
    time.sleep(SLEEP_SEC)

df_mols = pd.DataFrame(mol_records) if mol_records else pd.DataFrame(
    columns=["molecule_chembl_id", "pref_name", "max_phase"]
)
df_mols["max_phase"] = pd.to_numeric(df_mols["max_phase"], errors="coerce")
log.info(f"Moleculas recuperadas: {len(df_mols)}")

# ---------------------------------------------------------------------------
# Paso 4: Integrar y filtrar
# ---------------------------------------------------------------------------
log.info("\n--- Paso 4: Integrando y filtrando fase >= 3 ---")

if df_mechs.empty or df_targets.empty:
    log.error("Sin datos de mecanismos o targets. Revisar pasos 1-2.")
    # Exportar vacios para no bloquear pipeline
    pd.DataFrame().to_csv("results/tables/drug_targets/05_chembl_drugs.tsv",  sep="\t", index=False)
    pd.DataFrame().to_csv("results/tables/drug_targets/05_chembl_summary.tsv", sep="\t", index=False)
    sys.exit(0)

df_combined = (
    df_mechs
    .merge(df_targets, on="target_chembl_id", how="left")
    .merge(df_mols, on="molecule_chembl_id", how="left")
    .merge(gene_meta, on="uniprot_id", how="left")
)

df_phase3 = df_combined[
    df_combined["max_phase"].notna() & (df_combined["max_phase"] >= MIN_PHASE)
].copy()

def phase_label(p):
    if p == 4:   return "Approved"
    elif p == 3: return "Phase III"
    else:         return f"Phase {int(p)}"

df_phase3["phase_label"] = df_phase3["max_phase"].apply(phase_label)

df_out = (
    df_phase3
    .sort_values("max_phase", ascending=False)
    .drop_duplicates(subset=["uniprot_id", "molecule_chembl_id"])
    .reset_index(drop=True)
)

out_cols = [
    "symbol_org", "gene_symbol", "uniprot_id",
    "target_chembl_id", "target_pref_name",
    "molecule_chembl_id", "pref_name",
    "max_phase", "phase_label", "first_approval",
    "mechanism_of_action", "action_type",
    "molecule_type", "oral", "indication_class",
    "mw_freebase", "alogp",
    "logFC_TVsS", "direction", "adj.P.Val_TVsS",
]
out_cols = [c for c in out_cols if c in df_out.columns]
df_out = df_out[out_cols]

log.info(f"Interacciones totales (todas fases): {len(df_combined)}")
log.info(f"Interacciones fase >= {MIN_PHASE}:   {len(df_phase3)}")
log.info(f"Tras dedup (gen+mol):                {len(df_out)}")
log.info(f"Farmacos unicos fase >= {MIN_PHASE}: {df_out['molecule_chembl_id'].nunique()}")
log.info(f"Genes con farmaco fase >= {MIN_PHASE}: {df_out['uniprot_id'].nunique()}")

# Resumen por farmaco
summary = (
    df_out.groupby(
        ["molecule_chembl_id", "pref_name", "max_phase", "phase_label",
         "first_approval", "molecule_type", "indication_class"],
        dropna=False
    )
    .agg(
        n_targets=(   "uniprot_id",   "nunique"),
        target_genes=("symbol_org",   lambda x: "|".join(sorted(x.dropna().unique()))),
        n_up=(        "direction",    lambda x: (x == "up").sum()),
        n_down=(      "direction",    lambda x: (x == "down").sum()),
        mechanisms=(  "mechanism_of_action", lambda x: "|".join(sorted(x.dropna().unique()))),
    )
    .reset_index()
    .sort_values(["max_phase", "n_targets"], ascending=[False, False])
)

log.info(f"\n--- Por fase ---")
for phase, grp in summary.groupby("phase_label"):
    log.info(f"  {phase:<12} {len(grp):4d} farmacos")

log.info(f"\nTop 10 por n_targets (fase >= 3):")
for _, r in summary.head(10).iterrows():
    name = r["pref_name"] or r["molecule_chembl_id"]
    log.info(
        f"  {name:<35} fase={r['max_phase']}  "
        f"targets={r['n_targets']}  up={r['n_up']}  down={r['n_down']}"
    )

# Exportar
df_out.to_csv( "results/tables/drug_targets/05_chembl_drugs.tsv",   sep="\t", index=False)
summary.to_csv("results/tables/drug_targets/05_chembl_summary.tsv", sep="\t", index=False)

log.info(f"\nExportado: 05_chembl_drugs.tsv   ({len(df_out)} filas)")
log.info(f"Exportado: 05_chembl_summary.tsv ({len(summary)} farmacos)")
log.info(f"\nSiguiente: scripts/06_query_opentargets.py")
log.info(f"Fin: {datetime.now()}")
