#!/usr/bin/env python3
"""
11_clinicaltrials_pubmed.py
===========================
Para cada uno de los Top 20 candidatos, busca evidencia clínica y bibliográfica:

1. ClinicalTrials.gov API v2 — ensayos clínicos (drug + HNSCC)
2. PubMed via NCBI E-utilities — conteo de publicaciones (drug + HNSCC)

Outputs:
  results/tables/evidence/11_clinical_evidence.tsv   — tabla principal por candidato
  results/tables/evidence/11_trials_detail.tsv       — detalle de ensayos encontrados
  results/figures/evidence/11_evidence_bubble.pdf    — bubble chart evidencia
  results/figures/evidence/11_trials_phase_bar.pdf   — fases de ensayos por candidato

Uso:
  conda activate omics-py
  cd ~/bioinfo/projects/hnscc_drug_repurposing
  python scripts/11_clinicaltrials_pubmed.py
"""

import os
import sys
import time
import logging
import datetime
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

# ── Configuración ────────────────────────────────────────────────────────────
SCRIPT_DIR  = Path(__file__).resolve().parent
PROJ_DIR    = SCRIPT_DIR.parent
INPUT_TOP20 = PROJ_DIR / "results/tables/10_top20_candidates.tsv"
OUT_DIR     = PROJ_DIR / "results/tables/evidence"
FIG_DIR     = PROJ_DIR / "results/figures/evidence"
LOG_DIR     = PROJ_DIR / "logs"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_FILE  = LOG_DIR / f"11_clinicaltrials_pubmed_{TIMESTAMP}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

# ClinicalTrials.gov API v2 base URL
CT_BASE = "https://clinicaltrials.gov/api/v2/studies"

# NCBI E-utilities
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_EMAIL   = os.environ.get("NCBI_EMAIL", "jcarvajal@fucsalud.edu.co")

# HNSCC search terms
HNSCC_CT_TERMS = [
    "head and neck cancer",
    "head and neck squamous",
    "HNSCC",
    "squamous cell carcinoma head neck",
    "oral cancer",
    "oropharyngeal cancer",
]
HNSCC_PUBMED_TERM = ("(head neck squamous[MeSH Major Topic] OR HNSCC[Title/Abstract] "
                      "OR \"head and neck cancer\"[Title/Abstract] "
                      "OR \"squamous cell carcinoma\"[Title/Abstract])")

REQUEST_DELAY = 0.35   # segundos entre llamadas (NCBI permite ~3/s)

# ── Helpers ───────────────────────────────────────────────────────────────────

def get_json(url: str, params: dict, retries: int = 3, timeout: int = 30) -> dict | None:
    """GET con reintentos y delay."""
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except Exception as exc:
            log.warning(f"  Intento {attempt+1} fallido para {url}: {exc}")
            time.sleep(2 ** attempt)
    return None


def query_clinicaltrials(drug_name: str) -> dict:
    """
    Busca ensayos en ClinicalTrials.gov API v2 para 'drug_name' + HNSCC.
    Devuelve:
      total_trials: int
      hnscc_trials: int   (con término HNSCC en condition)
      by_phase: dict      {phase: count}
      trials: list[dict]  (detalle de ensayos encontrados)
    """
    all_trials = []
    page_token = None

    # Buscar por nombre de fármaco en intervention
    params_base = {
        "query.intr": drug_name,
        "pageSize": 100,
        "format": "json",
    }

    while True:
        params = params_base.copy()
        if page_token:
            params["pageToken"] = page_token

        data = get_json(CT_BASE, params)
        time.sleep(REQUEST_DELAY)

        if data is None:
            break

        studies = data.get("studies", [])
        for s in studies:
            proto = s.get("protocolSection", {})
            ident = proto.get("identificationModule", {})
            status = proto.get("statusModule", {})
            design = proto.get("designModule", {})
            cond   = proto.get("conditionsModule", {})
            inter  = proto.get("armsInterventionsModule", {})

            nct_id     = ident.get("nctId", "")
            title      = ident.get("briefTitle", "")
            phase_list = design.get("phases", [])
            phase      = ", ".join(phase_list) if phase_list else "N/A"
            overall_st = status.get("overallStatus", "")
            conditions = [c.lower() for c in cond.get("conditions", [])]

            # Comprobar si menciona HNSCC en condición (usa HNSCC_CT_TERMS)
            conditions_str = " ".join(conditions)
            is_hnscc = any(kw.lower() in conditions_str for kw in HNSCC_CT_TERMS)

            # Obtener nombre de intervención para confirmar el fármaco
            interventions = []
            for arm in inter.get("interventions", []):
                interventions.append(arm.get("name", "").lower())

            drug_confirmed = any(
                drug_name.lower().split()[0] in iv for iv in interventions
            ) if interventions else True   # si no hay detalle, asumir confirmado

            all_trials.append({
                "nct_id": nct_id,
                "title": title,
                "phase": phase,
                "status": overall_st,
                "is_hnscc": is_hnscc,
                "drug_confirmed": drug_confirmed,
            })

        next_token = data.get("nextPageToken")
        if not next_token:
            break
        page_token = next_token

    hnscc_trials = [t for t in all_trials if t["is_hnscc"] and t["drug_confirmed"]]
    by_phase = {}
    for t in hnscc_trials:
        p = t["phase"]
        by_phase[p] = by_phase.get(p, 0) + 1

    return {
        "total_trials": len(all_trials),
        "hnscc_trials": len(hnscc_trials),
        "by_phase": by_phase,
        "trials": hnscc_trials,
    }


SALT_SUFFIXES = [
    " hydrochloride", " hcl", " ditosylate", " mesylate", " tosylate",
    " sulfate", " fumarate", " maleate", " tartrate", " acetate",
    " sodium", " potassium", " calcium", " alfa", " beta",
    " monohydrate", " dihydrate", " anhydrous",
]


def simplify_drug_name(name: str) -> str:
    """Elimina sufijos de sales/formas farmacéuticas para mejorar búsqueda PubMed.
    Itera hasta eliminar todos los sufijos compuestos (ej: 'mesylate monohydrate').
    """
    s = name.lower().strip()
    changed = True
    while changed:
        changed = False
        for suffix in SALT_SUFFIXES:
            if s.endswith(suffix):
                s = s[: -len(suffix)].strip()
                changed = True
                break
    # Capitalizar primera letra
    return s.capitalize()


def query_pubmed(drug_name: str) -> dict:
    """
    Cuenta publicaciones en PubMed para 'drug_name' + HNSCC.
    Busca tanto con el nombre completo como con el nombre simplificado (sin sales).
    """
    base_name = simplify_drug_name(drug_name)
    # Usar el nombre simplificado si es distinto al original
    search_names = list({drug_name.title(), base_name})

    # Query 1: fármaco solo (para contexto)
    terms_drug = " OR ".join(f'"{n}"[Title/Abstract]' for n in search_names)
    q_drug = f"({terms_drug})"
    params_drug = {
        "db": "pubmed",
        "term": q_drug,
        "retmax": 0,
        "tool": "hnscc_drug_repurposing",
        "email": NCBI_EMAIL,
        "format": "json",
    }
    data_drug = get_json(NCBI_ESEARCH, params_drug)
    time.sleep(REQUEST_DELAY)
    total_drug = int(data_drug["esearchresult"]["count"]) if data_drug else 0

    # Query 2: fármaco + HNSCC (usando nombres simplificados)
    q_hnscc = f"({terms_drug}) AND {HNSCC_PUBMED_TERM}"
    params_hnscc = {
        "db": "pubmed",
        "term": q_hnscc,
        "retmax": 0,
        "tool": "hnscc_drug_repurposing",
        "email": NCBI_EMAIL,
        "format": "json",
    }
    data_hnscc = get_json(NCBI_ESEARCH, params_hnscc)
    time.sleep(REQUEST_DELAY)
    hnscc_papers = int(data_hnscc["esearchresult"]["count"]) if data_hnscc else 0

    return {
        "pubmed_total": total_drug,
        "pubmed_hnscc": hnscc_papers,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    log.info("=" * 60)
    log.info("Script 11: ClinicalTrials.gov + PubMed evidence")
    log.info("=" * 60)
    log.info(f"NCBI email: {NCBI_EMAIL}")

    # Cargar top 20
    if not INPUT_TOP20.exists():
        log.error(f"No encontrado: {INPUT_TOP20}")
        sys.exit(1)

    top20 = pd.read_csv(INPUT_TOP20, sep="\t")
    log.info(f"Top 20 candidatos cargados: {len(top20)} filas")
    log.info(f"Columnas: {list(top20.columns)}")

    # Identificar columna de nombre de fármaco
    name_col = None
    for c in ["drug_name", "drug_name_norm", "name", "compound"]:
        if c in top20.columns:
            name_col = c
            break
    if name_col is None:
        name_col = top20.columns[0]
    log.info(f"Columna de nombre de fármaco: '{name_col}'")

    drug_names = top20[name_col].tolist()
    log.info(f"Fármacos a consultar: {drug_names}")

    # ── Consultar ClinicalTrials y PubMed para cada candidato ────────────────
    results = []
    all_trial_details = []

    for i, drug in enumerate(drug_names):
        log.info(f"  [{i+1:02d}/{len(drug_names)}] {drug}")

        # ClinicalTrials
        ct = query_clinicaltrials(drug)
        log.info(f"    ClinicalTrials: {ct['total_trials']} totales, "
                 f"{ct['hnscc_trials']} en HNSCC")

        # PubMed
        pm = query_pubmed(drug)
        log.info(f"    PubMed: {pm['pubmed_total']} total, "
                 f"{pm['pubmed_hnscc']} HNSCC")

        # Fases de ensayos HNSCC
        phase_str = "; ".join(
            f"{ph}:{n}" for ph, n in sorted(ct["by_phase"].items())
        ) if ct["by_phase"] else "none"

        # Determinar si hay evidencia activa en HNSCC
        active_trials = sum(
            1 for t in ct["trials"]
            if t["status"] in {"RECRUITING", "ACTIVE_NOT_RECRUITING",
                                "NOT_YET_RECRUITING", "ENROLLING_BY_INVITATION"}
        )

        # Score de evidencia clínica (simple)
        evidence_score = (
            min(ct["hnscc_trials"] / 5.0, 1.0) * 0.5 +   # hasta 5 trials = max
            min(pm["pubmed_hnscc"] / 50.0, 1.0) * 0.3 +  # hasta 50 papers = max
            min(active_trials / 2.0, 1.0) * 0.2           # trials activos
        )

        row = top20.iloc[i].to_dict()
        row.update({
            "ct_total_trials": ct["total_trials"],
            "ct_hnscc_trials": ct["hnscc_trials"],
            "ct_active_trials": active_trials,
            "ct_phases_hnscc": phase_str,
            "pubmed_total": pm["pubmed_total"],
            "pubmed_hnscc": pm["pubmed_hnscc"],
            "evidence_score": round(evidence_score, 3),
        })
        results.append(row)

        # Guardar detalle de ensayos
        for t in ct["trials"]:
            t["drug_name"] = drug
            all_trial_details.append(t)

    # ── Exportar tablas ───────────────────────────────────────────────────────
    df_evidence = pd.DataFrame(results)
    df_evidence = df_evidence.sort_values("evidence_score", ascending=False)

    out_evidence = OUT_DIR / "11_clinical_evidence.tsv"
    df_evidence.to_csv(out_evidence, sep="\t", index=False)
    log.info(f"Tabla de evidencia guardada: {out_evidence}")

    if all_trial_details:
        df_trials = pd.DataFrame(all_trial_details)
        out_trials = OUT_DIR / "11_trials_detail.tsv"
        df_trials.to_csv(out_trials, sep="\t", index=False)
        log.info(f"Detalle de ensayos guardado: {out_trials}")
    else:
        log.warning("No se encontraron ensayos HNSCC para ningún candidato")

    # ── Figuras ───────────────────────────────────────────────────────────────
    # Color por clase
    class_colors = {"A": "#d62728", "B": "#ff7f0e", "C": "#2ca02c", "D": "#1f77b4"}
    class_col = next((c for c in ["drug_class", "repurposing_class", "class", "clase"] if c in df_evidence.columns), None)

    def get_color(row):
        if class_col and class_col in row:
            return class_colors.get(str(row[class_col]).upper(), "#aec7e8")
        return "#aec7e8"

    df_plot = df_evidence.head(20).copy()
    drug_labels = df_plot[name_col].str[:22].tolist()   # truncar etiquetas largas

    # --- Figura 1: Bubble chart (PubMed HNSCC vs Trials HNSCC) ---------------
    fig, ax = plt.subplots(figsize=(9, 7))
    colors = [get_color(r) for _, r in df_plot.iterrows()]
    sizes  = (df_plot["evidence_score"] * 500 + 50).values

    ax.scatter(
        df_plot["pubmed_hnscc"],
        df_plot["ct_hnscc_trials"],
        s=sizes,
        c=colors,
        alpha=0.75,
        edgecolors="white",
        linewidths=0.8,
    )

    for i, row in df_plot.iterrows():
        ax.annotate(
            row[name_col][:18],
            (row["pubmed_hnscc"], row["ct_hnscc_trials"]),
            fontsize=6.5,
            ha="left",
            va="bottom",
            xytext=(3, 3),
            textcoords="offset points",
        )

    ax.set_xlabel("Publicaciones PubMed (Drug + HNSCC)", fontsize=11)
    ax.set_ylabel("Ensayos ClinicalTrials (Drug + HNSCC)", fontsize=11)
    ax.set_title("Evidencia clínica y bibliográfica — Top 20 candidatos\n"
                 "(tamaño de burbuja = evidence score)", fontsize=12)

    legend_patches = [
        mpatches.Patch(color=v, label=f"Clase {k}") for k, v in class_colors.items()
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3, linestyle="--")

    plt.tight_layout()
    out_bubble = FIG_DIR / "11_evidence_bubble.pdf"
    plt.savefig(out_bubble, dpi=150)
    plt.close()
    log.info(f"Figura bubble: {out_bubble}")

    # --- Figura 2: Barras horizontales — trials HNSCC y PubMed HNSCC ---------
    fig, axes = plt.subplots(1, 2, figsize=(13, 7))

    y_pos = range(len(df_plot))
    bar_colors = [get_color(r) for _, r in df_plot.iterrows()]

    # Ensayos ClinicalTrials
    axes[0].barh(
        list(y_pos), df_plot["ct_hnscc_trials"].values,
        color=bar_colors, edgecolor="white", height=0.7
    )
    axes[0].set_yticks(list(y_pos))
    axes[0].set_yticklabels(drug_labels, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("N ensayos clínicos (HNSCC)", fontsize=10)
    axes[0].set_title("ClinicalTrials.gov", fontsize=11)
    axes[0].axvline(0, color="black", linewidth=0.5)
    for j, v in enumerate(df_plot["ct_hnscc_trials"].values):
        if v > 0:
            axes[0].text(v + 0.1, j, str(v), va="center", fontsize=7)

    # PubMed
    axes[1].barh(
        list(y_pos), df_plot["pubmed_hnscc"].values,
        color=bar_colors, edgecolor="white", height=0.7
    )
    axes[1].set_yticks(list(y_pos))
    axes[1].set_yticklabels(drug_labels, fontsize=8)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("N publicaciones PubMed (drug + HNSCC)", fontsize=10)
    axes[1].set_title("PubMed", fontsize=11)
    axes[1].axvline(0, color="black", linewidth=0.5)
    for j, v in enumerate(df_plot["pubmed_hnscc"].values):
        if v > 0:
            axes[1].text(v + 0.5, j, str(v), va="center", fontsize=7)

    legend_patches = [
        mpatches.Patch(color=v, label=f"Clase {k}") for k, v in class_colors.items()
    ]
    fig.legend(handles=legend_patches, fontsize=9, loc="lower center",
               ncol=4, bbox_to_anchor=(0.5, -0.02))

    plt.suptitle("Evidencia clínica + bibliográfica — Top 20 candidatos HNSCC",
                 fontsize=13, fontweight="bold")
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    out_bars = FIG_DIR / "11_trials_phase_bar.pdf"
    plt.savefig(out_bars, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"Figura barras: {out_bars}")

    # ── Resumen final ─────────────────────────────────────────────────────────
    log.info("\n" + "=" * 60)
    log.info("RESUMEN — Evidencia clínica y bibliográfica")
    log.info("=" * 60)
    log.info(f"Candidatos consultados: {len(results)}")
    log.info(f"Con ensayos HNSCC:      {(df_evidence['ct_hnscc_trials'] > 0).sum()}")
    log.info(f"Con trials activos:     {(df_evidence['ct_active_trials'] > 0).sum()}")
    log.info(f"Con papers HNSCC:       {(df_evidence['pubmed_hnscc'] > 0).sum()}")
    log.info("\nTop 10 por evidencia_score:")
    cols_show = [name_col, "ct_hnscc_trials", "ct_active_trials",
                 "pubmed_hnscc", "evidence_score"]
    cols_show = [c for c in cols_show if c in df_evidence.columns]
    log.info("\n" + df_evidence[cols_show].head(10).to_string(index=False))
    log.info(f"\nOutputs:\n  {out_evidence}\n  {out_bubble}\n  {out_bars}")
    log.info("Script 11 completado.")


if __name__ == "__main__":
    main()
