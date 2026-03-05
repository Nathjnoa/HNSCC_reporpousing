#!/usr/bin/env python3
"""
12_cosmic_overlap.py
====================
Cruza las proteínas DE de HNSCC con bases de datos de genes driver de cáncer:

1. COSMIC Cancer Gene Census (CGC) — si el archivo fue descargado manualmente en
   data/raw/cosmic/cancer_gene_census.csv  (requiere cuenta COSMIC, gratuita)
   Descarga: https://cancer.sanger.ac.uk/cosmic/download -> Cancer Gene Census

2. IntOGen Compendium of Mutational Drivers — descarga programática pública.
   Contiene genes driver con FDR < 0.05 en TCGA HNSC entre otros tumores.
   URL: https://www.intogen.org/api/

3. NCG (Network of Cancer Genes) — lista curada de oncogenes y tumor suppressors.
   URL: http://ncg.kcl.ac.uk/ (archivo público)

Outputs:
  results/tables/evidence/12_cosmic_overlap.tsv      — genes DE en CGC/IntOGen/NCG
  results/tables/evidence/12_cancer_drivers_annot.tsv — anotación completa por gen
  results/figures/evidence/12_venn_driver_sources.pdf — Venn de fuentes
  results/figures/evidence/12_driver_logfc_plot.pdf   — volcano con drivers resaltados

Descarga manual de COSMIC CGC (opcional pero recomendada):
  1. Crear cuenta en https://cancer.sanger.ac.uk/ (gratuita)
  2. Descargar "Cancer Gene Census" en formato CSV
  3. Copiar a: data/raw/cosmic/cancer_gene_census.csv

Uso:
  conda activate omics-py
  cd ~/bioinfo/projects/hnscc_drug_repurposing
  python scripts/12_cosmic_overlap.py
"""

import os
import sys
import io
import logging
import datetime
import zipfile
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# ── Configuración ────────────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).resolve().parent
PROJ_DIR     = SCRIPT_DIR.parent
INPUT_DE     = PROJ_DIR / "results/tables/de_limma/01_TVsS_all_proteins.tsv"
INPUT_SIG    = PROJ_DIR / "results/tables/de_limma/02_TVsS_significant_with_ids.tsv"
INPUT_TOP20  = PROJ_DIR / "results/tables/10_top20_candidates.tsv"
COSMIC_FILE  = PROJ_DIR / "data/raw/cosmic/cancer_gene_census.csv"
OUT_DIR      = PROJ_DIR / "results/tables/evidence"
FIG_DIR      = PROJ_DIR / "results/figures/evidence"
LOG_DIR      = PROJ_DIR / "logs"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)

TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_FILE  = LOG_DIR / f"12_cosmic_overlap_{TIMESTAMP}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)

# IntOGen drivers download URL
INTOGEN_DRIVERS_URL = "https://www.intogen.org/api/compendium/mutations/drivers/download?format=csv"
# Alternative: static compendium (may change)
INTOGEN_ALT_URL = "https://www.intogen.org/download?file=Compendium_Cancer_Genes.tsv.gz"

# NCG7 canonical cancer genes (public list embedded for robustness)
# Compiled from NCG7 (Repana et al. 2019, Genome Biology)
NCG7_ONCOGENES = {
    "ABL1", "AKT1", "ALK", "APC", "AR", "ARAF", "ARFRP1", "ARID1A", "ATM",
    "BRAF", "BRCA1", "BRCA2", "BTK", "CBL", "CCND1", "CCND2", "CCND3",
    "CCNE1", "CDH1", "CDK4", "CDK6", "CDKN2A", "CDKN2B", "CDKN2C", "CEBPA",
    "CSF1R", "CTNNB1", "DDR2", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ESR1",
    "EZH2", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT3", "FOXL2", "GATA3",
    "GNA11", "GNAQ", "GNAS", "HNF1A", "HRAS", "IDH1", "IDH2", "IGF1R",
    "JAK1", "JAK2", "JAK3", "KDM5C", "KDM6A", "KIT", "KRAS", "MAP2K1",
    "MAP2K2", "MAP2K4", "MDM2", "MDM4", "MET", "MLH1", "MSH2", "MSH6",
    "MTOR", "MYC", "MYCL", "MYCN", "MYD88", "NF1", "NF2", "NFE2L2", "NOTCH1",
    "NOTCH2", "NRAS", "NTRK1", "NTRK2", "NTRK3", "PALB2", "PBRM1", "PDGFRA",
    "PDGFRB", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PMS2", "POLE", "PTEN",
    "PTPN11", "RAC1", "RAF1", "RB1", "RET", "RHOA", "ROS1", "RUNX1", "SF3B1",
    "SMAD2", "SMAD4", "SMARCA4", "SMARCB1", "SMO", "SRC", "STK11", "TERT",
    "TET2", "TP53", "TSC1", "TSC2", "U2AF1", "VHL", "WT1", "ZNRF3",
    # HNSCC-relevant adicionales
    "CASP8", "CSMD1", "CSMD3", "FAT1", "FBXW7", "KDM6A", "MLL2",
    "NOTCH3", "NOTCH4", "PIK3R2", "RB1", "SOX2", "TP63",
}

# HNSCC-specific drivers from literature (TCGA HNSCC 2015, Alexandrov 2013)
HNSCC_DRIVERS = {
    "TP53", "CDKN2A", "EGFR", "HRAS", "PIK3CA", "PTEN", "FAT1", "NOTCH1",
    "CASP8", "FBXW7", "NSD1", "KDM6A", "MLL2", "ZNF750", "TGFBR2",
    "SMAD4", "CTCF", "RB1", "FGFR1", "MET", "CCND1", "SOX2",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def download_intogen(timeout: int = 60) -> pd.DataFrame | None:
    """Descarga IntOGen Compendium of Cancer Drivers. Devuelve DataFrame o None."""
    log.info("Descargando IntOGen Compendium of Cancer Drivers...")
    for url in [INTOGEN_DRIVERS_URL, INTOGEN_ALT_URL]:
        try:
            r = requests.get(url, timeout=timeout, stream=True)
            r.raise_for_status()
            ctype = r.headers.get("Content-Type", "")

            if "gzip" in ctype or url.endswith(".gz"):
                import gzip
                content = gzip.decompress(r.content)
                df = pd.read_csv(io.BytesIO(content), sep="\t", low_memory=False)
            elif "zip" in ctype or url.endswith(".zip"):
                z = zipfile.ZipFile(io.BytesIO(r.content))
                name = z.namelist()[0]
                df = pd.read_csv(z.open(name), sep="\t", low_memory=False)
            else:
                sep = "\t" if "tsv" in url else ","
                df = pd.read_csv(io.StringIO(r.text), sep=sep, low_memory=False)

            log.info(f"  IntOGen: {len(df)} filas, columnas: {list(df.columns[:6])}")
            return df
        except Exception as exc:
            log.warning(f"  IntOGen URL {url} falló: {exc}")

    log.warning("  IntOGen: no se pudo descargar. Usando genes HNSCC embebidos.")
    return None


def load_cosmic(filepath: Path) -> pd.DataFrame | None:
    """Carga COSMIC Cancer Gene Census desde archivo local."""
    if not filepath.exists():
        log.warning(f"COSMIC CGC no encontrado en: {filepath}")
        log.warning("Para incluirlo:")
        log.warning("  1. Crear cuenta gratuita en https://cancer.sanger.ac.uk/")
        log.warning("  2. Descargar Cancer Gene Census (CSV)")
        log.warning(f"  3. Copiar a: {filepath}")
        return None

    df = pd.read_csv(filepath, low_memory=False)
    log.info(f"COSMIC CGC cargado: {len(df)} genes, columnas: {list(df.columns[:6])}")
    return df


def extract_cosmic_hnscc_genes(df: pd.DataFrame) -> set:
    """Extrae genes con evidencia en tumores de cabeza y cuello del CGC."""
    if df is None:
        return set()
    # Columna de tumor types varía con la versión
    tumor_cols = [c for c in df.columns if "tumour" in c.lower() or "tumor" in c.lower()]
    gene_col   = next((c for c in df.columns if "gene" in c.lower() and "symbol" in c.lower()), None)
    if gene_col is None:
        gene_col = df.columns[0]

    if not tumor_cols:
        log.warning("COSMIC: no se encontró columna de tumor types. Usando todos los genes CGC.")
        return set(df[gene_col].dropna().astype(str).str.upper())

    hnscc_genes = set()
    hnscc_kw = ["head", "neck", "oral", "pharyn", "laryn", "squamous"]
    for tc in tumor_cols:
        mask = df[tc].astype(str).str.lower().apply(
            lambda x: any(kw in x for kw in hnscc_kw)
        )
        hnscc_genes.update(df.loc[mask, gene_col].dropna().astype(str).str.upper())

    log.info(f"COSMIC genes con evidencia HNSCC: {len(hnscc_genes)}")
    return hnscc_genes


def extract_intogen_hnscc_genes(df: pd.DataFrame) -> set:
    """Extrae genes driver de HNSC/HNSCC del compendio IntOGen."""
    if df is None:
        return set()

    # Detectar columnas relevantes
    gene_col  = next((c for c in df.columns if c.upper() in {"GENE", "SYMBOL", "GENE_ID", "HUGO"}), None)
    tumor_col = next((c for c in df.columns if c.upper() in {"CANCER_TYPE", "COHORT", "TUMOR_TYPE", "TUMOR"}), None)

    if gene_col is None:
        gene_col = df.columns[0]
    if tumor_col is None:
        log.warning("IntOGen: columna de tumor type no identificada. Usando todos los genes.")
        return set(df[gene_col].dropna().astype(str).str.upper())

    hnscc_kw = ["hnsc", "hnscc", "head", "neck", "oral", "laryn", "pharyn"]
    mask = df[tumor_col].astype(str).str.lower().apply(
        lambda x: any(kw in x for kw in hnscc_kw)
    )
    genes = set(df.loc[mask, gene_col].dropna().astype(str).str.upper())
    log.info(f"IntOGen genes HNSCC-specific: {len(genes)}")

    if len(genes) == 0:
        log.info("IntOGen: 0 genes HNSCC-specific. Usando todos los genes del compendio.")
        genes = set(df[gene_col].dropna().astype(str).str.upper())
        log.info(f"IntOGen genes totales (pan-cancer): {len(genes)}")

    return genes


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    log.info("=" * 60)
    log.info("Script 12: COSMIC / IntOGen / NCG overlap con proteínas DE")
    log.info("=" * 60)

    # ── Cargar datos DE ────────────────────────────────────────────────────────
    df_sig = pd.read_csv(INPUT_SIG, sep="\t")
    log.info(f"Proteinas DE significativas: {len(df_sig)}")

    # Columna de gene symbol
    sym_col = next((c for c in df_sig.columns if c.lower() in
                    {"gene_symbol", "symbol", "symbol_org", "gene_name"}), None)
    if sym_col is None:
        sym_col = df_sig.columns[0]
    log.info(f"Columna gene symbol: '{sym_col}'")

    de_genes = set(df_sig[sym_col].dropna().astype(str).str.upper())
    log.info(f"Genes DE unicos: {len(de_genes)}")

    # Cargar logFC para volcano
    df_all = pd.read_csv(INPUT_DE, sep="\t")
    lfc_col  = next((c for c in df_all.columns
                     if c.lower() in {"logfc", "log2fc", "logfc_tvss"}), None)
    padj_col = next((c for c in df_all.columns
                     if c.lower() in {"adj.p.val_tvss", "adj.p.val", "adj_p_val", "padj", "fdr"}), None)
    # Fallback: buscar por patrón
    if lfc_col is None:
        lfc_col = next((c for c in df_all.columns if "logfc" in c.lower()), None)
    if padj_col is None:
        padj_col = next((c for c in df_all.columns if "adj.p" in c.lower()), None)

    # ── Cargar fuentes de cancer drivers ─────────────────────────────────────
    cosmic_raw  = load_cosmic(COSMIC_FILE)
    intogen_raw = download_intogen()

    cosmic_genes  = extract_cosmic_hnscc_genes(cosmic_raw)
    intogen_genes = extract_intogen_hnscc_genes(intogen_raw)
    ncg_genes     = NCG7_ONCOGENES
    hnscc_lit     = HNSCC_DRIVERS

    log.info(f"\nGenes en cada fuente:")
    log.info(f"  COSMIC CGC (HNSCC): {len(cosmic_genes)}")
    log.info(f"  IntOGen HNSCC:      {len(intogen_genes)}")
    log.info(f"  NCG7 cancer genes:  {len(ncg_genes)}")
    log.info(f"  HNSCC lit drivers:  {len(hnscc_lit)}")

    # ── Cruzar con genes DE ───────────────────────────────────────────────────
    cosmic_de  = de_genes & cosmic_genes
    intogen_de = de_genes & intogen_genes
    ncg_de     = de_genes & ncg_genes
    lit_de     = de_genes & hnscc_lit

    log.info(f"\nCruce con genes DE significativos ({len(de_genes)}):")
    log.info(f"  COSMIC ∩ DE:   {len(cosmic_de)}")
    log.info(f"  IntOGen ∩ DE:  {len(intogen_de)}")
    log.info(f"  NCG7 ∩ DE:     {len(ncg_de)}")
    log.info(f"  LitHNSCC ∩ DE: {len(lit_de)}")

    any_driver = cosmic_de | intogen_de | ncg_de | lit_de
    log.info(f"  Union (any source): {len(any_driver)}")

    # ── Construir tabla de anotación ──────────────────────────────────────────
    # Partir de proteinas DE significativas
    records = []
    for _, row in df_sig.iterrows():
        gene = str(row.get(sym_col, "")).upper()
        in_cosmic   = gene in cosmic_genes
        in_intogen  = gene in intogen_genes
        in_ncg      = gene in ncg_genes
        in_lit_hnscc = gene in hnscc_lit
        n_sources = sum([in_cosmic, in_intogen, in_ncg, in_lit_hnscc])

        r = row.to_dict()
        r.update({
            "gene_upper": gene,
            "in_cosmic_cgc": in_cosmic,
            "in_intogen": in_intogen,
            "in_ncg7": in_ncg,
            "in_hnscc_lit": in_lit_hnscc,
            "n_driver_sources": n_sources,
            "is_cancer_driver": n_sources >= 1,
        })
        records.append(r)

    df_annot = pd.DataFrame(records)
    df_annot = df_annot.sort_values("n_driver_sources", ascending=False)

    out_annot = OUT_DIR / "12_cancer_drivers_annot.tsv"
    df_annot.to_csv(out_annot, sep="\t", index=False)
    log.info(f"\nTabla de anotación guardada: {out_annot}")

    # Tabla de solapamiento (solo los drivers)
    df_overlap = df_annot[df_annot["is_cancer_driver"]].copy()
    out_overlap = OUT_DIR / "12_cosmic_overlap.tsv"
    df_overlap.to_csv(out_overlap, sep="\t", index=False)
    log.info(f"Tabla de overlap guardada: {out_overlap} ({len(df_overlap)} genes)")

    # ── Figuras ───────────────────────────────────────────────────────────────

    # ---- Figura 1: Barras de overlap por fuente --------------------------------
    sources = {
        "COSMIC CGC\n(HNSCC)": len(cosmic_de),
        "IntOGen\n(HNSCC)": len(intogen_de),
        "NCG7\n(pan-cancer)": len(ncg_de),
        "Lit. HNSCC\ndrivers": len(lit_de),
        "Union\n(cualquiera)": len(any_driver),
    }
    colors = ["#e41a1c", "#ff7f00", "#4daf4a", "#377eb8", "#984ea3"]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(list(sources.keys()), list(sources.values()),
                  color=colors, edgecolor="white", width=0.6)
    for bar, val in zip(bars, sources.values()):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(val), ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylabel("N proteínas DE en base de datos", fontsize=11)
    ax.set_title("Solapamiento proteínas DE (HNSCC) con genes driver de cáncer",
                 fontsize=12)
    ax.set_ylim(0, max(sources.values()) * 1.2)
    ax.axhline(0, color="black", linewidth=0.5)
    plt.tight_layout()
    out_bar = FIG_DIR / "12_driver_overlap_bar.pdf"
    plt.savefig(out_bar, dpi=150)
    plt.close()
    log.info(f"Figura barras: {out_bar}")

    # ---- Figura 2: Volcano con drivers resaltados ------------------------------
    if lfc_col and padj_col and lfc_col in df_all.columns and padj_col in df_all.columns:
        df_vol = df_all[[sym_col, lfc_col, padj_col]].dropna().copy()
        df_vol.columns = ["gene", "logFC", "padj"]
        df_vol["gene_upper"] = df_vol["gene"].str.upper()
        df_vol["neg_log10p"] = -np.log10(df_vol["padj"].clip(lower=1e-300))
        df_vol["is_sig"] = (
            (df_vol["padj"] < 0.05) & (df_vol["logFC"].abs() > 1)
        )
        df_vol["is_driver"] = df_vol["gene_upper"].isin(any_driver)

        fig, ax = plt.subplots(figsize=(9, 7))

        # Fondo: no significativos
        mask_ns = ~df_vol["is_sig"]
        ax.scatter(df_vol.loc[mask_ns, "logFC"],
                   df_vol.loc[mask_ns, "neg_log10p"],
                   s=4, c="#cccccc", alpha=0.4, rasterized=True)

        # Significativos no-driver
        mask_sig_nd = df_vol["is_sig"] & ~df_vol["is_driver"]
        ax.scatter(df_vol.loc[mask_sig_nd, "logFC"],
                   df_vol.loc[mask_sig_nd, "neg_log10p"],
                   s=12, c="#74c0fc", alpha=0.6, rasterized=True)

        # Drivers DE resaltados
        mask_driver = df_vol["is_sig"] & df_vol["is_driver"]
        ax.scatter(df_vol.loc[mask_driver, "logFC"],
                   df_vol.loc[mask_driver, "neg_log10p"],
                   s=50, c="#e03131", alpha=0.9, edgecolors="white",
                   linewidths=0.6, zorder=5)

        # Etiquetar drivers
        for _, gr in df_vol[mask_driver].iterrows():
            ax.annotate(gr["gene"], (gr["logFC"], gr["neg_log10p"]),
                        fontsize=6.5, ha="left", va="bottom",
                        xytext=(3, 2), textcoords="offset points",
                        color="#c92a2a")

        # Líneas de referencia
        ax.axvline(-1, color="gray", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.axvline(1,  color="gray", linewidth=0.8, linestyle="--", alpha=0.6)
        ax.axhline(-np.log10(0.05), color="gray", linewidth=0.8, linestyle="--", alpha=0.6)

        legend_elems = [
            mpatches.Patch(color="#cccccc", label="No significativo"),
            mpatches.Patch(color="#74c0fc", label="DE significativo"),
            mpatches.Patch(color="#e03131", label=f"Cancer driver DE ({mask_driver.sum()})"),
        ]
        ax.legend(handles=legend_elems, fontsize=9, loc="upper left")

        ax.set_xlabel("log2 Fold Change (Tumor vs Normal)", fontsize=11)
        ax.set_ylabel("-log10(adj.P.Val)", fontsize=11)
        ax.set_title("Volcano: proteínas DE con genes driver de cáncer resaltados",
                     fontsize=12)
        plt.tight_layout()
        out_vol = FIG_DIR / "12_driver_logfc_plot.pdf"
        plt.savefig(out_vol, dpi=150)
        plt.close()
        log.info(f"Figura volcano: {out_vol}")
    else:
        log.warning("Columnas logFC/padj no encontradas; omitiendo volcano plot")
        out_vol = "N/A"

    # ── Anotar top 20 con estado de cancer driver ─────────────────────────────
    top20 = pd.read_csv(INPUT_TOP20, sep="\t")
    driver_info = (
        df_annot.groupby(sym_col)
        .agg(
            in_cosmic_cgc=("in_cosmic_cgc", "any"),
            in_intogen=("in_intogen", "any"),
            in_ncg7=("in_ncg7", "any"),
            in_hnscc_lit=("in_hnscc_lit", "any"),
            n_driver_sources=("n_driver_sources", "max"),
        )
        .reset_index()
    )

    # Expandir de_genes a filas por gen para cruzar con top 20
    # La columna de_genes contiene lista separada por "/"
    if "de_genes" in top20.columns:
        def annotate_row(de_genes_str):
            if pd.isna(de_genes_str):
                return 0, False
            genes_up = [g.strip().upper() for g in str(de_genes_str).split("/")]
            n_drivers = sum(g in any_driver for g in genes_up)
            return n_drivers, n_drivers > 0

        top20[["n_target_drivers", "has_driver_target"]] = pd.DataFrame(
            top20["de_genes"].apply(annotate_row).tolist(), index=top20.index
        )
        out_top20_annot = OUT_DIR / "12_top20_with_drivers.tsv"
        top20.to_csv(out_top20_annot, sep="\t", index=False)
        log.info(f"Top 20 anotado con drivers: {out_top20_annot}")
        log.info(f"Top 20 con target en driver: {top20['has_driver_target'].sum()}/20")

    # ── Resumen final ─────────────────────────────────────────────────────────
    log.info("\n" + "=" * 60)
    log.info("RESUMEN — Solapamiento con genes driver de cáncer")
    log.info("=" * 60)
    log.info(f"Proteinas DE totales:          {len(de_genes)}")
    log.info(f"En COSMIC CGC (HNSCC):         {len(cosmic_de)}")
    log.info(f"En IntOGen (HNSCC/pan-cancer): {len(intogen_de)}")
    log.info(f"En NCG7:                       {len(ncg_de)}")
    log.info(f"En literatura HNSCC:           {len(lit_de)}")
    log.info(f"En cualquier fuente:           {len(any_driver)}")
    log.info(f"\nTop drivers DE (>= 2 fuentes):")
    extra_cols = [c for c in [lfc_col, padj_col] if c and c in df_annot.columns]
    top_drivers = df_annot[df_annot["n_driver_sources"] >= 2][[sym_col] + extra_cols + ["n_driver_sources"]].head(20)
    log.info("\n" + top_drivers.to_string(index=False))
    log.info(f"\nOutputs:\n  {out_overlap}\n  {out_annot}\n  {out_bar}\n  {out_vol}")
    log.info("Script 12 completado.")


if __name__ == "__main__":
    main()
