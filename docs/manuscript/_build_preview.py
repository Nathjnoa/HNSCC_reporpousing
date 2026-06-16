#!/usr/bin/env python3
"""Assemble a clean preview manuscript from the section drafts.

Strips internal note blocks, ⚠️ flags, markdown emphasis (plain output per author),
maps Introduction placeholder citekeys to real verified keys, injects key method
citations, and now embeds the main-body multipanel figures (Fig 1-6) with their
legends and the prioritized-candidate main table (Table 1, assembled from the
Tab4 + Tab5 TSVs). Supplementary items are NOT included. Produces
manuscript_PREVIEW.md for pandoc -> docx.

Citation handling: every injected token is a verified Better BibTeX key present in
references.bib (checked against the bib before adding new anchors), so the markers
docx scans cleanly in Zotero (REFERENCES_WORKFLOW.md, Path B).
"""
import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TAB_DIR = HERE / ".." / ".." / "results" / "tables" / "pub" / "main"
FIG_DIR = "../../results/figures/pub/main"  # relative to this file; pandoc run from HERE


def body_from(path):
    """Return content from the first level-2 heading onward (drops the note block)."""
    lines = path.read_text(encoding="utf-8").splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("## "):
            return "\n".join(lines[i:])
    raise ValueError(f"no '## ' heading in {path}")


def clean_common(text):
    text = re.sub(r"⚠️ TO-VERIFY\s*", "", text)
    # unwrap backticked citekeys: `[@key]` -> [@key]
    text = re.sub(r"`(\[@[^\]]+\])`", r"\1", text)
    # drop bold then italic emphasis (plain output)
    text = re.sub(r"\*\*([^*\n]+)\*\*", r"\1", text)
    text = re.sub(r"\*([^*\n]+)\*", r"\1", text)
    return text


# --- Introduction: map placeholder keys -> verified citekeys -------------------
INTRO_CITES = {
    "[@hnscc_epidemiology]": "[@johnson2020; @bray2024]",
    "[@hnscc_outcomes]": "[@johnson2020]",
    "[@keynote048]": "[@burtness2019]",
    "[@keynote689]": "[@uppaluri2025]",
    "[@hpv_subtypes]": "[@ang2010]",
    "[@proteomics_phenotype]": "[@aebersold2016]",
    "[@hnscc_proteomics]": "[@huang2021]",
    "[@repurposing_rationale]": "[@pushpakom2019; @tanoli2025]",
}

# --- Methods: inject method/tool citations at first single-token anchor ---------
METHOD_CITES = [
    ("proteoDA", "proteoDA [@ritchie2015]"),
    ("(GSEA)", "(GSEA) [@subramanian2005]"),
    ("clusterProfiler", "clusterProfiler [@wu2021]"),
    ("DGIdb", "DGIdb [@cannon2024]"),
    ("ChEMBL", "ChEMBL [@mendez2019]"),
    ("Open Targets", "Open Targets [@ochoa2021]"),
    ("L2S2", "L2S2 [@evangelista2025]"),
    ("STRING", "STRING [@szklarczyk2023]"),
    ("Hallmark", "Hallmark [@liberzon2015]"),
    ("Huang et al., 2021", "Huang et al., 2021 [@huang2021]"),
]

# --- Figure legends: inject the one author citation that appears in captions -----
FIGURE_CITES = [
    ("Huang et al., 2021", "Huang et al., 2021 [@huang2021]"),
]

# Main-body multipanel figures: file stem -> figure number (legend pulled from figures.md)
MAIN_FIGURES = [
    ("Fig1_workflow", 1),
    ("Fig2_multipanel", 2),
    ("Fig3_multipanel", 3),
    ("Fig4_multipanel", 4),
    ("Fig5_multipanel", 5),
    ("Fig6_multipanel", 6),
]


def inject_first(text, anchor, replacement):
    return text.replace(anchor, replacement, 1)


# ---- Figure legends keyed by number, harvested from figures.md -----------------
def figure_legends():
    text = (HERE / "figures.md").read_text(encoding="utf-8")
    # stop before the supplementary block
    text = text.split("## Supplementary figures")[0]
    legends = {}
    blocks = re.split(r"\n## (Figure \d+\.)", text)
    # blocks[0] is the header; then pairs (title, body)
    for i in range(1, len(blocks), 2):
        title = blocks[i].strip()
        body = blocks[i + 1].strip()
        num = int(re.search(r"Figure (\d+)\.", title).group(1))
        legend = f"{title} {body}"
        for anchor, repl in FIGURE_CITES:
            legend = inject_first(legend, anchor, repl)
        legends[num] = clean_common(legend)
    return legends


def figures_section():
    legends = figure_legends()
    out = ["# Figures\n"]
    for stem, num in MAIN_FIGURES:
        out.append(f"![]({FIG_DIR}/{stem}.png)\n")
        out.append(legends[num] + "\n")
    return "\n".join(out)


# ---- Table 1: assemble from Tab4 (EGFR control) + Tab5 (tiers) ------------------
PHASE_EN = {
    "Aprobado (Fase IV)": "Approved (Phase IV)",
    "Fase III": "Phase III",
    "Fase II": "Phase II",
    "Fase I": "Phase I",
}
SUBCLASS_EN = {
    "Anticuerpo/ADC anti-EGFR": "Anti-EGFR antibody/ADC",
    "Inhibidor multi-quinasa (EGFR+)": "Multikinase inhibitor (EGFR+)",
    "TKI EGFR 1ª-2ª generación": "EGFR TKI (1st-2nd gen)",
    "TKI EGFR 3ª generación": "EGFR TKI (3rd gen)",
}
# Module label per anchor hub, matching the Fig 4 curated labels exactly
# (MODULE_NAMES in scripts/17_pub_figures.R; keyed by hub to avoid GO-name ambiguity).
HUB_MODULE_LABEL = {
    "NDUFS2": "OXPHOS",
    "PSMA2": "Proteasome",
    "MMP9": "Immune response",
    "RPS11": "Translation",
    "CDA": "Amino acid metabolism",
    "DNMT1": "Chromatin remodeling",
    "PDE6D": "Amino acid metabolism",
    "PKLR": "Carbohydrate catabolism",
    "ALDH5A1": "Carboxylic acid metabolism",
    "IMPDH2": "Amino acid metabolism",
    "CA2": "Cell adhesion / membrane",
    "MAOA": "Amino acid metabolism",
    "MMP8": "Immune response",
}


def _read_tsv(name):
    with open(TAB_DIR / name, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def table1_section():
    egfr = _read_tsv("Tab4_EGFR_validation.tsv")
    tiers = _read_tsv("Tab5_novel_candidates_by_module.tsv")  # includes Doxycycline (N_TAB5=13)

    header = (
        "| Drug | Block / tier | Module / subclass | Anchor target | Max clinical phase "
        "| HNSCC-approved | Composite | TP | DV | Sources (/4) | LOD-stable | Robustness (n/6) |"
    )
    sep = "|" + "---|" * 12

    rows = []

    # EGFR control block: the 7 leaders (composite >= 0.689); collapse the long tail
    leaders = [r for r in egfr if float(r["Composite score"]) >= 0.689]
    leaders.sort(key=lambda r: -float(r["Composite score"]))
    n_more = len(egfr) - len(leaders)
    for r in leaders:
        rows.append(
            f"| {r['Fármaco']} | EGFR control | "
            f"{SUBCLASS_EN.get(r['Subclase'], r['Subclase'])} | {r['Target primario']} "
            f"| {PHASE_EN.get(r['Fase clínica máx.'], r['Fase clínica máx.'])} "
            f"| {'Yes' if r['Aprobado HNSCC'] == 'Sí' else 'No'} | {r['Composite score']} "
            f"| {r['TP']} | {r['DV']} | {r['N fuentes']} "
            f"| {'Yes' if r['LOD-stable'] == 'Sí' else 'No'} | {r['Robustez pesos']} |"
        )
    rows.append(
        f"| *plus {n_more} additional anti-EGFR agents* | EGFR control | anti-EGFR (mixed) "
        f"| EGFR | Approved-to-Phase I | mixed | 0.50-0.66 | 0.605 | varies | 2-3 "
        f"| mixed | mixed |"
    )

    tier_label = {"Hub central de red": "Hub-central", "Periférico diferencial": "Peripheral-differential"}
    for r in tiers:
        rows.append(
            f"| {r['Fármaco']} | {tier_label.get(r['Tier'], r['Tier'])} | "
            f"{HUB_MODULE_LABEL.get(r['Hub / diana'], r['Módulo funcional'])} | {r['Hub / diana']} "
            f"| {PHASE_EN.get(r['Fase clínica máx.'], r['Fase clínica máx.'])} | No "
            f"| {r['Composite score']} | {r['TP']} | {r['DV']} | {r['N fuentes']} "
            f"| {'Yes' if r['LOD-stable'] == 'Sí' else 'No'} | {r['Robustez pesos']} |"
        )

    legend = clean_common(
        body_from(HERE / "tables.md")
        .split("## Supplementary tables")[0]
        .split("## Table 1.")[1]
        .strip()
    )
    legend = "Table 1. " + legend

    return "\n".join(["# Tables\n", header, sep, *rows, "", legend])


# ---- assemble ------------------------------------------------------------------
intro = clean_common(body_from(HERE / "introduction.md"))
for k, v in INTRO_CITES.items():
    intro = intro.replace(k, v)

methods = clean_common(body_from(HERE / "methods.md"))
for anchor, repl in METHOD_CITES:
    methods = inject_first(methods, anchor, repl)

results = clean_common(body_from(HERE / "results.md"))

discussion = clean_common(body_from(HERE / "discussion.md"))

header = (
    "# Network-anchored proteomic drug repurposing in HNSCC (working draft, preview)\n\n"
    "Partial preview for pipeline testing: Introduction, Methods, Results, Discussion, main-body "
    "figures (Fig 1-6) and the prioritized-candidate main table (Table 1). Abstract, Title and "
    "complete Methods/Results citations are pending (WRITING_PLAN.md steps 6-8). Citations "
    "shown are PubMed-verified; the bibliography is generated by pandoc / scanned in Zotero.\n\n"
)

out = "\n\n".join([
    header.rstrip(),
    intro,
    methods,
    results,
    discussion,
    figures_section(),
    table1_section(),
    "# References\n",
])
(HERE / "manuscript_PREVIEW.md").write_text(out + "\n", encoding="utf-8")
print("wrote manuscript_PREVIEW.md")
