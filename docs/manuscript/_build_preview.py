#!/usr/bin/env python3
"""Assemble a clean preview manuscript (Intro+Methods+Results) from the section drafts.

Strips internal note blocks, ⚠️ flags, markdown emphasis (plain output per author),
maps Introduction placeholder citekeys to real verified keys, and injects key method
citations. Produces manuscript_PREVIEW.md for pandoc → docx. Partial preview: Discussion,
Abstract, Title and full Methods/Results citations are pending (steps 5–8).
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


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
    "[@hnscc_epidemiology]": "[@johnson2020; @sung2021]",
    "[@hnscc_outcomes]": "[@johnson2020]",
    "[@keynote048]": "[@burtness2019]",
    "[@keynote689]": "[@imamura2025]",
    "[@hpv_subtypes]": "[@ang2010]",
    "[@proteomics_phenotype]": "[@huang2021]",
    "[@hnscc_proteomics]": "[@huang2021]",
    "[@repurposing_rationale]": "[@pushpakom2019]",
}

# --- Methods: inject method/tool citations at first single-token anchor ---------
METHOD_CITES = [
    ("proteoDA", "proteoDA [@ritchie2015]"),
    ("(GSEA)", "(GSEA) [@subramanian2005]"),
    ("clusterProfiler", "clusterProfiler [@wu2021]"),
    ("DGIdb", "DGIdb [@freshour2021]"),
    ("ChEMBL", "ChEMBL [@mendez2019]"),
    ("Open Targets", "Open Targets [@ochoa2021]"),
    ("L2S2", "L2S2 [@evangelista2025]"),
    ("STRING", "STRING [@szklarczyk2023]"),
    ("Hallmark", "Hallmark [@liberzon2015]"),
    ("Huang et al., 2021", "Huang et al., 2021 [@huang2021]"),
]


def inject_first(text, anchor, replacement):
    """Replace only the first occurrence of anchor with replacement."""
    return text.replace(anchor, replacement, 1)


intro = clean_common(body_from(HERE / "introduction.md"))
for k, v in INTRO_CITES.items():
    intro = intro.replace(k, v)

methods = clean_common(body_from(HERE / "methods.md"))
for anchor, repl in METHOD_CITES:
    methods = inject_first(methods, anchor, repl)

results = clean_common(body_from(HERE / "results.md"))

header = (
    "# Network-anchored proteomic drug repurposing in HNSCC (working draft, preview)\n\n"
    "Partial preview for pipeline testing: Introduction, Methods and Results only. "
    "Discussion, Abstract, Title and complete Methods/Results citations are pending "
    "(WRITING_PLAN.md steps 5–8). Citations shown are PubMed-verified; the bibliography "
    "is generated automatically by pandoc.\n\n"
)

out = "\n\n".join([header.rstrip(), intro, methods, results, "# References\n"])
(HERE / "manuscript_PREVIEW.md").write_text(out + "\n", encoding="utf-8")
print("wrote manuscript_PREVIEW.md")
