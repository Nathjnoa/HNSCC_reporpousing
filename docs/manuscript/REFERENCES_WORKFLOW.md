# References workflow — Markdown → Word with **live Zotero fields**

> Author decision (2026-06-14): live Zotero fields in Word, via Zotero **RTF/ODF Scan**.
> This file fixes the mechanism so it is not re-litigated. Related: `WRITING_PLAN.md` step 8.

## Chosen path: RTF/ODF Scan (live fields)

1. **Drafting:** prose carries readable `[@key]` placeholders flagged `⚠️ TO-VERIFY`. At assembly
   (step 9) these become RTF/ODF-Scan markers `{Author, Year}` (e.g., `{Burtness, 2019}`).
2. **Convert:** `pandoc manuscript.md -o manuscript.odt` — to **ODT, no citeproc** (Zotero scans
   ODF/RTF, not docx).
3. **Scan (author, in Zotero):** Tools → RTF/ODF Scan → input `manuscript.odt`. Zotero matches each
   `{Author, Year}` marker against the library (interactive dialog for ambiguous ones) and writes a
   new ODT with **live Zotero citation fields**.
4. **Word:** open the scanned ODT in Word/LibreOffice with the Zotero plugin → live citations
   (refresh, restyle, auto-bibliography). Save as `.docx`. Co-authors can edit citations in Word.

## Prerequisites

- Every cited item must exist in the author's Zotero library → **step 8** compiles the
  PubMed-verified citation list (author, year, title, **PMID/DOI**) for import. Markers match only
  real library items.
- Recommend **Better BibTeX** with pinned citekeys for stability (scan still matches by
  Author–Year + dialog).

## Implication for assembly (step 9)

- Output is **`.odt` with `{Author, Year}` markers**, not `.docx` directly.
- `⚠️ TO-VERIFY [@key]` placeholders persist during drafting; materialize to real `{Author, Year}`
  only once verified (rule `06-clinical-evidence-rigor`).

## Open / to confirm before step 9

- RTF/ODF Scan exists and is standard, but exact menu steps may differ between Zotero 6 and 7.
  **Verify current steps against Zotero docs (WebFetch) before relying on it at step 9** — do not
  give from-memory menu instructions.
