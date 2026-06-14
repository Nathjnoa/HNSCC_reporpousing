# References workflow — Markdown → Word with **live Zotero fields**

> Author decision (2026-06-14): live Zotero fields in Word, via Zotero **RTF/ODF Scan**.
> This file fixes the mechanism so it is not re-litigated. Related: `WRITING_PLAN.md` step 8.

## Chosen path: Zotero scan markers → live fields

Drafting carries readable `[@key]` placeholders flagged `⚠️ TO-VERIFY`. At assembly these become
scan markers `{Author, Year}` (e.g., `{Burtness, 2019}`).

**VERIFIED 2026-06-14 (Zotero docs + forums):** a scanned **`.odt` opens ONLY in LibreOffice** —
Word cannot open it (errors / flattens citations to plain text). The built-in Zotero
*Tools → RTF/ODF Scan* handles only RTF and ODF. **DOCX scanning requires the add-on
"ODF/DOCX Scan for Zotero"** (Juris-M: github.com/Juris-M/zotero-odf-scan-plugin).

### Path B — Word-native (chosen, since author uses Word)

**VERIFIED 2026-06-14 from the add-on dialog:** the ODF/DOCX Scan plugin does **not** accept plain
`{Author, Year}` markers. It accepts only (a) **Scannable Cite** `{ | cite | locator | suffix |
zu:USERID:ITEMKEY}` (needs the item URI/key) or (b) **pandoc syntax** `[@citekey]` which matches
the Zotero **Citation Key** field (provided by **Better BibTeX**). We use pandoc syntax.

1. Install **Better BibTeX** AND the **ODF/DOCX Scan** add-on in Zotero.
2. Import `references.bib` (BBT keeps citekeys like `burtness2019`). Verify item Citation Keys match
   the `[@key]` used in the document; if BBT regenerated them, rebuild the doc with the new keys.
3. `pandoc manuscript.md -f markdown-citations --reference-doc=reference.docx -o manuscript_markers.docx`
   (the `-citations` extension keeps `[@key]` as literal text; no citeproc).
4. Zotero → Tools → ODF/DOCX Scan → **"Pandoc → Zotero citations"** → input/output `.docx`.
5. Open the converted `.docx` in Word → live Zotero fields → set style, insert bibliography.

### Path A — LibreOffice (fallback / quick test, no add-on)
1. `pandoc manuscript_markers.md -o manuscript_markers.odt`.
2. Zotero → Tools → RTF/ODF Scan → ODF → produces an `.odt` with live fields.
3. **Open in LibreOffice (NOT Word).** Set Document Preferences, style, Insert Bibliography.
4. To reach Word: set Zotero "Store citations as" → **Bookmarks**, then Save As `.docx`
   (bookmarks are more fragile than reference marks if heavily edited).

## Prerequisites

- Every cited item must exist in the author's Zotero library → **step 8** compiles the
  PubMed-verified citation list (author, year, title, **PMID/DOI**) for import. Markers match only
  real library items.
- Recommend **Better BibTeX** with pinned citekeys for stability (scan still matches by
  Author–Year + dialog).

## Implication for assembly (step 9)

- Output is a **`.docx` with `{Author, Year}` markers** (Path B, Word-native) built with
  `reference.docx` styling; the author scans it with the ODF/DOCX Scan add-on.
- `⚠️ TO-VERIFY [@key]` placeholders persist during drafting; materialize to real `{Author, Year}`
  only once verified (rule `06-clinical-evidence-rigor`).

## Test status — ✅ VALIDATED (2026-06-14)

- Path B works end-to-end in **Word**: pandoc `[@citekey]` docx + Better BibTeX + ODF/DOCX Scan in
  **"Pandoc → Zotero citations"** mode produced a `.docx` with live Zotero fields that opens in Word.
- Failed attempts (removed): built-in RTF/ODF Scan → `.odt` (Word can't open ODT fields); plain
  `{Author, Year}` markers (the plugin only accepts Scannable Cite or pandoc syntax, not bare
  author-year).
- Working artifact kept: `manuscript_PREVIEW_pandoc_markers.docx` (pandoc `[@key]` tokens, input to
  the scan). Build: `pandoc <md> -f markdown-citations --reference-doc=reference.docx -o <out>.docx`.

## Fragile point (the one thing that can break Path B)

`[@citekey]` matches a Zotero item **only if the document key equals that item's Better BibTeX
Citation Key exactly**. Keys matched this time because BBT preserved the imported `.bib` keys. Risk
going forward: BBT's citekey format is configurable, so a newly added reference may get a different
key (e.g., `burtnessPembrolizumab2019` instead of `burtness2019`), and that citation would silently
fail to convert. Mitigation at step 8: either (a) pin citekeys in BBT to the values in
`references.bib`, or (b) the author shares the actual BBT keys / exports the library, and the
document is rebuilt with those exact keys before scanning.
