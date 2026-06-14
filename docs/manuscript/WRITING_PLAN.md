# Plan de redacción — Manuscrito HNSCC Drug Repurposing

> **Cómo retomar en una sesión nueva:** abre Claude Code en
> `~/bioinfo/projects/hnscc_drug_repurposing` y di: *"Continuemos la redacción del manuscrito
> siguiendo `docs/manuscript/WRITING_PLAN.md`; empieza por el STORYBOARD (paso 0)."*
> Cada sección se redacta en prosa completa y pasa por aprobación antes de avanzar.

## Contexto

Auditoría pre-manuscrito cerrada (commit `938d133`): figuras Fig1–6 + 3 supp y 6 tablas fuente
(`results/tables/pub/`). **Display-item scheme del manuscrito (2026-06-14):** 1 tabla main =
shortlist priorizado (Tab4+Tab5 fundidas → Table 1); conteos cohorte/DE in-text; resto → Tables
S1–S5 (ver `tables.md` / `README.md`).
**finalizadas**, números actualizados (fix M2 → cetuximab #1), borradores de texto auditados con citas
verificadas. **No existe aún borrador de artículo.** Material reutilizable: tesis en español
(`docs/thesis/ARTICULO_2.0_REVJC.md`), `docs/METHODS.md` (production-ready), `docs/REFERENCES.md`,
`docs/manuscript/AUDIT_TEXT_DRAFTS.md`, y `docs/manuscript/README.md` (ángulo central + notas de
figure legends).

**Objetivo:** manuscrito IMRAD en inglés, submission-ready, ~3500–4500 palabras. Trabajo de
redacción/organización, **no analítico** (figuras/tablas no se regeneran).

**Decisiones del autor:** (1) borrador **agnóstico** de revista, se fija antes de formatear; (2)
**sección por sección con checkpoint** de aprobación; (3) **prosa completa en inglés lista para
editar**.

**Ángulo central (`manuscript/README.md`, no se fuerza titular):** EGFR como validación interna +
priorización hub-céntrica de dos niveles (`composite = 0.60·TargetPriority + 0.40·DrugViability`).
EGFR encabeza y recupera fármacos aprobados → coherencia del método; Tier 1 hub-central
(OXPHOS/metformina, proteasoma); Tier 2 periférico (epigenético, antimetabolitos, repurposing
clásico). Metformina = **uno de varios candidatos**, enmarcada como vulnerabilidad metabólica (B2).

## Orden de redacción (anclado en evidencia, minimiza re-trabajo)

0. **Storyboard / spine** (`docs/manuscript/STORYBOARD.md`, nuevo): claim central en 1 párrafo, mapa
   figura/tabla → mensaje, outline con presupuesto de palabras. **Checkpoint** antes de prosa.
1. **Leyendas de figuras y tablas** (`figures.md`, `tables.md`): notas detalladas ya en
   `manuscript/README.md` ("Nota FigX"). Autocontenidas, anclan Results.
2. **Results**: números ya existen en tablas. Integra framing M4 ("corroboración", no "validación") +
   sensibilidad B1.
3. **Methods**: `docs/METHODS.md` ya al día (scoring v3); condensar/traducir + integrar M1/M4/M6 +
   software.
4. **Introduction**: condensar/traducir intro de tesis (~2.4k→600–800 palabras) → brecha → objetivo.
5. **Discussion**: estructura de tesis + reencuadres auditados (B1, B2, M1, M6); pilares, limitaciones,
   perspectivas.
6. **Abstract** (~250 palabras): se destila del paper terminado.
7. **Title + keywords + highlights.**
8. **References** (`references.bib`): desde `REFERENCES.md` + DOIs de `AUDIT_TEXT_DRAFTS.md`; **cada
   cita clínica/de literatura se verifica vía PubMed al usarse** (regla `06-clinical-evidence-rigor`).
9. **Ensamblaje** `manuscript.md` (IMRAD completo) + `supplementary.md`.
10. **Fase final (post-aprobación):** fijar revista → formateo a su plantilla + cover letter.

## Mapa de reutilización (qué archivo alimenta qué sección)

- **Intro** ← `docs/thesis/ARTICULO_2.0_REVJC.md` (≈L201–314) — condensar + traducir.
- **Methods** ← `docs/METHODS.md` + `AUDIT_TEXT_DRAFTS.md` (M1/M4/M6) + sección Software.
- **Results** ← `results/tables/pub/` (Table 1 = Tab4+Tab5; S1–S5 = resto) + notas de `manuscript/README.md` +
  `AUDIT_TEXT_DRAFTS.md` (B1/M4) + `results/audit/`.
- **Discussion** ← tesis (≈L701–1158) + `AUDIT_TEXT_DRAFTS.md` (B1/B2/M1/M6).
- **Legends** ← `manuscript/README.md` ("Nota Fig2–6") + estructura de cada `TabN`.
- **References** ← `docs/REFERENCES.md` + DOIs auditados (Pujalte-Martin 2024, Curry 2018, Crist 2022,
  Kemnade 2023, Evangelista 2025).

## Guardrails (no negociables)

- **Números:** todo dato debe coincidir con tablas/figuras **post-auditoría** (cetuximab composite
  0.73 #1; CPTAC r=0.789/86.3%; TCGA r=0.601/76.2%). Verificar contra `results/tables/pub/` y
  `results/audit/`.
- **Citas:** afirmación clínica/de literatura con cita inline verificada vía PubMed; registrar
  criterios nuevos en `docs/CLINICAL_CRITERIA.md`.
- **Metformina:** marco B2 (vulnerabilidad metabólica + precedente clínico HNSCC; **no** "inhibir
  sobreexpresión"). No forzar como titular único.
- **Lenguaje:** "corroboration/replication of target directionality", no "validation" (M4).
- **Sesgo de confirmación:** dejar hablar al ranking; exclusiones solo por criterio uniforme.
- **No regenerar figuras/tablas** — finalizadas; si un número no cuadra, se reporta, no se re-corre.

## Archivos a crear (en `docs/manuscript/`)

`STORYBOARD.md` · `figures.md` · `tables.md` · `manuscript.md` · `supplementary.md` ·
`references.bib` · (fase final) `cover_letter.md`.

## Verificación (cómo sabemos que está listo)

- Consistencia numérica vs `results/tables/pub/*.tsv` y `results/audit/*`.
- Cada referencia verificada en PubMed (DOI/PMID); `CLINICAL_CRITERIA.md` actualizado.
- Cuerpo 3500–4500 palabras (excl. refs); abstract ≤250.
- Callouts de figura = leyendas; pilares (EGFR/Proteasome/Epigenetic/OXPHOS) consistentes.
- Cada sección aprobada por el autor antes de avanzar.

## Estado de progreso (actualizar al avanzar)

- [x] 0. STORYBOARD ✅ (`STORYBOARD.md` — claim central en 5 movimientos, encuadre clínico/traslacional, brechas hedged; aprobado 2026-06-14)
- [x] 1. Leyendas (`figures.md` Fig1–6+S1–3, `tables.md` Table1+S1–5) ✅ esquema 1-main-table propagado a README/WRITING_PLAN/STORYBOARD/AUDIT; aprobado 2026-06-14
- [x] 2. Results (`results.md`, R1–R6 ~1.4k palabras; R3 detalle básico-riguroso, números finos a Methods) ✅ aprobado 2026-06-14
- [~] 3. Methods (`methods.md`, ~1.0k palabras; números finos de red, M1/m4/M6, "corroboration") — redactado, en revisión del autor
- [~] 4. Introduction (`introduction.md`, ~720 palabras; brecha hedged, citas placeholder ⚠️ TO-VERIFY) — redactado, en revisión del autor
- [ ] 5. Discussion
- [ ] 6. Abstract
- [ ] 7. Title + keywords
- [ ] 8. references.bib
- [ ] 9. Ensamblaje manuscript.md + supplementary.md
- [ ] 10. Revista + formateo + cover letter
