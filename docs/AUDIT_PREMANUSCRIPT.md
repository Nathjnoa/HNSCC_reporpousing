# Auditoría pre-manuscrito — HNSCC Drug Repurposing

**Fecha:** 2026-06-13 · **Auditor:** revisión independiente tipo PI externo (sin sesgo de
confirmación) · **Calibración:** revista especializada Q1/Q2 (oncología / bioinformática
traslacional) · **Alcance:** lógica analítica, scripts, resultados, figuras y tablas del pipeline
completo (01–18).

Metodología de la auditoría: (1) tres exploraciones independientes del código/resultados/figuras;
(2) **siete re-verificaciones empíricas** ejecutadas sobre los datos reales (`scripts/audit_rechecks.R`
→ `results/audit/`); (3) verificación de literatura en PubMed para los claims clínicos. Cada hallazgo
se etiqueta **✅ VERIFICADO empíricamente** o **⚠️ Asumido / pendiente** conforme a la regla
`06-clinical-evidence-rigor`.

---

## 0. Estado de implementación (2026-06-13)

Fixes **aplicados** en esta sesión y pipeline re-corrido (10→18):

- **M2** (✅ aplicado): `scripts/10_prioritization_scoring.R` renormaliza DrugViability cuando falta
  firma L2S2. Efecto verificado: **Cetuximab/Nimotuzumab/Panitumumab pasan a rank 1** (composite
  0.730 > gefitinib 0.717); Tab3/Tab4/Tab5 y Fig5 actualizadas.
- **M5** (✅ aplicado): Fig3A re-etiquetada "Highest clinical phase (any indication)".
- **M3** (✅ aplicado, ajustado): el plot de Fig4C vuelve a top-1 por módulo (legibilidad); el
  **top-3 se conserva como tabla suplementaria** (`network/17_module_enrichment_top3.tsv`, 70 filas).
  La circularidad se atiende por texto (anotación descriptiva, no "validación").
- **m1** (✅ aplicado): `set.seed(2026)` antes de DESeq2 en `scripts/16_external_validation.R`.
- **Borradores de texto** (B1/B2/M1/M4/M6): en `docs/manuscript/AUDIT_TEXT_DRAFTS.md` (con citas
  PubMed verificadas) — para integrar en la fase de escritura.

**Pendiente de decisión del autor (no son fixes de código):** resolver el framing de **B2**
(marco metabólico de Metformina) y declarar **M1/M4/M6/B1** en el manuscrito usando los borradores.

> Nota: tras re-correr, la corroboración externa se mantiene (TCGA 11/14 concordantes, 12/14
> FDR<0.05; CPTAC sin cambios). El núcleo del análisis es estable a los fixes.

---

## 1. Veredicto

**El material está sustancialmente listo para escribir, condicionado a resolver UN punto bloqueante
de framing (B2) y a incorporar 6 ajustes mayores de honestidad/transparencia.** Ninguno de los
hallazgos invalida el análisis ni obliga a re-correr el pipeline. El núcleo científico —recuperación
de EGFR, supresión robusta de OXPHOS, hubs de proteasoma, validación externa en dos cohortes— **es
robusto y sobrevivió a todas las pruebas de estrés empíricas** que diseñé para tumbarlo.

El hallazgo originalmente sospechado como bloqueante (**contaminación tisular**, B1) fue
**empíricamente acotado a "mayor"**: existe, pero no contamina las conclusiones priorizadas.

---

## 2. Tabla de hallazgos (severity-ranked)

| ID | Severidad | Hallazgo | Estado empírico |
| -- | --------- | -------- | --------------- |
| **B2** | 🔴 Bloqueante (framing) | El candidato titular (Metformina) contradice el racional direccional del propio score: NDUFS2/Complejo I está **fuertemente DOWN**, pero el score "premia inhibir dianas UP". Metformina se ranquea alto solo por centralidad de red. | ✅ VERIFICADO + literatura |
| **M1** | 🟠 Mayor | Sesgo direccional de π no declarado; en tensión directa con B2. | ✅ VERIFICADO |
| **M2** | 🟠 Mayor | L2S2 `NA→0` penaliza anticuerpos aprobados: cetuximab/panitumumab caen de **rank 1 → rank 9**. | ✅ VERIFICADO |
| **M3** | 🟠 Mayor | Circularidad en Fig4C: nombres de módulo derivados de GO, "validados" con GO. | ✅ VERIFICADO |
| **M4** | 🟠 Mayor | "Validación" sobre-vendida: corrobora dianas, no el ranking; 3–4 anclas validan débil/discordante. | ✅ VERIFICADO |
| **B1** | 🟠 Mayor (era 🔴) | Composición tisular (músculo/mucosa) domina las proteínas down extremas, pero **NO** altera hubs ni OXPHOS. | ✅ VERIFICADO (defused) |
| **M5** | 🟠 Mayor | "Approved" en Fig3A = cualquier indicación, no HNSCC. | ✅ VERIFICADO |
| **M6** | 🟠 Mayor | El composite está dominado por centralidad (Spearman ≈0.78–0.81); el aporte drug-level es modesto. | ✅ VERIFICADO |
| **m1** | 🟡 Menor | `set.seed` ausente antes de DESeq2 (script 16); fechas de query API no registradas. | revisión código |
| **m2** | 🟡 Menor | Sesgo de imputación MinProb no testeado formalmente. | parcial (ver B1) |
| **m3** | 🟡 Menor | Transparencia de etiquetado en volcano/heatmap; nº de gene-sets testeados en GSEA. | revisión figuras |
| **m4** | 🟡 Menor | Clases A/B/C/D heurísticas + overrides manuales; pesos sin base farmacológica. | revisión código |
| **m5** | 🟡 Menor | Consistencia de logFC de EGFR entre figuras (dato=2.23 vs etiquetas previas). | verificar |
| **m6** | 🟡 Menor | Claim de metilación (CDKN2A) ⚠️ NO-VERIFICADO en `CLINICAL_CRITERIA.md:13`; script 19 removido. | evitar huérfanos |

---

## 3. Hallazgos detallados

### 🔴 B2 — Incoherencia entre el racional del score y el candidato titular (Metformina)

**Claim auditado.** El manuscrito posiciona Metformina (Complejo I OXPHOS) como candidato novedoso
top, junto a EGFR.

**Evidencia empírica (✅ VERIFICADO).**
- Complejo I está **fuertemente subexpresado** en el proteoma tumoral: NDUFS2 logFC=**−3.05**
  (FDR 3.3e-4), NDUFS3 −3.11, NDUFA9 −1.60 (`results/tables/de_limma/01_TVsS_all_proteins.tsv`).
- La supresión de OXPHOS **es real, no artefacto de composición ni de n pequeño**: replica en CPTAC
  (NDUFS2 logFC=−0.64, FDR=1e-18) y TCGA (−0.60, FDR=8e-18) (`results/audit/recheck5_fig6c_per_gene.tsv`),
  y la firma Hallmark OXPHOS se **refuerza** al remover contaminantes (NES −2.21 → −2.32,
  `results/audit/recheck2_gsea_sensitivity.tsv`).
- El score `scoring_v3` declara premiar **inhibir dianas UP** (`scripts/10_prioritization_scoring.R:159-165`;
  `docs/METHODS.md:~87`). Para NDUFS2, `primary_pi` es bajo (diana down) → Metformina **no** gana por
  ese racional, sino **solo por centralidad de red** (NDUFS2 es hub de alto grado, deg=51).

**El problema (perspectiva PI).** Tal como está, la narrativa se autocontradice: el método dice
"inhibimos lo sobreexpresado", pero su segundo titular es un **inhibidor de un complejo suprimido**,
rescatado por topología. Un revisor preguntará: *"si Complejo I ya está bajo, ¿por qué inhibirlo más?"*

**Resolución (literatura — según PubMed).** El candidato **es defendible**, pero con otro mecanismo:
- Metformina es inhibidor (suave) de Complejo I y OXPHOS es "una vulnerabilidad al ser blanco en
  células cancerosas"; sin embargo los ensayos clínicos de metformina "no fueron como se esperaba"
  (Pujalte-Martin et al. 2024, *Mol Oncol*, [DOI](https://doi.org/10.1002/1878-0261.13583)).
- Existe precedente clínico **específico en HNSCC**: metformina inhibe Complejo I/OXPHOS y se asoció
  a mejores desenlaces, con efecto diferencial HPV−/HPV+ (Curry et al. 2018, *Front Oncol*,
  [DOI](https://doi.org/10.3389/fonc.2018.00436)); su efecto antitumoral en HNSCC es en parte
  **inmunometabólico y AMPK-independiente** (NK/CXCL1) (Crist et al. 2022, *J Immunother Cancer*,
  [DOI](https://doi.org/10.1136/jitc-2022-005632)); ensayo Fase I/II como radiosensibilizador
  (Kemnade et al. 2023, *Oral Oncol*, [DOI](https://doi.org/10.1016/j.oraloncology.2023.106536)).

**Recomendación.** Antes de escribir, decidir y declarar el marco para dianas metabólicas:
1. **No** presentar OXPHOS/Complejo I como "sobreexpresado a inhibir" — está **down** (dato propio +
   2 cohortes). Enmarcar como **vulnerabilidad/dependencia metabólica** y precedente clínico, no como
   lógica direccional.
2. Reconciliar con M1: declarar explícitamente que el score usa centralidad como ancla
   **independiente de la dirección**, y que para dianas down la justificación es topológica +
   evidencia clínica externa, no el π direccional.
3. Citar el matiz honesto de que la monoterapia con metformina ha sido decepcionante; el valor está
   en combinación (quimio-radio/IO).

**Riesgo si no se corrige:** rechazo en revisión por inconsistencia lógica del titular.

---

### 🟠 M1 — Sesgo direccional de π no declarado

`π_directional` (min-max con signo, `10_prioritization_scoring.R:151-165`) premia inhibir
sobreexpresadas y **penaliza estructuralmente las subexpresadas** (tumor-supresores / dianas
loss-of-function). Es decisión de diseño legítima pero no está declarada, y entra en tensión con B2
(Metformina la evade vía centralidad). **Recomendación:** declarar el sesgo en Methods/Limitations y
reconciliarlo con el ancla de centralidad (ver B2). **Riesgo:** revisor lo lee como sesgo oculto.

---

### 🟠 M2 — `NA→0` de L2S2 penaliza anticuerpos aprobados (✅ VERIFICADO con números)

`s_reversal` asigna 0 cuando no hay firma L2S2 (`10_prioritization_scoring.R:260-265`). Como los
anticuerpos no se ensayan en líneas celulares (artefacto del ensayo, no biología), esto **deprime
fármacos aprobados en HNSCC**. Recomputando DrugViability con NA **excluido** (renormalizando
reg+breadth):

| Fármaco | composite (NA→0) | rank | composite (NA excluido) | rank |
| ------- | ---------------- | ---- | ----------------------- | ---- |
| Cetuximab / Panitumumab / Nimotuzumab | 0.583 | **9** | 0.730 | **1** |

(`results/audit/recheck7_l2s2_na_handling.tsv`; 149 fármacos tienen `cmap NA`.) Es decir, la
elección de imputación **mueve al cetuximab del top-1 al rank 9**. **Recomendación:** tratar `NA` como
faltante (renormalizar) o, como mínimo, reportar el sesgo y mostrar el ranking alternativo en
suplementario. **Riesgo:** un revisor clínico notará que el método sub-ranquea el anti-EGFR aprobado.

---

### 🟠 M3 — Circularidad en Fig4C (nombres de módulo ↔ GO)

Los nombres de módulo se derivan del término GO BP top-1 (`09_string_network.R`,
`network/17_module_enrichment_top3.tsv`) y Fig4C "valida" esos nombres con el mismo enriquecimiento
GO. **Recomendación:** quitar el framing de "validación"; presentar el enriquecimiento GO como
**anotación descriptiva** de módulos derivados de topología (Louvain), mostrando top-3 términos con
FDR. **Riesgo:** lectura de razonamiento circular.

---

### 🟠 M4 — "Validación" vs "corroboración"; anclas débiles (✅ VERIFICADO por gen)

Las cohortes externas confirman la **direccionalidad de las dianas**, no la eficacia de los fármacos
ni el **orden del ranking**. Además, la validación por gen (`results/audit/recheck5_fig6c_per_gene.tsv`)
no es uniforme:

- **Fuertes (concordantes + FDR<0.05 en ambas):** EGFR, NDUFS2, DNMT1, PSMA2, ALDH5A1, MAOA, MMP8,
  IMPDH2, ALDH2.
- **Débiles/discordantes:** PDE6D (pentoxifilina; CPTAC logFC≈0, p=0.91), CA2 (acetazolamida;
  discordante en TCGA), CDA (cedazuridina; discordante TCGA), RPS11 (MT-3724; discordante TCGA).

**Recomendación:** (1) cambiar "validación" → "corroboración/replicación de dianas" en texto;
(2) en Fig6C/Tab6 reportar n, logFC y FDR **por gen** y marcar las discordantes; (3) no contar
"14/14"; (4) atenuar el peso narrativo de pentoxifilina/acetazolamida/cedazuridina/MT-3724.
**Riesgo:** sobre-venta de la validación.

---

### 🟠 B1 — Composición tisular: real pero acotada (era 🔴, empíricamente degradada)

**Evidencia (✅ VERIFICADO).** Los 20 proteínas más down son músculo esquelético/cardíaco
(SLC25A4, MYLPF, MYL2, CASQ1, FHL1, MYH2/3/7, MYL1/3, PYGM, TCAP, LDB3, ATP2A1) + mucosa respiratoria
(BPIFB1, BPIFA1) — patrón de normal adyacente de lengua/laringe. Cuantificación
(`results/audit/recheck1_tissue_markers.tsv`):
- Marcadores tisulares: **12.8%** del set down (43/337) vs **0.6%** del up — pero **80% del top-20
  down** y **50% del top-50 down**. La contaminación domina los **extremos**, no el conjunto.

**Pruebas de estrés (✅ VERIFICADO — el análisis resiste):**
- **Red** (`recheck3`): removiendo los 38 genes contaminantes (498→445 nodos), **todos los hubs
  priorizados conservan estatus idéntico** (NDUFS2 deg 51→51, NDUFS3 52→52, EGFR 24→23, PSMA2 16→16).
  El músculo formó su propio módulo ("Muscle contraction", 29 genes) **sin inflar** a los hubs clave.
  **Ningún candidato priorizado está anclado a un contaminante.**
- **GSEA** (`recheck2`): al excluir contaminantes, **OXPHOS se refuerza** (NES −2.21→−2.32) y
  E2F/G2M/EMT/IFN-γ se mantienen. Solo **Myogenesis cae fuerte** (−2.00→−1.56) — única firma con
  contribución tisular sustancial; Adipogenesis/Heme cambian poco.

**Recomendación:** (1) reportar la composición tisular como **limitación honesta** y mostrar el
análisis de sensibilidad (ya hecho) como evidencia de robustez; (2) reformular Myogenesis (y atenuar
Adipogenesis/Heme) como **composición tisular del normal adyacente**, no biología tumoral;
(3) en volcano/heatmap, considerar marcar o excluir los contaminantes extremos para no destacarlos
como hallazgo. **Riesgo si se ignora:** revisor desconfía de todo el eje down; al mostrar la
sensibilidad, se neutraliza.

---

### 🟠 M5 — "Approved" en Fig3A engaña

345 "aprobados" mezcla off-label (metformina, acetazolamida…) con los 5 realmente aprobados en HNSCC
(`results/tables/drug_targets/08_multi_source_candidates.tsv`; `17e_fig3_multipanel.R`).
**Recomendación:** re-etiquetar a "estado de desarrollo clínico (cualquier indicación)" y separar
visualmente los aprobados-en-HNSCC. **Riesgo:** lectura de "345 fármacos listos para HNSCC".

---

### 🟠 M6 — El composite está dominado por la centralidad (✅ VERIFICADO)

Ablación (`results/audit/recheck6_baseline_ablation.tsv`): Spearman del composite vs **centralidad
sola = 0.781**, vs **π solo = 0.776**, vs TargetPriority = 0.811. (El solapamiento de top-14 quedó
inflado por empates —drogas que comparten diana comparten centralidad— así que el Spearman es la
medida fiable.) Es decir, el ranking es en gran medida un re-ordenamiento de la centralidad de red.
No es fatal, pero un revisor Q1/Q2 pedirá demostrar **valor agregado** sobre un baseline simple.
**Recomendación:** presentar esta comparación honestamente (la centralidad es el driver declarado), y
enmarcar el aporte drug-level como **desempate clínico/farmacológico**, no como discriminación
principal. Opcional: incluir el baseline centrality-only como suplementario.

---

### 🟡 Menores

- **m1** — Añadir `set.seed()` antes de DESeq2 en `16_external_validation.R`; registrar fechas de
  query API (scripts 04–07, 09, 16) y versiones de paquetes para reproducibilidad de snapshot.
- **m2** — Test formal de sesgo MinProb (p.ej. `missing_pct ~ dirección`); se solapa con B1.
- **m3** — Declarar criterio de etiquetado en volcano/heatmap y nº de gene-sets testeados en GSEA.
- **m4** — Documentar la curación de clases A/B/C/D y los overrides manuales
  (`08_integrate_drug_targets.R:367-442`); enmarcar 1.0/0.75/0.5/0.25 como **ordinal**, sin pretender
  base farmacológica de los valores exactos.
- **m5** — Verificar consistencia del logFC de EGFR entre dato (2.23) y etiquetas de figuras.
- **m6** — Asegurar que el claim de metilación (CDKN2A, `CLINICAL_CRITERIA.md:13`, ⚠️ NO-VERIFICADO)
  **no entre** al manuscrito, dado que el script 19 fue removido.

---

## 4. Lo que es robusto (confirmado bajo estrés)

- **Recuperación de EGFR** como control positivo: sólida (up en proteoma propio +2.23, CPTAC y TCGA;
  LOD-estable; top en las 6 configuraciones de pesos).
- **Supresión de OXPHOS**: real y replicada en 2 cohortes; se refuerza al limpiar contaminantes.
- **Hubs priorizados**: invariantes a la remoción de genes contaminantes.
- **DE pareado (limma + duplicateCorrelation, HPV como covariable)**: estadísticamente apropiado para
  10 pares.
- **Robustez (LOD + sensibilidad de pesos)**: bien diseñada; la eliminación de la permutación está
  justificada para un score de priorización.
- **Reproducibilidad local**: scripts deterministas, config centralizada, semillas en GSEA/Louvain.

---

## 5. Apéndice — resultados de los re-checks empíricos

Script: `scripts/audit_rechecks.R` · Salidas: `results/audit/` · Log: `results/audit/rechecks_log.txt`

1. **Composición tisular:** down=337, up=329; marcadores tisulares 12.8% down vs 0.6% up; 80% del
   top-20 down. → `recheck1_tissue_markers.tsv`
2. **Sensibilidad GSEA:** OXPHOS −2.21→−2.32 (refuerza); Myogenesis −2.00→−1.56 (cae); resto estable.
   → `recheck2_gsea_sensitivity.tsv`
3. **Sensibilidad de red:** hubs invariantes (NDUFS2/NDUFS3/EGFR/PSMA2 sin cambio). →
   `recheck3_network_hub_sensitivity.tsv`
4. **Literatura Metformina/Complejo I:** ver B2 (PubMed; 4 DOIs).
5. **Validación por gen:** CPTAC 11/12 concordantes (10 FDR<0.05); TCGA 11/14 (11 FDR<0.05);
   discordantes PDE6D, CA2, CDA, RPS11. → `recheck5_fig6c_per_gene.tsv`
6. **Ablación baseline:** Spearman composite~centralidad 0.781, ~π 0.776. →
   `recheck6_baseline_ablation.tsv`
7. **NA-handling L2S2:** cetuximab/panitumumab rank 9→1 al excluir NA. →
   `recheck7_l2s2_na_handling.tsv`

---

## 6. Checklist antes de escribir

Estado: ✅ resuelto · ✍️ texto redactado en `AUDIT_TEXT_DRAFTS.md`, pendiente integrar al manuscrito
al escribir · ⬜ abierto.

- ✅✍️ **B2** — **Revisado** (autor decide no delimitar titular; metformina = uno de varios candidatos
  de reposicionamiento). Framing aceptado: priorización por **centralidad (STRING, independiente del
  DE)**, no por inversión de firma; vía alternativa probable + precedente clínico; candidato de
  combinación. Texto listo (B2). Pendiente: integrar al escribir.
- ✅✍️ **M1** — **Confirmado como decisión de diseño** (no se cambia el método). π aplica un *prior*
  pro-inhibición de sobreexpresión (paradigma dominante en repurposing); la centralidad mantiene
  abiertas las dependencias independientes de dirección (p.ej. metabólicas). Declarar el prior + sus
  dos limitaciones (sintético-letal/metabólica; supresores LoF río abajo). Texto listo (M1).
- ✅ **M2** — Aplicado (renormalización L2S2; cetuximab → rank 1). Pipeline re-corrido.
- ✅✍️ **M4** — Fig6C/Tab6 ya reportan por gen; reescribir "validación"→"corroboración de dianas" y
  marcar anclas débiles (PDE6D, CA2, CDA, RPS11). Texto listo (M4).
- ✅✍️ **B1** — Análisis de sensibilidad hecho (hubs/OXPHOS robustos). Texto listo (B1: limitación +
  reformular Myogenesis). Pendiente: integrar.
- ✅ **M3** — Aplicado (Fig4C plot top-1 legible; top-3 en tabla supp). Reframe textual en M4/leyenda.
- ✅ **M5** — Aplicado (Fig3A re-etiquetada).
- ✅✍️ **M6** — Cuantificado (ρ≈0.78). Texto listo (M6: centralidad domina; drug-level = desempate).
- ✅ **m1** — Aplicado (`set.seed` DESeq2).
- ✅✍️ **m2** — MinProb: sesgo cuantificado (down 11.2% vs up 6.8% missing, p=2.5e-7; solo 1% >50%
  imputado; se solapa con B1). Texto listo (m2).
- ✅✍️ **m3** — Transparencia: 50 Hallmarks testeados → 11 sig; criterio de labels volcano/heatmap.
  Texto listo (m3).
- ✅✍️ **m4** — Clases A/B/C/D documentadas (ordinal) + 3 bloques de override. Texto listo (m4).
- ✅ **m5** — Resuelto: EGFR log2FC=2.23 único y consistente (no había inconsistencia real).
- ✅ **m6** — `CLINICAL_CRITERIA.md` fila 13 marcada FUERA DE ALCANCE (script 19 removido); no entra
  al manuscrito.

> Fixes de código aplicados: M2, M3, M5, m1 (pipeline 10→18 re-corrido). Borradores de texto:
> `docs/manuscript/AUDIT_TEXT_DRAFTS.md` (B1, B2, M1, M4, M6).
