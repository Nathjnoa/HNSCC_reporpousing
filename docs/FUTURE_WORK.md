# Guía de Revisión Bibliográfica — HNSCC Drug Repurposing

El pipeline computacional (scripts 01–15 + 17) está completo. El panel definitivo
son los **23 candidatos LOD-stable** (`results/tables/15_lod_stability.tsv`).
Los 12 candidatos restantes del pool top-35 no-LOD-stable van a Supplementary Table.

Esta guía organiza la búsqueda bibliográfica por grupo mecanístico.

---

## Panel definitivo: 23 candidatos LOD-stable

| Clase | Fármacos |
|-------|----------|
| **A — Aprobado HNSCC** | Cetuximab |
| **B — Aprobado otro cáncer** | Afatinib, Amivantamab, Azacitidine, Crizotinib, Dacomitinib, Decitabine, Forodesine, Gefitinib, Lapatinib, Lazertinib, Mobocertinib, Necitumumab, Neratinib, Olmutinib, Osimertinib, Panitumumab, Vandetanib |
| **C — No oncológico aprobado** | Digoxin, Metformin, Mitapivat, Tranylcypromine, Valproic Acid |

---

## Grupos mecanísticos y preguntas clave

### 1. Inhibidores EGFR/HER (15 candidatos)
**Fármacos**: Cetuximab, Afatinib, Amivantamab, Dacomitinib, Gefitinib, Lapatinib,
Lazertinib, Mobocertinib, Necitumumab, Neratinib, Olmutinib, Osimertinib, Panitumumab,
Vandetanib, Crizotinib

**Preguntas a responder:**
- ¿Qué inhibidores EGFR más allá de Cetuximab/Afatinib tienen evidencia clínica en HNSCC?
- ¿Hay datos de resistencia a Cetuximab y sensibilización con inhibidores de segunda/tercera generación?
- ¿EGFR amplification/mutation rates en HNSCC HPV+ vs HPV-? (relevante para estratificación)
- Cetuximab es control positivo de clase A — sirve como ancla para calibrar la narrativa

**Búsqueda sugerida:**
```
("HNSCC" OR "head and neck squamous") AND ("EGFR inhibitor" OR "cetuximab" OR "afatinib") AND ("resistance" OR "clinical trial")
```

---

### 2. OXPHOS / Metabolismo (1 candidato principal)
**Fármacos**: Metformin (Complex I), Mitapivat (piruvato quinasa PKLR)

**Preguntas a responder:**
- ¿Metformina tiene evidencia clínica o epidemiológica en HNSCC?
- ¿NDUFA13 u otras subunidades del Complejo I tienen valor pronóstico en HNSCC?
- ¿Cuál es la racionalidad de inhibir el Complejo I en cáncer HPV+?
- Mitapivat (activador PKLR): PKLR está DE en nuestra proteómica — ¿hay datos en cáncer?

**Búsqueda sugerida:**
```
("metformin" OR "OXPHOS" OR "Complex I") AND ("HNSCC" OR "head and neck") AND ("repurposing" OR "antitumor")
```

---

### 3. Glucósidos cardíacos (1 candidato)
**Fármacos**: Digoxin (Na+/K+-ATPase, ATP1A1)

**Preguntas a responder:**
- ¿Digoxina tiene evidencia antitumoral en HNSCC u otros carcinomas escamosos?
- ¿ATP1A1 está documentado como target en cáncer?
- Epidemiología: ¿pacientes con digoxina tienen menor incidencia de HNSCC?

**Búsqueda sugerida:**
```
("digoxin" OR "cardiac glycoside") AND ("cancer" OR "carcinoma") AND ("antitumor" OR "repurposing")
```

---

### 4. Inhibidores DNMT / Epigenética (3 candidatos)
**Fármacos**: Decitabine, Azacitidine (DNMT1), Valproic Acid (HDAC/ALDH5A1)

**Preguntas a responder:**
- ¿Hay evidencia de hipermetilación en HNSCC que justifique inhibidores DNMT?
- ¿Decitabine o Azacitidine tienen datos en tumores sólidos (más allá de hematológicos)?
- Ácido Valproico: evidencia como agente epigenético en HNSCC o sinérgico con otros agentes

**Búsqueda sugerida:**
```
("decitabine" OR "azacitidine" OR "valproic acid") AND ("HNSCC" OR "head and neck") AND ("epigenetic" OR "methylation")
```

---

### 5. Inhibidores de enzimas específicas (2 candidatos)
**Fármacos**: Tranylcypromine (MAO-A/LSD1), Forodesine (purina nucleósido fosforilasa, PNP)

**Preguntas a responder:**
- Tranylcypromine: ¿LSD1 (KDM1A) está sobreexpresado en HNSCC? ¿Hay datos en carcinomas?
- Forodesine: mecanismo en cánceres sólidos — PNP es un target metabólico, ¿relevante en HNSCC?

**Búsqueda sugerida:**
```
("tranylcypromine" OR "LSD1 inhibitor" OR "KDM1A") AND ("cancer" OR "carcinoma")
("forodesine" OR "PNP inhibitor") AND ("cancer")
```

---

## Estrategia general de búsqueda

**Base de datos primaria**: PubMed
**Período**: 2015–2026
**Idioma**: inglés

**Búsqueda general de drug repurposing en HNSCC:**
```
("head and neck" OR "HNSCC") AND ("drug repurposing" OR "drug repositioning")
AND ("proteomic" OR "transcriptomic" OR "computational") AND (2015:2026[dp])
```

**Para cada candidato individual:**
```
("<nombre del fármaco>") AND ("HNSCC" OR "head and neck squamous" OR "oral cancer")
```

---

## Prioridad de revisión

| Prioridad | Fármacos | Razón |
|-----------|----------|-------|
| Alta | Cetuximab, Gefitinib, Afatinib, Metformin | Candidatos top con mayor evidencia clínica previa |
| Media | Decitabine, Azacitidine, Digoxin, Valproic Acid | Mecanismo relevante pero evidencia en HNSCC escasa |
| Baja | Tranylcypromine, Forodesine, Mitapivat | Candidatos novedosos — la novedad es el argumento |

---

## Output esperado de la revisión

Para cada candidato, documentar:
1. Evidencia preclínica en HNSCC (líneas celulares, modelos in vivo)
2. Trials clínicos activos o completados en HNSCC
3. Mecanismo de acción y conexión con la proteómica (¿cuál proteína DE es el target?)
4. Grado de novedad (¿ya reportado en HNSCC computacional?)

Consolidar en tabla para sección de Resultados/Discusión del manuscrito.

---

*Actualizado: 2026-04-05 — Pipeline computacional completo (v3); panel definitivo 23 candidatos LOD-stable*
