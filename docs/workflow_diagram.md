# Diagrama de flujo — Reposicionamiento de fármacos en HNSCC

> **¿Qué hace este proyecto?**
> A partir de muestras de tejido tumoral de pacientes con cáncer de cabeza y cuello (HNSCC),
> identificamos qué proteínas están alteradas y usamos bases de datos farmacológicas para proponer
> qué fármacos ya aprobados podrían funcionar contra este cáncer —
> estrategia llamada **reposicionamiento de fármacos**.

---

```mermaid
flowchart TD

    DATA[("🔬 Punto de partida
    20 muestras · 10 Tumor / 10 Normal
    3 352 proteínas cuantificadas")]

    PREP["📂 Fase 1 — Preparación y control de calidad
    Verificar datos e identificar proteínas alteradas
    520 proteínas significativas · 248 aumentadas · 272 disminuidas"]

    BIO["🧬 Fase 2 — Biología del tumor
    Rutas celulares sobreactivadas o apagadas
    Metabolismo · Ciclo celular · ECM · Inmunidad"]

    DB1[("DGIdb
    2 252 fármacos")]

    DB2[("ChEMBL
    55 aprobados")]

    DB3[("Open Targets
    66 candidatos")]

    DB4[("L2S2
    1,044 drugs FDA")]

    INT["💊 Fase 3 — Integración de fármacos
    4 bases de datos consultadas en paralelo
    2 421 fármacos únicos · 187 en ≥ 2 fuentes"]

    NET["🕸️ Fase 4A — Red de proteínas
    403 proteínas · 2 001 conexiones
    Hubs: Complejo I mitocondrial → Metformina"]

    SCO["Fase 4B — Puntuación multi-criterio
    6 criterios · 177 candidatos evaluados
    Expresión · Clínica · L2S2 · Red · Rutas"]

    VAL["🏥 Fase 5 — Validación clínica
    11/20 candidatos con ensayos activos en HNSCC
    EGFR: gen driver principal · logFC = 4.33"]

    TOP["🏆 Resultado — Top 20 candidatos a reposicionamiento
    Ranking: 60 % puntuación · 40 % evidencia clínica
    🥇 Erlotinib · 🥈 Cetuximab · 🥉 Metformina"]

    OUT[("📄 Reporte final
    Excel 5 hojas · 14 figuras
    Niveles de evidencia 1–4")]

    DATA --> PREP
    PREP --> BIO
    PREP --> DB1 & DB2 & DB3 & DB4
    DB1 & DB2 & DB3 & DB4 --> INT
    BIO --> SCO
    INT --> NET --> SCO --> VAL --> TOP --> OUT

    classDef input  fill:#0d5c8c,color:#ffffff,stroke:#063d5e,font-weight:bold
    classDef p1     fill:#dbeafe,color:#1e3a5f,stroke:#3b82f6
    classDef p2     fill:#ede9fe,color:#3b1f7a,stroke:#7c3aed
    classDef dbase  fill:#fef9c3,color:#78350f,stroke:#f59e0b
    classDef p3     fill:#f59e0b,color:#ffffff,stroke:#b45309,font-weight:bold
    classDef p4     fill:#dcfce7,color:#14532d,stroke:#16a34a
    classDef p5     fill:#fee2e2,color:#7f1d1d,stroke:#dc2626
    classDef result fill:#14532d,color:#ffffff,stroke:#052e16,font-weight:bold
    classDef output fill:#22c55e,color:#ffffff,stroke:#15803d

    class DATA input
    class PREP p1
    class BIO p2
    class DB1,DB2,DB3,DB4 dbase
    class INT p3
    class NET,SCO p4
    class VAL p5
    class TOP result
    class OUT output
```

---

## Descripción detallada de cada fase

| Fase | Paso | Qué se hace | Resultado clave |
| ------ | ------ | ------------- | ----------------- |
| 📂 Preparación | ① Control de calidad | Verificar integridad, distribución y validez estadística de los datos proteómicos | 520 proteínas significativamente alteradas (248 ↑ · 272 ↓) |
| | ② Traducción de IDs | Convertir códigos internos UniProt a nombres de genes reconocibles | 3 263 de 3 352 proteínas mapeadas (97.3 %) |
| 🧬 Biología | ③ Análisis de rutas | Identificar qué funciones celulares están sobreactivadas o apagadas en el tumor usando GO, KEGG, Reactome y Hallmarks | Metabolismo mitocondrial · ciclo celular · matriz extracelular · inmunidad |
| 💊 Fármacos | ④ Consulta a bases de datos | Buscar en 4 bases de datos independientes qué fármacos conocidos actúan sobre las proteínas alteradas | 2 421 fármacos únicos · 187 respaldados por ≥ 2 fuentes |
| 🕸️ Red | ⑤ Red de proteínas | Construir mapa de conexiones entre las proteínas tumorales para identificar las más importantes (hubs) | 403 proteínas · 2 001 conexiones · Hubs = Complejo I mitocondrial |
| | ⑥ Puntuación integrada | Asignar un puntaje a cada candidato combinando 6 criterios independientes | 177 candidatos evaluados con ranking objetivo |
| 🏥 Validación | Ensayos clínicos | Buscar en ClinicalTrials.gov si los candidatos ya tienen estudios en HNSCC | 11 de 20 con ensayos · 8 activos en 2026 |
| | Genes driver | Cruzar los targets con oncogenes conocidos de cáncer de cabeza y cuello | EGFR: gen driver principal (logFC = 4.33) |
| 🏆 Resultado | ⑦ Ranking final | Combinar puntuación multi-criterio (60 %) con evidencia clínica (40 %) | Top 20 candidatos · #1 Erlotinib · #2 Cetuximab · #3 Metformina |

---

## Leyenda de colores

| Color | Fase |
|-------|------|
| 🔵 Azul oscuro | Datos de entrada |
| 🔵 Azul claro | Fase 1 — Preparación |
| 🟣 Morado | Fase 2 — Biología |
| 🟡 Amarillo | Fase 3 — Bases de datos de fármacos |
| 🟠 Naranja | Integración de fuentes |
| 🟢 Verde claro | Fase 4 — Red y puntuación |
| 🔴 Rojo claro | Fase 5 — Validación clínica |
| 🟢 Verde oscuro | Resultado final |

---

*Pipeline completado: 2026-03-04 · 17 scripts · R + Python*
