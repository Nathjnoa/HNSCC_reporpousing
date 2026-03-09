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
    3 352 proteínas cuantificadas
    6 pares HPV+ · 4 pares HPV−")]

    PREP["📂 Fase 1 — Preparación y control de calidad
    Modelo limma HPV-ajustado · diseño pareado
    666 proteínas significativas · 329 ↑ · 337 ↓"]

    BIO["🧬 Fase 2 — Biología del tumor
    GSEA con ranking π-estadístico
    Metabolismo oxidativo · Ciclo celular · ECM · Inmunidad"]

    DB1[("DGIdb
    2 697 fármacos")]

    DB2[("ChEMBL
    113 aprobados/fase III")]

    DB3[("Open Targets
    173 candidatos
    score HNSCC ≥ 0.2")]

    DB4[("L2S2
    1 044 drugs FDA")]

    INT["💊 Fase 3 — Integración de fármacos
    4 bases de datos con snapshot de fecha
    Candidatos en ≥ 2 fuentes priorizados"]

    NET["🕸️ Fase 4A — Red de proteínas
    498 proteínas · 2 698 conexiones
    22 módulos Louvain · 74 hubs por módulo"]

    SCO["Fase 4B — Puntuación multi-criterio
    6 criterios · π-stat · Clínica · L2S2 · Red · Rutas · Evidencia
    Diversidad de módulos penaliza blancos redundantes"]

    VAL["🏥 Fase 5 — Validación clínica
    ClinicalTrials.gov · PubMed · COSMIC/NCG7
    Candidatos con ensayos activos en HNSCC"]

    SENS["🔁 Fase 6 — Análisis de sensibilidad
    6 configuraciones de pesos · drop-one-database
    Permutation test n=1 000"]

    TOP["🏆 Resultado — Top 20 candidatos a reposicionamiento
    Ranking: 60 % puntuación · 40 % evidencia clínica
    🥇 Gefitinib · 🥈 Metformina · 🥉 Mavacamten
    9 candidatos altamente estables en análisis de sensibilidad"]

    OUT[("📄 Reporte final
    Excel 5 hojas · 13 figuras de publicación
    Niveles de evidencia 1–4")]

    DATA --> PREP
    PREP --> BIO
    PREP --> DB1 & DB2 & DB3 & DB4
    DB1 & DB2 & DB3 & DB4 --> INT
    BIO --> SCO
    INT --> NET --> SCO --> VAL --> SENS --> TOP --> OUT

    classDef input  fill:#0d5c8c,color:#ffffff,stroke:#063d5e,font-weight:bold
    classDef p1     fill:#dbeafe,color:#1e3a5f,stroke:#3b82f6
    classDef p2     fill:#ede9fe,color:#3b1f7a,stroke:#7c3aed
    classDef dbase  fill:#fef9c3,color:#78350f,stroke:#f59e0b
    classDef p3     fill:#f59e0b,color:#ffffff,stroke:#b45309,font-weight:bold
    classDef p4     fill:#dcfce7,color:#14532d,stroke:#16a34a
    classDef p5     fill:#fee2e2,color:#7f1d1d,stroke:#dc2626
    classDef sens   fill:#f0fdf4,color:#14532d,stroke:#86efac
    classDef result fill:#14532d,color:#ffffff,stroke:#052e16,font-weight:bold
    classDef output fill:#22c55e,color:#ffffff,stroke:#15803d

    class DATA input
    class PREP p1
    class BIO p2
    class DB1,DB2,DB3,DB4 dbase
    class INT p3
    class NET,SCO p4
    class VAL p5
    class SENS sens
    class TOP result
    class OUT output
```

---

## Descripción detallada de cada fase

| Fase | Paso | Qué se hace | Resultado clave |
| ------ | ------ | ------------- | ----------------- |
| 📂 Preparación | ① Control de calidad | Modelo limma HPV-ajustado con `duplicateCorrelation()` para diseño pareado tumor/normal | 666 proteínas significativas (|logFC|>1, FDR<0.05) · 329 ↑ · 337 ↓ |
| | ② Traducción de IDs | Convertir códigos internos UniProt a nombres de genes reconocibles | 3 262 de 3 352 proteínas mapeadas (97.3 %) |
| 🧬 Biología | ③ Análisis de rutas | GSEA con ranking π-estadístico (`sign(logFC) × |logFC| × −log10(FDR)`). ORA para GO, KEGG, Reactome | Metabolismo oxidativo · ciclo celular · matriz extracelular · inmunidad |
| 💊 Fármacos | ④ Consulta a bases de datos | Buscar en 4 bases de datos independientes qué fármacos conocidos actúan sobre las proteínas alteradas. Open Targets filtrado por score HNSCC ≥ 0.2 | Candidatos únicos en ≥ 2 fuentes priorizados |
| 🕸️ Red | ⑤ Red de proteínas | Red STRING ≥ 700. Detección de módulos Louvain (22 módulos). Hubs = top 10 % betweenness dentro de cada módulo | 498 proteínas · 2 698 conexiones · 74 hubs modulares |
| | ⑥ Puntuación integrada | 6 criterios: π-stat (0.25), clínica (0.20), rutas (0.15), red (0.15), evidencia (0.15), L2S2 (0.10) | Top 20 candidatos con ranking objetivo |
| 🏥 Validación | ⑦ Evidencia clínica | ClinicalTrials.gov + PubMed por candidato. Cruce con oncogenes COSMIC/NCG7 | Candidatos con ensayos activos en HNSCC identificados |
| 🔁 Sensibilidad | ⑧ Robustez del ranking | 6 configuraciones de pesos + drop-one-database + permutation test (n=1 000) | 9 candidatos altamente estables en todas las configuraciones |
| 🏆 Resultado | ⑨ Ranking final | Score combinado 60 % multi-criterio + 40 % evidencia clínica | Top 20 · #1 Gefitinib · #2 Metformina · #3 Mavacamten |

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
| 🟢 Verde muy claro | Fase 6 — Análisis de sensibilidad |
| 🟢 Verde oscuro | Resultado final |

---

*Pipeline completado: 2026-03-08 · 16 scripts (01–15, 17) · R + Python*
*Correcciones metodológicas aplicadas: 2026-03-08 (commits a46ca45–2a4e572)*
