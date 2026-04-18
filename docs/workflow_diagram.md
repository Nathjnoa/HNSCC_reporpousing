# Diagrama de flujo — Reposicionamiento de fármacos en HNSCC

> **¿Qué hace este proyecto?**
> A partir de muestras de tejido tumoral de pacientes con cáncer de cabeza y cuello (HNSCC),
> identificamos qué proteínas están alteradas y usamos bases de datos farmacológicas para proponer
> qué fármacos ya aprobados podrían funcionar contra este cáncer —
> estrategia llamada **reposicionamiento de fármacos**.

---

```mermaid
flowchart TD

    DATA[("Punto de partida
    20 muestras · 10 Tumor / 10 Normal
    3 352 proteínas cuantificadas
    6 pares HPV+ · 4 pares HPV−")]

    PREP["Fase 1 — Preparación
    Modelo limma HPV-ajustado · diseño pareado
    666 proteínas significativas · 329 ↑ · 337 ↓"]

    DB1[("DGIdb v5")]
    DB2[("ChEMBL 34")]
    DB3[("Open Targets")]
    DB4[("L2S2")]

    INT["Fase 2 — Integración de fármacos
    4 bases de datos independientes
    Candidatos en ≥ 2 fuentes priorizados"]

    NET["Fase 3A — Red de proteínas
    STRING v12 · score ≥ 700
    Módulos Louvain · hubs por módulo"]

    SCO["Fase 3B — Puntuación multi-criterio
    5 dimensiones · π-stat · Clínica · L2S2 · Red · Rutas
    Pool top 35 para análisis LOD"]

    SENS["Fase 4 — Análisis de sensibilidad
    6 configuraciones de pesos
    LOD (leave-one-database) · n=1 000 permutaciones
    32 candidatos LOD-stable = panel final"]

    OUT[("Outputs de publicación
    7 figuras · 5 tablas main · 1 tabla supp")]

    DATA --> PREP
    PREP --> DB1 & DB2 & DB3 & DB4
    DB1 & DB2 & DB3 & DB4 --> INT
    INT --> NET --> SCO --> SENS --> OUT

    classDef input  fill:#0d5c8c,color:#ffffff,stroke:#063d5e,font-weight:bold
    classDef p1     fill:#dbeafe,color:#1e3a5f,stroke:#3b82f6
    classDef dbase  fill:#fef9c3,color:#78350f,stroke:#f59e0b
    classDef p3     fill:#f59e0b,color:#ffffff,stroke:#b45309,font-weight:bold
    classDef p4     fill:#dcfce7,color:#14532d,stroke:#16a34a
    classDef sens   fill:#f0fdf4,color:#14532d,stroke:#86efac
    classDef output fill:#22c55e,color:#ffffff,stroke:#15803d

    class DATA input
    class PREP p1
    class DB1,DB2,DB3,DB4 dbase
    class INT p3
    class NET,SCO p4
    class SENS sens
    class OUT output
```

---

## Descripción de cada fase

| Fase | Scripts | Qué se hace | Resultado clave |
| --- | --- | --- | --- |
| Preparación | 01, 02 | Modelo limma HPV-ajustado con `duplicateCorrelation()`. Mapeo UniProt → Entrez/Symbol | 666 proteínas DE (\|logFC\|>1, FDR<0.05) · 97.3 % mapeadas |
| Consulta BD | 04–07 | DGIdb, ChEMBL (fase ≥ 3), Open Targets, L2S2 (reversión transcriptómica) | Candidatos en ≥ 2 fuentes priorizados |
| Integración | 08 | Unificar 4 fuentes. Clasificar A/B/C/D. Multi-source candidates | Tabla maestra · 187 candidatos multi-fuente |
| Red PPI | 09 | STRING v12 REST API. Módulos Louvain. Hubs = top 10 % betweenness por módulo | Red + métricas de centralidad · hubs druggables |
| Scoring | 10 | 5 dimensiones ponderadas (ver `config/analysis_params.yaml`) | Pool top 35 candidatos |
| Sensibilidad | 15 | 6 configs de pesos + LOD + permutation test (n=1 000) | 32 candidatos LOD-stable = panel definitivo |
| Publicación | 17, 18 | Figuras y tablas publication-ready (PDF + PNG 300 DPI) | 7 figuras · 6 tablas |

---

Pipeline: 12 scripts (01–02, 04–10, 15, 17–18) · R + Python · Completado: 2026-03-08
