# References — HNSCC Drug Repurposing Pipeline

*Evidencia bibliográfica por fase del pipeline. Recuperado del commit 5963496.*
*Nota: la cita Parisi 2022 fue removida intencionalmente (commit e7ab397).*

---

## Background — Drug Repurposing & HNSCC

**Drug repurposing rationale:**

- Pushpakom S et al. (2019) Drug repurposing: progress, challenges and recommendations. *Nature Reviews Drug Discovery* 18:41–58. [PMID 30310233](https://pubmed.ncbi.nlm.nih.gov/30310233/) — Revisión sistemática del estado del arte del reposicionamiento, incluyendo estrategias computacionales basadas en datos ómicos.

**HNSCC epidemiology & clinical context:**

- Johnson DE et al. (2020) Head and neck squamous cell carcinoma. *Nature Reviews Disease Primers* 6:92. [PMID 33243986](https://pubmed.ncbi.nlm.nih.gov/33243986/) — Revisión clínica comprensiva sobre epidemiología, biología molecular y tratamiento del HNSCC.
- Sung H et al. (2021) Global Cancer Statistics 2020. *CA: A Cancer Journal for Clinicians* 71:209–249. [PMID 33538338](https://pubmed.ncbi.nlm.nih.gov/33538338/) — Datos epidemiológicos globales que ubican el HNSCC entre los 10 cánceres más frecuentes.

**Similar multi-modal computational repurposing pipelines:**

- Cheng F et al. (2019) [PMID 31596908](https://pubmed.ncbi.nlm.nih.gov/31596908/)
- Ye H et al. (2018) [PMID 29856687](https://pubmed.ncbi.nlm.nih.gov/29856687/)

---

## Data — DIA Proteomics & Differential Expression (Script 01)

**DIA mass spectrometry:**

- Ludwig C et al. (2018) Data-independent acquisition-based SWATH-MS for quantitative proteomics: a tutorial. *Molecular Systems Biology* 14:e8126. [PMID 29848443](https://pubmed.ncbi.nlm.nih.gov/29848443/) — Tutorial de referencia que establece DIA/SWATH-MS como estándar para proteómica cuantitativa con alta reproducibilidad y cobertura.
- Guo T et al. (2015) Rapid mass spectrometric conversion of tissue biopsy samples into permanent quantitative digital proteome maps. *Nature Medicine* 21:407–413. [PMID 25774927](https://pubmed.ncbi.nlm.nih.gov/25774927/) — Aplicabilidad de DIA-MS en biopsias tumorales para perfiles proteómicos cuantitativos reproducibles.

**limma / empirical Bayes:**

- Ritchie ME et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43:e47. [PMID 25605792](https://pubmed.ncbi.nlm.nih.gov/25605792/) — Artículo de referencia del paquete limma (>50,000 citas); valida la estimación empírica de Bayes para estudios con pocos replicados.
- Smyth GK (2004) Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. *Statistical Applications in Genetics and Molecular Biology* 3:Article3. [PMID 16646809](https://pubmed.ncbi.nlm.nih.gov/16646809/) — Artículo fundacional que introduce el moderado t-statistic y la base estadística de limma.

**FDR correction:**

- Benjamini Y & Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B* 57:289–300. — Artículo seminal de la corrección FDR, método estándar para control de errores en análisis ómicos.

---

## Phase 1: Identifier Mapping (Script 02)

- Wu T et al. (2021) clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2:100141. [PMID 34557778](https://pubmed.ncbi.nlm.nih.gov/34557778/) — Documenta la función `bitr()` para conversión de IDs entre sistemas de identificación de genes.
- Carlson M (2019) org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2. Bioconductor. — Base de datos de anotación del genoma humano; referencia estándar para mapeo de identificadores en Bioconductor.

---

## Phase 2: Functional Enrichment Analysis (Script 03)

**GSEA methodology:**

- Subramanian A et al. (2005) Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545–15550. [PMID 16199517](https://pubmed.ncbi.nlm.nih.gov/16199517/) — Artículo fundacional del GSEA (~50,000 citas); establece la metodología de ranking continuo y permutación.
- Mootha VK et al. (2003) PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. *Nature Genetics* 34:267–273. [PMID 12808457](https://pubmed.ncbi.nlm.nih.gov/12808457/) — Primera publicación del enfoque GSEA; introduce el concepto de detectar cambios coordinados en sets génicos.

**MSigDB Hallmarks:**

- Liberzon A et al. (2015) The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Systems* 1:417–425. [PMID 26771021](https://pubmed.ncbi.nlm.nih.gov/26771021/) — Define la colección Hallmarks como sets depurados de redundancia para capturar procesos biológicos bien definidos en cáncer.

**clusterProfiler (ORA + GSEA):**

- Wu T et al. (2021) clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2:100141. [PMID 34557778](https://pubmed.ncbi.nlm.nih.gov/34557778/) — Paquete R utilizado para ORA y GSEA.

**Background universe for ORA:**

- Huang DW et al. (2009) Systematic and integrative analysis of large gene lists using DAVID bioinformatics resources. *Nature Protocols* 4:44–57. [PMID 19131956](https://pubmed.ncbi.nlm.nih.gov/19131956/) — Discute la importancia del universo de fondo adecuado para análisis de enriquecimiento; justifica el uso de proteínas cuantificadas como background.

---

## Phase 3: Drug-Target Database Queries (Scripts 04–07)

**DGIdb (Script 04):**

- Freshour SL et al. (2021) Integration of the Drug-Gene Interaction Database (DGIdb 4.0) with open crowdsource efforts. *Nucleic Acids Research* 49:D1144–D1151. [PMID 33237278](https://pubmed.ncbi.nlm.nih.gov/33237278/) — Describe la arquitectura y cobertura de DGIdb v4/v5, incluyendo la integración de >30 fuentes de evidencia.

**ChEMBL (Script 05):**

- Mendez D et al. (2019) ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Research* 47:D930–D940. [PMID 30398643](https://pubmed.ncbi.nlm.nih.gov/30398643/) — Describe la infraestructura y cobertura de ChEMBL, incluyendo el sistema de fases clínicas que justifica el filtro max_phase ≥ 3.

**Open Targets (Script 06):**

- Ochoa D et al. (2021) Open Targets Platform: supporting systematic drug–target identification and prioritisation. *Nucleic Acids Research* 49:D1302–D1310. [PMID 33634751](https://pubmed.ncbi.nlm.nih.gov/33634751/) — Describe la integración de evidencia multi-fuente en Open Targets y su uso para priorización de targets, incluyendo la asociación a indicaciones oncológicas.

**CMap / signatureSearch (Script 07):**

- Lamb J et al. (2006) The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. *Science* 313:1929–1935. [PMID 17008526](https://pubmed.ncbi.nlm.nih.gov/17008526/) — Artículo fundacional de CMap; introduce el concepto de reversión de firma transcriptómica como estrategia de reposicionamiento.
- Duan Y et al. (2022) signatureSearch: environment for gene expression signature searching and functional interpretation. *Nucleic Acids Research* 50:e49. [PMID 35089344](https://pubmed.ncbi.nlm.nih.gov/35089344/) — Describe el paquete R signatureSearch utilizado para el análisis CMap2, incluyendo la función `gess_cmap()`.

---

## Phase 3: Drug Integration & Classification A–D (Script 08)

- Pushpakom S et al. (2019) [PMID 30310233](https://pubmed.ncbi.nlm.nih.gov/30310233/) — Propone clasificar candidatos de reposicionamiento según su nivel de evidencia clínica existente; fundamenta el esquema de clases A–D.

---

## Phase 4: PPI Network (Script 09)

- Szklarczyk D et al. (2023) The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any of 12,535 organisms. *Nucleic Acids Research* 51:D638–D646. [PMID 36370105](https://pubmed.ncbi.nlm.nih.gov/36370105/) — Referencia actual de STRING v12; describe el sistema de puntuación combinada y los umbrales de confianza (700 = alta confianza).
- Jeong H et al. (2001) Lethality and centrality in protein networks. *Nature* 411:41–42. [PMID 11333967](https://pubmed.ncbi.nlm.nih.gov/11333967/) — Demuestra que los hubs en redes PPI son letales cuando se eliminan; establece la base teórica para usar hubs como dianas terapéuticas.
- Barabási AL et al. (2011) Network medicine: a network-based approach to human disease. *Nature Reviews Genetics* 12:56–68. [PMID 21164525](https://pubmed.ncbi.nlm.nih.gov/21164525/) — Revisión fundacional de la medicina de redes; justifica el uso de centralidad y módulos de red para identificar dianas farmacológicas.

---

## Phase 4: Multi-Criteria Scoring (Script 10)

- Napolitano F et al. (2013) Drug repositioning: a machine-learning approach through data integration. *Journal of Cheminformatics* 5:30. [PMID 23830665](https://pubmed.ncbi.nlm.nih.gov/23830665/) — Establece el principio de integrar múltiples dimensiones de evidencia en un score compuesto; valida el enfoque de normalización y ponderación.
- Lotfi Shahreza M et al. (2018) A review of network-based approaches to drug repositioning. *Briefings in Bioinformatics* 19:878–892. [PMID 28334169](https://pubmed.ncbi.nlm.nih.gov/28334169/) — Sustenta la inclusión de centralidad de red como criterio de scoring; muestra que candidatos respaldados por topología de red tienen mayor tasa de validación.

---

## Phase 5: Evidence Validation (Scripts 11–12)

**ClinicalTrials.gov + PubMed (Script 11):**

- Zarin DA et al. (2011) The ClinicalTrials.gov results database — update and key issues. *New England Journal of Medicine* 364:852–860. [PMID 21366476](https://pubmed.ncbi.nlm.nih.gov/21366476/) — Describe ClinicalTrials.gov como repositorio de referencia para evidencia clínica de fármacos; justifica su uso como fuente de validación externa.
- Sayers EW et al. (2022) Database resources of the National Center for Biotechnology Information. *Nucleic Acids Research* 50:D20–D26. [PMID 34850941](https://pubmed.ncbi.nlm.nih.gov/34850941/) — Referencia oficial de las E-utilities de NCBI/PubMed utilizadas para la búsqueda bibliográfica automatizada.

**COSMIC CGC & NCG7 (Script 12):**

- Sondka Z et al. (2018) The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers. *Nature Reviews Cancer* 18:696–705. [PMID 30293088](https://pubmed.ncbi.nlm.nih.gov/30293088/) — Describe el Cancer Gene Census como catálogo más riguroso de genes driver en cáncer.
- Repana D et al. (2019) The Network of Cancer Genes (NCG): a comprehensive catalogue of known and candidate cancer genes from cancer sequencing screens. *Genome Biology* 20:1. [PMID 30606230](https://pubmed.ncbi.nlm.nih.gov/30606230/) — Describe NCG7 como recurso complementario para clasificación de genes driver.

---

## Phase 5: Final Integration (Script 13)

- Mullard A (2021) Drug repurposing programmes get lift off. *Nature Reviews Drug Discovery* 20:86–89. [PMID 33473228](https://pubmed.ncbi.nlm.nih.gov/33473228/) — Describe cómo la integración de evidencia computacional y clínica en un score único mejora la reproducibilidad y la comunicación de prioridades de reposicionamiento.

---

## Sensitivity Analysis (Script 15)

- Saltelli A et al. (2019) Why so many published sensitivity analyses are false: a systematic review of sensitivity analysis practices. *Environmental Modelling & Software* 114:29–39. — Establece principios para análisis de sensibilidad rigurosos en modelos computacionales, incluyendo variación sistemática de pesos en scores multi-criterio.
- Iqbal SA et al. (2016) Reproducible research practices and transparency across the biomedical literature. *PLOS Biology* 14:e1002333. [PMID 26726926](https://pubmed.ncbi.nlm.nih.gov/26726926/) — Argumenta que el análisis de sensibilidad es componente esencial de la reproducibilidad en investigación computacional biomédica.

---

## Publication Figures (Script 17)

- Okabe M & Ito K (2008) Color Universal Design (CUD): How to make figures and presentations that are friendly to colorblind people. *Journal of Japanese Society of Color Education.* — Define la paleta Okabe-Ito como estándar para visualización científica accesible (recomendada por Nature, Science, Cell).

---

*Para instrucciones de ejecución, ver `docs/RUNBOOK.md`. Para outputs generados, ver `docs/OUTPUTS.md`.*
