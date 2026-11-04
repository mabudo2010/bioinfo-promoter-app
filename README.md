# 🧬 Bioinfo Promoter Analyser Pro

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](http://mabudo.divasecret.co.mz)


**Bioinfo Promoter Analyser Pro** is a microbiological data visualization and computational analysis tool focused on the physical mapping of genomic promoters, evaluation of mutational impact, and phenotypic analysis of clinical isolates (e.g., *Candida glabrata*).

The application enables massive cross-referencing of mutational profiles (VCFs) against reference genomes (FASTA/GFF) to pinpoint loss/gain of function at Transcription Factor (TF) binding sites, calculating statistical resistance significance to drugs such as Fluconazole.

### 🌐 App Access
The tool is cloud-hosted and can be accessed via the clinical domain:
**[http://mabudo.divasecret.co.mz](http://mabudo.divasecret.co.mz)**

*(Note: The application is currently under development and requires restricted authentication).*

### 🛠️ How to Use
The application automatically pulls the reference genome (`.fasta`) and annotations (`.gff`) directly from this GitHub repository. To perform custom runs, researchers should interact with the application's sidebar:
1. **Original VCF File:** Upload the Clinical Variants file (`.vcf` or `.vcf.gz`). The system is optimized (via iterative chunk processing) to read massive files seamlessly.
2. **Metadata File (Optional but recommended):** Upload the clinical matrix (`.csv` or `.xlsx`) containing the isolate IDs and MIC50 resistance readouts (e.g., Fluconazole).

### ✨ Main Features
* **Pan & Core Regulomes:** Parametric identification classifying if a given TF regulates a fixed set of genes across the cohort (Core Regulon) or variably (Accessory Regulon).
* **In Silico Mutagenesis:** Disruption testing of biological footprints grounded on YEASTRACT empirical evidence.
* **Statistical Association:** Automated calculations for Fisher's Exact Test, Chi-Square, Student's t-Test, and Mann-Whitney U for robust Genotype-Phenotype validation.

---
**Autor / Author:** Martins dos Santos Abudo Mupueleque
