# 🧬 Bioinfo Promoter Analyser Pro

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](http://mabudo.divasecret.co.mz)

**Bioinfo Promoter Analyser Pro** é uma ferramenta de visualização de dados microbiológicos e análise computacional, focada no mapeamento físico de promotores genómicos, avaliação de impacto mutacional e análise fenotípica de isolados clínicos (ex: *Candida glabrata*).

*(Scroll down for English version).*

---

## 🇵🇹 Versão em Português

A aplicação permite o cruzamento massivo de perfis mutacionais (VCFs) com genomas de referência (FASTA/GFF) para identificar perdas e ganhos de função nos locais de ligação de Fatores de Transcrição (TFs), calculando a significância estatística de resistência associada a fármacos como o Fluconazol.

### 🌐 Acesso à Aplicação
A ferramenta está alojada na nuvem e pode ser acedida através do domínio clínico:
**[http://mabudo.divasecret.co.mz](http://mabudo.divasecret.co.mz)**

*(Nota: A aplicação encontra-se atualmente em fase de desenvolvimento e requer autenticação restrita).*

### 🛠️ Como Utilizar
A aplicação lê automaticamente o genoma de referência (`.fasta`) e as anotações (`.gff`) diretamente deste repositório GitHub. Para realizar análises customizadas, os investigadores devem interagir com a barra lateral da aplicação:
1. **Ficheiro VCF Original:** Carregar o ficheiro de Variantes Clínicas (`.vcf` ou `.vcf.gz`). O sistema está otimizado (via processamento em *chunks* iterativos) para ler ficheiros maciços diretamente na RAM de forma eficiente.
2. **Ficheiro de Metadados (Opcional, mas recomendado):** Carregar a matriz clínica (`.csv` ou `.xlsx`) contendo a coluna correspondente ao ID dos isolados e aos valores de resistência MIC50 (ex: Fluconazol).

### ✨ Funcionalidades Principais
* **Reguloma Pan & Core:** Identificação paramétrica se um determinado Fator de Transcrição regula um conjunto fixo de genes na população (Core Regulon) ou de forma variável (Accessory Regulon).
* **Mutagénese in Silico:** Teste de disrupção de Footprints biológicos baseados nas evidências empíricas da base de dados YEASTRACT.
* **Associação Estatística:** Cálculo automatizado do Teste Exato de Fisher, Qui-Quadrado, Teste T e Mann-Whitney para validação da correlação Genótipo-Fenótipo.

---

## 🇬🇧 English Version

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
