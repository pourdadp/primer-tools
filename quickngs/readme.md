
# 🧬 QuickNGS – Smart NGS & Sanger Analysis Suite

A smart, browser‑based platform for sequencing data analysis.  
Just upload your files—QuickNGS figures out the rest.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Flask](https://img.shields.io/badge/Flask-2.3+-lightgrey.svg)](https://flask.palletsprojects.com/)

---

## 🎯 What It Does

| You Upload | QuickNGS Detects | Analysis Performed |
|-----------|------------------|-------------------|
| 2 FASTQ files (R1 & R2) | **NGS Data** | Reference‑Based Pipeline (FastQC → BWA → Samtools → FreeBayes → SnpEff) |
| AB1 / FASTA / seq files | **Sanger Reads** | De Novo Assembly (Greedy / Semi‑Global / De Bruijn) |
| AB1 + Reference FASTA | **Sanger + Reference** | Guided Assembly (BWA alignment + consensus) |

---

## ✨ Features

### Core
- 🔍 **Automatic Input Detection** – No need to select mode manually
- 🧬 **NGS Pipeline** – Full clinical report with variant annotation
- 🧬 **Sanger Assembly** – Three algorithms + intelligent clustering
- 🧬 **Guided Assembly** – Align to reference and call consensus
- ✂️ **Quality Trimming** – Phred‑based trimming for AB1 files
- 🧬 **Translation** – 5 methods, 6 genetic codes, CAI scoring

### User Experience
- 📱 **Mobile‑First Design** – Works on phone, tablet, desktop
- 💡 **Built‑in Help** – Tooltip explanations for every parameter
- 📊 **Visual Results** – Contig graphs, variant tables, exportable
- 💾 **Run History** – All analyses saved, with delete option
- 📝 **Smart Project Naming** – Auto‑generated, editable

### Technical
- 💾 **Auto Storage Management** – Uses largest available disk
- 🛡️ **User‑Friendly Errors** – No technical jargon
- 📂 **Docker Ready** – Single‑command deployment
- 🔧 **Extensible** – Modular design (assembly, translation separate)

---

## 🚀 Quick Start

### Prerequisites
- Python 3.8+
- Bioinformatics tools: bwa, samtools, freebayes, snpeff, fastqc, trimmomatic, bcftools

### Option 1: Direct (WSL2/Linux)
```bash
# Install tools
sudo apt update
sudo apt install bwa samtools freebayes snpeff fastqc trimmomatic bcftools -y

# Install Python packages
pip install -r requirements.txt

# Run
python3 app.py
# Open http://127.0.0.1:5002
```

Option 2: Docker

```bash
docker-compose up -d

# Open http://localhost:5002
```

---

📂 Project Structure

```
quickngs/
├── app.py                 # Flask web application
├── assembly.py            # Assembly engine (Greedy, Semi‑Global, De Bruijn, Clustering, Merging, Guided)
├── translation.py         # Translation & Codon Usage (CAI) module
├── requirements.txt       # Python dependencies
├── Dockerfile             # Docker image definition
├── docker-compose.yml     # Docker Compose configuration
├── templates/             # HTML templates
│   ├── index.html         # Smart upload page
│   ├── results.html       # Progress page (NGS)
│   ├── results_final.html # Final report (NGS)
│   ├── assemble_result.html # Assembly results (Sanger)
│   └── history.html       # Run history
├── test_data/             # Sample data for testing
└── TODO.md                # Future improvements
```

---

🧬 Supported File Types

Extension Type Used For
.fastq, .fq, .gz NGS reads Reference‑Based NGS
.ab1 Sanger trace De Novo / Guided Assembly
.fasta, .fa, .fna FASTA sequences Sanger assembly
.seq, .txt Plain text Sanger assembly

---

🔧 Tools Used

Tool Purpose
BWA Read alignment
Samtools SAM/BAM processing
FreeBayes Variant calling
SnpEff Variant annotation
FastQC Quality control
Trimmomatic Adapter trimming
BCFtools Consensus calling
Biopython AB1/FASTA parsing, translation

---
👨‍🔬 Author

Pourdad Panahi – Biotechnologist & Computational Biologist
18+ years wet‑lab experience (cell culture, real‑time PCR, virus cultivation, ELISA)
Building digital tools for the life sciences.

· GitHub: github.com/pourdadp
· Portfolio: pourdadp.github.io

---

📄 License

MIT License – see LICENSE file
