
# 🧬 Primer Tools – Complete Molecular Biology Suite

A comprehensive toolkit for molecular biology laboratories: primer management, PCR protocol design, targeted sequencing panel optimization, and NGS data analysis.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Flask](https://img.shields.io/badge/Flask-2.0-lightgrey.svg)](https://flask.palletsprojects.com/)
[![Bootstrap](https://img.shields.io/badge/Bootstrap-5-purple.svg)](https://getbootstrap.com/)
[![DeepSeek](https://img.shields.io/badge/Built_with-DeepSeek_AI-6C5CE7.svg)](https://deepseek.com)

---

## 📖 Overview

**Primer Tools** is a suite of three integrated web applications designed by a molecular biologist with 18+ years of laboratory experience. These tools address real-world needs in diagnostic and research labs:

| Tool | Purpose | Key Users |
|------|---------|-----------|
| 🧬 **Primer Database Manager (PDM)** | Centralized primer/probe inventory & PCR protocol design | Lab technicians, Researchers |
| 🧬 **Sequencing Panel Optimizer (SPO)** | Optimize primer panels for full-length targeted sequencing | Bioinformaticians, NGS specialists |
| 🧬 **QuickNGS** | One‑click NGS analysis pipeline (FASTQ → VCF → clinical report) | Researchers, Clinicians |

All three tools work together seamlessly: **manage primers in PDM → optimize panels in SPO → analyze NGS data with QuickNGS**.

---

## 🧬 Primer Database Manager (PDM)

### What it does
A full-featured LIMS-style application for managing primers, probes, and PCR programs with multi-user support.

### ✨ Highlights
- 👥 **User roles**: Admin, Editor, Viewer with secure authentication
- 🧪 **Primer CRUD** with detailed metadata (Gene, Organism, Tm, Amplicon length, etc.)
- 🔬 **Probe management** with Reporter/Quencher tracking
- 🧪 **PCR Program Designer** with cycle groups, temperature steps, and Real-Time read detection
- 📋 **Multiplex Panels** – Group multiple primer pairs into single reactions
- 🔒 **Editing Lock** – Prevents simultaneous edits; Admin can break locks
- 🏷️ **Custom Fields & Aliases** – Flexible metadata for any primer
- 💾 **One-click Backup & Restore** with auto-cleanup (30 days)
- 🖨️ **Print-optimized** pages (A4)
- 📊 **Audit Log** – Track all user actions

### 📂 Location
`/primer_database_manager/`

### 🚀 Quick Start
```bash
cd primer_database_manager
pip install -r requirements.txt
python app.py
# Open http://127.0.0.1:5001
# Login: admin / admin123
```

---

🧬 Sequencing Panel Optimizer (SPO)

What it does

Select the best combination of existing primers to generate overlapping amplicons for full-length targeted sequencing (e.g., viral genomes, gene panels).

✨ Highlights

· 🧬 Load primers manually or import directly from PDM
· 🔍 IUPAC support – Handles degenerate bases (R, Y, S, W, K, M, B, D, H, V, N) and Inosine (I)
· 🧮 Semi-global alignment – Detects binding with indels and 5′ tails
· ⚠️ Multi-binding detection – Warns if a primer binds to multiple sites
· 🚫 Both-strand rejection – Automatically rejects ambiguous primers
· ✅ Validation – Amplicon length, Tm compatibility, dimer formation checks
· 📊 Greedy optimization – Finds the optimal primer set for full coverage
· 📈 Visual binding map – Graphical output (PNG) of primer binding sites
· 📥 Downloadable JSON reports – Binding results, valid pairs, optimal panel
· 📜 History tracking – Logs all runs with status and parameters

📂 Location

/sequencing_panel_optimizer/

🚀 Quick Start

```bash
cd sequencing_panel_optimizer
pip install -r requirements.txt
python app.py
# Open http://127.0.0.1:5000
```

---

🧬 QuickNGS

What it does

A mobile‑first web pipeline that takes raw FASTQ files and produces a variant report (VCF + annotated HTML) – no command line required.

✨ Highlights

· 📱 Mobile‑first interface – touch‑friendly, collapsible cards, FAB button
· 📂 FASTQ upload – paired‑end reads with instant visual feedback
· 🧬 Reference selection – built‑in human/mouse genomes or custom FASTA
· ⚙️ Adjustable parameters – quality, depth, variant quality, adapter trimming
· 🎯 Optional target region – analyze a specific genomic interval
· ⏳ Live progress – step‑by‑step status with polling
· 📊 Coverage plot – interactive (Plotly) coverage chart in the final report
· 📥 Downloadable outputs – VCF, annotated variants, and report

📂 Location

/quickngs/

🚀 Quick Start

```bash
cd quickngs
pip install -r requirements.txt
python app.py
# Open http://127.0.0.1:5002
```

ℹ️ The current version includes a mock pipeline that simulates the NGS workflow. The full bioinformatics tools (BWA, FreeBayes, SnpEff) will be integrated via Docker in the next release.

---

🔄 Workflow Integration

```
┌──────────────────────────────────────────────────────────────────┐
│                       PRIMER TOOLS SUITE                          │
├───────────────────┬───────────────────┬───────────────────────────┤
│   PRIMER DATABASE │ SEQUENCING PANEL  │      QUICKNGS             │
│   MANAGER (PDM)   │ OPTIMIZER (SPO)   │                            │
├───────────────────┼───────────────────┼───────────────────────────┤
│  • Store primers  │  • Import primers │  • Upload FASTQ            │
│  • Design PCR     │    from PDM       │  • Select reference        │
│  • Create panels  │  • Optimize       │  • Analyze NGS data        │
│  • Manage stock   │    coverage       │  • Generate report         │
│  • Backup/restore │  • Visualize      │  • Export VCF              │
└───────────────────┴───────────────────┴───────────────────────────┘
         │                    ▲                       ▲
         │   Import primers   │                       │
         └────────────────────┘                       │
                   │                                  │
                   └──────────────────────────────────┘
                        Load designed panel to analyze its NGS data
```

---

📂 Repository Structure

```
primer-tools/
├── README.md                           # ← You are here
├── LICENSE
├── primer_database_manager/            # 🧬 PDM Application
│   ├── app.py
│   ├── database.py
│   ├── auth.py
│   ├── config.py
│   ├── requirements.txt
│   ├── install.bat
│   ├── run.bat
│   ├── templates/
│   ├── static/
│   └── backups/
├── sequencing_panel_optimizer/         # 🧬 SPO Application
│   ├── app.py
│   ├── requirements.txt
│   ├── install.bat
│   ├── run.bat
│   ├── templates/
│   ├── scripts/
│   └── results/
└── quickngs/                           # 🧬 QuickNGS Application
    ├── app.py
    ├── requirements.txt
    ├── config.yaml
    ├── .gitignore
    ├── templates/
    │   ├── base.html
    │   ├── index.html
    │   └── results.html
    └── static/
        └── style.css
```

---

⚙️ Installation

Prerequisites

· Python 3.8 or higher
· pip (Python package manager)

Option 1: Windows (Recommended)

```bash
# For PDM
cd primer_database_manager
install.bat
run.bat

# For SPO
cd sequencing_panel_optimizer
install.bat
run.bat

# For QuickNGS
cd quickngs
pip install -r requirements.txt
python app.py
```

💡 For users in Iran/China: The install.bat scripts for PDM and SPO use Tsinghua mirror for faster downloads. QuickNGS dependencies are lightweight and can be installed directly.

Option 2: Manual (All Platforms)

```bash
# PDM
cd primer_database_manager
pip install flask werkzeug
python app.py

# SPO
cd sequencing_panel_optimizer
pip install flask pyyaml matplotlib tqdm
python app.py

# QuickNGS
cd quickngs
pip install flask biopython pyyaml plotly
python app.py
```

---

🛠️ Technologies Used

Technology PDM SPO QuickNGS Purpose
Python ✅ ✅ ✅ Backend logic
Flask ✅ ✅ ✅ Web framework
SQLite ✅ ❌ ❌ Data storage
Bootstrap 5 ✅ ✅ ✅ Responsive UI
Chart.js ✅ ❌ ❌ PCR profile visualization
Matplotlib ❌ ✅ ❌ Binding map generation
Plotly ❌ ❌ ✅ Interactive coverage plots
Biopython ❌ ❌ ✅ FASTA/FASTQ handling
PyYAML ❌ ❌ ✅ Configuration files
Jinja2 ✅ ✅ ✅ Template rendering
Werkzeug ✅ ❌ ❌ Password hashing

---

👤 Default Login

Application URL Username Password
PDM http://127.0.0.1:5001 admin admin123
SPO http://127.0.0.1:5000 (no login required) —
QuickNGS http://127.0.0.1:5002 (no login required) —

⚠️ Important: Change the default admin password in PDM immediately after first login.

---

👨‍💻 Developer

Pourdad Panahi – Biotechnologist & Computational Biologist

· 🧪 18+ years of molecular biology laboratory experience
· 💻 Leveraging modern AI-assisted development to build scientific software
· 🧬 Domain expertise in PCR, primer design, and sequencing workflows
· 🎯 Seeking Bioinformatics / Computational Biology positions in Europe (Germany / Netherlands)

https://img.shields.io/badge/GitHub-pourdadp-black?logo=github
https://img.shields.io/badge/Portfolio-Website-blue

---

💡 From Wet Lab to Web App

As a molecular biologist who understands exactly what the lab needs, I used DeepSeek AI to accelerate the development process — translating years of domain knowledge into functional tools in record time. The AI assisted with code generation while I provided:

· Scientific requirements and validation
· Algorithm design (semi-global alignment, greedy coverage optimization)
· UI/UX decisions based on real lab workflows
· Testing against authentic experimental scenarios

The future belongs to scientists who can direct AI to build what they envision.

Built in just 10 days with DeepSeek AI as a collaborative development partner.

---

📄 License

This project is open-source and available under the MIT License.

```
MIT License

Copyright (c) 2025 Pourdad Panahi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

🙏 Acknowledgments

· Designed with molecular biologists and diagnostic labs in mind
· Inspired by real-world needs for primer inventory management, sequencing panel optimization, and NGS analysis
· Developed with the assistance of DeepSeek AI
· Built by a scientist, for scientists

---

📧 Contact & Support

If you encounter any issues or have feature requests:

· 📂 Open an issue on GitHub
· 📧 Contact the developer directly

---

<p align="center">
  <b>Made with 🧬 by a scientist who codes — powered by DeepSeek AI</b>
</p>
