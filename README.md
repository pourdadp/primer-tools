# primer-tools
🧬 Primer Tools Suite

A comprehensive, open‑source toolkit for PCR primer panel optimization and primer database management. Designed for molecular biologists and bioinformaticians, this suite helps you select existing primers for full‑length targeted sequencing and manage your primer inventory with an advanced web interface.

📦 Repository Structure

```
primer-tools/
├── sequencing_panel_optimizer/   # Optimize primer panels for targeted sequencing
│   ├── app.py
│   ├── requirements.txt
│   ├── install.bat
│   ├── run.bat
│   ├── templates/
│   └── scripts/
└── primer_database_manager/      # Manage primer inventory, probes, and PCR programs
    ├── app.py
    ├── requirements.txt
    ├── install.bat
    ├── run.bat
    ├── templates/
    └── static/
```

🧬 1. Sequencing Panel Optimizer (SPO)

Select the best combination of existing primers to generate overlapping amplicons for full‑length targeted sequencing.

✨ Key Features

· Sequence Analysis: Input your target DNA sequence.
· Primer Loading: Manually paste primers or load them directly from the Primer Database Manager.
· Automatic Direction Detection: Detects if primers bind in Forward/Reverse orientation.
· IUPAC & Inosine Support: Fully supports degenerate primers (e.g., R, Y, N, I).
· Filtering & Validation:
  · Rejects primers that bind to multiple sites or in both directions.
  · Validates amplicon length (user‑defined min/max).
  · Checks primer Tm compatibility and dimer formation.
· Coverage Optimization: Uses a greedy algorithm to find the optimal set of primers covering the whole target.
· Visualization: Generates a graphical map of primer binding sites (PNG).
· Outputs: Downloadable JSON reports (binding results, valid pairs, optimal panel).
· History Tracking: Logs all runs with status and parameters.

🗄️ 2. Primer Database Manager (PDM)

A comprehensive web‑based management system for primers, probes, and PCR programs.

✨ Key Features

👥 User Management

· Multi‑user roles: Admin, Editor, Viewer.
· Admin can create/disable users and change roles.
· Secure password hashing (Werkzeug).
· Password Reset: Users request a reset; Admin receives a notification and sets a new password (no email required).

🧬 Primer Management

· Full CRUD (Create, Read, Update, Delete) for primers.
· Store detailed metadata: Gene, Organism, Strain/Serotype, Tm (estimated/experimental), Amplicon length, Binding region, and more.
· Editing Lock: Prevents simultaneous edits (one user at a time). Admin can break the lock.
· Custom Fields: Users can add dynamic key‑value fields to any primer.

🔬 Probe Management

· Store probes with separate fields for Reporter and Quencher.
· Supports probe type, modifications, and notes.

🧪 PCR Program Designer

· Store multiple PCR programs per primer pair.
· Default Program selection.
· Define steps: Denaturation, Annealing, Extension, Reverse Transcription, Melt Curve, Hold.
· Mark steps for Real‑Time data acquisition (Read Step).
· Schematic Visualization: Generates a graphical timeline of the PCR protocol.
· Linear Text Output: Displays the program as a step‑by‑step text protocol.

📋 Multiplex Panels

· Group multiple primer pairs (and probes) into single reaction panels.

💾 Backup & Restore

· One‑click backup of the entire SQLite database.
· Auto‑cleanup of backups older than 30 days.
· Restore any previous backup directly from the interface.

🖨️ Print Compatibility

· All pages (except login) are optimized for A4 printing.
· Print button available; non‑essential elements are hidden automatically.

🔗 Integration Between Tools

The Sequencing Panel Optimizer can directly import primers from the Primer Database Manager.

· In SPO, click the "📂 Load from Database" button.
· SPO will fetch active primers (name and sequence) from the shared primers.db file.
· This creates a seamless workflow: Manage → Optimize → Sequence.

⚙️ Installation & Setup

Prerequisites

· Python 3.8 or higher (download from python.org).
· Important: During Python installation, check "Add Python to PATH".

Step 1: Clone or Download

```bash
git clone https://github.com/pourdadp/primer-tools.git
cd primer-tools
```

Step 2: Install Dependencies

Each tool has its own install.bat (Windows) or you can use pip manually.

For Sequencing Panel Optimizer:

```bash
cd sequencing_panel_optimizer
install.bat
```

Or manually:

```bash
pip install -r requirements.txt
```

For Primer Database Manager:

```bash
cd ../primer_database_manager
install.bat
```

Or manually:

```bash
pip install -r requirements.txt
```

💡 The install.bat script uses a Chinese mirror (Tsinghua) for faster downloads, especially within Iran/China. If you prefer the default PyPI, just run pip install -r requirements.txt.

Step 3: Run the Applications

Each tool runs on a different port to avoid conflicts.

· SPO (Port 5000):
  ```bash
  cd sequencing_panel_optimizer
  run.bat
  ```
  Access: http://127.0.0.1:5000
· PDM (Port 5001):
  ```bash
  cd primer_database_manager
  run.bat
  ```
  Access: http://127.0.0.1:5001

Default Admin Login (PDM):

· Username: admin
· Password: admin123

🔧 Manual Installation (for Linux/macOS)

If you are on Linux/macOS, simply use pip:

```bash
# For SPO
pip install flask pyyaml matplotlib tqdm
python app.py

# For PDM
pip install flask werkzeug
python app.py
```

📸 Screenshots (Placeholder)

(Add screenshots of the Dashboard, Primer Detail, PCR Program Editor, and SPO Results here)

🛠️ Technologies Used

· Backend: Python, Flask
· Database: SQLite (embedded)
· Frontend: Bootstrap 5, Jinja2
· Visualization: Matplotlib (for SPO charts)
· Authentication: Werkzeug (password hashing)

📜 License

This project is open‑source and available under the MIT License.

👨‍💻 Developer

Pourdad Panahi
GitHub: https://github.com/pourdadp
Project Page: https://github.com/pourdadp/pourdadp.github.io

📧 Contact & Support

If you encounter any issues or have feature requests, please open an issue on GitHub or contact the developer directly.

🙏 Acknowledgments

Designed with molecular biologists and diagnostic labs in mind. Inspired by real‑world needs for PCR optimization and primer inventory management.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
