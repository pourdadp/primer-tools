<!-- quickngs/CHANGELOG.md -->
# Changelog – QuickNGS

## v2.0 (2025-06-27)
### Added
- Smart detection of input type (NGS vs Sanger)
- Sanger De Novo assembly (Greedy OLC, Semi‑Global Alignment, De Bruijn Graph)
- Sanger Reference‑Guided assembly (BWA + consensus)
- Quality trimming for AB1 files (Phred scores, configurable threshold/window)
- Read direction control (Auto / Forward / Reverse)
- Intelligent read clustering and cluster merging (scaffolding)
- Tooltip help icons for technical parameters
- About, Help, and History pages
- Footer with version badge and navigation links
- Dockerfile for containerized execution

### Changed
- Unified smart upload form replaces separate tabs
- Improved error messages (user‑friendly, no technical jargon)
- Automatic storage management (finds largest available partition)

### Fixed
- NGS file naming mismatch between frontend and backend
- Flask/Werkzeug compatibility issue (flask>=2.3.0)
- Missing base.html template

## v1.0 (2025-06-20)
### Added
- Reference‑Based NGS pipeline (FastQC → Trimmomatic → BWA → Samtools → FreeBayes → SnpEff)
- Real reference genome download (hg38, hg19, mm10)
- Live progress polling with step badges
- Run history with status tracking
- File download endpoint
- Responsive mobile‑first design with Bootstrap 5
