<!-- quickngs/CHANGELOG.md -->
# Changelog – QuickNGS

All notable changes to QuickNGS will be documented in this file.

## v2.0 (2025-06-28)

### Added
- **Smart Input Detection** – Automatically distinguishes NGS (FASTQ) from Sanger (AB1/FASTA/seq)
- **Sanger De Novo Assembly** – Three algorithms: Greedy OLC, Semi‑Global Alignment, De Bruijn Graph
- **Sanger Guided Assembly** – Align to reference with BWA, generate consensus with bcftools
- **Quality Trimming for AB1** – Phred‑based trimming with configurable threshold/window
- **Intelligent Clustering & Merging** – Reads grouped by overlap, assembled separately, then merged if possible
- **Translation Module** – 5 methods + 6 genetic codes (including Mycoplasma) with Biopython CAI
- **Read Direction Control** – Auto / Forward / Reverse per Sanger read
- **Tooltip Help Icons** – Clickable `?` icons for every technical parameter
- **Visual Contig Graphs** – Canvas‑based bars showing coverage and cluster merging
- **About & Help Pages** – Built‑in documentation
- **Delete Run Button** – In History page with confirmation
- **Footer with Version Badge** – v2.0 with navigation links
- **Dockerfile & docker‑compose.yml** – Containerized deployment

### Changed
- **Unified Smart Upload Form** – One page for both NGS and Sanger
- **Separated Translation Module** – `translation.py` extracted from `assembly.py`
- **Improved Trimmomatic Errors** – Clear messages for missing adapter files
- **merge_clusters Now Checks Reverse Complement** – Better scaffolding
- **/start_ngs Supports Custom Reference Upload** – User can provide their own FASTA
- **Automatic Storage Management** – Finds largest available partition

### Fixed
- **NGS file naming mismatch** – Frontend now sends `fastq_r1`/`fastq_r2` to `/start_ngs`
- **Flask/Werkzeug compatibility** – `flask>=2.3.0` in requirements.txt
- **Missing `base.html` template** – Added with footer signature
- **`cluster_reads` variable collision** – Renamed internal variable
- **`parse_uploaded_files` typo** – Fixed missing underscore
- **bcftools redundancy** – Removed duplicate check in `/assemble_guided`
- **Custom reference in NGS** – Now properly saved and indexed
- **favicon.ico 404** – Removed broken link, using emoji fallback

---

## v1.0 (2025-06-20)

### Added
- **Reference‑Based NGS Pipeline** – FastQC → Trimmomatic → BWA → Samtools → FreeBayes → SnpEff
- **Real Reference Genome Download** – hg38, hg19, mm10 with automatic indexing
- **Live Progress Polling** – Step‑by‑step status with animated progress bar
- **Run History** – All analyses saved with status tracking
- **File Download Endpoint** – Download VCF, BAM, and report files
- **Responsive Mobile‑First Design** – Bootstrap 5 interface
