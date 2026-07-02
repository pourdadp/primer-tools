#!/bin/bash
# =======================================================
# QuickNGS – One-Command Installer
# Installs all system tools and Python packages.
# Run: bash install.sh
# Powered by Pourdad Panahi – Built with DeepSeek AI
# =======================================================

set -e  # Exit on first error

echo "======================================="
echo "  🧬 QuickNGS Installer"
echo "======================================="
echo ""

# ---------- 1. Update package list ----------
echo "[1/5] Updating package list..."
sudo apt update -qq

# ---------- 2. Install bioinformatics tools ----------
echo "[2/5] Installing bioinformatics tools..."
sudo apt install -y -qq bwa samtools fastqc trimmomatic freebayes bcftools snpeff muscle 2>/dev/null || {
    echo "⚠️  Some tools could not be installed via apt."
    echo "    You may need to install them manually:"
    echo "    - snpeff:   conda install -c bioconda snpeff"
    echo "    - muscle:   conda install -c bioconda muscle"
    echo ""
}

# ---------- 3. Install system dependencies for WeasyPrint ----------
echo "[3/5] Installing WeasyPrint system dependencies..."
sudo apt install -y -qq libpango-1.0-0 libpangocairo-1.0-0 libffi-dev libcairo2-dev libgdk-pixbuf2.0-0 shared-mime-info 2>/dev/null || {
    echo "⚠️  Some WeasyPrint dependencies could not be installed."
    echo "    PDF report may not work without them."
}

# ---------- 4. Install Python packages ----------
echo "[4/5] Installing Python packages..."
pip install -r requirements.txt --quiet || pip3 install -r requirements.txt --quiet

# ---------- 5. Verify installations ----------
echo "[5/5] Verifying installations..."
echo ""

MISSING=""
for tool in bwa samtools fastqc trimmomatic freebayes bcftools snpeff muscle; do
    if command -v $tool &>/dev/null; then
        echo "  ✅ $tool"
    else
        echo "  ❌ $tool — not found"
        MISSING="$MISSING $tool"
    fi
done

# Check Python packages
python3 -c "import flask, Bio, yaml, plotly" 2>/dev/null && echo "  ✅ Python packages (flask, biopython, pyyaml, plotly)" || echo "  ❌ Some Python packages missing"

# WeasyPrint is optional
python3 -c "import weasyprint" 2>/dev/null && echo "  ✅ weasyprint" || echo "  ⚠️ weasyprint not installed (PDF report won't work)"

echo ""
if [ -n "$MISSING" ]; then
    echo "⚠️  Some tools are missing: $MISSING"
    echo "    QuickNGS will still run, but some features may be limited."
    echo "    Install missing tools with: sudo apt install $MISSING"
else
    echo "✅ All system tools installed successfully!"
fi

echo ""
echo "======================================="
echo "  🚀 Installation complete!"
echo "  Start QuickNGS: python3 app.py"
echo "  Open: http://127.0.0.1:5002"
echo "======================================="
