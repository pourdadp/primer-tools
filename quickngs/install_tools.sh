#!/bin/bash
# quickngs/install_tools.sh
# Install all required system tools for QuickNGS on Ubuntu/Debian.
# Run: bash install_tools.sh

echo "Installing bioinformatics tools for QuickNGS..."

sudo apt update
sudo apt install -y bwa samtools freebayes fastqc trimmomatic bcftools snpeff

echo ""
echo "✅ All tools installed!"
echo ""
echo "Next steps:"
echo "  pip install -r requirements.txt"
echo "  python3 app.py"
