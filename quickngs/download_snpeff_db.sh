#!/bin/bash
# quickngs/download_snpeff_db.sh
# Run this once to download the GRCh38.99 SnpEff database into the volume.
# Usage: docker-compose run quickngs bash download_snpeff_db.sh

echo "Downloading SnpEff database GRCh38.99..."
java -jar /usr/share/java/snpeff.jar download -v GRCh38.99
echo "Database downloaded to $SNPEFF_DATA"
