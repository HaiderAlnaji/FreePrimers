#!/usr/bin/env bash
#
# Builds local BLAST nucleotide databases for specificity checking.
# Run this once at setup time (and periodically to refresh miRBase),
# not on every request — that's the whole point of the local backend.
#
# Usage:
#   PRIMERFORGE_BLAST_DB_DIR=/data/blastdb ./scripts/build_databases.sh
#
# Requires: BLAST+ (makeblastdb, on PATH or set
# PRIMERFORGE_MAKEBLASTDB_PATH), curl, gunzip.

set -euo pipefail

DB_DIR="${PRIMERFORGE_BLAST_DB_DIR:-/data/blastdb}"
MAKEBLASTDB="${PRIMERFORGE_MAKEBLASTDB_PATH:-makeblastdb}"
# Verified live as of June 2026 by directly fetching the URL.
# miRBase has changed its layout before; if this 404s again, check
# https://www.mirbase.org/download/ for the current location.
MIRBASE_RELEASE_URL="https://www.mirbase.org/download/mature.fa"

mkdir -p "$DB_DIR"
cd "$DB_DIR"

echo "==> Building '${DB_DIR}/mirbase_mature' from miRBase mature.fa"
echo "    (check https://www.mirbase.org/download/ for the current release URL"
echo "     if this fails — miRBase occasionally changes its download paths)"

curl -fsSL "$MIRBASE_RELEASE_URL" -o mature.fa.tmp
mv mature.fa.tmp mature.fa

# miRBase FASTA headers are RNA (U); BLAST nucleotide DBs expect DNA
# convention (T) for blastn against a DNA-style query. Translating U->T
# in the reference avoids spurious near-misses from base-letter
# mismatches that aren't real biological mismatches.
awk '/^>/{print; next} {gsub(/U/,"T"); gsub(/u/,"t"); print}' mature.fa > mirbase_mature.fa

"$MAKEBLASTDB" \
  -in mirbase_mature.fa \
  -dbtype nucl \
  -out mirbase_mature \
  -title "miRBase mature miRNAs (U->T converted)" \
  -parse_seqids

echo "==> Done. Database 'mirbase_mature' is ready in ${DB_DIR}."
echo ""
echo "To add a transcriptome or genome-region database for standard"
echo "qPCR/PCR primer specificity, download the relevant FASTA (e.g. a"
echo "RefSeq transcript set) and run:"
echo ""
echo "  ${MAKEBLASTDB} -in your_file.fa -dbtype nucl -out your_db_name -parse_seqids"
echo ""
echo "Then reference 'your_db_name' as the database parameter in"
echo "POST /specificity requests."
