#!/bin/bash
# Prepare canonical chromosome hg38 reference for Stomata validation testing
# Combines chr1-22, chrX, chrY, chrM into single multi-FASTA and compresses

set -euo pipefail

# Paths
INDIV_CHR_DIR="path/to/genomes/GRCh38.p14/indiv_chr"
OUTPUT_DIR="path/to/genomes/GRCh38.p14"
OUTPUT_FASTA="${OUTPUT_DIR}/hg38.p14.canonical.fa"
OUTPUT_GZIP="${OUTPUT_FASTA}.gz"

# Canonical chromosomes in order
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)

echo "Combining canonical chromosomes into ${OUTPUT_FASTA}..."

# Remove old files if they exist
rm -f "${OUTPUT_FASTA}" "${OUTPUT_GZIP}"

# Concatenate all chromosomes
for chr in "${CHROMOSOMES[@]}"; do
    chr_file="${INDIV_CHR_DIR}/${chr}.fa"
    if [[ -f "${chr_file}" ]]; then
        echo "  Adding ${chr}..."
        cat "${chr_file}" >> "${OUTPUT_FASTA}"
    else
        echo "  WARNING: ${chr_file} not found, skipping"
    fi
done

echo "Compressing with gzip..."
gzip -c "${OUTPUT_FASTA}" > "${OUTPUT_GZIP}"

echo "Creating FASTA index..."
samtools faidx "${OUTPUT_GZIP}"

# Report statistics
echo ""
echo "=== Statistics ==="
echo "Uncompressed size: $(du -h ${OUTPUT_FASTA} | cut -f1)"
echo "Compressed size:   $(du -h ${OUTPUT_GZIP} | cut -f1)"
echo "Number of sequences: $(grep -c '^>' ${OUTPUT_FASTA})"
echo ""
echo "Output files:"
echo "  ${OUTPUT_FASTA}"
echo "  ${OUTPUT_GZIP}"
echo "  ${OUTPUT_GZIP}.fai"
echo ""
echo "Done!"
