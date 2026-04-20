#!/bin/bash
# Demonstration: why batch processing is fast
# This simulates searching 5 spacers by loading genome once in memory

set -euo pipefail

CYSTIDIA_BIN="./build/src/cystidia"
GENOME="path/to/genomes/GRCh38.p14/hg38.p14.canonical.fa"
THRESHOLD=4

declare -a SPACERS=(
    "GAGTCCGAGCAGAAGAAGAA"
    "GGAATCCCTTCTGCAGCACC"
    "GACCCCCTCCACCCCGCCTC"
    "GGTGAGTGAGTGTGTGCGTG"
    "GAGTGAGGCTCCCCTGACCC"
)

echo "Batch Search Demo"
echo "================="
echo ""
echo "Current approach: Load genome separately for each spacer"
echo "Better approach: Load once, search N times (not yet implemented)"
echo ""

# Simulate current approach (5x genome loading)
echo "Simulated time for 5 spacers (current approach):"
echo "  Genome loading: 53.7s × 5 = 268.5s"
echo "  Search time:     4.2s × 5 =  21.0s"
echo "  Total:                      289.5s (~4.8 minutes)"
echo ""

# Proposed batch approach
echo "Proposed time with batch loading (not yet implemented):"
echo "  Genome loading: 53.7s × 1 = 53.7s"
echo "  Search time:     4.2s × 5 = 21.0s"
echo "  Total:                      74.7s (~1.2 minutes)"
echo ""
echo "Speedup: 3.9x faster for 5 spacers"
echo "         5.0x faster for 10 spacers"
echo "        10.0x faster for 50 spacers"
echo ""
echo "This requires adding batch mode to cystidia CLI (future enhancement)"
