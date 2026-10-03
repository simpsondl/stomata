#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STOMATA="${STOMATA:-$ROOT/build/src/stomata}"
"$STOMATA" --genome "$ROOT/tests/data/synthetic/test_small.fa" --spacer-file "$ROOT/tests/data/synthetic/test_spacers.txt" --threshold 1 --summary --summary-format tsv
