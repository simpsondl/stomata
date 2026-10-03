#!/usr/bin/env bash
# Combine supplied UCSC-named chromosome FASTAs; never guess/download an assembly.
set -euo pipefail
if [[ $# != 2 ]]; then
  echo "Usage: $0 CHROMOSOME_DIRECTORY OUTPUT.fa" >&2
  exit 2
fi
python3 - "$1" "$2" <<'PY'
import sys
from pathlib import Path
source, destination = map(Path, sys.argv[1:])
chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
inputs = [source / f"{name}.fa" for name in chromosomes]
missing = [str(p) for p in inputs if not p.is_file()]
if missing:
    raise SystemExit("Missing chromosomes: " + ", ".join(missing))
with destination.open("xb") as output:
    for path in inputs:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                output.write(chunk)
        output.write(b"\n")
print(destination)
PY
