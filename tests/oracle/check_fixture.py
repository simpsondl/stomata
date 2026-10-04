#!/usr/bin/env python3
"""Compare public CLI output with the independent DP oracle on the bundled fixture."""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from fuzz_oracle_vs_stomata import compare_files

ROOT = Path(__file__).resolve().parents[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", type=Path, default=ROOT / "build/src/stomata")
    ap.add_argument("--cpu-only", action="store_true")
    ap.add_argument("--output-dir", type=Path, default=ROOT / "tests/oracle/results")
    args = ap.parse_args()
    binary = str(args.binary.resolve())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    genome = ROOT / "tests/regression/oracle_genome.fa"
    spacers = ROOT / "tests/regression/oracle_fixture.tsv"
    version = subprocess.check_output([binary, "--version"], text=True)
    (args.output_dir / "version.txt").write_text(version)
    failed = False
    for mode in ("hamming", "levenshtein"):
        for pam in ("NGG", ""):
            label = f"{mode}-{pam or 'any'}"
            expected = args.output_dir / f"{label}-oracle.tsv"
            actual = args.output_dir / f"{label}-stomata.tsv"
            subprocess.run([sys.executable, str(Path(__file__).with_name("oracle_lev_search.py")),
                            "--spacers", str(spacers), "--genome", str(genome),
                            "--threshold", "3", "--pam", pam, "--distance-mode", mode,
                            "--out", str(expected)], check=True)
            cmd = [binary, "--spacer-file", str(spacers), "--genome", str(genome),
                   "--threshold", "3", "--distance-mode", mode, "--output", str(actual)]
            if pam:
                cmd += ["--pam", pam]
            if args.cpu_only:
                cmd += ["--cpu-only"]
            subprocess.run(cmd, check=True)
            diff = compare_files(actual, expected)
            (args.output_dir / f"{label}.diff.txt").write_text(diff)
            print(f"{label}: {'FAIL' if diff else 'PASS'}", flush=True)
            failed |= bool(diff)
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
